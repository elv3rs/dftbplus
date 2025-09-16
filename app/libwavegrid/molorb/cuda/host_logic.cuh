/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#pragma once
#include "kernel.cuh"
#include "device_params.cuh"

#include <assert.h>
#include <algorithm>
#include <cstdio>

// amount of shared memory set aside for nEig accumulators
constexpr float SHARED_MEM_FACTOR = 0.95f;
// max output array share of free global memory
constexpr float GLOBAL_MEM_FACTOR = 0.80f;
// Threads per block, multiple of warp size 32
constexpr int block_size = 256;

struct GpuLaunchConfig {
    int    deviceId;            // CUDA device ID
    int    z_count;             // Number of Z-slices for this GPU
    int    z_start;             // Starting Z-slice for this GPU
    int    z_per_batch;         // Number of Z-slices to process per kernel launch
    int    nEig_per_pass;       // Number of eigenstates to accumulate in shared memory
    size_t shared_mem_for_pass; // Amount of shared memory per block for the accumulators
    // Kernel template parameters
    const bool isRealInput;
    const bool calcAtomicDensity;
    const bool calcTotalChrg;
    const bool useRadialLut;

    /** 
     * @brief Splits work among GPUs using Z-slices, and calculate gpu launch parameters based on available memory.
     *
     * @param deviceId      The CUDA GPU device ID.
     * @param numGpus       Total number of GPUs used for splitting the work.
     * @param grid          Pointer to the GridParams structure containing grid dimensions.
     * @param calc          Pointer to the CalculationParams structure containing calculation flags.
     * @param useRadialLut  Flag indicating whether to use radial lookup table for sto evaluation.
     * @return A GpuLaunchConfig structure with the calculated parameters for the specified GPU.
     */
    GpuLaunchConfig(int deviceId, int numGpus, const GridParams* grid, const CalculationParams* calc, bool useRadialLut):
        deviceId(deviceId),
        isRealInput(calc->isRealInput),
        calcAtomicDensity(calc->calcAtomicDensity),
        calcTotalChrg(calc->calcTotalChrg),
        useRadialLut(useRadialLut) {
        // Evenly split Z-slices among GPUs
        z_count = grid->nPointsZ / numGpus;
        z_start = deviceId * z_count;
        // Handle uneven Z-slice count: Last GPU takes remaining slices
        if (deviceId == numGpus - 1) z_count = grid->nPointsZ - z_start;

        // Query available memory sizes with safety margins
        cudaDeviceProp prop;
        size_t         free_global_mem, total_global_mem;
        CHECK_CUDA(cudaMemGetInfo(&free_global_mem, &total_global_mem));
        CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));
        size_t available_shared = prop.sharedMemPerBlock * SHARED_MEM_FACTOR;
        size_t available_global = static_cast<size_t>(free_global_mem * GLOBAL_MEM_FACTOR);

        size_t accumulator_number_size = calc->isRealInput ? sizeof(double) : sizeof(complexd);
        size_t output_number_size      = calc->isRealOutput ? sizeof(double) : sizeof(complexd);

        // Determine number of eigenstates that fit into shared memory
        nEig_per_pass       = available_shared / (block_size * accumulator_number_size);
        nEig_per_pass       = std::min(calc->nEigIn, std::max(1, nEig_per_pass));  // clamp to [1, nEigIn]
        shared_mem_for_pass = (size_t)nEig_per_pass * block_size * accumulator_number_size;

        // Determine max Z-slices that can fit in available (global) memory
        size_t bytes_per_slice = (size_t)grid->nPointsX * grid->nPointsY * calc->nEigOut * output_number_size;
        z_per_batch     = std::min(z_count, (int)available_global / (int)bytes_per_slice);
        z_per_batch     = std::max(1, z_per_batch);  // at least 1

        // Debug output
        if (deviceId == 0 && debug) {
            printf("\n--- GPU %d Configuration ---\n", deviceId);
            printf("  Z-slice workload: %d (from index %d to %d)\n", z_count, z_start,
                z_start + z_count - 1);
            printf("  Block size: %d threads, %zub shared mem per block, %d eigs of %d per pass\n", block_size,
                shared_mem_for_pass, nEig_per_pass, calc->nEigIn);
            size_t total_size_valueOut =
                (size_t)grid->nPointsX * grid->nPointsY * grid->nPointsZ * calc->nEigOut * sizeof(double);
            if (!calc->isRealOutput) total_size_valueOut *= 2;
            printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n", free_global_mem / 1e9,
                grid->nPointsX, grid->nPointsY, grid->nPointsZ, calc->nEigOut, total_size_valueOut / 1e9);
            printf("  Processing Z-slices in batches of %d\n", z_per_batch);
        }
    }


};





/**
 * @brief A simple GPU timer using CUDA events.
 *
 * Wraps CUDA event calls to provide a simple way to time GPU execution.
 * Use elapsed() to retrieve the accumulated time without stopping the timer.
 * The timer is initially stopped unless startNow=true is passed to the constructor.
 */
class GpuTimer {
public:

    /**
     * @brief Constructor. Creates the CUDA events, optionally starts the timer.
     * @param startNow If true, starts the timer immediately.
     */
    explicit GpuTimer(bool startNow = false) : _accumulated_ms(0.0f), _running(false) {
        CHECK_CUDA(cudaEventCreate(&_startEvent));
        CHECK_CUDA(cudaEventCreate(&_stopEvent));
        if (startNow) start();
        
    }

    ~GpuTimer() {
        cudaEventDestroy(_startEvent);
        cudaEventDestroy(_stopEvent);
    }
   
    /**
     * @brief Starts the timer. If already running, does nothing.
     */
    void start() {
        if (!_running) {
            CHECK_CUDA(cudaEventRecord(_startEvent));
            _running = true;
        }
    }
    /**
     * @brief Returns the total elapsed time, without stopping the timer.
     * If the timer is running, it records the current time and updates the accumulated time.
     * @return Elapsed time in milliseconds.
     */
    float elapsed_ms() {
        if (_running) {
            CHECK_CUDA(cudaEventRecord(_stopEvent));
            CHECK_CUDA(cudaEventSynchronize(_stopEvent));
            float elapsed;
            CHECK_CUDA(cudaEventElapsedTime(&elapsed, _startEvent, _stopEvent));
            _accumulated_ms += elapsed;
        }
        return _accumulated_ms;
    }

    /**
     * @brief Stops the timer and returns the total elapsed time.
     * If the timer is not running, simply returns the accumulated time.
     * @return Total elapsed time in milliseconds.
     */
    float stop() {
        _accumulated_ms = elapsed_ms();
        _running = false;
        return _accumulated_ms;
    }

    // Disallow copy and assign
    GpuTimer(const GpuTimer&) = delete;
    GpuTimer& operator=(const GpuTimer&) = delete;

private:
    cudaEvent_t _startEvent{};
    cudaEvent_t _stopEvent{};
    float _accumulated_ms;
    bool _running;
};


// Function declarations
void copyD2H(void* d_src_ptr, void* h_dest_ptr, int nPointsX, int nPointsY, int nPointsZ, int z_per_batch,
    int z_offset_global, const CalculationParams* calc);
