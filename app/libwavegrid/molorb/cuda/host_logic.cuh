/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#pragma once
#include "kernel.cuh"
#include "device_params.cuh"

// amount of shared memory set aside for nEig accumulators
constexpr float SHARED_MEM_FACTOR = 0.95f;
// max output array share of free global memory
constexpr float GLOBAL_MEM_FACTOR = 0.80f;
// Threads per block, multiple of warp size 32
constexpr int block_size = 256;

struct GpuLaunchConfig {
    int    deviceId;
    int    z_start;              // Starting Z-slice for this GPU
    int    z_count;              // Number of Z-slices for this GPU
    int    z_per_batch;          // Number of Z-slices to process per kernel launch
    int    nEig_per_pass;        // Number of eigenstates to accumulate in shared memory
    size_t shared_mem_for_pass;  // Amount of shared memory per block for the accumulators
};


/* 
 * Wraps CUDA event calls to provide a simple way to time GPU execution.
 * Use elapsed() to retrieve the accumulated time without stopping the timer.
 * The timer is initially stopped unless startNow=true is passed to the constructor.
 */
class GpuTimer {
public:
    explicit GpuTimer(bool startNow = false) : _accumulated_ms(0.0f), _running(false) {
        CHECK_CUDA(cudaEventCreate(&_startEvent));
        CHECK_CUDA(cudaEventCreate(&_stopEvent));
        if (startNow) start();
        
    }

    ~GpuTimer() {
        cudaEventDestroy(_startEvent);
        cudaEventDestroy(_stopEvent);
    }

    void start() {
        if (!_running) {
            CHECK_CUDA(cudaEventRecord(_startEvent));
            _running = true;
        }
    }

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
GpuLaunchConfig setup_gpu_config(int deviceId, int numGpus, const GridParams* grid, const CalculationParams* calc);
void dispatchKernel(const DeviceKernelParams* params, bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg,
    bool useRadialLut, int grid_size, size_t shared_mem_for_pass);
void copyD2H(void* d_src_ptr, void* h_dest_ptr, int nPointsX, int nPointsY, int nPointsZ, int z_per_batch,
    int z_offset_global, const CalculationParams* calc);
