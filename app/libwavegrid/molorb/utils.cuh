/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/

#ifndef UTILS_CUH_
#define UTILS_CUH_
#include <cuda_runtime.h>
#include <string>
#include <stdexcept>

// Helper macro for robust CUDA calls
#define CHECK_CUDA(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while (0)


// Helper macros for column-major (Fortran-style) index calculations
#define IDX2F(i, j, lda) ((j) * (size_t)(lda) + (i))
#define IDX3F(i, j, k, lda, ldb) (((k) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))
#define IDX4F(i, j, k, l, lda, ldb, ldc) ((((l) * (size_t)(ldc) + (k)) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))


/*
 * A simple RAII wrapper for device memory.
 * Use get() to retrieve the raw pointer. 
 */
template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;

    explicit DeviceBuffer(size_t count) {
        allocate(count);
    }

    // Allocate and copy from host
    DeviceBuffer(const T* host_ptr, size_t count) {
        allocate(count);
        copy_to_device(host_ptr, count);
    }

    // Destructor automatically frees the memory
    ~DeviceBuffer() {
        deallocate();
    }

    void deallocate() {
        if (_devicePtr) {
            CHECK_CUDA(cudaFree(_devicePtr));
            _devicePtr = nullptr;
            _count = 0;
        }
    }
    
    // Assign a new size and copy from host
    void assign(const T* host_ptr, size_t count) {
        deallocate();
        allocate(count);
        copy_to_device(host_ptr, count);
    }

    // Disable copy semantics
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    // Enable move semantics
    DeviceBuffer(DeviceBuffer&& other) noexcept : _devicePtr(other._devicePtr), _count(other._count) {
        other._devicePtr = nullptr;
        other._count = 0;
    }
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (_devicePtr) cudaFree(_devicePtr);
            _devicePtr = other._devicePtr;
            _count = other._count;
            other._devicePtr = nullptr;
            other._count = 0;
        }
        return *this;
    }
    
    
    T* get() { return _devicePtr; }
    const T* get() const { return _devicePtr; }
    size_t size() const { return _count; }

    void copy_to_host(T* host_ptr, size_t count_to_copy) const {
        if (count_to_copy > _count) {
             fprintf(stderr, "Error: trying to copy more elements than buffer contains.\n");
             exit(EXIT_FAILURE);
        }
        CHECK_CUDA(cudaMemcpy(host_ptr, _devicePtr, count_to_copy * sizeof(T), cudaMemcpyDeviceToHost));
    }

    void copy_to_device(const T* host_ptr, size_t count_to_copy) {
        if (!_devicePtr) {
            fprintf(stderr, "Error: device pointer is null. Cannot copy to device.\n");
            exit(EXIT_FAILURE);
        }
        if (count_to_copy > _count) {
            fprintf(stderr, "Error: trying to copy more elements than buffer contains.\n");
            exit(EXIT_FAILURE);
        }
        CHECK_CUDA(cudaMemcpy(_devicePtr, host_ptr, count_to_copy * sizeof(T), cudaMemcpyHostToDevice));
    }

private:
    void allocate(size_t count) {
        _count = count;
        if (count > 0) {
            CHECK_CUDA(cudaMalloc(&_devicePtr, count * sizeof(T)));
        } else {
            _devicePtr = nullptr;
        }
    }
    T* _devicePtr = nullptr;
    size_t _count = 0;
};



/*
 * We implement the LUT as a 2D texture with a single float channel.
 * The radial Functions are stored row-wise.
 * We assume identical radial grids for all STOs.
 * GPUs have dedicated hardware for texture access & interpolation.
 */
class GpuLutTexture {
public:
    GpuLutTexture(const double* lutData, int nPoints, int nStos) {
        // Convert the Fortran passed doubles to floats
        size_t totalValues = (size_t)nStos * nPoints;
        std::vector<float> lutFloats(totalValues);
        for (size_t i = 0; i < totalValues; ++i) 
            lutFloats[i] = static_cast<float>(lutData[i]);
        

        // Allocate memory on device
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        CHECK_CUDA(cudaMallocArray(&_lutArray, &channelDesc, nPoints, nStos, 0));

        // Copy data to array
        CHECK_CUDA(cudaMemcpy2DToArray(
            _lutArray,                  // dst array
            0, 0,                       // no offset in dst
            lutFloats.data(),           // src pointer
            nPoints * sizeof(float),    // src pitch (for alignment, bytes to next row)
            nPoints * sizeof(float),    // width in bytes
            nStos,                      // height (number of cached stos)
            cudaMemcpyHostToDevice
        ));

        // Prepare texture object properties
        cudaResourceDesc resDesc{};
        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = _lutArray;

        cudaTextureDesc texDesc{};
        // OOB access clamped to edge values
        // (Should not occur if sto_cutoffs are set correctly)
        texDesc.addressMode[0] = cudaAddressModeClamp; 
        texDesc.addressMode[1] = cudaAddressModeClamp;
        // Enable linear interpolation
        texDesc.filterMode = cudaFilterModeLinear;
        // Do not normalize the lut values
        texDesc.readMode = cudaReadModeElementType;
        // Access using texel coords, i.e. [0, N-1]
        // Imagine a pixel, add 0.5f to get the center.
        texDesc.normalizedCoords = 0;
        
        // Create texture object
        CHECK_CUDA(cudaCreateTextureObject(&_textureObject, &resDesc, &texDesc, nullptr));
    }

    ~GpuLutTexture() {
        if (_textureObject) cudaDestroyTextureObject(_textureObject);
        if (_lutArray) cudaFreeArray(_lutArray);
    }

    // Disable copy
    GpuLutTexture(const GpuLutTexture&) = delete;
    GpuLutTexture& operator=(const GpuLutTexture&) = delete;

    cudaTextureObject_t get() const { return _textureObject; }

private:
    cudaArray_t _lutArray = nullptr;
    cudaTextureObject_t _textureObject = 0;
};






#endif // UTILS_CUH_
