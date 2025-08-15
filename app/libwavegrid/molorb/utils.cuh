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


// Helper macros for column-major (Fortran-style) indexing.
#define IDX2F(i, j, lda) ((j) * (size_t)(lda) + (i))
#define IDX3F(i, j, k, lda, ldb) (((k) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))
#define IDX4F(i, j, k, l, lda, ldb, ldc) ((((l) * (size_t)(ldc) + (k)) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))


// A simple RAII wrapper for device memory.
// Use get() to retrieve the raw pointer.
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

#endif // UTILS_CUH_
