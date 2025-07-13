#ifndef UTILS_CUH_
#define UTILS_CUH_

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





// In a new header, e.g., "cuda_utils.cuh"
#include <cuda_runtime.h>
#include <stdexcept>




// A simple RAII wrapper for device memory.
// Use get() to retrieve the raw pointer.
//
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
        CHECK_CUDA(cudaMemcpy(d_ptr_, host_ptr, count * sizeof(T), cudaMemcpyHostToDevice));
    }

    // Destructor automatically frees the memory
    ~DeviceBuffer() {
        if (d_ptr_) {
            cudaFree(d_ptr_);
        }
    }

    // Disable copy semantics
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    // Enable move semantics
    DeviceBuffer(DeviceBuffer&& other) noexcept : d_ptr_(other.d_ptr_), count_(other.count_) {
        other.d_ptr_ = nullptr;
        other.count_ = 0;
    }
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (d_ptr_) cudaFree(d_ptr_);
            d_ptr_ = other.d_ptr_;
            count_ = other.count_;
            other.d_ptr_ = nullptr;
            other.count_ = 0;
        }
        return *this;
    }
    
    
    T* get() { return d_ptr_; }
    const T* get() const { return d_ptr_; }
    size_t size() const { return count_; }

    void copy_to_host(T* host_ptr, size_t count_to_copy) const {
        if (count_to_copy > count_) {
             fprintf(stderr, "Error: trying to copy more elements than buffer contains.\n");
             exit(EXIT_FAILURE);
        }
        CHECK_CUDA(cudaMemcpy(host_ptr, d_ptr_, count_to_copy * sizeof(T), cudaMemcpyDeviceToHost));
    }

private:
    void allocate(size_t count) {
        count_ = count;
        if (count > 0) {
            CHECK_CUDA(cudaMalloc(&d_ptr_, count * sizeof(T)));
        }
    }
    T* d_ptr_ = nullptr;
    size_t count_ = 0;
};

#endif // UTILS_CUH_
