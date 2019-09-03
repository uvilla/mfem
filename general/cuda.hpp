// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef MFEM_CUDA_HPP
#define MFEM_CUDA_HPP

#include "../config/config.hpp"
#include "error.hpp"

#ifdef MFEM_USE_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#endif

// CUDA block size used by MFEM.
#define MFEM_CUDA_BLOCKS 256

#ifdef MFEM_USE_CUDA
#define MFEM_DEVICE __device__
#define MFEM_HOST_DEVICE __host__ __device__
// Define a CUDA error check macro, MFEM_GPU_CHECK(x), where x returns/is of
// type 'cudaError_t'. This macro evaluates 'x' and raises an error if the
// result is not cudaSuccess.
#define MFEM_GPU_CHECK(x) \
   do \
   { \
      cudaError_t err = (x); \
      if (err != cudaSuccess) \
      { \
         mfem_cuda_error(err, #x, _MFEM_FUNC_NAME, __FILE__, __LINE__); \
      } \
   } \
   while (0)
#else
#define MFEM_DEVICE
#define MFEM_HOST_DEVICE
#endif // MFEM_USE_CUDA

#if defined(__CUDA_ARCH__)
#define INNER_LAMBDA [=] __device__
#else
#define INNER_LAMBDA [=]
#endif

//kernel helper...
template<typename LAMBDA>
#if defined(__CUDA_ARCH__)
__device__
#endif
inline void team_loop(int Nx, LAMBDA &&body)
{
#if defined(__CUDA_ARCH__)
  for(int tz=threadIdx.z; tz<1; tz+=blockDim.z)
  for(int ty=threadIdx.y; ty<1; ty+=blockDim.y)
  for(int tx=threadIdx.x; tx<Nx; tx+=blockDim.x)
  {
    body(tx);
  }
#else
  for(int tx=0; tx<Nx; tx++)
  {
    body(tx);
  }
#endif
}

template<typename LAMBDA>
#if defined(__CUDA_ARCH__)
__device__
#endif
inline void team_loop(int Nx, int Ny, LAMBDA &&body)
{
#if defined(__CUDA_ARCH__)
  for(int tz=threadIdx.z; tz<1; tz+=blockDim.z)
  for(int ty=threadIdx.y; ty<Ny; ty+=blockDim.y)
  for(int tx=threadIdx.x; tx<Nx; tx+=blockDim.x)
  {
    body(tx,ty);
  }
#else
  for(int ty=0; ty<Ny; ty++)
  for(int tx=0; tx<Nx; tx++)
  {
    body(tx,ty);
  }
#endif
}

template<typename LAMBDA>
#if defined(__CUDA_ARCH__)
__device__
#endif
inline void team_loop(int Nx, int Ny, int Nz, LAMBDA &&body)
{
#if defined(__CUDA_ARCH__)
  for(int tz=threadIdx.z; tz<Nz; tz+=blockDim.z)
  for(int ty=threadIdx.y; ty<Ny; ty+=blockDim.y)
  for(int tx=threadIdx.x; tx<Nx; tx+=blockDim.x)
  {
    body(tx, ty, tz);
  }
#else
  for(int tz=0; tz<Nz; tz++)
  for(int ty=0; ty<Ny; ty++)
  for(int tx=0; tx<Nx; tx++)
  {
    body(tx, ty, tz);
  }
#endif
}


// Define the MFEM inner threading macros
#if defined(MFEM_USE_CUDA) && defined(__CUDA_ARCH__)
#define MFEM_SHARED __shared__
#define MFEM_SYNC_THREAD __syncthreads()
#define MFEM_THREAD_ID(k) threadIdx.k
#define MFEM_THREAD_SIZE(k) blockDim.k
#define MFEM_FOREACH_THREAD(i,k,N) for(int i=threadIdx.k; i<N; i+=blockDim.k)
#else
#define MFEM_SHARED
#define MFEM_SYNC_THREAD
#define MFEM_THREAD_ID(k) 0
#define MFEM_THREAD_SIZE(k) 1
#define MFEM_FOREACH_THREAD(i,k,N) for(int i=0; i<N; i++)
#endif

namespace mfem
{

#ifdef MFEM_USE_CUDA
// Function used by the macro MFEM_GPU_CHECK.
void mfem_cuda_error(cudaError_t err, const char *expr, const char *func,
                     const char *file, int line);
#endif

/// Allocates device memory
void* CuMemAlloc(void **d_ptr, size_t bytes);

/// Frees device memory
void* CuMemFree(void *d_ptr);

/// Copies memory from Host to Device
void* CuMemcpyHtoD(void *d_dst, const void *h_src, size_t bytes);

/// Copies memory from Host to Device
void* CuMemcpyHtoDAsync(void *d_dst, const void *h_src, size_t bytes);

/// Copies memory from Device to Device
void* CuMemcpyDtoD(void *d_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Device
void* CuMemcpyDtoDAsync(void *d_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Host
void* CuMemcpyDtoH(void *h_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Host
void* CuMemcpyDtoHAsync(void *h_dst, const void *d_src, size_t bytes);

} // namespace mfem

#endif // MFEM_CUDA_HPP
