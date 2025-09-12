/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#pragma once
#include "device_params.cuh"

template <bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg, bool useRadialLut>
__global__ void evaluateKernel(const DeviceKernelParams p);
