// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FFT_H
#define FFT_H

#include "ComMod.h"

#include <vector>

void fft(const int np, const std::vector<std::vector<double>>& temporal_values, fcType& gt);

// `delay` (default 0) shifts the evaluation time to com_mod.time - delay, e.g. for
// per-element active-stress activation delay.
void ifft(const ComMod& com_mod, const fcType& gt, Vector<double>& Y, Vector<double>& dY,
    double delay = 0.0);

void igbc(const ComMod& com_mod, const MBType& gm, Array<double>& Y, Array<double>& dY);

#endif

