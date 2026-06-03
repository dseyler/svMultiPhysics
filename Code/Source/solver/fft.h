// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FFT_H
#define FFT_H

#include "ComMod.h"

#include <vector>
#include <limits>

void fft(const int np, const std::vector<std::vector<double>>& temporal_values, fcType& gt);

// time_override (default NaN -> use com_mod.time) evaluates the curve at a shifted
// time, e.g. com_mod.time - delay for per-element active-stress activation delay.
void ifft(const ComMod& com_mod, const fcType& gt, Vector<double>& Y, Vector<double>& dY,
    double time_override = std::numeric_limits<double>::quiet_NaN());

void igbc(const ComMod& com_mod, const MBType& gm, Array<double>& Y, Array<double>& dY);

#endif

