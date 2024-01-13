#pragma once
#include "Common/Image.h"

#ifndef MARKER_EPSILON
#define MARKER_EPSILON 0.2f//0.0001f
#endif // MARKER_EPSILON


LibTIM::Image<LibTIM::TLabel> markers(const LibTIM::Image<LibTIM::U8>& source, float sigma = 50);