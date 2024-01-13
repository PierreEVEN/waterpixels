#pragma once
#include "Common/Image.h"

LibTIM::Image<LibTIM::TLabel> markers(const LibTIM::Image<LibTIM::RGB>& source, float sigma = 50);
