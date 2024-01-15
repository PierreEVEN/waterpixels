#pragma once

#include "Common/Types.h"
#include "Common/Image.h"
#include "glm/vec3.hpp"
#include "glm/mat3x3.hpp"

namespace WP
{
	glm::vec3 rgbToCIELAB(LibTIM::RGB pixelRGB);

	std::vector<glm::vec2> makePoints(uint32_t width, uint32_t height, float sigma);

	LibTIM::Image<LibTIM::U8> sobelFilter(LibTIM::Image<LibTIM::U8> image);

	LibTIM::Image<LibTIM::U8> applyFilter(LibTIM::Image<LibTIM::U8> image, std::vector<std::vector<float>> kernel);

	LibTIM::Image<LibTIM::U8> gaussianFilter3x3(LibTIM::Image<LibTIM::U8> image);

	LibTIM::Image<LibTIM::U8> spatialRegularization(LibTIM::Image<LibTIM::U8> image, float sigma, float k);

	LibTIM::Image<LibTIM::U8> labelToBinaryImage(const LibTIM::Image<LibTIM::TLabel>& image);
	LibTIM::Image<LibTIM::TLabel> imageToBinaryLabel(const LibTIM::Image<LibTIM::U8>& image);
}
