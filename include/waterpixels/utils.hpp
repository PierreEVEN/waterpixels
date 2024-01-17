#pragma once

#include <libtim/Common/Types.h>
#include <libtim/Common/Image.h>
#include <glm/glm.hpp>

namespace WP
{
	glm::vec3 rgbToCIELAB(LibTIM::RGB pixelRGB);

	std::vector<glm::ivec2> makePoints(uint32_t width, uint32_t height, float sigma);

	LibTIM::Image<LibTIM::U8> sobelFilter(LibTIM::Image<LibTIM::U8> image);

	LibTIM::Image<LibTIM::U8> applyFilter(LibTIM::Image<LibTIM::U8> image, std::vector<std::vector<float>> kernel);

	LibTIM::Image<LibTIM::U8> gaussianFilter(LibTIM::Image<LibTIM::U8> image, int radius, float sigma);
	
	LibTIM::Image<LibTIM::U8> labelToImage(const LibTIM::Image<LibTIM::TLabel>& image);
	LibTIM::Image<LibTIM::U8> labelToBinaryImage(const LibTIM::Image<LibTIM::TLabel>& image);
	LibTIM::Image<LibTIM::TLabel> imageToBinaryLabel(const LibTIM::Image<LibTIM::U8>& image);
}
