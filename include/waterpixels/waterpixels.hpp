#pragma once
#include <libtim/Common/Image.h>
#include <libtim/Common/Types.h>
#include "config.hpp"

#include "glm/vec2.hpp"

#ifndef WP_MARKER_EPSILON
#define WP_MARKER_EPSILON 0.0001
#endif // WP_MARKER_EPSILON


#ifndef PREFER_CELL_CENTER
#define PREFER_CELL_CENTER true
#endif // PREFER_CELL_CENTER

#ifndef USE_LINF_REG_DISTANCE
#define USE_LINF_REG_DISTANCE false
#endif // USE_LINF_REG_DISTANCE

namespace WP
{
	class VoronoiGraph;
	LibTIM::Image<LibTIM::U8> rgbImageIntensity(const LibTIM::Image<LibTIM::RGB>& image);
	LibTIM::Image<LibTIM::TLabel> makeWatershedMarkers(const LibTIM::Image<LibTIM::U8>& source, const VoronoiGraph& voronoiCells, float sigma, float cellScale);
	
	/*
	* Spatial regularization according to a grid with cells of length sigma
	@param image: a grayscale image
	@param sigma: the size of a cell in the grid
	@param k: the regularization parameter (if equals 0 then no regularization)
	*/
	[[nodiscard]] LibTIM::Image<LibTIM::U8> spatialRegularization(const LibTIM::Image<LibTIM::U8>& image, const VoronoiGraph& voronoiCells, float sigma, float k);
	[[nodiscard]] LibTIM::Image<LibTIM::TLabel> waterpixel(const LibTIM::Image<LibTIM::U8>& grayScaleImage, const std::vector<glm::ivec2> cellCenters, float sigma, float k, float cellScale);
}
