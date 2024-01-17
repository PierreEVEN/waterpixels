#pragma once
#include <libtim/Common/Image.h>
#include <libtim/Common/Types.h>

#ifndef WP_MARKER_EPSILON
#define WP_MARKER_EPSILON 0.0001
#endif // WP_MARKER_EPSILON


#ifndef PREFER_CELL_CENTER
#define PREFER_CELL_CENTER true
#endif // PREFER_CELL_CENTER

namespace WP
{
	LibTIM::Image<LibTIM::U8> rgbImageIntensity(const LibTIM::Image<LibTIM::RGB>& image);
	LibTIM::Image<LibTIM::TLabel> makeWatershedMarkers(const LibTIM::Image<LibTIM::U8>& source, float sigma, float cellScale = 0.8f);
	
	/*
	* Spatial regularization according to a grid with cells of length sigma
	@param image: a grayscale image
	@param sigma: the size of a cell in the grid
	@param k: the regularization parameter (if equals 0 then no regularization)
	*/
	LibTIM::Image<LibTIM::U8> spatialRegularization(LibTIM::Image<LibTIM::U8> image, float sigma = 50, float k = 5);
	LibTIM::Image<LibTIM::TLabel> waterpixel(const LibTIM::Image<LibTIM::U8>& grayScaleImage, float sigma = 50, float k = 5, float cellScale = 0.8f);
}
