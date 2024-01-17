#include <iostream>
#include <filesystem>

#include "Algorithms/ConnectedComponents.hxx"
#include "Algorithms/Watershed.hxx"
#include "Common/FlatSE.h"
#include "waterpixels/utils.hpp"
#include "waterpixels/waterpixels.hpp"

int main(int argc, char** argv)
{
	// Read parameters
	if (argc < 5)
	{
		std::cerr <<
			"Wrong usage : waterpixels <input> <output> <sigma> <k>\n\timage : pgm P6 input image path\n\toutput : ppm output image path\n\tsigma : default = 50\n\tk : default = 5"
			<< std::endl;
		return -1;
	}
	const auto sigma = static_cast<float>(atof(argv[3]));
	const auto k = static_cast<float>(atof(argv[4]));
	const float cellScale = argc > 5 ? static_cast<float>(atof(argv[5])) : 0.8f;
	const int blurRadius = argc > 6 ? static_cast<float>(atof(argv[6])) : 5;

	// Load image
	if (!std::filesystem::exists(argv[1]))
	{
		std::cerr << "Failed to find input file '" << argv[1] << "'." << std::endl;
		return -1;
	}
	auto image = LibTIM::Image<LibTIM::RGB>();
	LibTIM::Image<LibTIM::RGB>::load(argv[1], image);

	auto bluredImageIntensity = WP::gaussianFilter(WP::rgbImageIntensity(image), blurRadius, 1);

	// Waterpixels algorithm
	const auto markers = WP::waterpixel(bluredImageIntensity, sigma, k, cellScale);

	// Save image
	WP::labelToBinaryImage(markers).save(argv[2]);


	// return 0;

	// Move to the derivative space
	auto sobelImg = WP::sobelFilter(bluredImageIntensity);
	sobelImg.save("images/sobel.ppm");

	// This will serve as guide to the watershed algorithm
	auto regularizedSobelImg = WP::spatialRegularization(sobelImg, sigma, k);
	regularizedSobelImg.save("images/spatialRegularization.ppm");

	// Generate watershed origins
	const auto watershedSources = WP::makeWatershedMarkers(sobelImg, sigma, cellScale);
	WP::labelToBinaryImage(watershedSources).save("images/sources.ppm");
	auto sources = LibTIM::Image<LibTIM::RGB>(image.getSizeX(), image.getSizeY());
	for (size_t x = 0; x < sources.getSizeX(); ++x)
		for (size_t y = 0; y < sources.getSizeY(); ++y)
		{
			sources(x, y) = LibTIM::RGB({0, 0, 0});
			if (x % static_cast<int>(sigma) == 0 || y % static_cast<int>(sigma) == 0)
				sources(x, y)[0] = 255;
			if (watershedSources(x, y))
				sources(x, y)[1] = 255;
		}
	for (const auto& point : WP::makePoints(sobelImg.getSizeX(), sobelImg.getSizeY(), sigma))
		if (sources.isPosValid(point.x, point.y))
			sources(point.x, point.y)[2] = 255;
	sources.save("images/sources.pgm");

	auto result = LibTIM::Image<LibTIM::RGB>(image.getSizeX(), image.getSizeY());
	for (size_t x = 0; x < result.getSizeX(); ++x)
		for (size_t y = 0; y < result.getSizeY(); ++y)
		{
			result(x, y) = image(x, y);
			if (markers(x, y))
				result(x, y) = LibTIM::RGB({255, 255, 255});
		}
	result.save("images/final.pgm");
}
