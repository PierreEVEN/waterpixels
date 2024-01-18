#include <iostream>
#include <filesystem>

#include <libtim/Algorithms/Morphology.h>
#include <libtim/Common/FlatSE.h>
#include <waterpixels/utils.hpp>
#include <waterpixels/waterpixels.hpp>
#include <waterpixels/config.hpp>

#ifndef OUTPUT_DEBUG
#define OUTPUT_DEBUG false
#endif

#if OUTPUT_DEBUG
#include <Algorithms/ConnectedComponents.h>
#endif

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
	const float cellScale = argc > 5 ? static_cast<float>(atof(argv[5])) : 2 / 3.f;
	const int blurRadius = argc > 6 ? static_cast<float>(atof(argv[6])) : 5;

	/****** 1) Load the image ******/
	if (!std::filesystem::exists(argv[1]))
	{
		std::cerr << "Failed to find input file '" << argv[1] << "'." << std::endl;
		return -1;
	}
	auto image = LibTIM::Image<LibTIM::RGB>();
	{
		MEASURE_DURATION(loading, "Load image");
		LibTIM::Image<LibTIM::RGB>::load(argv[1], image);
	}

	/****** 2) Pre filter image to remove details ******/
	LibTIM::FlatSE filter;
	filter.make2DEuclidianBall(blurRadius);
	LibTIM::Image<LibTIM::U8> preFilteredImage;
	{
		MEASURE_DURATION(prefiltering, "Image pre-filtering");
		if (blurRadius >= 1)
		{
			preFilteredImage = closing(opening(WP::rgbImageIntensity(image), filter), filter);
		}
		else
			preFilteredImage = WP::rgbImageIntensity(image);
	}

	/****** 3) Choose cell centers ******/
	std::vector<glm::ivec2> cellCenters;
	{
		MEASURE_DURATION(centers, "Generate cell centers");
		cellCenters = WP::makeRectGrid2D(preFilteredImage.getSizeX(), preFilteredImage.getSizeY(), sigma);
	}

	/****** 4) Execute waterpixel algorithm ******/
	LibTIM::Image<LibTIM::TLabel> markers;
	{
		MEASURE_DURATION(watershed, "Waterpixel algorithm");
		markers = WP::waterpixel(preFilteredImage, cellCenters, sigma, k, cellScale);
	}

	filter.make2DN4();
	const auto markerDelimitation = morphologicalGradient(markers, filter);
	
	// Save image
	WP::labelToBinaryImage(markerDelimitation).save(argv[2]);

#if OUTPUT_DEBUG

	const WP::VoronoiGraph voronoi(preFilteredImage.getSizeX(), preFilteredImage.getSizeY(), cellCenters);

	// Move to the derivative space
	filter.make2DN4();
	auto gradient = morphologicalGradient(preFilteredImage, filter);
	gradient.save("images/imageGradient.ppm");

	// This will serve as guide to the watershed algorithm
	auto regularizedSobelImg = spatialRegularization(gradient, voronoi, sigma, k);
	regularizedSobelImg.save("images/spatialRegularizationGradient.ppm");

	// Generate watershed origins
	const auto watershedSources = makeWatershedMarkers(gradient, voronoi, sigma, cellScale);

	auto gridDebugImage = voronoi.debugVisualization();

	for (size_t x = 0; x < gridDebugImage.getSizeX(); ++x)
		for (size_t y = 0; y < gridDebugImage.getSizeY(); ++y)
			if (watershedSources(x, y))
				gridDebugImage(x, y)[2] = 255;
	gridDebugImage.save("images/watershedSources.pgm");

	// Return the delimitation instead of just connected components

	auto result = LibTIM::Image<LibTIM::RGB>(image.getSizeX(), image.getSizeY());
	for (size_t x = 0; x < result.getSizeX(); ++x)
		for (size_t y = 0; y < result.getSizeY(); ++y)
		{
			result(x, y) = image(x, y);
			if (markerDelimitation(x, y))
				result(x, y) = LibTIM::RGB({255, 255, 255});
		}
	result.save("images/combined.pgm");
#endif
}
