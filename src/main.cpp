#include <iostream>
#include <filesystem>

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
	const auto sigma = static_cast<float>(atof(argv[2]));
	const auto k = static_cast<float>(atof(argv[3]));

	// Load image
	if (!std::filesystem::exists(argv[1]))
	{
		std::cerr << "Failed to find input file '" << argv[1] << "'." << std::endl;
		return -1;
	}
	auto image = LibTIM::Image<LibTIM::RGB>();
	LibTIM::Image<LibTIM::RGB>::load(argv[1], image);

	// Waterpixels algorithm
	const auto markers = WP::waterpixel(WP::rgbImageIntensity(image), 50, 5);

	// Save image
	WP::labelToBinaryImage(markers).save(argv[4]);
}
