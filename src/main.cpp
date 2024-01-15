#include "waterpixels.hpp"

int main(int, char**)
{
	constexpr float sigma = 50;
	constexpr float k = 5.f;
	auto image = LibTIM::Image<LibTIM::RGB>();
	LibTIM::Image<LibTIM::RGB>::load("Images/landscape.ppm", image);

	// Waterpixels algorithm
	const auto markers = WP::waterpixel(WP::rgbImageIntensity(image), sigma, k);
	
	LibTIM::Image<LibTIM::RGB> stackedResults(image.getSizeX(), image.getSizeY());
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
		{
			stackedResults(x, y)[0] = markers(x, y) ? 255 : 0;
			if (x % static_cast<int>(sigma) == 0 || y % static_cast<int>(sigma) == 0)
				stackedResults(x, y)[2] = 255;
		}

	stackedResults.save("Images/watershed.pgm");
}
