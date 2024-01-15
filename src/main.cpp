#include "utils.hpp"
#include "waterpixels.hpp"

int main(int argc, char** argv)
{
	auto image = LibTIM::Image<LibTIM::RGB>();
	LibTIM::Image<LibTIM::RGB>::load("Images/landscape.ppm", image);

	LibTIM::Image<LibTIM::U8> grayScaleImage{
		static_cast<LibTIM::TSize>(image.getSizeX()), static_cast<LibTIM::TSize>(image.getSizeY())
	};
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
			grayScaleImage(x, y) = static_cast<LibTIM::U8>(WP::rgbToCIELAB(image(x, y)).r / 100.f * 255.f);


	int sigma = 50;
	float k = 5.f;
	const auto markers = WP::waterpixel(grayScaleImage, static_cast<float>(sigma), k);
	
	LibTIM::Image<LibTIM::RGB> stackedResults(image.getSizeX(), image.getSizeY());
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
		{
			stackedResults(x, y)[0] = markers(x, y) ? 255 : 0;
			if (x % sigma == 0 || y % sigma == 0)
				stackedResults(x, y)[2] = 255;
		}

	stackedResults.save("Images/watershed.pgm");
}
