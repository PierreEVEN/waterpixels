#include <cstdint>

#include <libtim/Algorithms/Watershed.h>

#include <filesystem>

#include "markers.hpp"
#include "utils.hpp"
#include "Algorithms/ConnectedComponents.h"
#include "Algorithms/Misc.h"
#include "Algorithms/RegionGrowing.h"

int main(int argc, char** argv)
{
	auto image = LibTIM::Image<LibTIM::RGB>();
	LibTIM::Image<LibTIM::RGB>::load("Images/landscape.ppm", image);

	LibTIM::Image<LibTIM::U8> grayScaleImage{
		static_cast<LibTIM::TSize>(image.getSizeX()), static_cast<LibTIM::TSize>(image.getSizeY())
	};
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
			grayScaleImage(x, y) = static_cast<LibTIM::U8>(rgbToCIELAB(image(x, y)).r / 100.f * 255.f);

	/*
	grayScaleImage.save("Images/landscape_CIELAB.ppm");

	auto sobelImg = sobelFilter(grayScaleImage);
	sobelImg.save("Images/imgSobel.ppm");
	*/

	int sigma = 50;
	LibTIM::Image<LibTIM::TLabel> ms = markers(image, sigma);


	LibTIM::Image<LibTIM::U8> minimumPoints(image.getSizeX(), image.getSizeY());
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
			minimumPoints(x, y) = ms(x, y) ? 255 : 0;


	LibTIM::FlatSE connex;
	connex.make2DN4();
	LibTIM::Image<LibTIM::TLabel> minLabeled = LibTIM::labelConnectedComponents(minimumPoints, connex);

	LibTIM::FlatSE flatSe;
	flatSe.make2DN8();
	LibTIM::watershedMeyer<uint8_t>(grayScaleImage, minLabeled, flatSe);

	LibTIM::Image<LibTIM::U8> finalMap = LibTIM::morphologicalGradient(minLabeled, flatSe);

	LibTIM::Image<LibTIM::RGB> stackedResults(image.getSizeX(), image.getSizeY());
	for (int x = 0; x < image.getSizeX(); x++)
		for (int y = 0; y < image.getSizeY(); y++)
		{
			stackedResults(x, y)[0] = finalMap(x, y) ? 255 : 0;
			stackedResults(x, y)[1] = ms(x, y) ? 255 : 0;
			if (x % sigma == 0 || y % sigma == 0)
				stackedResults(x, y)[2] = 255;
		}

	stackedResults.save("Images/watershed.ppm");
}
