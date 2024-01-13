#include <cstdint>

#include <libtim/Algorithms/Watershed.h>

#include <filesystem>

#include "markers.hpp"
#include "utils.hpp"

int main(int argc, char** argv)
{
	auto image = LibTIM::Image<LibTIM::RGB>();

	const auto val = rgbToCIELAB(LibTIM::RGB({58, 59, 60}));
	std::cout << val.x << "," << val.y << "," << val.z << std::endl;

	std::cout << std::filesystem::current_path() << std::endl;

	LibTIM::Image<LibTIM::RGB>::load("Images/landscape.ppm", image);

	int sigma = 20;

	LibTIM::Image<LibTIM::TLabel> ms = markers(image, sigma);

	auto out = LibTIM::Image<LibTIM::U8>(image.getSizeX(), image.getSizeY());
	for (int x = 0; x < image.getSizeX(); x++)
	{
		for (int y = 0; y < image.getSizeY(); y++)
		{
			if (x % sigma == 0 || y% sigma == 0)
				out(x, y) = 128;

			if (ms(x, y))
				out(x, y) = 255;
		}
	}
	out.save("Images/imgMarkers.ppm");
	return 0;

	//// stocke les dimensions de l'image dans dx et dy
	int dx = image.getSizeX();
	int dy = image.getSizeY();

	//LibTIM::Image<LibTIM::RGB> imgInCIELAB{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};
	LibTIM::Image<LibTIM::U8> imgInCIELAB{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };


	for (int y = 0; y < dy; y++)
		for (int x = 0; x < dx; x++)
		{
			glm::vec3 pixelInLAB = rgbToCIELAB(image(x, y)); // Here the value are in [0, 1]
			imgInCIELAB(x, y) = static_cast<LibTIM::U8>(pixelInLAB.r / 100.f * 255.f);
		}
	imgInCIELAB.save("Images/landscape_CIELAB.ppm");

	auto sobelImg = sobelFilter(imgInCIELAB);
	sobelImg.save("Images/imgSobel.ppm");

	auto markers = LibTIM::Image<LibTIM::TLabel>(image.getSizeX(), image.getSizeY());


	LibTIM::FlatSE flatSe;

	//LibTIM::watershedMeyer<uint8_t>(image, markers, flatSe, false);
}
