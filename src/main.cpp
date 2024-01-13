#include <cstdint>

#include <libtim/Algorithms/Watershed.h>

#include <filesystem>

#include "utils.hpp"

int main(int argc, char** argv)
{
	auto image = LibTIM::Image<LibTIM::RGB>();

	const auto val = rgbToCIELAB(LibTIM::RGB({58, 59, 60}));
	std::cout << val.x << "," << val.y << "," << val.z << std::endl;

	std::cout << std::filesystem::current_path() << std::endl;

	LibTIM::Image<LibTIM::RGB>::load("Images/Peyto_Lake_Panorama.ppm", image);


	// stocke les dimensions de l'image dans dx et dy
	int dx = image.getSizeX();
	int dy = image.getSizeY();

	LibTIM::Image<LibTIM::RGB> imgInCIELAB{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};
	for (int y = 0; y < dy; y++)
		for (int x = 0; x < dx; x++)
		{
			glm::vec3 pixelInLAB = rgbToCIELAB(image(x, y)); // Here the value are in [0, 1]
			imgInCIELAB(x, y) = LibTIM::RGB({
				static_cast<LibTIM::U8>(pixelInLAB.z / 100.f * 255.f), static_cast<LibTIM::U8>(pixelInLAB.z / 100.f * 255.f),
				static_cast<LibTIM::U8>(pixelInLAB.z / 100.f * 255.f)
			});
		}
	imgInCIELAB.save("Images/imgInCIELAB.ppm");


	auto markers = LibTIM::Image<LibTIM::TLabel>(image.getSizeX(), image.getSizeY());


	LibTIM::FlatSE flatSe;

	//LibTIM::watershedMeyer<uint8_t>(image, markers, flatSe, false);
}
