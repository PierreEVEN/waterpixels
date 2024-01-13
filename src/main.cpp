
#include <cstdint>

#include <libtim/Algorithms/Watershed.h>

#include <filesystem>

int main(int argc, char** argv) {
	auto image = LibTIM::Image<LibTIM::RGB>();

	std::cout << std::filesystem::current_path() << std::endl;

	LibTIM::Image<LibTIM::RGB>::load("Images/Peyto_Lake_Panorama.ppm", image);

	auto markers = LibTIM::Image<LibTIM::TLabel>(image.getSizeX(), image.getSizeY());

	LibTIM::FlatSE flatSe;

	//LibTIM::watershedMeyer<uint8_t>(image, markers, flatSe, false);
}
