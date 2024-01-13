
#include <cstdint>

#include <libtim/Algorithms/Watershed.h>

#include <filesystem>

int main(int argc, char** argv) {
	auto image = LibTIM::Image<uint8_t>();

	std::cout << std::filesystem::current_path() << std::endl;

	LibTIM::Image<uint8_t>::load("Images/Peyto_Lake_Panorama.jpg", image);

	auto markers = LibTIM::Image<LibTIM::TLabel>(image);

	LibTIM::FlatSE flatSe;

	//LibTIM::watershedMeyer<uint8_t>(image, markers, flatSe, false);
}
