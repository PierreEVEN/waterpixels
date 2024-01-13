#include "markers.hpp"


static std::vector<> makePoints(uint32_t width, uint32_t height, float sigma)
{
	
}

LibTIM::Image<LibTIM::TLabel> markers(const LibTIM::Image<LibTIM::RGB>& source, float sigma)
{
	auto markers = LibTIM::Image<LibTIM::TLabel>(source.getSizeX(), source.getSizeY());
	markers.fill(0);



}
