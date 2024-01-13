#include "markers.hpp"

#include <glm/glm.hpp>

static std::vector<glm::vec2> makePoints(uint32_t width, uint32_t height, float sigma)
{
	std::vector<glm::vec2> points;
	points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));
	for (uint32_t x = 0; x < static_cast<uint32_t>(width / sigma); ++x)
	{
		for (uint32_t y = 0; y < static_cast<uint32_t>(height / sigma); ++y)
		{
			points.emplace_back(glm::vec2{x * sigma, height * sigma});
		}
	}

	return points;
}

static std::vector<glm::vec2> getMinimums(const LibTIM::Image<LibTIM::RGB>& source, float sigma,
                                          const glm::vec2& center)
{
	const uint32_t exploredWidth = static_cast<uint32_t>(source.getSizeX() / sigma);
	const uint32_t exploredHeight = static_cast<uint32_t>(source.getSizeY() / sigma);

	float min = FLT_MAX;

	for (uint32_t x = std::max(0u, static_cast<uint32_t>(center.x - exploredWidth / 2));
		x <= std::min(static_cast<uint32_t>(source.getSizeX()), static_cast<uint32_t>(center.x + exploredWidth / 2));
		++x)
	{
		for (uint32_t y = std::max(0u, static_cast<uint32_t>(center.x - exploredWidth / 2));
			y <= std::min(static_cast<uint32_t>(source.getSizeX()), static_cast<uint32_t>(center.x + exploredWidth / 2));
			++y)
		{
			source(x, y) = LibTIM::RGB({ 1, 2, 3 });
		}		
	}


		std::vector<glm::vec2> minimums(exploredWidth * exploredHeight);
}


LibTIM::Image<LibTIM::TLabel> markers(const LibTIM::Image<LibTIM::RGB>& source, float sigma)
{
	auto markers = LibTIM::Image<LibTIM::TLabel>(source.getSizeX(), source.getSizeY());
	markers.fill(0);

	const auto centers = makePoints(source.getSizeX(), source.getSizeY(), sigma);


	return markers;
}
