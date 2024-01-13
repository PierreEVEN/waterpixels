#include "markers.hpp"

#include <glm/glm.hpp>

#include "utils.hpp"

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


static void getLocalMinimums(const LibTIM::Image<LibTIM::RGB>& source, LibTIM::Image<LibTIM::TLabel>& minimums,
                             float sigma,
                             const glm::vec2& center)
{
	const uint32_t exploredWidth = static_cast<uint32_t>(source.getSizeX() / sigma);
	const uint32_t exploredHeight = static_cast<uint32_t>(source.getSizeY() / sigma);

	const int32_t minX = std::max(0u, static_cast<uint32_t>(center.x - exploredWidth / 2));
	const int32_t maxX = std::min(static_cast<uint32_t>(source.getSizeX() - 1),
	                              static_cast<uint32_t>(center.x + exploredWidth / 2));
	const int32_t minY = std::max(0u, static_cast<uint32_t>(center.y - exploredHeight / 2));
	const int32_t maxY = std::min(static_cast<uint32_t>(source.getSizeY() - 1),
	                              static_cast<uint32_t>(center.y + exploredHeight / 2));

	// Get Min
	float min = FLT_MAX;
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			float lValue = rgbToCIELAB(source(x, y)).x;
			if (lValue < min)
				min = lValue;
		}
	}

	// Get Min Points
	if (minimums.getSizeX() != exploredWidth || minimums.getSizeY() != exploredHeight)
		minimums.setSize(exploredWidth, exploredHeight, 1);
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			float lValue = rgbToCIELAB(source(x, y)).x;
			if (std::abs(lValue - min) <= MARKER_EPSILON)
			{
				minimums(static_cast<LibTIM::TCoord>(x - center.x), static_cast<LibTIM::TCoord>(y - center.y)) = 1;
			}
			else
				minimums(static_cast<LibTIM::TCoord>(x - center.x), static_cast<LibTIM::TCoord>(y - center.y)) = 0;
		}
	}

	// Find the closest connected primitive
	std::vector<std::pair<uint32_t, float>> componentDistances;
	const auto getComponentDistanceAndUpdateValue = [&](const glm::vec2& origin)
	{
		std::vector<glm::ivec2> pointsToUpdate = {origin};

		uint32_t newValue = componentDistances.size() + 2;

		float minDistance = FLT_MAX;

		while (!pointsToUpdate.empty())
		{
			const auto coord = glm::ivec2(pointsToUpdate.back());

			minimums(static_cast<LibTIM::TCoord>(coord.x - center.x),
			         static_cast<LibTIM::TCoord>(coord.y - center.y)) = newValue;

			float distance = std::sqrt(std::pow(coord.x - center.x, 2) + std::pow(coord.y - center.y, 2));
			if (distance < minDistance)
				minDistance = distance;

			pointsToUpdate.pop_back();

			for (int32_t x = -1; x <= 1; ++x)
				for (int32_t y = -1; y <= 1; ++y)
					if (!(x == 0 && y == 0))
					{
						const auto newCoord = coord + glm::ivec2(x, y);
						if (minimums(static_cast<LibTIM::TCoord>(newCoord.x - center.x),
						             static_cast<LibTIM::TCoord>(newCoord.y - center.y)) == 1)
							pointsToUpdate.emplace_back(newCoord);
					}
		}
		componentDistances.push_back({newValue, minDistance});
	};


	for (uint32_t x = minX; x <= maxX; ++x)
		for (uint32_t y = minY; y <= maxY; ++y)
			if (minimums(x, y) == 1)
				getComponentDistanceAndUpdateValue(glm::ivec2{x, y});

	uint32_t minComponent = componentDistances[0].first;
	float minCDistance = FLT_MAX;
	for (const auto p : componentDistances)
		if (p.second < minCDistance)
			minComponent = p.first;

	// Erase undesired components
	if (minimums.getSizeX() != exploredWidth || minimums.getSizeY() != exploredHeight)
		minimums.setSize(exploredWidth, exploredHeight, 1);
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			if (minimums(static_cast<LibTIM::TCoord>(x - center.x), static_cast<LibTIM::TCoord>(y - center.y)) ==
				minComponent)
				minimums(static_cast<LibTIM::TCoord>(x - center.x), static_cast<LibTIM::TCoord>(y - center.y)) = 1;
			else
				minimums(static_cast<LibTIM::TCoord>(x - center.x), static_cast<LibTIM::TCoord>(y - center.y)) = 0;
		}
	}
}


LibTIM::Image<LibTIM::TLabel> markers(const LibTIM::Image<LibTIM::RGB>& source, float sigma)
{
	auto markers = LibTIM::Image<LibTIM::TLabel>(source.getSizeX(), source.getSizeY());
	markers.fill(0);

	const auto centers = makePoints(source.getSizeX(), source.getSizeY(), sigma);

	LibTIM::Image<LibTIM::TLabel> minimums(1, 1);
	for (const auto& center : centers)
	{
		getLocalMinimums(source, minimums, sigma, center);
	}


	return markers;
}
