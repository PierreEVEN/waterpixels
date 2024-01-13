#include "markers.hpp"

#include <glm/glm.hpp>

#include "utils.hpp"

static std::vector<glm::vec2> makePoints(uint32_t width, uint32_t height, float sigma)
{
	std::vector<glm::vec2> points;
	points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));
	for (uint32_t x = 0; x < static_cast<uint32_t>(width / sigma) + 1; ++x)
	{
		for (uint32_t y = 0; y < static_cast<uint32_t>(height / sigma) + 1; ++y)
		{
			points.emplace_back(glm::vec2{x * sigma + sigma / 2, y * sigma + sigma / 2});
		}
	}
	return points;
}


static void getLocalMinimums(const LibTIM::Image<LibTIM::RGB>& source, LibTIM::Image<LibTIM::TLabel>& markers,
                             float sigma,
                             const glm::vec2& center)
{
	const int32_t exploredWidth = static_cast<int32_t>(sigma);
	const int32_t exploredHeight = static_cast<int32_t>(sigma);

	uint32_t imageWidth = exploredWidth + 1;
	uint32_t imageHeight = exploredWidth + 1;
	if (markers.getSizeX() != imageWidth || markers.getSizeY() != imageHeight)
		markers.setSize(imageWidth, imageHeight, 1);
	markers.fill(0);

	const int32_t minX = std::max(0, static_cast<int32_t>(center.x - std::floor(exploredWidth / 2 - 1)));
	const int32_t maxX = std::min(source.getSizeX() - 1,
	                              static_cast<int32_t>(center.x + std::floor(exploredWidth / 2)));
	const int32_t minY = std::max(0, static_cast<int32_t>(center.y - std::floor(exploredHeight / 2 - 1)));
	const int32_t maxY = std::min(source.getSizeY() - 1,
	                              static_cast<int32_t>(center.y + std::floor(exploredHeight / 2)));

	const auto getMarkers = [&](uint32_t globalX, uint32_t globalY)
	{
		return markers(globalX - center.x + exploredWidth / 2 - 1, globalY - center.y + exploredHeight / 2 - 1);
	};

	const auto setMarkers = [&](int32_t globalX, int32_t globalY, unsigned long value)
	{
		markers(globalX - center.x + exploredWidth / 2 - 1, globalY - center.y + exploredHeight / 2 - 1) = value;
	};

	// Get Min
	float minChunkValue = FLT_MAX;
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			float lValue = rgbToCIELAB(source(x, y)).x;
			if (lValue < minChunkValue)
				minChunkValue = lValue;
		}
	}

	// Get Min Points
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			const float lValue = rgbToCIELAB(source(x, y)).x;
			if (std::abs(lValue - minChunkValue) <= MARKER_EPSILON)
			{
				setMarkers(x, y, 1);
			}
			else
				setMarkers(x, y, 0);
		}
	}

	// Find the closest connected primitive
	std::vector<std::pair<uint32_t, uint32_t>> componentDistances;
	const auto getComponentDistanceAndUpdateValue = [&](const glm::vec2& origin)
	{
		std::vector<glm::ivec2> pointsToUpdate = {origin};

		uint32_t newValue = componentDistances.size() + 2;

		uint32_t nbPoints = 0;

		float minDistance = FLT_MAX;

		while (!pointsToUpdate.empty())
		{
			const auto coord = glm::ivec2(pointsToUpdate.back());
			nbPoints++;
			setMarkers(coord.x, coord.y, newValue);

			const float distance = std::sqrt(std::pow(coord.x - center.x, 2) + std::pow(coord.y - center.y, 2));
			if (distance < minDistance)
				minDistance = distance;

			pointsToUpdate.pop_back();

			for (int32_t x = -1; x <= 1; ++x)
				for (int32_t y = -1; y <= 1; ++y)
					if (!(x == 0 && y == 0))
					{
						const auto newCoord = glm::ivec2(std::clamp(coord.x + x, minX, maxX),
						                                 std::clamp(coord.y + y, minY, maxY));
						if (getMarkers(newCoord.x, newCoord.y) == 1)
							pointsToUpdate.emplace_back(newCoord);
					}
		}
		componentDistances.emplace_back(std::pair{newValue, nbPoints});
	};

	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{			
			if (getMarkers(x, y) == 1) {
				getComponentDistanceAndUpdateValue(glm::ivec2{ x, y });
			}
		}
	}

	if (componentDistances.empty())
		throw std::runtime_error("Failed to find minimum distance");

	uint32_t minComponent = componentDistances[0].first;
	float maxSize = 0;
	for (const auto p : componentDistances)
		if (p.second > maxSize)
		{
			maxSize = p.second;
			minComponent = p.first;
		}

	// Erase undesired components
	for (int32_t x = minX; x <= maxX; ++x)
	{
		for (int32_t y = minY; y <= maxY; ++y)
		{
			if (getMarkers(x, y) == minComponent)
			{
				setMarkers(x, y, 1);
			}
			else
				setMarkers(x, y, 0);
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

		// Copy each local minimum to global markers
		for (int x = 0; x < minimums.getSizeX(); ++x)
		{
			if (x + center.x - minimums.getSizeX() / 2 + 1 >= source.getSizeX())
				break;
			for (int y = 0; y < minimums.getSizeY(); ++y)
			{
				if (!minimums(x, y))
					continue;
				if (y + center.y - minimums.getSizeY() / 2 + 1 >= source.getSizeY())
					break;

				markers(x + center.x - minimums.getSizeX() / 2 + 1, y + center.y - minimums.getSizeY() / 2 + 1) =
					minimums(x, y);
			}
		}
	}

	return markers;
}
