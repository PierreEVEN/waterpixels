#include "waterpixels/waterpixels.hpp"

#include "waterpixels/utils.hpp"

#include <libtim/Algorithms/ConnectedComponents.hxx>
#include <libtim/Algorithms/Watershed.hxx>
#include <libtim/Common/FlatSE.h>

#include <glm/glm.hpp>

namespace WP
{
	LibTIM::Image<LibTIM::U8> rgbImageIntensity(const LibTIM::Image<LibTIM::RGB>& image)
	{
		LibTIM::Image<LibTIM::U8> lImage(image.getSizeX(), image.getSizeY());
		for (int x = 0; x < image.getSizeX(); x++)
			for (int y = 0; y < image.getSizeY(); y++)
				lImage(x, y) = static_cast<LibTIM::U8>(rgbToCIELAB(image(x, y)).r / 100.f * 255.f);
		return lImage;
	}

	static void getLocalMinimums(const LibTIM::Image<LibTIM::U8>& source, LibTIM::Image<LibTIM::TLabel>& markers,
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
				float pVal = source(x, y);
				if (pVal < minChunkValue)
					minChunkValue = pVal;
			}
		}

		// Get Min Points
		for (int32_t x = minX; x <= maxX; ++x)
		{
			for (int32_t y = minY; y <= maxY; ++y)
			{
				const float pVal = source(x, y);
				if (std::abs(pVal - minChunkValue) <= WP_MARKER_EPSILON)
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
				if (getMarkers(x, y) == 1)
				{
					getComponentDistanceAndUpdateValue(glm::ivec2{x, y});
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
	
	LibTIM::Image<LibTIM::U8> spatialRegularization(LibTIM::Image<LibTIM::U8> image, float sigma, float k)
	{
		int dx = image.getSizeX();
		int dy = image.getSizeY();

		LibTIM::Image<LibTIM::U8> regularizedImg{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};

		for (int y = 0; y < dy; y++)
		{
			for (int x = 0; x < dx; x++)
			{
				glm::vec2 closestCenter{
					int(static_cast<float>(x) / sigma) * sigma + sigma / 2.f,
					int(static_cast<float>(y) / sigma) * sigma + sigma / 2.f
				};
				glm::vec2 currentPixel{x, y};
				float d = std::max(std::abs(currentPixel.x - closestCenter.x),
				                   std::abs(currentPixel.y - closestCenter.y));

				regularizedImg(x, y) = static_cast<LibTIM::U8>(image(x, y) + k * ((2.f * d) / sigma));
			}
		}

		return regularizedImg;
	}

	LibTIM::Image<LibTIM::TLabel> makeWatershedMarkers(const LibTIM::Image<LibTIM::U8>& source, float sigma)
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


	LibTIM::Image<LibTIM::TLabel> waterpixel(const LibTIM::Image<LibTIM::U8>& grayScaleImage, float sigma, float k)
	{
		// Move to the derivative space
		const auto sobelImg = sobelFilter(grayScaleImage);

		// This will serve as guide to the watershed algorithm
		auto regularizedSobelImg = spatialRegularization(sobelImg, sigma, k);

		// Generate watershed origins
		auto watershedSources = makeWatershedMarkers(sobelImg, sigma);

		LibTIM::FlatSE connectivity;
		connectivity.make2DN4();
		
		LibTIM::Image<LibTIM::TLabel> minLabeled = labelConnectedComponents(watershedSources, connectivity);
		LibTIM::watershedMeyer<uint8_t>(regularizedSobelImg, minLabeled, connectivity);

		return imageToBinaryLabel(morphologicalGradient(minLabeled, connectivity));
	}
}
