#include "waterpixels/waterpixels.hpp"

#include "waterpixels/utils.hpp"

#include <libtim/Algorithms/ConnectedComponents.hxx>
#include <libtim/Algorithms/Watershed.hxx>
#include <libtim/Common/FlatSE.h>

#include <glm/glm.hpp>

namespace WP
{
	template <typename Lambda_T>
	void iterateConnectedComponent(LibTIM::Image<LibTIM::TLabel> background, int ccOffset, const glm::ivec2& start,
	                               Lambda_T l)
	{
		const auto centerClamped = glm::ivec2{
			std::clamp(start.x, 0, background.getSizeX() - 1), std::clamp(start.y, 0, background.getSizeY() - 1)
		};
		const unsigned long componentValue = background(centerClamped.x, centerClamped.y);
		std::vector testedPoints = {centerClamped};
		while (!testedPoints.empty())
		{
			const auto point = testedPoints.back();
			testedPoints.pop_back();

			l(point);

			const auto testPosition = [&](const glm::ivec2& pos)
			{
				if (!background.isPosValid(pos.x, pos.y) || background(pos.x, pos.y) != componentValue)
					return;
				testedPoints.emplace_back(pos);
				background(pos.x, pos.y) += ccOffset;
			};

			testPosition(point + glm::ivec2{1, 0});
			testPosition(point + glm::ivec2{-1, 0});
			testPosition(point + glm::ivec2{0, 1});
			testPosition(point + glm::ivec2{0, -1});
		}
	}


	LibTIM::Image<LibTIM::U8> rgbImageIntensity(const LibTIM::Image<LibTIM::RGB>& image)
	{
		LibTIM::Image<LibTIM::U8> lImage(image.getSizeX(), image.getSizeY());
		for (int x = 0; x < image.getSizeX(); x++)
			for (int y = 0; y < image.getSizeY(); y++)
				lImage(x, y) = static_cast<LibTIM::U8>(rgbToCIELAB(image(x, y)).r / 100.f * 255.f);
		return lImage;
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
		LibTIM::Image<LibTIM::TLabel> markers(source.getSizeX(), source.getSizeY());
		LibTIM::Image<LibTIM::TLabel> voronoiMap(source.getSizeX(), source.getSizeY());
		LibTIM::Image<LibTIM::TLabel> ccMarkerMap(source.getSizeX(), source.getSizeY());
		markers.fill(0);

		// 0) Create the centers
		const auto centers = makePoints(source.getSizeX(), source.getSizeY(), sigma);

		// 1) Generate the voronoi diagram
		// @TODO : accelerate using a hash map
		const auto getClosestCenter = [&](const glm::ivec2& pos) -> unsigned long
		{
			float closestDistance = FLT_MAX;
			size_t closestIndex = 0;
			for (size_t i = 0; i < centers.size(); ++i)
			{
				const auto delta = centers[i] - pos;
				const float distance = static_cast<float>(std::sqrt(delta.x * delta.x + delta.y * delta.y));
				if (distance < closestDistance)
				{
					closestDistance = distance;
					closestIndex = i;
				}
			}
			return static_cast<unsigned int>(closestIndex);
		};

		for (LibTIM::TSize x = 0; x < source.getSizeX(); ++x)
			for (LibTIM::TSize y = 0; y < source.getSizeY(); ++y)
				voronoiMap(x, y) = getClosestCenter({x, y}) + 1;

		for (const auto& center : centers)
		{
			// 2) Search the local minimum in each voronoi cell
			float minValue = FLT_MAX;
			iterateConnectedComponent(voronoiMap, centers.size(), center, [&](const glm::ivec2& point)
			{
				if (const float val = source(point.x, point.y); val < minValue)
					minValue = val;
			});

			iterateConnectedComponent(voronoiMap, -centers.size(), center, [&](const glm::ivec2& point)
			{
				if (const float val = source(point.x, point.y); std::abs(val - minValue) < WP_MARKER_EPSILON)
					markers(point.x, point.y) = 1;
			});

			// 3) Keep only largest cell
			size_t newCCIndex = 2;
			std::vector<std::pair<size_t, size_t>> componentSizes;
			iterateConnectedComponent(voronoiMap, centers.size(), center, [&](const glm::ivec2& point)
			{
					if (markers(point.x, point.y) == 1) {
						size_t CCSize = 0;
						iterateConnectedComponent(markers, newCCIndex, point, [&](const glm::ivec2 newPoint)
							{
								CCSize++;
							});
						componentSizes.emplace_back(std::pair{ newCCIndex + 1, CCSize });
						newCCIndex++;
					}
			});

			auto maxCCSize = 0;
			auto maxCC = componentSizes[0].first;
			for (const auto& p: componentSizes)
			{
				if (p.second > maxCCSize) {
					maxCC = p.first;
					maxCCSize = p.second;
				}
			}
			iterateConnectedComponent(voronoiMap, centers.size(), center, [&](const glm::ivec2& point)
				{
					if (markers(point.x, point.y) == maxCC)
						markers(point.x, point.y) = 1;
					else
						markers(point.x, point.y) = 0;
				});
		}


		labelToBinaryImage(markers).save("images/test.ppm");
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
