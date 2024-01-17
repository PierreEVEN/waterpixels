#include "waterpixels/waterpixels.hpp"

#include <unordered_set>

#include "waterpixels/utils.hpp"

#include <libtim/Algorithms/ConnectedComponents.hxx>
#include <libtim/Algorithms/Watershed.hxx>
#include <libtim/Common/FlatSE.h>

#include <glm/glm.hpp>

template <>
struct std::hash<glm::ivec2>
{
	std::size_t operator()(const glm::ivec2& k) const
	{
		using std::size_t;
		using std::hash;
		using std::string;

		// Compute individual hash values for first,
		// second and third and combine them using XOR
		// and bit shifting:

		return ((hash<int>()(k.x)
			^ (hash<int>()(k.y) << 1)) >> 1);
	}
};


namespace WP
{
	glm::ivec2 clampPointToImage(LibTIM::Image<LibTIM::TLabel> image, glm::ivec2 point)
	{
		return glm::ivec2{
			std::clamp(point.x, 0, image.getSizeX() - 1), std::clamp(point.y, 0, image.getSizeY() - 1)
		};
	}

	template <typename Lambda_T>
	void iterateConnectedComponent(LibTIM::Image<LibTIM::TLabel> background, int ccOffset, const glm::ivec2& start,
	                               Lambda_T callback)
	{
		const auto centerClamped = clampPointToImage(background, start);
		const unsigned long componentValue = background(centerClamped.x, centerClamped.y);
		std::vector testedPoints = {centerClamped};
		background(centerClamped.x, centerClamped.y) += ccOffset;
		while (!testedPoints.empty())
		{
			const auto point = testedPoints.back();
			testedPoints.pop_back();

			callback(point);

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
				const auto delta = currentPixel - closestCenter;

				// Compute L1 distance from point to center
				float d = std::max(std::abs(currentPixel.x - closestCenter.x),
				                   std::abs(currentPixel.y - closestCenter.y));

				regularizedImg(x, y) = static_cast<LibTIM::U8>(std::min(
					static_cast<int>(image(x, y) + k * (2.f * d / sigma)), 255));
			}
		}

		return regularizedImg;
	}

	LibTIM::Image<LibTIM::TLabel> makeWatershedMarkers(const LibTIM::Image<LibTIM::U8>& source, float sigma,
	                                                   float cellScale)
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
			std::vector<glm::ivec2> cellPoints;
			iterateConnectedComponent(voronoiMap, centers.size(), center, [&](const glm::ivec2& point)
			{
				cellPoints.emplace_back(point);
			});

			// 3) Apply homothety on the cell points;
			std::unordered_set<glm::ivec2> movedPoints;
			for (const auto& point : cellPoints)
			{
				const auto pointToCenter = glm::vec2(center - point);
				const auto distance = length(pointToCenter);
				const auto newPos = center + glm::ivec2(normalize(-pointToCenter) * distance * cellScale);
				if (source.isPosValid(newPos.x, newPos.y))
					movedPoints.insert(newPos);
			}

			cellPoints.clear();
			for (const auto& point : movedPoints)
				cellPoints.emplace_back(point);

			// 2) Search the local minimum in each voronoi cell
			float minValue = FLT_MAX;
			for (const auto& point : cellPoints)
				if (const float val = source(point.x, point.y); val < minValue)
					minValue = val;

			for (const auto& point : cellPoints)
				if (const float val = source(point.x, point.y); std::abs(val - minValue) < WP_MARKER_EPSILON)
					markers(point.x, point.y) = 1;

			// 3) Keep only largest cell
#if PREFER_CELL_CENTER
			float selectedValue = FLT_MAX;
#else
			int selectedValue = 0;
#endif
			size_t selectedIndex = 0;
			std::vector<std::vector<glm::ivec2>> componentSizes;
			for (const auto& point : cellPoints)
			{
				if (markers(point.x, point.y) == 1)
				{
#if PREFER_CELL_CENTER
					/* VERSION WITH CLOSEST Connected Component */
					float minDistance = FLT_MAX;
					std::vector<glm::ivec2> componentElements;
					iterateConnectedComponent(markers, 1, point, [&](const glm::ivec2& pt)
					{
						componentElements.emplace_back(pt);
						const auto delta = pt - center;
						const auto distance = std::sqrt(delta.x * delta.x + delta.y * delta.y);
						if (distance < minDistance)
							minDistance = distance;
					});
					if (minDistance < selectedValue)
					{
						selectedValue = minDistance;
						selectedIndex = componentSizes.size();
					}
#else
					/* VERSION WITH LARGEST Connected Component */
					std::vector<glm::ivec2> componentElements;
					iterateConnectedComponent(markers, 1, point, [&](const glm::ivec2& pt)
					{
						componentElements.emplace_back(pt);
					});
					if (componentElements.size() > selectedValue)
					{
						selectedValue = componentElements.size();
						selectedIndex = componentSizes.size();
					}
#endif
					componentSizes.emplace_back(componentElements);
				}
			}
			if (componentSizes.empty())
				componentSizes = {{clampPointToImage(source, center)}};

			for (const auto& point : cellPoints)
				markers(point.x, point.y) = 0;
			for (const auto& ccPoint : componentSizes[selectedIndex])
				markers(ccPoint.x, ccPoint.y) = 1;
		}

		labelToBinaryImage(markers).save("images/test.ppm");
		return markers;
	}


	LibTIM::Image<LibTIM::TLabel> waterpixel(const LibTIM::Image<LibTIM::U8>& grayScaleImage, float sigma, float k, float cellScale)
	{
		// Move to the derivative space
		const auto sobelImg = sobelFilter(grayScaleImage);

		// This will serve as guide to the watershed algorithm
		auto regularizedSobelImg = spatialRegularization(sobelImg, sigma, k);

		// Generate watershed origins
		auto watershedSources = makeWatershedMarkers(sobelImg, sigma, cellScale);

		LibTIM::FlatSE connectivity;
		connectivity.make2DN4();

		LibTIM::Image<LibTIM::TLabel> minLabeled = labelConnectedComponents(watershedSources, connectivity);
		LibTIM::watershedMeyer<uint8_t>(regularizedSobelImg, minLabeled, connectivity);

		return imageToBinaryLabel(morphologicalGradient(minLabeled, connectivity));
	}
}
