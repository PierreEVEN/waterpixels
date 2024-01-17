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
	void iterateConnectedComponent(LibTIM::Image<LibTIM::TLabel>& background, int ccOffset, const glm::ivec2& start,
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
		#pragma omp parallel for collapse(2)
		for (int x = 0; x < image.getSizeX(); x++)
			for (int y = 0; y < image.getSizeY(); y++)
				lImage(x, y) = static_cast<LibTIM::U8>(rgbToCIELAB(image(x, y)).r / 100.f * 255.f);
		return lImage;
	}

	LibTIM::Image<LibTIM::U8> spatialRegularization(const LibTIM::Image<LibTIM::U8>& source,
	                                                const VoronoiGraph& voronoiCells, float sigma, float k)
	{
		LibTIM::Image<LibTIM::U8> result(source.getSizeX(), source.getSizeY());
		for (const auto& cell : voronoiCells.cells())
		{
			const auto& center = cell.first;
			for (const auto& point : cell.second)
			{
				// Compute Linf distance from point to center
				//const float d = std::max(std::abs(point.x - center.x),std::abs(point.y - center.y));
				const auto delta = point - center;
				const float d = std::sqrt(delta.x * delta.x + delta.y * delta.y);

				result(point.x, point.y) = static_cast<LibTIM::U8>(std::min(
					static_cast<int>(source(point.x, point.y) + k * (2.f * d / sigma)), 255));
			}
		}
		return result;
	}

	LibTIM::Image<LibTIM::TLabel> makeWatershedMarkers(const LibTIM::Image<LibTIM::U8>& source,
	                                                   const VoronoiGraph& voronoiCells, float sigma,
	                                                   float cellScale)
	{
		LibTIM::Image<LibTIM::TLabel> markers(source.getSizeX(), source.getSizeY());
		LibTIM::Image<LibTIM::TLabel> voronoiMap(source.getSizeX(), source.getSizeY());
		LibTIM::Image<LibTIM::TLabel> ccMarkerMap(source.getSizeX(), source.getSizeY());
		markers.fill(0);

		size_t cellIndex = 1;
		MEASURE_AVERAGE_DURATION(cellMarkerAvg, "Generate watershed markers for one cell");
		MEASURE_CUMULATIVE_DURATION(cellHomotTot, "Apply homothety for one cell");
		MEASURE_CUMULATIVE_DURATION(searchAllMin, "Search all pixels with minimum value in cell");
		MEASURE_CUMULATIVE_DURATION(cellFilterBestTarget,
		                            "Filter and keep only the best connected component for each cell");
		MEASURE_CUMULATIVE_DURATION(iterateSubCellComponents,
			"Local cell component iteration");

		const auto& cells = voronoiCells.cells();

		std::mutex lastMutex;
		auto last = std::chrono::steady_clock::now();

#pragma omp parallel for
		for (int64_t ci = 0; ci < cells.size(); ++ci)
		{
			if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - last).count() >= 1000)
			{
				std::lock_guard m(lastMutex);
				last = std::chrono::steady_clock::now();
				std::cout << "Watershed markers generation : " << static_cast<float>(ci) / cells.size() * 100 << "% ..." << std::endl;
			}


			const auto& cell = cells[ci];
			MEASURE_ADD_CUMULATOR(cellMarkerAvg);
			const auto& center = cell.first;

			auto cellPoints = cell.second;

			// 3) Apply homothety on the cell points;
			std::unordered_set<glm::ivec2> movedPoints;
			{
				MEASURE_ADD_CUMULATOR(cellHomotTot);
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
			}

			// 2) Search the local minimum in each voronoi cell
			float minValue = FLT_MAX;
			{
				MEASURE_ADD_CUMULATOR(searchAllMin);
				for (const auto& point : cellPoints)
					if (const float val = source(point.x, point.y); val < minValue)
						minValue = val;

				for (const auto& point : cellPoints)
					if (const float val = source(point.x, point.y); std::abs(val - minValue) < WP_MARKER_EPSILON)
						markers(point.x, point.y) = 1;
			}

			// 3) Keep only largest cell
#if PREFER_CELL_CENTER
			float selectedValue = FLT_MAX;
#else
			int selectedValue = 0;
#endif
			size_t selectedIndex = 0;
			std::vector<std::vector<glm::ivec2>> componentSizes;
			{
				MEASURE_ADD_CUMULATOR(cellFilterBestTarget);
				for (const auto& point : cellPoints)
				{
					if (markers(point.x, point.y) == 1)
					{
#if PREFER_CELL_CENTER
						/* VERSION WITH CLOSEST Connected Component */
						float minDistance = FLT_MAX;
						std::vector<glm::ivec2> componentElements;
						{
							MEASURE_ADD_CUMULATOR(iterateSubCellComponents);
							iterateConnectedComponent(markers, 1, point, [&](const glm::ivec2& pt)
								{
									componentElements.emplace_back(pt);
									const auto delta = pt - center;
									const auto distance = std::sqrt(delta.x * delta.x + delta.y * delta.y);
									if (distance < minDistance)
										minDistance = distance;
								});
						}
						if (minDistance < selectedValue)
						{
							selectedValue = minDistance;
							selectedIndex = componentSizes.size();
						}
#else
						/* VERSION WITH LARGEST Connected Component */
						std::vector<glm::ivec2> componentElements;
						{
							MEASURE_ADD_CUMULATOR(iterateSubCellComponents);
							iterateConnectedComponent(markers, 1, point, [&](const glm::ivec2& pt)
								{
									componentElements.emplace_back(pt);
								});
						}
						if (componentElements.size() > selectedValue)
						{
							selectedValue = componentElements.size();
							selectedIndex = componentSizes.size();
						}
#endif
						componentSizes.emplace_back(componentElements);
					}
				}
			}

			if (componentSizes.empty())
				componentSizes = {{clampPointToImage(source, center)}};
			for (const auto& point : cellPoints)
				markers(point.x, point.y) = 0;
			for (const auto& ccPoint : componentSizes[selectedIndex])
				markers(ccPoint.x, ccPoint.y) = cellIndex;
			cellIndex++;
		}

		return markers;
	}


	LibTIM::Image<LibTIM::TLabel> waterpixel(const LibTIM::Image<LibTIM::U8>& grayScaleImage,
	                                         const std::vector<glm::ivec2> cellCenters, float sigma, float k,
	                                         float cellScale)
	{
		// Generate a voronoi graph from source points
		VoronoiGraph voronoi;
		{
			MEASURE_DURATION(voronoiCells, "Generate voronoi cells");
			voronoi = VoronoiGraph(grayScaleImage.getSizeX(), grayScaleImage.getSizeY(), cellCenters);
		}

		// Move to the derivative space
		LibTIM::Image<LibTIM::U8> gradient;
		{
			LibTIM::FlatSE filter;
			filter.make2DN4();
			MEASURE_DURATION(grad, "Compute image gradient");
			gradient = morphologicalGradient(grayScaleImage, filter);
		}

		// This will serve as guide to the watershed algorithm
		LibTIM::Image<LibTIM::U8> gradientWithRegularization;
		{
			MEASURE_DURATION(gradReg, "Add spatial regularization");
			gradientWithRegularization = spatialRegularization(gradient, voronoi, sigma, k);
		}

		// Generate watershed origins by finding the lowest connected component for each voronoi cell
		LibTIM::Image<LibTIM::TLabel> watershedSources;
		{
			MEASURE_DURATION(watMark, "Generate watershed markers");
			watershedSources = makeWatershedMarkers(gradient, voronoi, sigma, cellScale);
		}

		// Finally run watershed-meyer algorithm on markers
		MEASURE_DURATION(watMark, "Run watershed-meyer algorithm");
		LibTIM::FlatSE connectivity;
		connectivity.make2DN4();
		LibTIM::watershedMeyer<uint8_t>(gradientWithRegularization, watershedSources, connectivity);

		return watershedSources;
	}
}
