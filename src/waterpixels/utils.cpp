#include "waterpixels/utils.hpp"

#include <corecrt_math_defines.h>
#include <vector>

#include <glm/glm.hpp>
#include <libtim/Common/Types.h>

#include "Algorithms/Morphology.hxx"

namespace WP
{
	std::vector<glm::ivec2> makeRectGrid2D(uint32_t width, uint32_t height, float sigma)
	{
		std::vector<glm::ivec2> points;
		points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));

		for (uint32_t x = sigma / 2; x < width + sigma; x += sigma)
			for (uint32_t y = sigma / 2; y < height + sigma; y += sigma)
				points.emplace_back(glm::ivec2{x, y});
		return points;
	}

	std::vector<glm::ivec2> makeHexGrid2D(uint32_t width, uint32_t height, float sigma)
	{
		std::vector<glm::ivec2> points;
		points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));

		for (uint32_t x = sigma / 2; x < width + sigma; x += sigma)
			for (uint32_t y = x % 2 == 0 ? sigma / 2 : 0; y < height + sigma; y += sigma)
				points.emplace_back(glm::ivec2{x, y});
		return points;
	}

	std::vector<glm::ivec2> makeRandGrid2D(uint32_t width, uint32_t height, float sigma)
	{
		std::vector<glm::ivec2> points;
		points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));

		const auto randFloatInRange = [](float min, float max)
		{
			return rand() / static_cast<float>(RAND_MAX) * (max - min) - min;
		};

		for (uint32_t x = sigma / 2; x < width + sigma; x += sigma)
			for (uint32_t y = sigma / 2; y < height + sigma; y += sigma)
				points.emplace_back(glm::ivec2{
					x + randFloatInRange(-sigma / 2, sigma / 2), y + randFloatInRange(-sigma / 2, sigma / 2)
				});
		return points;
	}

	glm::vec3 rgbToCIELAB(LibTIM::RGB& pixelRGB)
	{
		const glm::vec3 normalizedRGB(static_cast<float>(pixelRGB[0]), static_cast<float>(pixelRGB[1]),
		                              static_cast<float>(pixelRGB[2]));

		const glm::mat3 rgbToXYZ(
			0.618f, 0.299f, 0.0f, // First column
			0.177f, 0.587f, 0.056f,
			0.205f, 0.114f, 0.944f);


		const glm::vec3 pixelXYZ = rgbToXYZ * normalizedRGB;

		const glm::vec3 XYZn = rgbToXYZ * glm::vec3(255.f);

		float L;

		if (pixelXYZ.y / XYZn.y > 0.008856)
			L = 116.f * std::pow(pixelXYZ.y / XYZn.y, 1.f / 3.f) - 16;
		else
			L = 903.3f * pixelXYZ.y / XYZn.y;

		auto f = [](float t)
		{
			if (t > 0.008856)
				return std::pow(t, 1.f / 3.f);
			
			return 7.7787f * t + 16.f / 116.f;
		};

		float a = 500.f * (f(pixelXYZ.x / XYZn.x) - f(pixelXYZ.y / XYZn.y));
		float b = 200.f * (f(pixelXYZ.y / XYZn.y) - f(pixelXYZ.z / XYZn.z));

		return {L, a, b};
	}

	LibTIM::Image<LibTIM::U8> sobelFilter(LibTIM::Image<LibTIM::U8> image)
	{
		int dx = image.getSizeX();
		int dy = image.getSizeY();

		LibTIM::Image<LibTIM::U8> imgSobelX{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};
		LibTIM::Image<LibTIM::U8> imgSobelY{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};
		LibTIM::Image<LibTIM::U8> resImg{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};

		float kernelX[3][3] = {{-1 - 2, -1}, {0, 0, 0}, {1, 2, 1}};
		float kernelY[3][3] = {{1, 0, -1}, {2, 0, -2}, {1, 0, -1}};
		#pragma omp parallel for collapse(2)
		for (int y = 0; y < dy; y++)
		{
			for (int x = 0; x < dx; x++)
			{
				int startY = y - 1;
				int startX = x - 1;
				int endY = y + 1;
				int endX = x + 1;

				float resX = 0.0f;
				float resY = 0.0f;
				for (int i = startY; i <= endY; i++)
				{
					if (i < 0 || i >= dy)
						continue;
					for (int j = startX; j <= endX; j++)
					{
						if (j < 0 || j >= dx)
							continue;
						float p = image(j, i);
						resX += kernelX[j - startX][i - startY] * p;
						resY += kernelY[j - startX][i - startY] * p;
					}
				}
				imgSobelX(x, y) = static_cast<LibTIM::U8>(resX);
				imgSobelY(x, y) = static_cast<LibTIM::U8>(resY);

				resImg(x, y) = static_cast<LibTIM::U8>(std::sqrt(resX * resX + resY * resY));
			}
		}
		return resImg;
	}

	LibTIM::Image<LibTIM::U8> labelToImage(const LibTIM::Image<LibTIM::TLabel>& image)
	{
		LibTIM::Image<LibTIM::TLabel> binaryImage(image.getSizeX(), image.getSizeY());
		#pragma omp parallel for collapse(2)
		for (int x = 0; x < binaryImage.getSizeX(); x++)
			for (int y = 0; y < binaryImage.getSizeY(); y++)
				binaryImage(x, y) = image(x, y);
		return binaryImage;
	}

	LibTIM::Image<LibTIM::U8> labelToBinaryImage(const LibTIM::Image<LibTIM::TLabel>& image)
	{
		LibTIM::Image<LibTIM::TLabel> binaryImage(image.getSizeX(), image.getSizeY());
		#pragma omp parallel for collapse(2)
		for (int x = 0; x < binaryImage.getSizeX(); x++)
			for (int y = 0; y < binaryImage.getSizeY(); y++)
				binaryImage(x, y) = image(x, y) ? 255 : 0;
		return binaryImage;
	}

	LibTIM::Image<LibTIM::TLabel> imageToBinaryLabel(const LibTIM::Image<LibTIM::U8>& image)
	{
		LibTIM::Image<LibTIM::U8> binaryLabels(image.getSizeX(), image.getSizeY());
		#pragma omp parallel for collapse(2)
		for (int x = 0; x < image.getSizeX(); x++)
			for (int y = 0; y < image.getSizeY(); y++)
				binaryLabels(x, y) = image(x, y) ? 1 : 0;
		return binaryLabels;
	}

	VoronoiGraph::VoronoiGraph() : width(0), height(0)
	{
	}

	VoronoiGraph::VoronoiGraph(size_t _width, size_t _height, const std::vector<glm::ivec2>& centers) :
		voronoiCells(centers.size()), width(_width), height(_height)
	{
		for (size_t i = 0; i < centers.size(); ++i)
			voronoiCells[i] = {centers[i], {}};

		// @TODO : accelerate using a grid map
		const auto getClosestCenter = [&](const glm::ivec2& pos) -> unsigned long
		{
			float closestDistance = FLT_MAX;
			size_t closestIndex = 0;
			for (size_t i = 0; i < centers.size(); ++i)
			{
				const auto delta = centers[i] - pos;
				const auto distance = static_cast<float>(std::sqrt(delta.x * delta.x + delta.y * delta.y));
				if (distance < closestDistance)
				{
					closestDistance = distance;
					closestIndex = i;
				}
			}
			return static_cast<unsigned long>(closestIndex);
		};

		for (int64_t x = 0; x < width; ++x)
			for (int64_t y = 0; y < height; ++y)
				voronoiCells[getClosestCenter({x, y})].second.emplace_back(glm::ivec2{x, y});
	}

	LibTIM::Image<LibTIM::RGB> VoronoiGraph::debugVisualization() const
	{
		LibTIM::Image<LibTIM::RGB> result(width, height);
		result.fill(LibTIM::RGB({0, 0, 0}));
		LibTIM::Image<LibTIM::U8> stamp(width, height);
		int ccId = 1;
		for (const auto& cc : voronoiCells)
		{
			for (const auto& pixel : cc.second)
				stamp(pixel.x, pixel.y) = ccId;

			if (result.isPosValid(cc.first.x, cc.first.y))
				result(cc.first.x, cc.first.y)[0] = 255;

			ccId++;
		}

		LibTIM::FlatSE connectivity;
		connectivity.make2DN4();
		const auto borders = morphologicalGradient(stamp, connectivity);

		#pragma omp parallel for collapse(2)
		for (int64_t x = 0; x < width; ++x)
			for (int64_t y = 0; y < height; ++y)
				result(x, y)[1] = borders(x, y) ? 255 : 0;
		return result;
	}

	static int profilerIndent = 0;
	static std::mutex profilerIndentMutex;

	Profiler::~Profiler()
	{
		if (print)
			Profiler::printDuration(description, getDuration());
		std::lock_guard m(profilerIndentMutex);
		profilerIndent--;
	}

	double Profiler::getDuration()
	{
		print = false;
		return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count() /
			1000.0;
	}

	void Profiler::printDuration(const std::string& description, double duration)
	{
		std::cout << Profiler::makeIndent() << duration << "ms elapsed : " << description << std::endl;
	}

	std::string Profiler::makeIndent()
	{
		std::string str = " ";
		std::lock_guard m(profilerIndentMutex);
		for (int i = 0; i < profilerIndent; ++i)
			str = "--" + str;
		return str;
	}

	Profiler::Profiler(std::string _description) : start(std::chrono::steady_clock::now()),
	                                               description(std::move(_description))
	{
		std::lock_guard m(profilerIndentMutex);
		profilerIndent++;
	}

	ProfilerCumulator::~ProfilerCumulator()
	{
		if (average)
			std::cout << Profiler::makeIndent() << total / calls << "ms elapsed in average : " << description <<
				" (" << calls <<
				" calls)" << std::endl;
		else
			std::cout << Profiler::makeIndent() << total << "ms elapsed in total : " << description << " (" <<
				calls << " calls)" << std::endl;
		std::lock_guard m(profilerIndentMutex);
		profilerIndent--;
	}

	ProfilerCumulator::Instance::Instance(ProfilerCumulator& _parent)
		: start(std::chrono::steady_clock::now()), parent(_parent)
	{
		std::lock_guard m(profilerIndentMutex);
		profilerIndent++;
	}

	ProfilerCumulator::Instance::~Instance()
	{
		parent.append(
			std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).
			count() / 1000.0);
		std::lock_guard m(profilerIndentMutex);
		profilerIndent--;
	}

	ProfilerCumulator::ProfilerCumulator(std::string _description, bool _average) :
		description(std::move(_description)), average(_average)
	{
		std::lock_guard m(profilerIndentMutex);
		profilerIndent++;
	}
}
