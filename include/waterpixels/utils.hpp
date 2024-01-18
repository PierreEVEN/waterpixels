#pragma once

#include <chrono>
#include <mutex>
#include <libtim/Common/Types.h>
#include <libtim/Common/Image.h>
#include <glm/glm.hpp>
#include "config.hpp"

#ifndef ENABLE_PROFILER
#define ENABLE_PROFILER true
#endif // ENABLE_PROFILER

#if ENABLE_PROFILER
#define MEASURE_DURATION(name, description) WP::Profiler name(description)
#define MEASURE_CUMULATIVE_DURATION(name, description) WP::ProfilerCumulator name(description, false)
#define MEASURE_AVERAGE_DURATION(name, description) WP::ProfilerCumulator name(description, true)
#define MEASURE_ADD_CUMULATOR(name) WP::ProfilerCumulator::Instance name##_inst(name)
#else
#define MEASURE_DURATION(name, description)
#define MEASURE_CUMULATIVE_DURATION(name, description) 
#define MEASURE_AVERAGE_DURATION(name, description) 
#define MEASURE_ADD_CUMULATOR(name) 
#endif

namespace WP
{
	/******** 2D GRID GENERATION ********/

	// Generate a rectangular 2D grid (each point spaced with a distance of sigma)
	std::vector<glm::ivec2> makeRectGrid2D(uint32_t width, uint32_t height, float sigma);
	// Generate a hexagonal 2D grid (each point spaced with a distance of sigma)
	std::vector<glm::ivec2> makeHexGrid2D(uint32_t width, uint32_t height, float sigma);
	// Generate a rectangular 2D grid where each point has a random position in its cell
	std::vector<glm::ivec2> makeRandGrid2D(uint32_t width, uint32_t height, float sigma);


	/******** IMAGE UTILITIES ********/

	glm::vec3 rgbToCIELAB(LibTIM::RGB& pixelRGB);

	// if image(x) > 0 then label(x) == 1
	LibTIM::Image<LibTIM::U8> labelToImage(const LibTIM::Image<LibTIM::TLabel>& image);
	// if label(x) > 0 then image(x) == 255
	LibTIM::Image<LibTIM::U8> labelToBinaryImage(const LibTIM::Image<LibTIM::TLabel>& image);
	// image(x) == label(x)
	LibTIM::Image<LibTIM::TLabel> imageToBinaryLabel(const LibTIM::Image<LibTIM::U8>& image);

	/******** MORPHOLOGICAL OPERATIONS ********/

	// Apply a 3x3 sobel filter
	LibTIM::Image<LibTIM::U8> sobelFilter(LibTIM::Image<LibTIM::U8> image);

	/******** VORONOI GRAPHS ********/

	class VoronoiGraph
	{
	public:
		VoronoiGraph();
		VoronoiGraph(size_t width, size_t height, const std::vector<glm::ivec2>& centers);

		// Get each voronoi cell of the graph. For each cell, first = cell center / second = list of cell points
		[[nodiscard]] const std::vector<std::pair<glm::ivec2, std::vector<glm::ivec2>>>& cells() const
		{
			return voronoiCells;
		}

		// Red channel = centers, Green channel = cell delimitation
		[[nodiscard]] LibTIM::Image<LibTIM::RGB> debugVisualization() const;
	private:
		std::vector<std::pair<glm::ivec2, std::vector<glm::ivec2>>> voronoiCells;
		size_t width;
		size_t height;
	};

	/******** PROFILING TOOLS ********/

	class ProfilerCumulator
	{
	public:
		ProfilerCumulator(std::string name, bool average);
		~ProfilerCumulator();

		void append(double time)
		{
			std::lock_guard l(m);
			total += time;
			calls++;
		}

		class Instance
		{
		public:
			Instance(ProfilerCumulator& _parent);

			~Instance();

		private:
			const std::chrono::steady_clock::time_point start;
			ProfilerCumulator& parent;
		};

	private:
		const std::string description;
		double total = 0;
		int calls = 0;
		bool average;
		std::mutex m;
	};

	class Profiler
	{
	public:
		Profiler(std::string name);
		~Profiler();
		double getDuration();
		static void printDuration(const std::string& description, double duration);
		static std::string makeIndent();
	private:
		bool print = true;
		const std::chrono::steady_clock::time_point start;
		const std::string description;
	};
}
