#include "waterpixels/utils.hpp"

#include <corecrt_math_defines.h>
#include <vector>

#include <glm/glm.hpp>
#include <libtim/Common/Types.h>

namespace WP
{
	std::vector<glm::ivec2> makePoints(uint32_t width, uint32_t height, float sigma)
	{
		std::vector<glm::ivec2> points;
		points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));
		for (uint32_t x = 0; x < static_cast<uint32_t>(width / sigma) + 1; ++x)
		{
			for (uint32_t y = 0; y < static_cast<uint32_t>(height / sigma) + 1; ++y)
			{
				points.emplace_back(glm::ivec2{x * sigma + sigma / 2, y * sigma + sigma / 2});
			}
		}
		return points;
	}

	glm::vec3 rgbToCIELAB(LibTIM::RGB pixelRGB)
	{
		const glm::vec3 normalizedRGB(static_cast<float>(pixelRGB[0]), static_cast<float>(pixelRGB[1]),
		                              static_cast<float>(pixelRGB[2]));

		glm::mat3 rgbToXYZ(
			0.618f, 0.299f, 0.0f, // First column
			0.177f, 0.587f, 0.056f,
			0.205f, 0.114f, 0.944f);


		glm::vec3 pixelXYZ = rgbToXYZ * normalizedRGB;

		glm::vec3 XYZn = rgbToXYZ * glm::vec3(255.f);

		float L;

		if (pixelXYZ.y / XYZn.y > 0.008856)
			L = 116.f * std::pow(pixelXYZ.y / XYZn.y, 1.f / 3.f) - 16;
		else
			L = 903.3f * pixelXYZ.y / XYZn.y;

		auto f = [](float t)
		{
			if (t > 0.008856)
				return std::pow(t, 1.f / 3.f);
			else
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

	LibTIM::Image<LibTIM::U8> applyFilter(LibTIM::Image<LibTIM::U8> image, std::vector<std::vector<float>> kernel)
	{
		int dx = image.getSizeX();
		int dy = image.getSizeY();

		LibTIM::Image<LibTIM::U8> resImg{static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy)};

		for (int y = 0; y < dy; y++)
		{
			for (int x = 0; x < dx; x++)
			{
				int startY = y - int(kernel[0].size() / 2);
				int startX = x - int(kernel.size() / 2);
				int endY = y + int(kernel[0].size() / 2);
				int endX = x + int(kernel.size() / 2);

				float res = 0;

				for (int i = startY; i <= endY; i++)
				{
					if (i < 0 || i >= dy)
						continue;
					for (int j = startX; j <= endX; j++)
					{
						if (j < 0 || j >= dx)
							continue;
						LibTIM::U8 p = image(j, i);
						res += kernel[j - startX][i - startY] * p;
					}
				}
				resImg(x, y) = std::clamp(static_cast<int>(res), 0, 255);
			}
		}
		return resImg;
	}

	LibTIM::Image<LibTIM::U8> gaussianFilter(LibTIM::Image<LibTIM::U8> image, int radius, float sigma)
	{
		double s = 2.0 * sigma * sigma;
		float sum = 0;
		std::vector<std::vector<float>> filter;
		for (int64_t x = -radius; x <= radius; ++x)
		{
			std::vector<float> values;
			for (int64_t y = -radius; y <= radius; ++y)
			{
				double r = sqrt(x * x + y * y);
				double kernel = (exp(-(r * r) / s)) / (M_PI * s);
				values.emplace_back(kernel);
				sum += kernel;
			}
			filter.emplace_back(values);
		}
		for (auto& row : filter)
			for (auto& val : row)
				val /= sum;

		return applyFilter(image, filter);
	}
	
	LibTIM::Image<LibTIM::U8> labelToImage(const LibTIM::Image<LibTIM::TLabel>& image)
	{
		LibTIM::Image<LibTIM::TLabel> binaryImage(image.getSizeX(), image.getSizeY());
		for (int x = 0; x < binaryImage.getSizeX(); x++)
			for (int y = 0; y < binaryImage.getSizeY(); y++)
				binaryImage(x, y) = image(x, y) ;
		return binaryImage;
	}

	LibTIM::Image<LibTIM::U8> labelToBinaryImage(const LibTIM::Image<LibTIM::TLabel>& image)
	{
		LibTIM::Image<LibTIM::TLabel> binaryImage(image.getSizeX(), image.getSizeY());
		for (int x = 0; x < binaryImage.getSizeX(); x++)
			for (int y = 0; y < binaryImage.getSizeY(); y++)
				binaryImage(x, y) = image(x, y) ? 255 : 0;
		return binaryImage;
	}

	LibTIM::Image<LibTIM::TLabel> imageToBinaryLabel(const LibTIM::Image<LibTIM::U8>& image)
	{
		LibTIM::Image<LibTIM::U8> binaryLabels(image.getSizeX(), image.getSizeY());
		for (int x = 0; x < image.getSizeX(); x++)
			for (int y = 0; y < image.getSizeY(); y++)
				binaryLabels(x, y) = image(x, y) ? 1 : 0;
		return binaryLabels;
	}
}
