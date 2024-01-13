#include "utils.hpp"

#include "Common/Types.h"

std::vector<glm::vec2> makePoints(uint32_t width, uint32_t height, float sigma)
{
	std::vector<glm::vec2> points;
	points.reserve(static_cast<uint32_t>((width / sigma) * (height / sigma)));
	for (uint32_t x = 0; x < static_cast<uint32_t>(width / sigma) + 1; ++x)
	{
		for (uint32_t y = 0; y < static_cast<uint32_t>(height / sigma) + 1; ++y)
		{
			points.emplace_back(glm::vec2{ x * sigma + sigma / 2, y * sigma + sigma / 2 });
		}
	}
	return points;
}

glm::vec3 rgbToCIELAB(LibTIM::RGB pixelRGB)
{
	const glm::vec3 normalizedRGB(static_cast<float>(pixelRGB[0]), static_cast<float>(pixelRGB[1]), static_cast<float>(pixelRGB[2]));

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

	return { L, a, b };
}

LibTIM::Image<LibTIM::U8> sobelFilter(LibTIM::Image<LibTIM::U8> image)
{
	int dx = image.getSizeX();
	int dy = image.getSizeY();

	LibTIM::Image<LibTIM::U8> imgSobelX{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };
	LibTIM::Image<LibTIM::U8> imgSobelY{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };
	LibTIM::Image<LibTIM::U8> resImg{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };

	float kernelX[3][3] = { {-1 - 2, -1}, {0, 0, 0}, {1, 2, 1} };
	float kernelY[3][3] = { {1, 0, -1}, {2, 0, -2}, {1, 0, -1} };
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

	LibTIM::Image<LibTIM::U8> resImg{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };

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

LibTIM::Image<LibTIM::U8> gaussianFilter3x3(LibTIM::Image<LibTIM::U8> image)
{
	std::vector<std::vector<float>> gaussian = { {1.f / 16.f, 2.f / 16.f, 1.f / 16.f}, {2.f / 16.f, 4 / 16.f, 2 / 16.f}, {1.f / 16.f, 2.f / 16.f, 1.f / 16.f} };
	return applyFilter(image, gaussian);
}

/*
* Spatial regularization according to a grid with cells of length sigma
@param image: a grayscale image
@param sigma: the size of a cell in the grid
@param k: the regularization parameter (if equals 0 then no regularization)
*/
LibTIM::Image<LibTIM::U8> spatialRegularization(LibTIM::Image<LibTIM::U8> image, float sigma, float k)
{
	int dx = image.getSizeX();
	int dy = image.getSizeY();

	LibTIM::Image<LibTIM::U8> regularizedImg{ static_cast<LibTIM::TSize>(dx), static_cast<LibTIM::TSize>(dy) };

	for (int y = 0; y < dy; y++)
	{
		for (int x = 0; x < dx; x++)
		{
			glm::vec2 closestCenter{ int(static_cast<float>(x) / sigma) * sigma + sigma / 2.f, int(static_cast<float>(y) / sigma) * sigma + sigma / 2.f };
			glm::vec2 currentPixel{ x, y };
			float d = std::max(std::abs(currentPixel.x - closestCenter.x), std::abs(currentPixel.y - closestCenter.y));

			regularizedImg(x, y) = static_cast<LibTIM::U8>(image(x, y) + k * ((2.f * d) / sigma));
		}
	}

	return regularizedImg;
}
