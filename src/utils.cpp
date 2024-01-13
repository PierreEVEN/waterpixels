#include "utils.hpp"

#include <glm/gtc/matrix_transform.hpp>

#include "Common/Types.h"

glm::vec3 rgbToCIELAB(LibTIM::RGB pixelRGB)
{
	const glm::vec3 normalizedRGB(static_cast<float>(pixelRGB[0]), static_cast<float>(pixelRGB[1]), static_cast<float>(pixelRGB[2]));

	constexpr glm::mat3 rgbToXYZ(
		0.618f, 0.299f, 0.0f, // First column
		0.177f, 0.587f, 0.056f,
		0.205f, 0.114f, 0.944f);

	glm::vec3 pixelXYZ = rgbToXYZ * normalizedRGB;

	constexpr  glm::vec3 XYZn = rgbToXYZ * glm::vec3(255.f);

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
			imgSobelX(x, y) = resX;
			imgSobelY(x, y) = resY;

			resImg(x, y) = std::sqrt(resX * resX + resY * resY);
		}
	}
	return resImg;
}