#include "utils.hpp"

#include "Common/Types.h"

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
