
#include "samplers.h"
#include "../util/rand.h"

constexpr bool IMPORTANCE_SAMPLING = true;

namespace Samplers {

Vec2 Rect::sample(RNG &rng) const {
	//A3T1 - step 2 - supersampling

    // Return a point selected uniformly at random from the rectangle [0,size.x)x[0,size.y)
    // Useful function: rng.unit()

	Vec2 randomPoint = Vec2(rng.unit() * size.x, rng.unit() * size.y);

    return randomPoint;
}

float Rect::pdf(Vec2 at) const {
	if (at.x < 0.0f || at.x > size.x || at.y < 0.0f || at.y > size.y) return 0.0f;
	return 1.0f / (size.x * size.y);
}

Vec2 Circle::sample(RNG &rng) const {
	//A3EC - bokeh - circle sampling

    // Return a point selected uniformly at random from a circle defined by its
	// center and radius.
    // Useful function: rng.unit()

    return Vec2{};
}

float Circle::pdf(Vec2 at) const {
	//A3EC - bokeh - circle pdf

	// Return the pdf of sampling the point 'at' for a circle defined by its
	// center and radius.

    return 1.f;
}

Vec3 Point::sample(RNG &rng) const {
	return point;
}

float Point::pdf(Vec3 at) const {
	return at == point ? 1.0f : 0.0f;
}

Vec3 Triangle::sample(RNG &rng) const {
	float u = std::sqrt(rng.unit());
	float v = rng.unit();
	float a = u * (1.0f - v);
	float b = u * v;
	return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

float Triangle::pdf(Vec3 at) const {
	float a = 0.5f * cross(v1 - v0, v2 - v0).norm();
	float u = 0.5f * cross(at - v1, at - v2).norm() / a;
	float v = 0.5f * cross(at - v2, at - v0).norm() / a;
	float w = 1.0f - u - v;
	if (u < 0.0f || v < 0.0f || w < 0.0f) return 0.0f;
	if (u > 1.0f || v > 1.0f || w > 1.0f) return 0.0f;
	return 1.0f / a;
}

Vec3 Hemisphere::Uniform::sample(RNG &rng) const {

	float Xi1 = rng.unit();
	float Xi2 = rng.unit();

	float theta = std::acos(Xi1);
	float phi = 2.0f * PI_F * Xi2;

	float xs = std::sin(theta) * std::cos(phi);
	float ys = std::cos(theta);
	float zs = std::sin(theta) * std::sin(phi);

	return Vec3(xs, ys, zs);
}

float Hemisphere::Uniform::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return 1.0f / (2.0f * PI_F);
}

Vec3 Hemisphere::Cosine::sample(RNG &rng) const {

	float phi = rng.unit() * 2.0f * PI_F;
	float cos_t = std::sqrt(rng.unit());

	float sin_t = std::sqrt(1 - cos_t * cos_t);
	float x = std::cos(phi) * sin_t;
	float z = std::sin(phi) * sin_t;
	float y = cos_t;

	return Vec3(x, y, z);
}

float Hemisphere::Cosine::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return dir.y / PI_F;
}

Vec3 Sphere::Uniform::sample(RNG &rng) const {
	//A3T7 - sphere sampler

    // Generate a uniformly random point on the unit sphere.
    // Tip: start with Hemisphere::Uniform

	// Sample full sphere by sampling hemisphere twice (once per hemisphere)
	float xi1 = rng.unit();
	float xi2 = rng.unit();

	float theta = std::acos(1.0f - 2.0f * xi1); // maps [0,1] to [0,pi]
	float phi = 2.0f * PI_F * xi2;

	return Vec3(
		std::sin(theta) * std::cos(phi),
		std::cos(theta),
		std::sin(theta) * std::sin(phi)
	);
}

float Sphere::Uniform::pdf(Vec3 dir) const {
	return 1.0f / (4.0f * PI_F);
}

Sphere::Image::Image(const HDR_Image& image) {
    //A3T7 - image sampler init

    // Set up importance sampling data structures for a spherical environment map image.
    // You may make use of the _pdf, _cdf, and total members, or create your own.

	const auto [_w, _h] = image.dimension();
	w = _w;
	h = _h;

	// Build PDF weighted by luminance * sin(theta) for each pixel
	_pdf.resize(w * h);
	for (uint32_t y = 0; y < h; y++) {
		// Pixel center in [0,1], then map to theta in [0, pi]
		// Note: y=0 is bottom of image, which maps to theta=pi
		float v = (float(y) + 0.5f) / float(h);
		float theta = PI_F * (1.0f - v); // flip: y=0 -> theta=pi, y=h -> theta=0
		float sin_theta = std::sin(theta);
		for (uint32_t x = 0; x < w; x++) {
			float luma = image.at(x, y).luma();
			_pdf[y * w + x] = luma * sin_theta;
		}
	}

	// Compute total and normalize PDF
	float total = 0.0f;
	for (float v : _pdf) total += v;
	if (total > 0.0f) {
		for (float& v : _pdf) v /= total;
	}

	// Build CDF
	_cdf.resize(w * h);
	_cdf[0] = _pdf[0];
	for (uint32_t i = 1; i < w * h; i++) {
		_cdf[i] = _cdf[i - 1] + _pdf[i];
	}
}

Vec3 Sphere::Image::sample(RNG &rng) const {
	if (!IMPORTANCE_SAMPLING) {
		Sphere::Uniform uniform;
		return uniform.sample(rng);
	}
	else {
		// Inversion sampling via binary search on CDF
		float xi = rng.unit();
		auto it = std::upper_bound(_cdf.begin(), _cdf.end(), xi);
		uint32_t idx = uint32_t(std::clamp(
			int(std::distance(_cdf.begin(), it)),
			0, int(w * h) - 1
		));

		uint32_t px = idx % w;
		uint32_t py = idx / w;

		// Convert pixel to (phi, theta)
		float phi = (float(px) + 0.5f) / float(w) * 2.0f * PI_F;
		float v = (float(py) + 0.5f) / float(h);
		float theta = PI_F * (1.0f - v); // flip y-axis

		return Vec3(
			std::sin(theta) * std::cos(phi),
			std::cos(theta),
			std::sin(theta) * std::sin(phi)
		);
	}
}

float Sphere::Image::pdf(Vec3 dir) const {
	if (!IMPORTANCE_SAMPLING) {
		Sphere::Uniform uniform;
		return uniform.pdf(dir);
	}
	else {
		// Convert direction to (phi, theta)
		float theta = std::acos(std::clamp(dir.y, -1.0f, 1.0f));
		float phi = std::atan2(dir.z, dir.x);
		if (phi < 0.0f) phi += 2.0f * PI_F;

		// Convert to pixel coordinates
		float v = 1.0f - theta / PI_F; // flip y-axis
		uint32_t px = uint32_t(std::clamp(phi / (2.0f * PI_F) * float(w), 0.0f, float(w) - 1.0f));
		uint32_t py = uint32_t(std::clamp(v * float(h), 0.0f, float(h) - 1.0f));

		float p = _pdf[py * w + px];

		// Apply Jacobian: (w*h) / (2*pi^2 * sin(theta))
		float sin_theta = std::sin(theta);
		if (sin_theta < 1e-6f) return 0.0f;

		return p * float(w * h) / (2.0f * PI_F * PI_F * sin_theta);
	}
}

} // namespace Samplers
