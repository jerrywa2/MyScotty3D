
#include "shape.h"
#include "../geometry/util.h"

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const {
	//A3T2 - sphere hit

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

	Vec3 o = ray.point;
	Vec3 d = ray.dir;

	float denom = 2 * dot(d, d);
	float radicand = 4 * pow(dot(o, d), 2) - 4 * dot(d, d) * (dot(o, o) - pow(radius, 2));
	float t1, t2, t;
	bool hitPresent = false;
	if (denom != 0.0f && radicand >= 0.0f)
	{
		t1 = ( -2 * dot(o, d) - sqrt(radicand) ) / denom;
		t2 = ( -2 * dot(o, d) + sqrt(radicand) ) / denom;
		hitPresent = (t1 >= ray.dist_bounds.x && t1 <= ray.dist_bounds.y) || (t2 >= ray.dist_bounds.x && t2 <= ray.dist_bounds.y);
		if (t1 >= ray.dist_bounds.x && t1 <= ray.dist_bounds.y) { t = t1; }
		else if (t2 >= ray.dist_bounds.x && t2 <= ray.dist_bounds.y) { t = t2; }
	}
	Vec3 hitPoint = o + t * d;

    PT::Trace ret;
    ret.origin = ray.point;
    ret.hit = hitPresent;       // was there an intersection?
    ret.distance = hitPresent ? t : 0.0f;   // at what distance did the intersection occur?
	ret.position = hitPresent ? hitPoint : Vec3{}; // where was the intersection?
	ret.normal = hitPresent ? hitPoint.unit() : Vec3{};   // what was the surface normal at the intersection?
	ret.uv = hitPresent ? uv(hitPoint) : Vec2{}; 	   // what was the uv coordinates at the intersection? (you may find Sphere::uv to be useful)
    return ret;
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}
