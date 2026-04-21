
#include "particles.h"

bool Particles::Particle::update(const PT::Aggregate &scene, Vec3 const &gravity, const float radius, const float dt) {

	//A4T4: particle update

	// Compute the trajectory of this particle for the next dt seconds.

	// (1) Build a ray representing the particle's path as if it travelled at constant velocity.

	// (2) Intersect the ray with the scene and account for collisions. Be careful when placing
	// collision points using the particle radius. Move the particle to its next position.

	// (3) Account for acceleration due to gravity after updating position.

	// (4) Repeat until the entire time step has been consumed.

	// (5) Decrease the particle's age and return 'false' if it should be removed.

	float remaining = dt;

	while (remaining > EPS_F) {

		// --- Build the ray for this sub-step ---
		float speed = velocity.norm();

		// If the particle is essentially stationary, just apply gravity and break
		if (speed < EPS_F) {
			position += velocity * remaining;
			velocity += gravity * remaining;
			break;
		}

		Vec3 dir = velocity / speed;           // unit direction
		//float max_dist = speed * remaining;    // how far we could travel

		Ray ray(position, velocity.unit(), Vec2{EPS_F, remaining * velocity.norm() + radius});

		PT::Trace hit = scene.hit(ray);

		if (hit.hit) {
			// Compute sphere-adjusted collision distance
			//float sin_angle = std::abs(dot(-dir, hit.normal));
			//sin_angle = std::max(sin_angle, EPS_F);

			float adjusted_dist = hit.distance - radius;

			//if (adjusted_dist < 0.0f) {
			//	// Sphere is already intersecting — don't move, just reflect outward
			//	velocity = velocity - 2.0f * dot(velocity, hit.normal) * hit.normal;
			//	remaining -= EPS_F;
			//	continue;
			//}

			// Time consumed reaching the collision point
			float time_used = adjusted_dist / speed;

			// Move particle to the collision point
			position += dir * adjusted_dist;

			velocity = velocity - 2.0f * dot(velocity, hit.normal) * hit.normal;
			velocity += gravity * time_used;
			remaining -= time_used;

			// Reflect velocity about the surface normal (elastic collision)
			

		}
		else {
			// No collision — travel the full remaining distance, apply gravity, done
			position += velocity * remaining;
			velocity += gravity * remaining;
			break;
		}
	}
	//printf("pos: %f %f %f, vel: %f %f %f\n", position.x, position.y, position.z, velocity.x, velocity.y, velocity.z);
	// Decrease age and return whether the particle should live on
	age -= dt;
	return age > 0.0f;
}

void Particles::advance(const PT::Aggregate& scene, const Mat4& to_world, float dt) {

	if(step_size < EPS_F) return;

	step_accum += dt;

	while(step_accum > step_size) {
		step(scene, to_world);
		step_accum -= step_size;
	}
}

void Particles::step(const PT::Aggregate& scene, const Mat4& to_world) {

	std::vector<Particle> next;
	next.reserve(particles.size());

	for(Particle& p : particles) {
		if(p.update(scene, gravity, radius, step_size)) {
			next.emplace_back(p);
		}
	}

	if(rate > 0.0f) {

		//helpful when emitting particles:
		float cos = std::cos(Radians(spread_angle) / 2.0f);

		//will emit particle i when i == time * rate
		//(i.e., will emit particle when time * rate hits an integer value.)
		//so need to figure out all integers in [current_step, current_step+1) * step_size * rate
		//compute the range:
		double begin_t = current_step * double(step_size) * double(rate);
		double end_t = (current_step + 1) * double(step_size) * double(rate);

		uint64_t begin_i = uint64_t(std::max(0.0, std::ceil(begin_t)));
		uint64_t end_i = uint64_t(std::max(0.0, std::ceil(end_t)));

		//iterate all integers in [begin, end):
		for (uint64_t i = begin_i; i < end_i; ++i) {
			//spawn particle 'i':

			float y = lerp(cos, 1.0f, rng.unit());
			float t = 2 * PI_F * rng.unit();
			float d = std::sqrt(1.0f - y * y);
			Vec3 dir = initial_velocity * Vec3(d * std::cos(t), y, d * std::sin(t));

			Particle p;
			p.position = to_world * Vec3(0.0f, 0.0f, 0.0f);
			p.velocity = to_world.rotate(dir);
			p.age = lifetime; //NOTE: could adjust lifetime based on index
			next.push_back(p);
		}
	}

	particles = std::move(next);
	current_step += 1;
}

void Particles::reset() {
	particles.clear();
	step_accum = 0.0f;
	current_step = 0;
	rng.seed(seed);
}

bool operator!=(const Particles& a, const Particles& b) {
	return a.gravity != b.gravity
	|| a.radius != b.radius
	|| a.initial_velocity != b.initial_velocity
	|| a.spread_angle != b.spread_angle
	|| a.lifetime != b.lifetime
	|| a.rate != b.rate
	|| a.step_size != b.step_size
	|| a.seed != b.seed;
}
