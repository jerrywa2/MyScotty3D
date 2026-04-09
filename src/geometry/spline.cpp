
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

    // Handle edge cases: empty spline, single knot, or time out of range
    if (!any()) return T();
    if (knots.size() == 1) return knots.begin()->second;
    if (time <= knots.begin()->first) return knots.begin()->second;
    if (time >= std::prev(knots.end())->first) return std::prev(knots.end())->second;

    // Find the knots surrounding the requested time: [prev, next)
    auto nextKnot = knots.upper_bound(time);
    auto prevKnot = std::prev(nextKnot);

    float timePrev = prevKnot->first, timeNext = nextKnot->first;
    T     posPrev = prevKnot->second, posNext = nextKnot->second;

    // Resolve the knot before prevKnot, mirroring if at the beginning
    float timeBeforePrev; T posBeforePrev;
    if (prevKnot == knots.begin()) {
        timeBeforePrev = timePrev - (timeNext - timePrev);
        posBeforePrev = posPrev - (posNext - posPrev);
    }
    else {
        auto beforePrevKnot = std::prev(prevKnot);
        timeBeforePrev = beforePrevKnot->first;
        posBeforePrev = beforePrevKnot->second;
    }

    // Resolve the knot after nextKnot, mirroring if at the end
    float timeAfterNext; T posAfterNext;
    auto afterNextKnot = std::next(nextKnot);
    if (afterNextKnot == knots.end()) {
        timeAfterNext = timeNext + (timeNext - timePrev);
        posAfterNext = posNext + (posNext - posPrev);
    }
    else {
        timeAfterNext = afterNextKnot->first;
        posAfterNext = afterNextKnot->second;
    }

    // Compute Catmull-Rom tangents scaled to the [prevKnot, nextKnot] interval
    float segmentDuration = timeNext - timePrev;
    T tangentAtPrev = (posNext - posBeforePrev) / (timeNext - timeBeforePrev) * segmentDuration;
    T tangentAtNext = (posAfterNext - posPrev) / (timeAfterNext - timePrev) * segmentDuration;

    // Normalize time to [0, 1] within the segment and evaluate
    float normalizedTime = (time - timePrev) / segmentDuration;
    return cubic_unit_spline(normalizedTime, posPrev, posNext, tangentAtPrev, tangentAtNext);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.
	
	float t_2 = time * time;
	float t_3 = time * time * time;

	float h00 = 2 * t_3 - 3 * t_2 + 1;
	float h10 = t_3 - 2 * t_2 + time;
	float h01 = -2 * t_3 + 3 * t_2;
	float h11 = t_3 - t_2;

	T spline = h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;

	return spline;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
