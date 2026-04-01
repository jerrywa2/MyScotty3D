
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
	nodes.clear();
	primitives = std::move(prims);

	// Construct a BVH from the given vector of primitives and maximum leaf
	// size configuration.

	if (primitives.empty()) return;

	constexpr size_t NUM_BUCKETS = 8;

	std::function<size_t(size_t, size_t)> build_node = [&](size_t start, size_t size) -> size_t {

		// Compute bbox over all primitives in this range
		BBox node_bbox;
		for (size_t i = start; i < start + size; i++)
			node_bbox.enclose(primitives[i].bbox());

		// Base case: make a leaf
		if (size <= max_leaf_size) {
			return new_node(node_bbox, start, size, 0, 0);
		}

		float best_cost = std::numeric_limits<float>::infinity();
		int   best_axis = -1;
		size_t best_split = 0;

		float parent_area = node_bbox.surface_area();

		for (int axis = 0; axis < 3; axis++) {

			float cmin = std::numeric_limits<float>::infinity();
			float cmax = -std::numeric_limits<float>::infinity();
			for (size_t i = start; i < start + size; i++) {
				float c = primitives[i].bbox().center()[axis];
				cmin = std::min(cmin, c);
				cmax = std::max(cmax, c);
			}
			if (cmin == cmax) continue;

			std::array<SAHBucketData, 8> buckets = {};

			for (size_t i = start; i < start + size; i++) {
				float c = primitives[i].bbox().center()[axis];
				size_t b = (size_t)(NUM_BUCKETS * (c - cmin) / (cmax - cmin));
				if (b == NUM_BUCKETS) b = NUM_BUCKETS - 1;
				buckets[b].bb.enclose(primitives[i].bbox());
				buckets[b].num_prims++;
			}

			for (size_t split = 1; split < NUM_BUCKETS; split++) {
				BBox left_bb, right_bb;
				size_t left_count = 0, right_count = 0;

				for (size_t b = 0; b < split; b++) {
					left_bb.enclose(buckets[b].bb);
					left_count += buckets[b].num_prims;
				}
				for (size_t b = split; b < NUM_BUCKETS; b++) {
					right_bb.enclose(buckets[b].bb);
					right_count += buckets[b].num_prims;
				}

				if (left_count == 0 || right_count == 0) continue;

				float cost = (left_bb.surface_area() * (float)left_count +
					right_bb.surface_area() * (float)right_count)
					/ parent_area;

				if (cost < best_cost) {
					best_cost = cost;
					best_axis = axis;
					best_split = split;
				}
			}
		}

		if (best_axis == -1) {
			return new_node(node_bbox, start, size, 0, 0);
		}

		float cmin = std::numeric_limits<float>::infinity();
		float cmax = -std::numeric_limits<float>::infinity();
		for (size_t i = start; i < start + size; i++) {
			float c = primitives[i].bbox().center()[best_axis];
			cmin = std::min(cmin, c);
			cmax = std::max(cmax, c);
		}

		auto mid_it = std::partition(
			primitives.begin() + start,
			primitives.begin() + start + size,
			[&](const Primitive& p) {
				float c = p.bbox().center()[best_axis];
				size_t b = (size_t)(NUM_BUCKETS * (c - cmin) / (cmax - cmin));
				if (b == NUM_BUCKETS) b = NUM_BUCKETS - 1;
				return b < best_split;
			}
		);

		size_t mid = (size_t)std::distance(primitives.begin(), mid_it);

		if (mid == start || mid == start + size)
			mid = start + size / 2;

		size_t node_idx = new_node(node_bbox, start, size, 0, 0);

		size_t left_idx = build_node(start, mid - start);
		size_t right_idx = build_node(mid, start + size - mid);

		nodes[node_idx].l = left_idx;
		nodes[node_idx].r = right_idx;

		return node_idx;
		};

	root_idx = build_node(0, primitives.size());
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
    //Trace ret;
    //for(const Primitive& prim : primitives) {
    //    Trace hit = prim.hit(ray);
    //    ret = Trace::min(ret, hit);
    //}
    //return ret;

	Trace ret;
	if (nodes.empty()) return ret;

	std::function<void(size_t)> traverse = [&](size_t node_idx) {
		const Node& node = nodes[node_idx];

		Vec2 times = ray.dist_bounds;
		if (!node.bbox.hit(ray, times)) return;

		if (node.is_leaf()) {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				Trace t = primitives[i].hit(ray);
				ret = Trace::min(ret, t);
			}
			return;
		}

		// Visit closer child first for better early termination
		Vec2 t_left = ray.dist_bounds, t_right = ray.dist_bounds;
		bool hit_left = nodes[node.l].bbox.hit(ray, t_left);
		bool hit_right = nodes[node.r].bbox.hit(ray, t_right);

		if (hit_left && hit_right) {
			// Visit closer child first
			size_t first = t_left.x <= t_right.x ? node.l : node.r;
			size_t second = t_left.x <= t_right.x ? node.r : node.l;
			float  second_t = t_left.x <= t_right.x ? t_right.x : t_left.x;
			traverse(first);
			// Skip second child if we already found something closer
			if (!ret.hit || second_t < ret.distance)
				traverse(second);
		}
		else if (hit_left) {
			traverse(node.l);
		}
		else if (hit_right) {
			traverse(node.r);
		}
		};

	traverse(root_idx);
	return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
