#include <unordered_set>
#include "skeleton.h"
#include "test.h"
#include <iostream>

void Skeleton::Bone::compute_rotation_axes(Vec3 *x_, Vec3 *y_, Vec3 *z_) const {
	assert(x_ && y_ && z_);
	auto &x = *x_;
	auto &y = *y_;
	auto &z = *z_;

	//y axis points in the direction of extent:
	y = extent.unit();
	//if extent is too short to normalize nicely, point along the skeleton's 'y' axis:
	if (!y.valid()) {
		y = Vec3{0.0f, 1.0f, 0.0f};
	}

	//x gets skeleton's 'x' axis projected to be orthogonal to 'y':
	x = Vec3{1.0f, 0.0f, 0.0f};
	x = (x - dot(x,y) * y).unit();
	if (!x.valid()) {
		//if y perfectly aligns with skeleton's 'x' axis, x, gets skeleton's z axis:
		x = Vec3{0.0f, 0.0f, 1.0f};
		x = (x - dot(x,y) * y).unit(); //(this should do nothing)
	}

	//z computed from x,y:
	z = cross(x,y);

	//x,z rotated by roll:
	float cr = std::cos(roll / 180.0f * PI_F);
	float sr = std::sin(roll / 180.0f * PI_F);
	// x = cr * x + sr * -z;
	// z = cross(x,y);
	std::tie(x, z) = std::make_pair(cr * x + sr * -z, cr * z + sr * x);
}

std::vector< Mat4 > Skeleton::bind_pose() const {
	//A4T2a: bone-to-skeleton transformations in the bind pose
	//(the bind pose does not rotate by Bone::pose)

	std::vector<Mat4> bind;
	bind.reserve(bones.size());

	// NOTE: bones is guaranteed to be ordered such that parents appear before children.
	for (auto const& bone : bones) {
		bool isRoot = (bone.parent == -1U);

		// A root bone translates from the skeleton base;
		// a child bone translates from its parent's tip (extent).
		Vec3 originInParentSpace = isRoot ? base : bones[bone.parent].extent;
		Mat4 localToParent = Mat4::translate(originInParentSpace);

		// Accumulate the parent's world transform if this isn't the root.
		Mat4 localToSkeleton = isRoot ? localToParent
			: bind[bone.parent] * localToParent;
		bind.emplace_back(localToSkeleton);
	}

	assert(bind.size() == bones.size());
	return bind;
}

std::vector< Mat4 > Skeleton::current_pose() const {
    //A4T2a: bone-to-skeleton transformations in the current pose

	//Similar to bind_pose(), but takes rotation from Bone::pose into account.
	// (and translation from Skeleton::base_offset!)

	//You'll probably want to write a loop similar to bind_pose().

	//Useful functions:
	//Bone::compute_rotation_axes() will tell you what axes (in local bone space) Bone::pose should rotate around.
	//Mat4::angle_axis(angle, axis) will produce a matrix that rotates angle (in degrees) around a given axis.

	std::vector<Mat4> pose;
	pose.reserve(bones.size());

	for (auto const& bone : bones) {
		bool isRoot = (bone.parent == -1U);

		// Build the local rotation from Bone::pose Euler angles,
		// applied around the bone's own local axes in Z -> Y -> X order.
		Vec3 axisX, axisY, axisZ;
		bone.compute_rotation_axes(&axisX, &axisY, &axisZ);
		Mat4 localRotation = Mat4::angle_axis(bone.pose.z, axisZ)
			* Mat4::angle_axis(bone.pose.y, axisY)
			* Mat4::angle_axis(bone.pose.x, axisX);

		// A root bone originates at the skeleton base plus the animated offset;
		// a child bone originates at its parent's tip (extent).
		Vec3 originInParentSpace = isRoot ? (base + base_offset)
			: bones[bone.parent].extent;
		Mat4 localToParent = Mat4::translate(originInParentSpace) * localRotation;

		// Accumulate the parent's world transform if this isn't the root.
		Mat4 localToSkeleton = isRoot ? localToParent
			: pose[bone.parent] * localToParent;
		pose.emplace_back(localToSkeleton);
	}

	assert(pose.size() == bones.size());
	return pose;

	return std::vector< Mat4 >(bones.size(), Mat4::I);

}

std::vector< Vec3 > Skeleton::gradient_in_current_pose() const {
    //A4T2b: IK gradient

    // Computes the gradient (partial derivative) of IK energy relative to each bone's Bone::pose, in the current pose.

	//The IK energy is the sum over all *enabled* handles of the squared distance from the tip of Handle::bone to Handle::target
	std::vector< Vec3 > gradient(bones.size(), Vec3{0.0f, 0.0f, 0.0f});

	//TODO: loop over handles and over bones in the chain leading to the handle, accumulating gradient contributions.
	//remember bone.compute_rotation_axes() -- should be useful here, too!

	auto pose = current_pose();
	for (auto const& handle : handles) {
		if (!handle.enabled) continue;

		// Compute the world-space position of this handle's bone tip and the displacement to target.
		Vec3 tip = (pose[handle.bone] * Vec4(bones[handle.bone].extent, 1.0f)).xyz();
		Vec3 diff = tip - handle.target;

		// Walk up the ancestor chain, accumulating each bone's contribution to the gradient.
		BoneIndex idx = handle.bone;
		while (idx != -1U) {
			// Get this bone's local rotation axes, then transform them into world space.
			// Note: x and y are fully transformed by this bone's pose (they live at the bone tip);
			//       z uses only the parent's transform because it is the twist axis at the bone root.
			Vec3 x, y, z;
			bones[idx].compute_rotation_axes(&x, &y, &z);
			Vec3 worldAxisX = pose[idx].rotate(x);
			Vec3 worldAxisY = pose[idx].rotate(y);
			Vec3 worldAxisZ = (bones[idx].parent == -1U) ? Mat4::I.rotate(z)
				: pose[bones[idx].parent].rotate(z);

			// r is the lever arm: the vector from this joint's origin to the handle tip.
			Vec3 jointPos = (pose[idx] * Vec4(0.0f, 0.0f, 0.0f, 1.0f)).xyz();
			Vec3 r = tip - jointPos;

			// Partial derivative of energy w.r.t. each Euler angle:
			// d(energy)/d(angle) = dot(diff, cross(worldAxis, r))
			gradient[idx].x += dot(diff, cross(worldAxisX, r));
			gradient[idx].y += dot(diff, cross(worldAxisY, r));
			gradient[idx].z += dot(diff, cross(worldAxisZ, r));

			idx = bones[idx].parent;
		}
	}

	assert(gradient.size() == bones.size());
	return gradient;
}

bool Skeleton::solve_ik(uint32_t steps) {
	//A4T2b - gradient descent
	//check which handles are enabled
	//run `steps` iterations
	
	//call gradient_in_current_pose() to compute d loss / d pose
	//add ...

	//if at a local minimum (e.g., gradient is near-zero), return 'true'.
	//if run through all steps, return `false`.

	const float stepSize = 0.1f;

	for (uint32_t step = 0; step < steps; step++) {
		std::vector<Vec3> gradient = gradient_in_current_pose();

		// Check for convergence: if no bone has a meaningful gradient, we're at a local minimum.
		bool atLocalMinimum = std::all_of(gradient.begin(), gradient.end(),
			[](Vec3 const& g) { return g.norm() <= 1e-4f; });
		if (atLocalMinimum) return true;

		// Descend along the negative gradient for each bone's pose.
		for (size_t i = 0; i < bones.size(); i++) {
			bones[i].pose -= stepSize * gradient[i];
		}
	}

	return false; // Step budget exhausted without converging.
}

Vec3 Skeleton::closest_point_on_line_segment(Vec3 const &a, Vec3 const &b, Vec3 const &p) {
	//A4T3: bone weight computation (closest point helper)

    // Return the closest point to 'p' on the line segment from a to b

	//Efficiency note: you can do this without any sqrt's! (no .unit() or .norm() is needed!)

	Vec3 ab = b - a;
	Vec3 ap = p - a;

	float ab_dot_ab = dot(ab, ab);

	// Degenerate segment (zero length): just return a
	if (ab_dot_ab == 0.0f) return a;

	// Project p onto the line, get parameter t
	float t = dot(ap, ab) / ab_dot_ab;

	// Clamp to [0,1] to stay on the segment
	t = std::clamp(t, 0.0f, 1.0f);

	return a + t * ab;

    return Vec3{};
}

void Skeleton::assign_bone_weights(Halfedge_Mesh *mesh_) const {
	assert(mesh_);
	auto &mesh = *mesh_;
	(void)mesh; //avoid complaints about unused mesh

	//A4T3: bone weight computation

	//visit every vertex and **set new values** in Vertex::bone_weights (don't append to old values)

	//be sure to use bone positions in the bind pose (not the current pose!)

	//you should fill in the helper closest_point_on_line_segment() before working on this function

	// Get bone positions in bind (rest) pose
	std::vector<Mat4> bind = bind_pose();

	for (auto v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		// Clear any old weights
		v->bone_weights.clear();

		// Compute raw (unnormalized) weight for each bone
		for (uint32_t j = 0; j < bones.size(); j++) {
			// Bone start and end in skeleton space
			Vec3 bone_start = (bind[j] * Vec4(0.0f, 0.0f, 0.0f, 1.0f)).xyz();
			Vec3 bone_end = (bind[j] * Vec4(bones[j].extent, 1.0f)).xyz();

			// Distance from vertex to closest point on bone
			Vec3 closest = closest_point_on_line_segment(bone_start, bone_end, v->position);
			float d = (v->position - closest).norm();

			float r = bones[j].radius;
			float w_hat = std::max(0.0f, r - d) / r;

			if (w_hat > 0.0f) {
				v->bone_weights.push_back({ j, w_hat });
			}
		}

		// Normalize weights to sum to 1
		float total = 0.0f;
		for (auto& bw : v->bone_weights) total += bw.weight;

		if (total > 0.0f) {
			for (auto& bw : v->bone_weights) bw.weight /= total;
		}
		else {
			// No bone influences this vertex — leave bone_weights empty
			v->bone_weights.clear();
		}
	}

}

Indexed_Mesh Skeleton::skin(Halfedge_Mesh const &mesh, std::vector< Mat4 > const &bind, std::vector< Mat4 > const &current) {
	assert(bind.size() == current.size());


	//A4T3: linear blend skinning

	//one approach you might take is to first compute the skinned positions (at every vertex) and normals (at every corner)
	// then generate faces in the style of Indexed_Mesh::from_halfedge_mesh

	std::vector<Mat4> M(bind.size());
	for (uint32_t j = 0; j < bind.size(); j++) {
		M[j] = current[j] * bind[j].inverse();
	}

	//---- step 1: figure out skinned positions ---

	std::unordered_map< Halfedge_Mesh::VertexCRef, Vec3 > skinned_positions;
	std::unordered_map< Halfedge_Mesh::HalfedgeCRef, Vec3 > skinned_normals;
	//reserve hash table space to (one hopes) avoid re-hashing:
	skinned_positions.reserve(mesh.vertices.size());
	skinned_normals.reserve(mesh.halfedges.size());

	//(you will probably want to precompute some bind-to-current transformation matrices here)

	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi) {
		Vec3 new_pos;
		Mat4 blended = Mat4{ Vec4{0,0,0,0}, Vec4{0,0,0,0},
							Vec4{0,0,0,0}, Vec4{0,0,0,0} };		//NOTE: vertices with empty bone_weights should remain in place.

		if (vi->bone_weights.empty()) {
			// No weights: leave vertex in place (identity transform)
			new_pos = vi->position;
			skinned_positions.emplace(vi, new_pos);

			auto h = vi->halfedge;
			do {
				skinned_normals.emplace(h, h->corner_normal);
				h = h->twin->next;
			} while (h != vi->halfedge);
			continue;
		}

		// Accumulate blended transform: sum of w_ij * M_j
		for (auto const& bw : vi->bone_weights) {
			blended += bw.weight * M[bw.bone];
		}

		new_pos = (blended * Vec4(vi->position, 1.0f)).xyz();
		skinned_positions.emplace(vi, new_pos);

		// Normal transform: inverse transpose of blended
		Mat4 normal_mat = blended.inverse().T();

		// Circulate halfedges (corners) at this vertex
		auto h = vi->halfedge;
		do {
			Vec3 n = (normal_mat * Vec4(h->corner_normal, 0.0f)).xyz();
			// Safely normalize
			if (n.norm() > 1e-6f) n = n.unit();
			skinned_normals.emplace(h, n);
			h = h->twin->next;
		} while (h != vi->halfedge);
	}

	//---- step 2: transform into an indexed mesh ---

	//Hint: you should be able to use the code from Indexed_Mesh::from_halfedge_mesh (SplitEdges version) pretty much verbatim, you'll just need to fill in the positions and normals.

	Indexed_Mesh result;
	std::vector<Indexed_Mesh::Vert> verts;
	std::vector<Indexed_Mesh::Index> idxs;

	for (auto fi = mesh.faces.begin(); fi != mesh.faces.end(); ++fi) {
		if (fi->boundary) continue;

		// Collect corners of this face
		std::vector<Halfedge_Mesh::HalfedgeCRef> face_he;
		auto h = fi->halfedge;
		do {
			face_he.push_back(h);
			h = h->next;
		} while (h != fi->halfedge);

		// Fan triangulation — emit one vertex per corner
		uint32_t base = uint32_t(verts.size());
		for (auto he : face_he) {
			Indexed_Mesh::Vert vert;
			vert.pos = skinned_positions.at(he->vertex);
			vert.norm = skinned_normals.at(he);
			vert.uv = he->corner_uv;
			vert.id = he->vertex->id;
			verts.push_back(vert);
		}

		// Triangulate the polygon as a fan from base vertex
		for (uint32_t i = 1; i + 1 < face_he.size(); i++) {
			idxs.push_back(base);
			idxs.push_back(base + i);
			idxs.push_back(base + i + 1);
		}
	}

	return Indexed_Mesh(std::move(verts), std::move(idxs));
}

void Skeleton::for_bones(const std::function<void(Bone&)>& f) {
	for (auto& bone : bones) {
		f(bone);
	}
}


void Skeleton::erase_bone(BoneIndex bone) {
	assert(bone < bones.size());
	//update indices in bones:
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (bones[b].parent == -1U) continue;
		if (bones[b].parent == bone) {
			assert(b > bone); //topological sort!
			//keep bone tips in the same place when deleting parent bone:
			bones[b].extent += bones[bone].extent;
			bones[b].parent = bones[bone].parent;
		} else if (bones[b].parent > bone) {
			assert(b > bones[b].parent); //topological sort!
			bones[b].parent -= 1;
		}
	}
	// erase the bone
	bones.erase(bones.begin() + bone);
	//update indices in handles (and erase any handles on this bone):
	for (uint32_t h = 0; h < handles.size(); /* later */) {
		if (handles[h].bone == bone) {
			erase_handle(h);
		} else if (handles[h].bone > bone) {
			handles[h].bone -= 1;
			++h;
		} else {
			++h;
		}
	}
}

void Skeleton::erase_handle(HandleIndex handle) {
	assert(handle < handles.size());

	//nothing internally refers to handles by index so can just delete:
	handles.erase(handles.begin() + handle);
}


Skeleton::BoneIndex Skeleton::add_bone(BoneIndex parent, Vec3 extent) {
	assert(parent == -1U || parent < bones.size());
	Bone bone;
	bone.extent = extent;
	bone.parent = parent;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &b : bones) {
		used.emplace(b.channel_id);
	}
	while (used.count(next_bone_channel_id)) ++next_bone_channel_id;
	bone.channel_id = next_bone_channel_id++;

	//all other parameters left as default.

	BoneIndex index = BoneIndex(bones.size());
	bones.emplace_back(bone);

	return index;
}

Skeleton::HandleIndex Skeleton::add_handle(BoneIndex bone, Vec3 target) {
	assert(bone < bones.size());
	Handle handle;
	handle.bone = bone;
	handle.target = target;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &h : handles) {
		used.emplace(h.channel_id);
	}
	while (used.count(next_handle_channel_id)) ++next_handle_channel_id;
	handle.channel_id = next_handle_channel_id++;

	HandleIndex index = HandleIndex(handles.size());
	handles.emplace_back(handle);

	return index;
}


Skeleton Skeleton::copy() {
	//turns out that there aren't any fancy pointer data structures to fix up here.
	return *this;
}

void Skeleton::make_valid() {
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (!(bones[b].parent == -1U || bones[b].parent < b)) {
			warn("bones[%u].parent is %u, which is not < %u; setting to -1.", b, bones[b].parent, b);
			bones[b].parent = -1U;
		}
	}
	if (bones.empty() && !handles.empty()) {
		warn("Have %u handles but no bones. Deleting handles.", uint32_t(handles.size()));
		handles.clear();
	}
	for (uint32_t h = 0; h < handles.size(); ++h) {
		if (handles[h].bone >= HandleIndex(bones.size())) {
			warn("handles[%u].bone is %u, which is not < bones.size(); setting to 0.", h, handles[h].bone);
			handles[h].bone = 0;
		}
	}
}

//-------------------------------------------------

Indexed_Mesh Skinned_Mesh::bind_mesh() const {
	return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges);
}

Indexed_Mesh Skinned_Mesh::posed_mesh() const {
	return Skeleton::skin(mesh, skeleton.bind_pose(), skeleton.current_pose());
}

Skinned_Mesh Skinned_Mesh::copy() {
	return Skinned_Mesh{mesh.copy(), skeleton.copy()};
}
