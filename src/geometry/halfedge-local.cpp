
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */


std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// Gather halfedges, vertices, and faces around the edge
	HalfedgeRef h0 = e->halfedge;
	HalfedgeRef h1 = h0->twin;

	HalfedgeRef h2 = h0->next;
	HalfedgeRef h3 = h2->next;
	HalfedgeRef h4 = h1->next;
	HalfedgeRef h5 = h4->next;

	VertexRef va = h2->vertex;
	VertexRef vb = h4->vertex;
	VertexRef vc = h3->vertex;
	VertexRef vd = h5->vertex;

	FaceRef f0 = h0->face;
	FaceRef f1 = h1->face;

	bool bo0 = f0->boundary;
	bool bo1 = f1->boundary;

	HalfedgeRef h0_prev = h0; while (h0_prev->next != h0) h0_prev = h0_prev->next;
	HalfedgeRef h1_prev = h1; while (h1_prev->next != h1) h1_prev = h1_prev->next;

	// Allocate the new vertex and the two halves of the split edge
	VertexRef m = emplace_vertex();
	m->position = e->center();

	EdgeRef ea = emplace_edge(false);
	EdgeRef eb = emplace_edge(false);

	HalfedgeRef a0 = emplace_halfedge();
	HalfedgeRef a1 = emplace_halfedge();
	HalfedgeRef b0 = emplace_halfedge();
	HalfedgeRef b1 = emplace_halfedge();

	ea->halfedge = a0;
	eb->halfedge = b1;

	// Conditionally allocate geometry for each non-boundary face
	EdgeRef ec, ed;
	FaceRef f2, f3;
	HalfedgeRef c0, c1, d0, d1;

	if (!f0->boundary) {
		ec = emplace_edge();
		c0 = emplace_halfedge();
		c1 = emplace_halfedge();
		f2 = emplace_face(false);
	}
	if (!f1->boundary) {
		ed = emplace_edge();
		d0 = emplace_halfedge();
		d1 = emplace_halfedge();
		f3 = emplace_face(false);
	}

	// Wire up the h0 side
	if (!f0->boundary) {
		a0->twin = a1; a0->next = h2; a0->vertex = m;  a0->edge = ea; a0->face = f2;
		b0->twin = b1; b0->next = c1; b0->vertex = vb; b0->edge = eb; b0->face = f0;
		c0->twin = c1; c0->next = a0; c0->vertex = vc; c0->edge = ec; c0->face = f2;
		c1->twin = c0; c1->next = h3; c1->vertex = m;  c1->edge = ec; c1->face = f0;

		ec->halfedge = c0;
		m->halfedge = a0;
		vc->halfedge = c0;
		vb->halfedge = b0;

		f2->halfedge = c0;
		f0->halfedge = h0_prev;

		h2->next = c0;
		h0_prev->next = b0;

		HalfedgeRef cur = a0; do { cur->face = f2; cur = cur->next; } while (cur != a0);
		cur = b0;             do { cur->face = f0; cur = cur->next; } while (cur != b0);
	}

	if (f0->boundary) {
		a0->twin = a1; a0->next = h2;  a0->vertex = m;  a0->edge = ea; a0->face = f0;
		b0->twin = b1; b0->next = a0;  b0->vertex = vb; b0->edge = eb; b0->face = f0;

		h0_prev->next = b0;

		ea->halfedge = a0;
		eb->halfedge = b0;
		f0->halfedge = a0;
		m->halfedge = a0;
	}

	// Wire up the h1 side
	if (!f1->boundary) {
		a1->twin = a0; a1->next = d1; a1->vertex = va; a1->edge = ea; a1->face = f3;
		b1->twin = b0; b1->next = h4; b1->vertex = m;  b1->edge = eb; b1->face = f1;
		d0->twin = d1; d0->next = b1; d0->vertex = vd; d0->edge = ed; d0->face = f1;
		d1->twin = d0; d1->next = h5; d1->vertex = m;  d1->edge = ed; d1->face = f3;

		ed->halfedge = d0;
		vd->halfedge = d0;
		va->halfedge = a1;

		f3->halfedge = d1;
		f1->halfedge = h4;

		h4->next = d0;
		h1_prev->next = a1;
		d0->next = b1;

		HalfedgeRef cur = a1; do { cur->face = f3; cur = cur->next; } while (cur != a1);
		cur = b1;             do { cur->face = f1; cur = cur->next; } while (cur != b1);
	}

	if (bo1) {
		a1->twin = a0; a1->next = b1;  a1->vertex = va; a1->edge = ea; a1->face = f1;
		b1->twin = b0; b1->next = h4;  b1->vertex = m;  b1->edge = eb; b1->face = f1;

		h1_prev->next = a1;

		ea->halfedge = a1;
		eb->halfedge = b1;
		f1->halfedge = b1;
		m->halfedge = b1;
	}

	// Remove the original edge and its halfedges
	erase_halfedge(h0);
	erase_halfedge(h1);
	erase_edge(e);

	return m;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face-
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	std::vector<VertexRef> vf; // Original vertices of face
	std::vector<VertexRef> vfm; // New vertices for extruded face
	std::vector<EdgeRef> edges; // New edges created by extrusion
	std::vector<EdgeRef> edges_n; // New edges parallel to original face edges
	std::vector<HalfedgeRef> hs; // New halfedges created by extrusion
	std::vector<HalfedgeRef> h_os;
	std::vector<HalfedgeRef> c_os; // New halfedges parallel to original face halfedges, but opposite direction
	std::vector<HalfedgeRef> c_cs; // New halfedges on the inside of the new face in the center

	auto get_cs = [](FaceRef f) {
		std::vector<HalfedgeRef> cs;
		HalfedgeRef start = f->halfedge;
		HalfedgeRef cur = start;
		do {
			cs.push_back(cur);
			cur = cur->next;
		} while (cur != start);
		return cs;
		};

	std::vector<HalfedgeRef> cs = get_cs(f);
	for (auto c : cs)
	{
		vf.push_back(c->vertex);
		VertexRef v = c->vertex;
				
		// Create edge from old vertex to new vertex
		EdgeRef e = emplace_edge();
		EdgeRef en = emplace_edge();
		HalfedgeRef h = emplace_halfedge();
		HalfedgeRef h_o = emplace_halfedge();
		HalfedgeRef c_o = emplace_halfedge();
		HalfedgeRef c_c = emplace_halfedge();
		//std::cout << "h: " << h->id << ", h_o: " << h_o->id << ", c_o: " << c_o->id << ", c_c: " << c_c->id << std::endl;
		e->halfedge = h;
		h->edge = e;
		en->halfedge = c_o;
		c_o->edge = en;
		c_c->edge = en;
		h->next = c;
		c->next = h_o;
		h_o->next = c_o;
		c_o->next = h;
		c_o->twin = c_c;
		c_c->twin = c_o;

		edges.push_back(e);
		edges_n.push_back(en);
		hs.push_back(h);
		h_os.push_back(h_o);
		c_os.push_back(c_o);
		c_cs.push_back(c_c);

		// Create new vertex for extruded face
		VertexRef vm = emplace_vertex();
		//std::cout << "Created vertex " << vm->id << std::endl;
		vfm.push_back(vm);
		vm->position = v->position;
		vm->bone_weights = v->bone_weights;
		vm->halfedge = h;
		h->vertex = vm;
		c_c->vertex = vm;

		// Create new face for face shared between c and h
		FaceRef fn = emplace_face();
		h->face = fn;
		c->face = fn;
		h_o->face = fn;
		c_o->face = fn;
		fn->halfedge = h;

		//std::cout << "-------------" << std::endl;
	}

	//std::cout << "Created " << vf.size() << " new vertices, " << std::endl;
	f->halfedge = c_cs[0];
	for (int i = 0; i < (int)hs.size(); i++)
	{
		int j = (i == 0) ? (int)hs.size() - 1 : i - 1; // previous index, with wraparound
		int k = (i + 1) % (int)hs.size(); // next index, with wraparound

		// Set up twins and vertex/edge pointers from halfedges from the last iteration
		HalfedgeRef t = h_os[j];
		hs[i]->twin = t;
		t->twin = hs[i];
		t->edge = edges[i];
		t->vertex = vf[i];
		c_os[j]->vertex = vfm[i];
		//std::cout << "Assigning vertex " << vfm[i]->id << " to halfedge " << c_os[j]->id << std::endl;

		// Set up new face
		c_cs[i]->face = f;
		c_cs[i]->next = c_cs[k];
		//std::cout << "Inner halfedge " << c_cs[i]->id << std::endl;
	}
	return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */

std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	if (e->on_boundary())
		return std::nullopt;

	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;

	// Invalid mesh check
	if (h->vertex->degree() <= 2 || t->vertex->degree() <= 2)
		return std::nullopt;

	HalfedgeRef h_next = h->next;
	HalfedgeRef h_next_next = h_next->next;
	HalfedgeRef t_next = t->next;
	HalfedgeRef t_next_next = t_next->next;

	HalfedgeRef prev = h;
	do {
		prev = prev->next;
	} while (prev->next->id != h->id);
	HalfedgeRef h_prev = prev;
	prev = t;
	do {
		prev = prev->next;
	} while (prev->next->id != t->id);
	HalfedgeRef t_prev = prev;

	FaceRef fh = h->face;
	FaceRef ft = t->face;

	VertexRef vh = t->vertex;
	VertexRef vt = h->vertex;
	VertexRef next_h_from_vertex = t_next->next->vertex;
	VertexRef next_t_from_vertex = h_next->next->vertex;
	
	h->next = h_next_next;
	h->vertex = next_h_from_vertex;
	t->next = t_next_next;
	t->vertex = next_t_from_vertex;

	h_prev->next = t_next;
	h_next->next = t;
	h_next->face = ft;
	t_prev->next = h_next;
	t_next->next = h;
	t_next->face = fh;

	if (vh->halfedge == t)
		vh->halfedge = h_next;
	if (vt->halfedge == h)
		vt->halfedge = t_next;

	if (fh->halfedge == h_next)
		fh->halfedge = h;
	if (ft->halfedge == t_next)
		ft->halfedge = t;

	return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	HalfedgeRef h = e->halfedge;
	//std::cout << "Collapsing h " << h->id << std::endl;
	HalfedgeRef t = h->twin;
	//std::cout << "Collapsing t " << t->id << std::endl;
	FaceRef     fh = h->face;
	FaceRef     ft = t->face;
	VertexRef   vh = h->vertex;
	VertexRef   vt = t->vertex;
	

	// --- Validity checks ---
	auto get_neighbors = [](VertexRef v) {
		std::vector<VertexRef> neighbors;
		HalfedgeRef start = v->halfedge;
		HalfedgeRef cur = start;
		do {
			neighbors.push_back(cur->twin->vertex);
			cur = cur->twin->next;
		} while (cur != start);
		return neighbors;
		};

	auto get_leaving = [](VertexRef v) {
		std::vector<HalfedgeRef> leaving;
		HalfedgeRef start = v->halfedge;
		HalfedgeRef cur = start;
		do {
			leaving.push_back(cur);
			cur = cur->twin->next;
		} while (cur != start);
		return leaving;
		};

	std::vector<VertexRef> nh = get_neighbors(vh);
	std::vector<VertexRef> nt = get_neighbors(vt);
	std::vector<HalfedgeRef> leaving_vh = get_leaving(vh);
	std::vector<HalfedgeRef> leaving_vt = get_leaving(vt);
	int shared = 0;
	for (auto a : nh) {
		if (a == vt) continue;
		for (auto b : nt) {
			if (b == vh) continue;
			if (a == b) shared++;
		}
	}
	if (!fh->boundary && !ft->boundary && shared > 2) return std::nullopt;
	if (fh->boundary && !ft->boundary && shared > 1) return std::nullopt;
	if (!fh->boundary && ft->boundary && shared > 1) return std::nullopt;
	if (fh->boundary && ft->boundary && shared > 0) return std::nullopt;
	if (fh->boundary && ft->boundary) return std::nullopt;

	// --- Collect halfedges ---
	HalfedgeRef h_next = h->next;
	HalfedgeRef t_next = t->next;
	HalfedgeRef h_prev = h->next;
	while (h_prev->next != h)
	{
		h_prev = h_prev->next;
	}
	HalfedgeRef t_prev = t->next;
	while (t_prev->next != t)
	{
		t_prev = t_prev->next;
	}
	//HalfedgeRef h_next_twin = h_next->twin;
	//HalfedgeRef h_prev_twin = h_prev->twin;
	//HalfedgeRef t_next_twin = t_next->twin;
	//HalfedgeRef t_prev_twin = t_prev->twin;


	// Assign halfedges leaving vt to vh
	VertexRef vn = emplace_vertex();
	//std::cout << "Collapsing edge " << e->id << " to vertex " << vn->id << std::endl;
	vn->position = e->center();
	interpolate_data({ vh, vt }, vn);
	for (auto l : leaving_vh)
	{
		if (l != h)
		{
			l->vertex = vn;
			vn->halfedge = l;
			//std::cout << "Reassigning halfedge " << l->id << " from vh " << vh->id << " to vn " << vn->id << std::endl;
		}
	}
	for (auto l : leaving_vt)
	{
		if (l != t)
		{
			l->vertex = vn;
			vn->halfedge = l;
			//std::cout << "Reassigning halfedge " << l->id << " from vt " << vt->id << " to vn " << vn->id << std::endl;
		}
	}

	if (shared == 0)
	{
		h_prev->next = h_next;
		t_prev->next = t_next;
	}
	else
	{
		// Find any shared neighbors of vh and vt, and if they exist, erase the face between them
		// as well as the edge and halfedges connecting them to vh and vt
		for (auto a : nh) {
			if (a == vt) continue;
			for (auto b : nt) {
				if (b == vh) continue;
				if (a == b)
				{
					//std::cout << "Shared neighbor " << a->id << std::endl;
					bool erased = false;
					for (auto l : leaving_vh)
					{
						// Find halfedge leaving vh that connects to neighbor, check if the halfedge is not boundary
						if (l->twin->vertex == a && !(l->face->boundary || l->twin->face->boundary))
						{
							// Two cases from here: h side and t side. We can check which side by checking if l->face is ft or l->twin->face is fh
							if (l->face == ft)
							{
								// We found the triangle face that will be erased, and it is on the t side
								HalfedgeRef l_twin_prev = l->twin->next;
								while (l_twin_prev->next != l->twin)
								{
									l_twin_prev = l_twin_prev->next;
								}
								l_twin_prev->next = l->next;
								l->next->next = l->twin->next;

								FaceRef f = l->twin->face;
								l->next->face = f;
								f->halfedge = l->next;
								//std::cout << "Reassigning halfedge " << l->next->id << " to face " << f->id << std::endl;

								vt->halfedge = l->next->twin;
								l->twin->vertex->halfedge = l->next;

								//std::cout << "Erasing face " << l->face->id << ", edge " << l->edge->id << ", halfedges " << l->id << " and " << l->twin->id << std::endl;
								erase_face(l->face);
								erase_edge(l->edge);
								erase_halfedge(l->twin);
								erase_halfedge(l);
								erased = true;
								break;
							}
							else if (l->twin->face == fh)
							{
								// We found the triangle face that will be erased, and it is on the h side
								HalfedgeRef l_prev = l->next;
								while (l_prev->next != l)
								{
									l_prev = l_prev->next;
								}
								l_prev->next = h->next;
								h->next->next = l->next;

								FaceRef f = l->face;
								h->next->face = f;
								f->halfedge = h->next;
								//std::cout << "Reassigning halfedge " << h->next->id << " to face " << f->id << std::endl;

								vt->halfedge = h->next;
								l->twin->vertex->halfedge = l->next;

								//std::cout << "Erasing face " << l->twin->face->id << ", edge " << l->edge->id << ", halfedges " << l->id << " and " << l->twin->id << std::endl;
								erase_face(l->twin->face);
								erase_edge(l->edge);
								erase_halfedge(l->twin);
								erase_halfedge(l);
								erased = true;
								break;
							}
						}
					}
					if (erased) continue;
					for (auto l : leaving_vt)
					{
						// Getting here would require the halfedge leaving vh to be a boundary edge
						// Find halfedge leaving vt that connects to neighbor, check if the halfedge is not boundary
						if (l->twin->vertex == b && !(l->face->boundary || l->twin->face->boundary))
						{
							// Two cases from here: h side and t side. We can check which side by checking if l->face is ft or l->twin->face is fh
							if (l->face == fh)
							{
								// We found the triangle face that will be erased, and it is on the h side
								HalfedgeRef l_twin_prev = l->twin->next;
								while (l_twin_prev->next != l->twin)
								{
									l_twin_prev = l_twin_prev->next;
								}
								l_twin_prev->next = l->next;
								l->next->next = l->twin->next;

								FaceRef f = l->twin->face;
								l->next->face = f;
								f->halfedge = l->next;
								//std::cout << "Reassigning halfedge " << l->next->id << " to face " << f->id << std::endl;

								vt->halfedge = l->twin->next;
								l->twin->vertex->halfedge = l->next;

								//std::cout << "Erasing face " << l->face->id << ", edge " << l->edge->id << ", halfedges " << l->id << " and " << l->twin->id << std::endl;
								erase_face(l->face);
								erase_edge(l->edge);
								erase_halfedge(l->twin);
								erase_halfedge(l);
								break;
							}
							if (l->twin->face == ft)
							{
								// We found the triangle face that will be erased, and it is on the t side
								HalfedgeRef l_prev = l->next;
								while (l_prev->next != l)
								{
									l_prev = l_prev->next;
								}
								l_prev->next = t->next;
								t->next->next = l->next;

								FaceRef f = l->face;
								t->next->face = f;
								f->halfedge = t->next;
								//std::cout << "Reassigning halfedge " << t->next->id << " to face " << f->id << std::endl;

								vt->halfedge = t->next;
								l->twin->vertex->halfedge = l->next;

								//std::cout << "Erasing face " << l->twin->face->id << ", edge " << l->edge->id << ", halfedges " << l->id << " and " << l->twin->id << std::endl;
								erase_face(l->twin->face);
								erase_edge(l->edge);
								erase_halfedge(l->twin);
								erase_halfedge(l);
								break;
							}
						}
					}
				}
			}
		}
	}
	if (fh->boundary)
	{
		fh->halfedge = h_next;
		h_prev->next = h_next;
	}
	if (ft->boundary)
	{
		ft->halfedge = t_next;
		t_prev->next = t_next;
	}
	erase_vertex(vh);
	erase_vertex(vt);
	erase_edge(e);
	erase_halfedge(h);
	erase_halfedge(t);

	return vt;

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)
	
    //return std::nullopt;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

	Vec3 center = face->center();
	auto get_halfedges = [](FaceRef f) {
		std::vector<HalfedgeRef> hs;
		HalfedgeRef start = f->halfedge;
		HalfedgeRef cur = start;
		do {
			hs.push_back(cur);
			cur = cur->next;
		} while (cur != start);
		return hs;
		};
	
	std::vector<HalfedgeRef> hs = get_halfedges(face);
	std::vector<VertexRef> vs;
	for (auto h : hs)
	{
		vs.push_back(h->vertex);
	}

	for (auto v : vs)
	{
		Vec3 start = v->position;
		Vec3 target = start - shrink * (start - center) + move;
		v->position = target;
	}
}


