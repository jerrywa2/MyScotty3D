
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
	// A2L2 (REQUIRED): split_edge
	
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	FaceRef     fh = h->face;
	FaceRef     ft = t->face;
	VertexRef   vh = h->vertex;
	VertexRef   vt = t->vertex;
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

	if (!fh->boundary || !ft->boundary)
	{
		// Split fh by adding a vertex at the midpoint of e and connecting it to the ccw-most vertex from h in fh
		VertexRef vm = emplace_vertex();
		//std::cout << "Splitting edge " << e->id << " with vertex " << vm->id << std::endl;
		vm->position = e->center();
		vm->halfedge = !fh->boundary ? h : t;

		EdgeRef em = emplace_edge();
		em->halfedge = t;
		h->vertex = vm;
		t->vertex = vm;
		t->edge = em;

		// New halfedge twin for h
		HalfedgeRef twinh = emplace_halfedge();
		//std::cout << "New halfedge twin for h " << twinh->id << std::endl;
		twinh->twin = h;
		h->twin = twinh;
		twinh->next = t;
		twinh->vertex = vt;
		vt->halfedge = twinh;
		twinh->edge = h->edge;
		twinh->face = ft;
		if (ft->boundary)
		{
			h->next->twin->next = twinh;
		}

		// New halfedge twin for t
		HalfedgeRef twint = emplace_halfedge();
		//std::cout << "New halfedge twin for t " << twint->id << std::endl;
		twint->twin = t;
		t->twin = twint;
		twint->next = h;
		twint->vertex = vh;
		vh->halfedge = twint;
		twint->edge = t->edge;
		twint->face = fh;
		if (fh->boundary)
		{
			t->next->twin->next = twint;
		}

		if (!fh->boundary)
		{
			h_prev->next = twint;

			// Create new face
			FaceRef fhm = emplace_face();
			//std::cout << "New face for h " << fhm->id << std::endl;
			fhm->halfedge = h;
			h->face = fhm;
			h->next->face = fhm;
			
			// Close new face with new edge and halfedges
			EdgeRef eh = emplace_edge();
			HalfedgeRef eh_h = emplace_halfedge();
			HalfedgeRef eh_t = emplace_halfedge();
			//std::cout << "New edge and halfedges for h " << eh->id << ", " << eh_h->id << ", " << eh_t->id << std::endl;
			eh->halfedge = eh_h;
			eh_h->twin = eh_t;
			eh_t->twin = eh_h;
			eh_h->next = h;
			eh_t->next = h->next->next;
			h->next->next = eh_h;
			eh_h->vertex = h->next->twin->vertex;
			eh_t->vertex = vm;
			eh_h->edge = eh;
			eh_t->edge = eh;
			eh_h->face = fhm;
			eh_t->face = fh;

			if (fh->halfedge == h || fh->halfedge == h->next)
			{
				fh->halfedge = eh_t;
			}

			// New next for twin t if new e created
			twint->next = eh_t;
		}
		
		if (!ft->boundary)
		{
			t_prev->next = twinh;

			// Create new face
			FaceRef ftm = emplace_face();
			//std::cout << "New face for t " << ftm->id << std::endl;
			ftm->halfedge = t;
			t->face = ftm;
			t->next->face = ftm;

			// Close new face with new edge and halfedges
			EdgeRef et = emplace_edge();
			HalfedgeRef et_h = emplace_halfedge();
			HalfedgeRef et_t = emplace_halfedge();
			et->halfedge = et_h;
			et_h->twin = et_t;
			et_t->twin = et_h;
			et_h->next = t;
			et_t->next = t->next->next;
			t->next->next = et_h;
			et_h->vertex = t->next->twin->vertex;
			et_t->vertex = vm;
			et_h->edge = et;
			et_t->edge = et;
			et_h->face = ftm;
			et_t->face = ft;

			if (ft->halfedge == t || ft->halfedge == t->next)
			{
				ft->halfedge = et_t;
			}

			// New next for twin t if new e created
			twinh->next = et_t;
		}

		return vm;

	}
    return std::nullopt;
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
	//A2L1: Flip Edge
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	FaceRef     fh = h->face;
	FaceRef     ft = t->face;
	VertexRef   v1;
	VertexRef   v2;

	if (fh->boundary || ft->boundary) return std::nullopt;

	// Change vertices
	if (h->vertex->halfedge == h) h->vertex->halfedge = t->next;
	if (t->vertex->halfedge == t) t->vertex->halfedge = h->next;
	v1 = t->next->twin->vertex;
	v2 = h->next->twin->vertex;
	if (v1 == v2) return std::nullopt; //Reject if degenerate edge would be created
	h->vertex = v1;
	t->vertex = v2;

	// Change faces
	h->next->face = ft;
	t->next->face = fh;
	if (fh->halfedge == h->next) fh->halfedge = h;
	if (ft->halfedge == t->next) ft->halfedge = t;

	// CHange next pointers
	HalfedgeRef h_next = h->next;
	HalfedgeRef t_next = t->next;
	HalfedgeRef h_prev = t->next->twin->next->twin;
	HalfedgeRef t_prev = h->next->twin->next->twin;

	h->next = h->next->next;
	t->next = t->next->next;
	h_next->next = t;
	t_next->next = h;

	h_prev->next = t_next;
	t_prev->next = h_next;
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
	HalfedgeRef h_next_twin = h_next->twin;
	HalfedgeRef h_prev_twin = h_prev->twin;
	HalfedgeRef t_next_twin = t_next->twin;
	HalfedgeRef t_prev_twin = t_prev->twin;


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


