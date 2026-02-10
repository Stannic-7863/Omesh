package mesh

import "core:slice"
import "core:image/netpbm"
import "core:math/linalg"
import "base:intrinsics"
import "base:runtime"
import "core:log"

Vec3f32 :: [3]f32

Half_Edge_Index :: distinct i32
Vertex_Index :: distinct i32
Face_Index :: distinct i32

Half_Edge :: struct {
	opposite: Half_Edge_Index,
	prev:     Half_Edge_Index,
	next:     Half_Edge_Index,
	vertex:   Vertex_Index,
	face:     Face_Index,
}

Face :: struct {
	edge: Half_Edge_Index,
}

Vertex :: struct {
	edge:     Half_Edge_Index,
	position: Vec3f32,
}

Lookup_Pair :: struct {
	source, target: Vertex_Index,
}

Mesh :: struct {
	faces:      [dynamic]Face,
	verts:      [dynamic]Vertex,
	edges:      [dynamic]Half_Edge,

	active_faces: [dynamic]Face_Index,
	active_verts: [dynamic]Vertex_Index,
	active_edges: [dynamic]Half_Edge_Index,

	free_faces: [dynamic]Face_Index,
	free_verts: [dynamic]Vertex_Index,
	free_edges: [dynamic]Half_Edge_Index,
	lookup:     map[Lookup_Pair]Half_Edge_Index, // Source -> Target
}

Face_Edge_Iterator :: struct {
	mesh:	 ^Mesh,
	step:	 i32,
	start:	 Half_Edge_Index,
	current: Half_Edge_Index,
}

Vertex_Edge_Iterator :: struct {
	mesh:    ^Mesh,
	step:	 i32,
	start:	 Half_Edge_Index,
	current: Half_Edge_Index,
}

Triangle_Emitter_Iterator :: struct {
	mesh:			^Mesh,
	walk_step:		i32,
	vertex_step:	i32,
	vertex_base:	i32,
	face_step:  	i32,
	face:			Face_Index,
	edge:			Half_Edge_Index,
	start:			Half_Edge_Index,
	normal:			Vec3f32,
}

Platonic_Solid :: enum {
	Tetrahedron, 	// dual of itself
	Cube, 			// dual of Octahedro
	Octahedron, 	// dual of Cube
	Dodecaheron, 	// dual of Icosahedron
	Icosahedron, 	// dual of Dodecahedron
}

Archimedean_Solid :: enum {
	Truncated_Tetrahedron,			// tT, 	truncate 		(tetrahedron)
	Cuboctahedron,					// aC, 	ambo 			(cube)
	Truncated_Cube,					// tC, 	truncate		(cube)
	Truncated_Octahedron,			// tO, 	truncate 		(octahedron)
	Rhombicuboctahedron,			// aaC,	ambo ambo		(cube)
	Truncated_Cuboctahedron,		// taC, truncate ambo 	(cube)
	Snub_Cube,						// sC, 	snub			(cube)
	Icosidodecahedron,				// aD,	ambo 			(dodecahedron)
	Truncated_dodecahedron,			// tD, 	truncate 		(tetrahedron)
	Truncated_icosahedron,			// tI, 	truncate 		(icosahedron)
	Rhombicosidodecahedron,			// aaD, ambo ambo 		(dodecahedron)
	Truncated_Icosidodecahedron,	// taD, truncate ambo 	(dodecahedron)
	Snub_dodecahedron,				// sD, 	snub 			(dodecahedron)
}

Catalan_Solid :: enum { // Conway operations generating the duals of Archimedean solids
    Triakis_Tetrahedron,          	// kT,  kis      		(tetrahedron)
    Rhombic_Dodecahedron,         	// jC,  join     		(cube)
    Triakis_Octahedron,           	// kO,  kis      		(octahedron)
    Tetrakis_Hexahedron,          	// kC,  kis      		(cube)
    Deltoidal_Icositetrahedron,   	// oC,  ortho    		(cube)
    Disdyakis_Dodecahedron,       	// mC,  meta     		(cube)
    Pentagonal_Icositetrahedron,  	// gC,  gyro     		(cube)
    Rhombic_Triacontahedron,      	// jD,  join     		(dodecahedron)
    Triakis_Icosahedron,          	// kI,  kis      		(icosahedron)
    Pentakis_Dodecahedron,        	// kD,  kis      		(dodecahedron)
    Deltoidal_Hexecontahedron,    	// oD,  ortho    		(dodecahedron)
    Disdyakis_Triacontahedron,    	// mD,  meta     		(dodecahedron)
    Pentagonal_Hexecontahedron,   	// gD,  gyro     		(dodecahedron)
}

// Ratios are calculated to ensure the resulting dual faces are congruent.
// Height is the distance to offset the centroid along the face normal.
CATALAN_TRI_TETRAHEDRON_KIS_HEIGHT  	:: 1.0 / 3.0
CATALAN_TETRA_HEXAHEDRON_KIS_HEIGHT 	:: 0.5
CATALAN_TRI_OCTAHEDRON_KIS_HEIGHT   	:: 0.414213
CATALAN_PENTA_DODECAHEDRON_KIS_HEIGHT 	:: 0.11135
CATALAN_TRI_ICOSAHEDRON_KIS_HEIGHT    	:: 0.15836

PHI :: 1.618033988749894
INV_PHI :: 0.618033988749894

Convay_Operation :: enum {
	Ambo,
	Bevel,					// Also called omnitruncation
	Dual,
	Expand,
	Gyro,
	Join,
	Kis,
	Meta,
	Ortho,
	Snub,
	Truncate,
	Needle,
	Zip, 					// Also called bitruncation
	Classical_Alternation, 	// Only works for meshes where vertices are 2 color-able. Basically Edge count of all faces must be even
	Classical_Snub, 		// same limitation as classical alternation
	Classical_Gyro, 		// same limitation as classical alternation
}

mesh_create :: proc(allocator: runtime.Allocator) -> Mesh {
	m := Mesh{}
	m.faces = make([dynamic]Face, allocator)
	m.verts = make([dynamic]Vertex, allocator)
	m.edges = make([dynamic]Half_Edge, allocator)

	m.free_faces = make([dynamic]Face_Index, allocator)
	m.free_verts = make([dynamic]Vertex_Index, allocator)
	m.free_edges = make([dynamic]Half_Edge_Index, allocator)

	m.active_faces = make([dynamic]Face_Index, allocator)
	m.active_verts = make([dynamic]Vertex_Index, allocator)
	m.active_edges = make([dynamic]Half_Edge_Index, allocator)

	m.lookup = make(map[Lookup_Pair]Half_Edge_Index, allocator)
	return m
}

mesh_destroy :: proc(mesh: Mesh) {
	delete(mesh.lookup)

	delete(mesh.faces)
	delete(mesh.verts)
	delete(mesh.edges)

	delete(mesh.free_faces)
	delete(mesh.free_verts)
	delete(mesh.free_edges)

	delete(mesh.active_faces)
	delete(mesh.active_verts)
	delete(mesh.active_edges)
}

meshes_destroy :: proc(meshes: ..Mesh) {
	for mesh in meshes {
		mesh_destroy(mesh)
	}
}

mesh_get_face :: proc (mesh: Mesh, index: Face_Index) -> Face {
	return mesh.faces[index]
}

mesh_get_vertex :: proc (mesh: Mesh, index: Vertex_Index) -> Vertex {
	return mesh.verts[index]
}

mesh_get_edge :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return mesh.edges[index]
}

mesh_get_face_ptr :: proc (mesh: Mesh, index: Face_Index) -> ^Face {
	return &mesh.faces[index]
}

mesh_get_vertex_ptr :: proc (mesh: Mesh, index: Vertex_Index) -> ^Vertex {
	return &mesh.verts[index]
}

mesh_get_edge_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return &mesh.edges[index]
}

mesh_get_edge_next :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return mesh.edges[mesh.edges[index].next]
}

mesh_get_edge_prev :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return mesh.edges[mesh.edges[index].prev]
}

mesh_get_edge_opposite :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return mesh.edges[mesh.edges[index].opposite]
}

mesh_get_edge_next_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return &mesh.edges[mesh.edges[index].next]
}

mesh_get_edge_prev_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return &mesh.edges[mesh.edges[index].prev]
}

mesh_get_edge_opposite_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return &mesh.edges[mesh.edges[index].opposite]
}

mesh_get_edge_source :: proc (mesh: Mesh, index: Half_Edge_Index) -> Vertex {
	return mesh.verts[mesh.edges[mesh.edges[index].prev].vertex]
}

mesh_get_edge_target :: proc (mesh: Mesh, index: Half_Edge_Index) -> Vertex {
	return mesh.verts[mesh.edges[index].vertex]
}

mesh_get_edge_source_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Vertex {
	return &mesh.verts[mesh.edges[mesh.edges[index].prev].vertex]
}

mesh_get_edge_target_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Vertex {
	return &mesh.verts[mesh.edges[index].vertex]
}

mesh_get_face_safe :: proc (mesh: Mesh, index: Face_Index) -> (face: Face, ok: bool) {
	for f in mesh.free_faces {
		if f == index {
			return {}, false
		}
	}

	return mesh.faces[index], true
}

mesh_get_vertex_safe :: proc (mesh: Mesh, index: Vertex_Index) -> (vertex: Vertex, ok: bool) {
	for v in mesh.free_verts {
		if v == index {
			return {}, false
		}
	}
	return mesh.verts[index], true
}

mesh_get_edge_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	for e in mesh.free_edges {
		if e == index {
			return {}, false
		}
	}
	return mesh.edges[index], true
}

mesh_get_face_ptr_safe :: proc (mesh: Mesh, index: Face_Index) -> (face: ^Face, ok: bool) {
	for f in mesh.free_faces {
		if f == index {
			return nil, false
		}
	}
	return &mesh.faces[index], true
}

mesh_get_vertex_ptr_safe :: proc (mesh: Mesh, index: Vertex_Index) -> (vertex: ^Vertex, ok: bool) {
	for v in mesh.free_verts {
		if v == index {
			return nil, false
		}
	}
	return &mesh.verts[index], true
}

mesh_get_edge_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	for e in mesh.free_edges {
		if e == index {
			return nil, false
		}
	}
	return &mesh.edges[index], true
}

mesh_get_edge_next_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := mesh_get_edge_safe(mesh, index) or_return
	return mesh_get_edge_safe(mesh, e.next)
}

mesh_get_edge_prev_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := mesh_get_edge_safe(mesh, index) or_return
	return mesh_get_edge_safe(mesh, e.prev)
}

mesh_get_edge_opposite_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := mesh_get_edge_safe(mesh, index) or_return
	return mesh_get_edge_safe(mesh, e.opposite)
}

mesh_get_edge_next_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := mesh_get_edge_ptr_safe(mesh, index) or_return
	return mesh_get_edge_ptr_safe(mesh, e.next)
}

mesh_get_edge_prev_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := mesh_get_edge_ptr_safe(mesh, index) or_return
	return mesh_get_edge_ptr_safe(mesh, e.prev)
}

mesh_get_edge_opposite_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := mesh_get_edge_ptr_safe(mesh, index) or_return
	return mesh_get_edge_ptr_safe(mesh, e.opposite)
}

mesh_get_edge_source_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: Vertex, ok: bool) {
	edge := mesh_get_edge_safe(mesh, index) or_return
	prev := mesh_get_edge_safe(mesh, edge.prev) or_return
	return mesh_get_vertex_safe(mesh, prev.vertex)
}

mesh_get_edge_target_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: Vertex, ok: bool) {
	edge := mesh_get_edge_safe(mesh, index) or_return
	return mesh_get_vertex_safe(mesh, edge.vertex)
}

mesh_get_edge_source_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: ^Vertex, ok: bool) {
	edge := mesh_get_edge_safe(mesh, index) or_return
	prev := mesh_get_edge_safe(mesh, edge.prev) or_return
	return mesh_get_vertex_ptr_safe(mesh, prev.vertex)
}

mesh_get_edge_target_ptr_safe :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: ^Vertex, ok: bool) {
	edge := mesh_get_edge_safe(mesh, index) or_return
	return mesh_get_vertex_ptr_safe(mesh, edge.vertex)
}

// TODO : Use the above procedures to make code more readable in certain places

mesh_create_face_edge_iterator :: proc(mesh: ^Mesh, face: Face_Index) -> Face_Edge_Iterator {
	return {
		current = mesh.faces[face].edge,
		start = mesh.faces[face].edge,
		mesh = mesh,
		step = 0,
	}
}

mesh_face_edge_forward_iter :: proc(iter: ^Face_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	prev := iter.current
	e := &iter.mesh.edges[iter.current]
	iter.current = e.next

	return e, prev, true
}

mesh_face_edge_backward_iter :: proc(iter: ^Face_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	prev := iter.current
	e := &iter.mesh.edges[iter.current]
	iter.current = e.prev

	return e, prev, true
}

mesh_create_vertex_edge_iterator :: proc(mesh: ^Mesh, vertex: Vertex_Index) -> Vertex_Edge_Iterator {
	edge := mesh.edges[mesh.verts[vertex].edge]

	if edge.vertex != vertex {return {}}

	start := mesh.verts[vertex].edge

	return {
		mesh = mesh,
		start = start,
		current = start,
		step = 0,
	}
}

mesh_vertex_incomming_edge_iter :: proc(iter: ^Vertex_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	e := &iter.mesh.edges[iter.current]
	prev := iter.current
	iter.current = iter.mesh.edges[e.next].opposite

	return e, prev, true
}

mesh_vertex_outgoing_edge_iter :: proc(iter: ^Vertex_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	e := &iter.mesh.edges[iter.current]
	e_r := &iter.mesh.edges[e.opposite]
	prev := e.opposite
	iter.current = iter.mesh.edges[e.next].opposite

	return e_r, prev, true
}

mesh_create_triangle_emitter_iter :: proc(mesh: ^Mesh) -> Triangle_Emitter_Iterator {
	face := mesh.faces[mesh.active_faces[0]]
	normal := mesh_calculate_face_normal(mesh, mesh.active_faces[0])
	return {
		mesh = mesh,
		face = mesh.active_faces[0],
		start = face.edge,
		edge = mesh.edges[face.edge].next,
		normal = normal
	}
}

mesh_triangle_emitter_indexed_flat_iter :: proc(iter: ^Triangle_Emitter_Iterator) -> (count: i32, positions: [3]Vec3f32, normal: Vec3f32, indices: [3]i32, ok: bool) {
	if iter.mesh.edges[iter.edge].next == iter.start { // loop till < n-1
		if iter.face_step < i32(len(iter.mesh.active_faces) - 1) {
			iter.walk_step = 0
			iter.face_step += 1
			iter.face = iter.mesh.active_faces[iter.face_step]
			iter.normal = mesh_calculate_face_normal(iter.mesh, iter.face)
			iter.start = iter.mesh.faces[iter.face].edge // 0
			iter.edge = iter.mesh.edges[iter.start].next // n
			iter.vertex_base = iter.vertex_step
		} else {
			return 0, 0, 0, 0, false
		}
	}

	edge := iter.mesh.edges[iter.edge]

	first := mesh_get_edge_target(iter.mesh^, iter.start) // 0
	n := mesh_get_edge_target(iter.mesh^, iter.edge) // n
	n_next := mesh_get_edge_target(iter.mesh^, edge.next) // n + 1

	if iter.walk_step == 0 {
		iter.vertex_step += 3
		iter.walk_step += 1
		iter.edge = edge.next
		return 3, {first.position, n.position, n_next.position}, iter.normal, {iter.vertex_base, iter.vertex_step - 2, iter.vertex_step - 1}, true
	}

	iter.walk_step += 1
	iter.vertex_step += 1
	iter.edge = edge.next

	return 1, {n_next.position, 0, 0}, iter.normal, {iter.vertex_base, iter.vertex_step - 2, iter.vertex_step - 1}, true
}

mesh_alloc_face :: proc(mesh: ^Mesh, face: Face) -> Face_Index {
	index, ok := pop_safe(&mesh.free_faces)
	defer append(&mesh.active_faces, index)
	if ok {
		mesh.faces[index] = face
		return index
	}
	index = Face_Index(len(mesh.faces))
	append(&mesh.faces, face)
	return index
}

mesh_alloc_half_edge :: proc(mesh: ^Mesh, half_edge: Half_Edge) -> Half_Edge_Index {
	index, ok := pop_safe(&mesh.free_edges)
	defer append(&mesh.active_edges, index)
	if ok {
		mesh.edges[index] = half_edge
		return index
	}
	index = Half_Edge_Index(len(mesh.edges))
	append(&mesh.edges, half_edge)
	return index
}

mesh_alloc_vertex :: proc(mesh: ^Mesh, vertex: Vertex) -> Vertex_Index {
	index, ok := pop_safe(&mesh.free_verts)
	defer append(&mesh.active_verts, index)
	if ok {
		mesh.verts[index] = vertex
		return index
	}
	index = Vertex_Index(len(mesh.verts))
	append(&mesh.verts, vertex)
	return index
}

mesh_alloc_vertices :: proc(mesh: ^Mesh, vertices: ..Vertex) {
	for v in vertices {
		mesh_alloc_vertex(mesh, v)
	}
}

mesh_free_half_edge :: proc(mesh: ^Mesh, edge: Half_Edge_Index) {
	if edge < 0 { return }
	append(&mesh.free_edges, edge)
	for i, j in mesh.active_edges {
		if i == edge {
			ordered_remove(&mesh.active_edges, j)
			return
		}
	}
}

mesh_free_vertex :: proc(mesh: ^Mesh, vertex: Vertex_Index) {
	if vertex < 0 { return }
	append(&mesh.free_verts, vertex)
	for i, j in mesh.active_verts {
		if i == vertex {
			ordered_remove(&mesh.active_verts, j)
			return
		}
	}
}

mesh_free_face :: proc(mesh: ^Mesh, face: Face_Index) {
	if face < 0 { return }
	append(&mesh.free_faces, face)
	for i, j in mesh.active_faces {
		if i == face {
			ordered_remove(&mesh.active_faces, j)
			return
		}
	}
}


mesh_dissolve_vertex_face_split :: proc(mesh: ^Mesh, vertex: Vertex_Index, temp_alloc := context.temp_allocator) -> (new_face: Face_Index) {
	// TODO: Should the indices be validated for being valid in the free list? or should the user be trusted?
	// Blender provides an option to dissolve vertex without face splits. TODO: Figure that out
	if vertex < 0 { return -1 }

	iter := mesh_create_vertex_edge_iterator(mesh, vertex)

	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
        if e.face == -1 {
         	// TODO: Add a way to dissolve boundary vertex?
            log.error("Cannot dissolve boundary vertex")
            return -1
        }
	}

	// If only two edges are incidence on the vertex, then delete the vertex and one pair of edges and reconnect the remaining edges
	if iter.step < 3 {
		v := mesh.verts[vertex]

		incomming_index := v.edge
		outgoing_index := mesh.edges[incomming_index].next

		incomming := &mesh.edges[incomming_index]
		outgoing := &mesh.edges[outgoing_index]
		incomming_op := &mesh.edges[incomming.opposite]
		outgoing_op := &mesh.edges[outgoing.opposite]

		source, target := incomming_op.vertex, outgoing_op.vertex

		outgoing.prev = incomming.prev
		outgoing_op.next = incomming_op.next
		outgoing_op.vertex = incomming_op.vertex

		mesh.edges[incomming.prev].next = outgoing_index
		mesh.edges[outgoing_op.next].prev = outgoing.opposite

		incomming_op_vertex := &mesh.verts[incomming_op.vertex]
		if incomming_op_vertex.edge == incomming.opposite {
			incomming_op_vertex.edge = outgoing.opposite
		}

		mesh_free_vertex(mesh, vertex)
		mesh_free_half_edge(mesh, incomming.opposite)
		mesh_free_half_edge(mesh, incomming_index)
		delete_key(&mesh.lookup, Lookup_Pair{source, vertex})
		delete_key(&mesh.lookup, Lookup_Pair{vertex, source})
		delete_key(&mesh.lookup, Lookup_Pair{target, vertex})
		delete_key(&mesh.lookup, Lookup_Pair{vertex, target})

		mesh.lookup[Lookup_Pair{source, target}] = outgoing_index
		mesh.lookup[Lookup_Pair{target, source}] = outgoing.opposite
		return
	}

	// Collect all the outgoing edges from the vertex. Split the face between target vertices of consecutive outgoing edges to get a rim.

	outgoing := make([dynamic]Half_Edge_Index, temp_alloc)

	iter = mesh_create_vertex_edge_iterator(mesh, vertex)
	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
		if mesh.verts[e.vertex].edge == i {
			mesh.verts[e.vertex].edge = mesh.edges[e.next].opposite
		}

		append(&outgoing, i)
	}

	for i := 0; i < len(outgoing); i += 1 {
		u := outgoing[i]
		v := outgoing[(i + 1) % len(outgoing)]

		if mesh.edges[u].next == v || mesh.edges[v].next == u { continue } // Skip adjacent vertices. This split function handles that but this is here to avoid logging the warnings

		mesh_split_face(mesh, mesh.edges[v].face, mesh.edges[u].vertex, mesh.edges[v].vertex)
	}

	face_edge := mesh.edges[mesh.edges[iter.start].opposite].next
	face := mesh_alloc_face(mesh, {-1})

	iter = mesh_create_vertex_edge_iterator(mesh, vertex)
	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
		if mesh.verts[e.vertex].edge == i {
			mesh.verts[e.vertex].edge = mesh.edges[e.next].opposite
		}

		mesh.edges[mesh.edges[i].next].prev = mesh.edges[mesh.edges[i].opposite].prev
		mesh.edges[mesh.edges[mesh.edges[i].opposite].prev].next = mesh.edges[i].next
		delete_key(&mesh.lookup, Lookup_Pair{vertex, e.vertex})
		delete_key(&mesh.lookup, Lookup_Pair{e.vertex, vertex})
		mesh_free_half_edge(mesh, i)
		mesh_free_half_edge(mesh, e.opposite)
		mesh_free_face(mesh, e.face)
	}

	s := mesh.edges[mesh.edges[mesh.verts[vertex].edge].opposite].next
	c := s

	for {
		mesh.edges[c].face = face
		c = mesh.edges[c].next
		if c == s {
			break
		}
	}

    mesh.faces[face].edge = face_edge
    mesh_free_vertex(mesh, vertex)

    return face
}

mesh_dissolve_half_edge :: proc(mesh: ^Mesh, edge: Half_Edge_Index) -> (kept_face: Face_Index) {
	if edge < 0 {return -1}

	// Wire edge.incomming-edge to the edge.opposite.outgoing edge and vice versa
	// Wire edge.outgoing-edge to the edge.opposite.incomming-edge and vice versa

	// This is blender's edge dissolve with the "dissolve vertex" option unselected

	e := &mesh.edges[edge]
	e_op := &mesh.edges[e.opposite]

	e_prev_index := e.prev
	e_op_prev_index := e_op.prev

	e_next_index := e.next
	e_op_next_index := e_op.next

	mesh.edges[e.next].prev = e_op_prev_index
	mesh.edges[e.prev].next = e_op_next_index

	mesh.edges[e_op.next].prev = e_prev_index
	mesh.edges[e_op.prev].next = e_next_index

	selected_edge_index := edge

	target_index := e.vertex
	source_index := e_op.vertex
	target, source := &mesh.verts[target_index], &mesh.verts[source_index]

	if target.edge == edge {
		target.edge = e_op.prev
	}

	if source.edge == e.opposite {
		source.edge = e.prev
	}

	if e.face == -1 {
		selected_edge_index = e.opposite
	}

	selected_edge := mesh.edges[selected_edge_index]
	to_be_deleted_index := selected_edge.opposite
	to_be_deleted := mesh.edges[selected_edge.opposite]

	if selected_edge.face != -1 {
		if mesh.faces[selected_edge.face].edge == selected_edge_index {
			mesh.faces[selected_edge.face].edge = selected_edge.next
		}

		iter := mesh_create_face_edge_iterator(mesh, selected_edge.face)
		for face_e in mesh_face_edge_forward_iter(&iter) {
			face_e.face = selected_edge.face
		}
	}

	delete_key(&mesh.lookup, Lookup_Pair{source_index, target_index})
	delete_key(&mesh.lookup, Lookup_Pair{target_index, source_index})
	mesh_free_half_edge(mesh, edge)
	mesh_free_half_edge(mesh, e.opposite)
	mesh_free_face(mesh, to_be_deleted.face)
	return selected_edge.face
}

mesh_dissolve_faces :: proc(mesh: ^Mesh, face_a: Face_Index, face_b: Face_Index) -> (kept_face: Face_Index) {
	// TODO: Make this take N-Faces instead of only two
	if face_a < 0 || face_b < 0 {
		return -1
	}

	common_edge := Half_Edge_Index(-1)

	iter := mesh_create_face_edge_iterator(mesh, face_a)
	for e, i in mesh_face_edge_forward_iter(&iter) {
		op := mesh.edges[e.opposite]
		if op.face == face_b {
			common_edge = i
			break
		}
	}

	return mesh_dissolve_half_edge(mesh, common_edge)
}

mesh_remove_face :: proc(mesh: ^Mesh, face: Face_Index) {
	iter := mesh_create_face_edge_iterator(mesh, face)
	for e in mesh_face_edge_forward_iter(&iter) {
		e.face = -1
	}

	mesh_free_face(mesh, face)
}

mesh_add_vertices :: proc(mesh: ^Mesh, positions: ..Vec3f32) {
	for position in positions {
		mesh_add_vertex(mesh, position)
	}
}

mesh_add_vertex :: proc(mesh: ^Mesh, position: Vec3f32) -> Vertex_Index {
	return mesh_alloc_vertex(mesh, {position = position, edge = -1})
}

mesh_add_faces :: proc(mesh: ^Mesh, faces: ..[]Vertex_Index) {
	for face in faces {
		mesh_add_face(mesh, face)
	}
}

mesh_add_face :: proc(mesh: ^Mesh, face: []Vertex_Index) -> Face_Index {
	n := len(face)

	first_edge_index := mesh_alloc_half_edge(mesh, {})
	face_index := mesh_alloc_face(mesh, Face{edge = first_edge_index})

	curr_edge_index := first_edge_index
	prev_edge_index := Half_Edge_Index(-1)

	for i := 0; i < n; i += 1 {
		curr_vert_index := face[i]
		next_vert_index := face[(i + 1) % n]

		curr_vert := &mesh.verts[curr_vert_index]
		next_vert := &mesh.verts[next_vert_index]

		lookup_pair := Lookup_Pair{curr_vert_index, next_vert_index}
		lookup_pair_op := Lookup_Pair{next_vert_index, curr_vert_index}

		op_index := mesh.lookup[lookup_pair_op] or_else -1
		mesh.lookup[lookup_pair] = curr_edge_index

		if op_index >= 0 {
			mesh.edges[op_index].opposite = curr_edge_index
		}

		if next_vert.edge < 0 {
			next_vert.edge = curr_edge_index
		}

		next_edge_index := Half_Edge_Index(-1)
		if i < n - 1 { next_edge_index = mesh_alloc_half_edge(mesh, {})
		} else { next_edge_index = first_edge_index }

		mesh.edges[curr_edge_index] = Half_Edge {
			face     = face_index,
			next     = next_edge_index,
			prev     = prev_edge_index,
			vertex   = next_vert_index,
			opposite = op_index,
		}

		prev_edge_index = curr_edge_index
		curr_edge_index = next_edge_index
	}

	mesh.edges[first_edge_index].prev = prev_edge_index
	return face_index
}

mesh_split_edges_all :: proc(mesh: ^Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.active_edges), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.active_edges), temp_alloc)
	copy(prev_edges[:], mesh.active_edges[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge(mesh, i, factor)
			lookup[mesh.edges[i].opposite] = {}
		}
	}
}

mesh_split_edges_twice_all :: proc(mesh: ^Mesh, factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	// Todo : Make a split varient for creating N-splits
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.active_edges), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.active_edges), temp_alloc)
	copy(prev_edges[:], mesh.active_edges[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge_twice(mesh, i, factor)
			lookup[mesh.edges[i].opposite] = {}
		}
	}
}

mesh_split_edge_twice :: proc(mesh: ^Mesh, half_edge_index: Half_Edge_Index, factor := f32(0.5)) {
	new_e_index := mesh_alloc_half_edge(mesh, {})
	new_e_op_index := mesh_alloc_half_edge(mesh, {})
	new_e1_index := mesh_alloc_half_edge(mesh, {})
	new_e1_op_index := mesh_alloc_half_edge(mesh, {})
	new_vertex_index := mesh_alloc_vertex(mesh, {})
	new_vertex1_index := mesh_alloc_vertex(mesh, {})

	e := &mesh.edges[half_edge_index]
	e_op := &mesh.edges[e.opposite]

	e_index := e_op.opposite
	e_op_index := e.opposite
	e_next_index := e.next
	e_op_next_index := e_op.next
	e_prev_index := e.prev
	e_op_prev_index := e_op.prev

	target_index, source_index := e.vertex, e_op.vertex
	target, source := &mesh.verts[target_index], &mesh.verts[source_index]

	new_vertex := &mesh.verts[new_vertex_index]
	new_vertex1 := &mesh.verts[new_vertex1_index]
	mid_point := (source.position + target.position) / 2
	new_vertex.position = target.position + (mid_point - target.position) * factor
	new_vertex1.position = source.position + (mid_point - source.position) * factor

	new_e := &mesh.edges[new_e_index]
	new_e_op := &mesh.edges[new_e_op_index]
	new_e1 := &mesh.edges[new_e1_index]
	new_e1_op := &mesh.edges[new_e1_op_index]

	new_vertex.edge = e_index
	new_vertex1.edge = new_e1_index
	source.edge = new_e1_op_index
	target.edge = new_e_index

	new_e^ = {
		prev     = e_index,
		next     = e_next_index,
		face     = e.face,
		vertex   = target_index,
		opposite = new_e_op_index,
	}

	new_e_op^ = {
		prev     = e_op.prev,
		next     = e_op_index,
		face     = e_op.face,
		vertex   = new_vertex_index,
		opposite = new_e_index,
	}

	new_e1^ = {
		prev = e_prev_index,
		next = e_index,
		face = e.face,
		opposite = new_e1_op_index,
		vertex = new_vertex1_index,
	}

	new_e1_op^ = {
		prev = e_op_index,
		next = e_op_next_index,
		face = e_op.face,
		vertex = source_index,
		opposite = new_e1_index,
	}

	mesh.edges[e_next_index].prev = new_e_index
	mesh.edges[e_op.prev].next = new_e_op_index

	mesh.edges[e_prev_index].next = new_e1_index
	mesh.edges[e_op_next_index].prev = new_e1_op_index

	e.next = new_e_index
	e_op.prev = new_e_op_index

	e.prev = new_e1_index
	e_op.next = new_e1_op_index

	e.vertex = new_vertex_index
	e_op.vertex = new_vertex1_index

	delete_key(&mesh.lookup, Lookup_Pair{source_index, target_index})
	delete_key(&mesh.lookup, Lookup_Pair{target_index, source_index})

	mesh.lookup[Lookup_Pair{source_index, new_vertex1_index}] 		= new_e1_index
	mesh.lookup[Lookup_Pair{new_vertex1_index, source_index}] 		= new_e1_op_index
	mesh.lookup[Lookup_Pair{new_vertex1_index, new_vertex_index}] 	= e_index
	mesh.lookup[Lookup_Pair{new_vertex_index, new_vertex1_index}] 	= e_op_index
	mesh.lookup[Lookup_Pair{new_vertex_index, target_index}] 		= new_e_index
	mesh.lookup[Lookup_Pair{target_index, new_vertex_index}] 		= new_e_op_index
}

// Factor defines where to place the vertex. Factor of 0 will place the vertex at the source, and factor of 1 will place the vertex at target
mesh_split_edge :: proc(mesh: ^Mesh, half_edge_index: Half_Edge_Index, factor := f32(0.5)) -> (Half_Edge_Index) {
	new_e_index := mesh_alloc_half_edge(mesh, {})
	new_e_op_index := mesh_alloc_half_edge(mesh, {})
	new_vertex_index := mesh_alloc_vertex(mesh, {})

	e := &mesh.edges[half_edge_index]
	e_op := &mesh.edges[e.opposite]

	e_index := e_op.opposite
	e_op_index := e.opposite
	e_next_index := e.next
	e_op_next_index := e_op.next

	e_next := &mesh.edges[e_next_index]
	e_op_next := &mesh.edges[e_op_next_index]

	target_index, source_index := e.vertex, e_op.vertex
	target, source := &mesh.verts[target_index], &mesh.verts[source_index]
	new_vertex := &mesh.verts[new_vertex_index]
	new_vertex.position = source.position + (target.position - source.position) * factor

	new_e := &mesh.edges[new_e_index]
	new_e_op := &mesh.edges[new_e_op_index]

	new_vertex.edge = e_index
	source.edge = e_op_index
	target.edge = new_e_index

	// We split edges such that the new edges are opposites of each other. This makes things quite simpler
	new_e^ = {
		prev     = e_index,
		next     = e_next_index,
		face     = e.face,
		vertex   = target_index,
		opposite = new_e_op_index,
	}

	new_e_op^ = {
		prev     = e_op.prev,
		next     = e_op_index,
		face     = e_op.face,
		vertex   = new_vertex_index,
		opposite = new_e_index,
	}

	mesh.edges[e_next_index].prev = new_e_index
	mesh.edges[e_op.prev].next = new_e_op_index

	e.next = new_e_index
	e_op.prev = new_e_op_index

	e.vertex = new_vertex_index

	delete_key(&mesh.lookup, Lookup_Pair{source_index, target_index})
	delete_key(&mesh.lookup, Lookup_Pair{target_index, source_index})
	mesh.lookup[Lookup_Pair{source_index, new_vertex_index}] = e_index
	mesh.lookup[Lookup_Pair{new_vertex_index, source_index}] = e_op_index
	mesh.lookup[Lookup_Pair{new_vertex_index, target_index}] = new_e_index
	mesh.lookup[Lookup_Pair{target_index, new_vertex_index}] = new_e_op_index

	return new_e_index
}

mesh_split_face :: proc(mesh: ^Mesh, face_index: Face_Index, a_index, b_index: Vertex_Index) -> (invalid: bool) {
	face := &mesh.faces[face_index]

	incomming_a_index, incomming_b_index := Half_Edge_Index(-1), Half_Edge_Index(-1)
	outgoing_a_index, outgoing_b_index := Half_Edge_Index(-1), Half_Edge_Index(-1)

	if a_index == b_index {
		log.warnf("Same vertex cannot be used to split a face")
		return true
	}

	{ 	// Validate inputs first
		found_a, found_b := false, false

		iter := mesh_create_face_edge_iterator(mesh, face_index)
		for e, i in mesh_face_edge_forward_iter(&iter) {
			if e.vertex == a_index {
				found_a = true
				incomming_a_index = i
			}
			if e.vertex == b_index {
				found_b = true
				incomming_b_index = i
			}
		}

		if !found_b || !found_a {
			log.warnf("Vertices do not belong to provided face. Face Index : %i, Face : %v, Vertex a : %i, Vertex b : %i", face_index, face, a_index, b_index)
			return true
		}
	}

	outgoing_a_index = mesh.edges[incomming_a_index].next
	outgoing_b_index = mesh.edges[incomming_b_index].next

	if mesh.edges[outgoing_a_index].next == outgoing_b_index || mesh.edges[outgoing_b_index].next == outgoing_a_index {
		log.warnf("Two adjacent vertices cannot be used to split a face. Vertex a : %i, Vertex b : %i", a_index, b_index)
		return true
	}

	a_to_b_index := mesh_alloc_half_edge(mesh, {})
	b_to_a_index := mesh_alloc_half_edge(mesh, {})
	a_to_b := &mesh.edges[a_to_b_index]
	b_to_a := &mesh.edges[b_to_a_index]

	incomming_a, incomming_b := &mesh.edges[incomming_a_index], &mesh.edges[incomming_b_index]
	outgoing_a, outgoing_b := &mesh.edges[outgoing_a_index], &mesh.edges[outgoing_b_index]

	incomming_a.next = a_to_b_index
	a_to_b.next = outgoing_b_index

	incomming_b.next = b_to_a_index
	b_to_a.next = outgoing_a_index

	outgoing_a.prev = b_to_a_index
	outgoing_b.prev = a_to_b_index

	b_to_a.prev = incomming_b_index
	a_to_b.prev = incomming_a_index

	a_to_b.vertex = b_index
	b_to_a.vertex = a_index

	a_to_b.opposite = b_to_a_index
	b_to_a.opposite = a_to_b_index

	a_to_b.face = face_index
	face.edge = a_to_b_index

	new_face_index := mesh_alloc_face(mesh, Face{edge = b_to_a_index})
	b_to_a.face = new_face_index

	{ 	// Set the correct face for half-edges
		iter := mesh_create_face_edge_iterator(mesh, face_index)
		for e in mesh_face_edge_forward_iter(&iter) {
			e.face = face_index
		}

		iter = mesh_create_face_edge_iterator(mesh, new_face_index)
		for e in mesh_face_edge_forward_iter(&iter) {
			e.face = new_face_index
		}
	}
	return false
}

mesh_add_boundaries :: proc(mesh: ^Mesh) {
	total_non_boundary_edges := len(mesh.active_edges)
	for i in mesh.active_edges {
		e := mesh.edges[i]
		if e.opposite != -1 {continue}

		op := Half_Edge {
			face     = -1,
			opposite = Half_Edge_Index(i),
			vertex   = mesh.edges[e.prev].vertex,
			next     = -1,
			prev     = -1,
		}

		op_index := mesh_alloc_half_edge(mesh, op)
		mesh.edges[i].opposite = op_index

		mesh.lookup[Lookup_Pair{op.vertex, e.vertex}] = op_index
	}

	for i in mesh.active_edges[total_non_boundary_edges:] {
		b_index := Half_Edge_Index(i)
		b := &mesh.edges[b_index]

		curr := b.opposite
		for {
			curr = mesh.edges[mesh.edges[curr].prev].opposite
			if mesh.edges[curr].face == -1 {
				b.next = curr
				mesh.edges[curr].prev = b_index
				break
			}
		}
	}
}

mesh_triangulate_face_from_centroid :: proc (mesh: ^Mesh, face: Face_Index, height := f32(0), temp_alloc := context.temp_allocator) -> Vertex_Index {
    collected_edges := make([dynamic]Half_Edge_Index, temp_alloc)
	centroid := Vec3f32{}
    iter := mesh_create_face_edge_iterator(mesh, face)
    for e, i in mesh_face_edge_forward_iter(&iter) {
        centroid += mesh.verts[e.vertex].position
        append(&collected_edges, i)
    }

	vertex_count := len(collected_edges)

	if vertex_count < 3 {
		return -1
	}

	normal := mesh_calculate_face_normal(mesh, face)
	centroid /= f32(vertex_count)
	centroid = centroid + normal * height

	centroid_vertex_index := mesh_add_vertex(mesh, centroid)

    for i := 0; i < vertex_count; i += 1 {
		edge_index := collected_edges[i]
		source := mesh.edges[mesh.edges[edge_index].opposite].vertex
		target := mesh.edges[edge_index].vertex

		edge_to_centroid_index := mesh_alloc_half_edge(mesh, {})
		edge_from_centroid_index := mesh_alloc_half_edge(mesh, {})
		new_face := mesh_alloc_face(mesh, {edge_index})

		mesh.verts[centroid_vertex_index].edge = edge_to_centroid_index

		edge_ptr := &mesh.edges[edge_index]

		edge_to_centroid := Half_Edge{face = new_face, next = edge_from_centroid_index, prev = edge_index, vertex = centroid_vertex_index}
		edge_from_centroid := Half_Edge{face = new_face, prev = edge_to_centroid_index, next = edge_index, vertex = source}

		edge_to_centroid_op := mesh.lookup[Lookup_Pair{centroid_vertex_index, target}] or_else -1
		edge_from_centroid_op := mesh.lookup[Lookup_Pair{source, centroid_vertex_index}] or_else -1

		edge_to_centroid.opposite = edge_to_centroid_op
		edge_from_centroid.opposite = edge_from_centroid_op

		if edge_from_centroid_op != -1 {
			mesh.edges[edge_from_centroid_op].opposite = edge_from_centroid_index
		}

		if edge_to_centroid_op != -1 {
			mesh.edges[edge_to_centroid_op].opposite = edge_to_centroid_index
		}

		edge_ptr.face = new_face
		edge_ptr.next = edge_to_centroid_index
		edge_ptr.prev = edge_from_centroid_index

		mesh.faces[new_face].edge = edge_to_centroid_index

		mesh.edges[edge_to_centroid_index] = edge_to_centroid
		mesh.edges[edge_from_centroid_index] = edge_from_centroid
		mesh.lookup[Lookup_Pair{target, centroid_vertex_index}] = edge_to_centroid_index
		mesh.lookup[Lookup_Pair{centroid_vertex_index, source}] = edge_from_centroid_index
	}

	mesh_free_face(mesh, face)
	return centroid_vertex_index
}

mesh_triangulate_face_from_vertex :: proc(mesh: ^Mesh, face: Face_Index, vertex: Vertex_Index) {}

mesh_calculate_face_normal :: proc(mesh: ^Mesh, face: Face_Index) -> Vec3f32 {
	normal := Vec3f32{}
	iter := mesh_create_face_edge_iterator(mesh, face)
	for e_c in mesh_face_edge_forward_iter(&iter) {
		v_c := mesh.verts[e_c.vertex].position
		v_n := mesh.verts[mesh.edges[e_c.next].vertex].position

		normal.x += (v_n.y - v_c.y) * (v_n.z + v_c.z)
		normal.y += (v_n.z - v_c.z) * (v_n.x + v_c.x)
		normal.z += (v_n.x - v_c.x) * (v_n.y + v_c.y)
	}
	return linalg.normalize0(normal)
}

mesh_normalize_onto_sphere :: proc(mesh: ^Mesh) {
	length := f32(0)
	for v in mesh.active_verts {
		length = max(length, linalg.length(mesh.verts[v].position))
	}

	for v in mesh.active_verts {
		mesh.verts[v].position /= length
	}
}

// Convay operations
mesh_convay_dual :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
    dual := mesh_create(mesh.edges.allocator)

    dual_verts := make([dynamic]Vertex_Index, len(mesh.faces), temp_alloc)

    for f in mesh.active_faces {
		centroid := Vec3f32{}
		iter := mesh_create_face_edge_iterator(mesh, f)
		for e, i in mesh_face_edge_backward_iter(&iter) {
			centroid += mesh.verts[e.vertex].position
		}
		centroid /= f32(iter.step)
		dual_verts[f] = mesh_add_vertex(&dual, centroid)
    }

	verts := make([dynamic]Vertex_Index, temp_alloc)
    for v in mesh.active_verts {
        iter := mesh_create_vertex_edge_iterator(mesh, v)
        for e in mesh_vertex_outgoing_edge_iter(&iter) {
            append(&verts, dual_verts[e.face])
        }
        slice.reverse(verts[:])
        mesh_add_face(&dual, verts[:])
		clear(&verts)
	}

	mesh_destroy(mesh^)
	mesh^ = dual
}

mesh_convay_kis :: proc(mesh: ^Mesh, kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	faces := make([dynamic]Face_Index, len(mesh.active_faces), temp_alloc)
	copy(faces[:], mesh.active_faces[:])

	for f in faces {
		mesh_triangulate_face_from_centroid(mesh, f, kis_height, temp_alloc)
	}
}

mesh_convay_ambo :: proc (mesh: ^Mesh, ambo_factor := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]Vertex_Index, len(mesh.active_verts), temp_alloc)
	copy(verts[:], mesh.active_verts[:])

	mesh_split_edges_all(mesh, ambo_factor, temp_alloc)

	for v in verts {
		mesh_dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

mesh_convay_truncate :: proc(mesh: ^Mesh, truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	verts :=  make([dynamic]Vertex_Index, len(mesh.active_verts), temp_alloc)
	copy(verts[:], mesh.active_verts[:])

	mesh_split_edges_twice_all(mesh, truncate_factor, temp_alloc)

	for v in verts {
		mesh_dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

mesh_convay_snub :: proc(mesh: ^Mesh, truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_gyro(mesh, truncate_factor, gyro_height, temp_alloc)
	mesh_convay_kis(mesh, kis_height, temp_alloc)
}

mesh_convay_gyro :: proc(mesh: ^Mesh, truncate_factor := f32(2.0/3.0), height := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]Vertex_Index, len(mesh.active_verts), temp_alloc)
	faces := make([dynamic]Face_Index, len(mesh.active_faces), temp_alloc)
	copy(verts[:], mesh.active_verts[:])
	copy(faces[:], mesh.active_faces[:])

	mesh_split_edges_twice_all(mesh, truncate_factor, temp_alloc)

	for f in faces {
		centroid_vert := mesh_triangulate_face_from_centroid(mesh, f, height, temp_alloc)

		to_dissolve := make([dynamic]Half_Edge_Index, temp_alloc)
		for v in verts {
			iter_v := mesh_create_vertex_edge_iterator(mesh, v)
			for e, i in mesh_vertex_outgoing_edge_iter(&iter_v) {
				if e.vertex == centroid_vert {
					append(&to_dissolve, i)
					break
				}
			}
		}

		for e in to_dissolve {
			mesh_dissolve_half_edge(mesh, mesh.edges[e].next)
			mesh_dissolve_half_edge(mesh, e)
		}

		delete(to_dissolve)
	}
}

mesh_convay_classical_alternation :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	lookup := make(map[Vertex_Index]bool, temp_alloc)
	queue := make([dynamic]Vertex_Index, temp_alloc)

	append(&queue, mesh.active_verts[0])
	lookup[queue[0]] = false

	for len(queue) > 0 {
		v := pop_front(&queue)

		iter := mesh_create_vertex_edge_iterator(mesh, v)
		for e in mesh_vertex_outgoing_edge_iter(&iter) {
			u := e.vertex
			if u not_in lookup {
				lookup[u] = !lookup[v]
				append(&queue, u)
			} else if lookup[u] == lookup[v] {
				log.warnf("Cannot alternate a polyhedron with non-even faces")
				return
			}
		}
	}

	for k, v in lookup {
		if v == true {
			mesh_dissolve_vertex_face_split(mesh, k, temp_alloc)
		}
	}

}

mesh_convay_classical_snub :: proc(mesh: ^Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	mesh_convay_ambo(mesh, ambo_factor, temp_alloc)
	mesh_convay_truncate(mesh, truncate_factor, temp_alloc)
	mesh_convay_classical_alternation(mesh, temp_alloc)
}

mesh_convay_classical_gyro :: proc(mesh: ^Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	mesh_convay_classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_bevel :: proc(mesh: ^Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	mesh_convay_ambo(mesh, ambo_factor, temp_alloc)
	mesh_convay_truncate(mesh, truncate_factor, temp_alloc)
}

mesh_convay_expand :: proc(mesh: ^Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_ambo(mesh, factor, temp_alloc)
	mesh_convay_ambo(mesh, factor, temp_alloc)
}

mesh_convay_join :: proc(mesh: ^Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_ambo(mesh, factor, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_ortho :: proc(mesh: ^Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_expand(mesh, factor, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_meta :: proc(mesh: ^Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_needle :: proc(mesh: ^Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_kis(mesh, height, temp_alloc)
}

mesh_convay_zip :: proc(mesh: ^Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	mesh_convay_kis(mesh, height, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_operation :: proc(mesh: ^Mesh, operation: Convay_Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	switch operation {
		case .Kis:						mesh_convay_kis(mesh, kis_height, temp_alloc)
		case .Zip:						mesh_convay_zip(mesh, kis_height, temp_alloc)
		case .Ambo:						mesh_convay_ambo(mesh, ambo_factor, temp_alloc)
		case .Dual:						mesh_convay_dual(mesh, temp_alloc)
		case .Snub:						mesh_convay_snub(mesh, truncate_factor, gyro_height, kis_height, temp_alloc)
		case .Join:						mesh_convay_join(mesh, ambo_factor, temp_alloc)
		case .Meta:						mesh_convay_meta(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Gyro:						mesh_convay_gyro(mesh, truncate_factor, gyro_height, temp_alloc)
		case .Ortho:					mesh_convay_ortho(mesh, ambo_factor, temp_alloc)
		case .Bevel:					mesh_convay_bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Needle:					mesh_convay_needle(mesh, kis_height, temp_alloc)
		case .Expand:					mesh_convay_expand(mesh, ambo_factor, temp_alloc)
		case .Truncate:					mesh_convay_truncate(mesh, truncate_factor, temp_alloc)
		case .Classical_Snub:			mesh_convay_classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Classical_Gyro:			mesh_convay_classical_gyro(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Classical_Alternation:	mesh_convay_classical_alternation(mesh, temp_alloc)
	}
}

mesh_convay_operations :: proc(mesh: ^Mesh, operations: ..Convay_Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	for operation in operations {
		mesh_convay_operation(mesh, operation, temp_alloc, ambo_factor, truncate_factor, gyro_height, kis_height)
	}
}

// Polygon generation
mesh_generate_tetrahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, {1, 1, 1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1})
	mesh_add_faces(&mesh, {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2})
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_cube :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1})
	mesh_add_faces(&mesh, {3, 2, 1, 0}, {3, 0, 5, 4}, {4, 5, 6, 7}, {7, 6, 1, 2}, {6, 5, 0, 1}, {2, 3, 4, 7})
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_octahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, { 1,  0,  0},{-1,  0,  0},{ 0,  1,  0},{ 0, -1,  0},{ 0,  0,  1},{ 0,  0, -1})
	mesh_add_faces(&mesh,{0, 4, 2},{2, 4, 1},{1, 4, 3},{3, 4, 0},{2, 5, 0},{1, 5, 2},{3, 5, 1},{0, 5, 3})
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_icosahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh,
		{-1,  PHI,  0}, { 1,  PHI,  0}, {-1, -PHI,  0}, { 1, -PHI,  0},
		{ 0, -1,  PHI}, { 0,  1,  PHI}, { 0, -1, -PHI}, { 0,  1, -PHI},
		{ PHI,  0, -1}, { PHI,  0,  1}, {-PHI,  0, -1}, {-PHI,  0,  1},
	)
	mesh_add_faces(&mesh,
		{5, 11, 0}, {1, 5, 0}, {7, 1, 0}, {10, 7, 0},
		{11, 10, 0}, {9, 5, 1}, {4, 11, 5}, {2, 10, 11},
	 	{6, 7, 10}, {8, 1, 7}, {4, 9, 3}, {2, 4, 3},
		{6, 2, 3}, {8, 6, 3}, {9, 8, 3}, {5, 9, 4},
		{11, 4, 2}, {10, 2, 6}, {7, 6, 8}, {1, 8, 9}
	)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_dodecahedron :: proc(allocator := context.allocator) -> Mesh {
    mesh := mesh_create(allocator)

    mesh_add_vertices(&mesh,
        { 1,  1,  1}, { 1,  1, -1}, { 1, -1,  1}, { 1, -1, -1}, // 0-3: Cube vertices
        {-1,  1,  1}, {-1,  1, -1}, {-1, -1,  1}, {-1, -1, -1}, // 4-7: Cube vertices
        { 0, INV_PHI,  PHI}, { 0, INV_PHI, -PHI}, { 0, -INV_PHI,  PHI}, { 0, -INV_PHI, -PHI}, // 8-11
        { INV_PHI,  PHI, 0}, { INV_PHI, -PHI, 0}, {-INV_PHI,  PHI, 0}, {-INV_PHI, -PHI, 0}, // 12-15
        { PHI, 0,  INV_PHI}, { PHI, 0, -INV_PHI}, {-PHI, 0,  INV_PHI}, {-PHI, 0, -INV_PHI}, // 16-19
    )

    mesh_add_faces(&mesh,
        {0, 16, 2, 10, 8},   {0, 8, 4, 14, 12},   {0, 12, 1, 17, 16},
        {3, 17, 1, 9, 11},   {3, 11, 7, 15, 13},  {3, 13, 2, 16, 17},
        {5, 9, 1, 12, 14},   {5, 14, 4, 18, 19},  {5, 19, 7, 11, 9},
        {6, 10, 2, 13, 15},  {6, 15, 7, 19, 18},  {6, 18, 4, 8, 10},
    )
    mesh_normalize_onto_sphere(&mesh)
    return mesh
}

mesh_generate_platonic_solid :: proc(solid : Platonic_Solid, allocator := context.allocator) -> Mesh {
	switch solid {
	case .Cube:			return mesh_generate_cube(allocator)
	case .Octahedron:	return mesh_generate_octahedron(allocator)
	case .Tetrahedron:	return mesh_generate_tetrahedron(allocator)
	case .Dodecaheron:	return mesh_generate_dodecahedron(allocator)
	case .Icosahedron:	return mesh_generate_icosahedron(allocator)
	}
	unreachable()
}

mesh_generate_all_platonic_solids :: proc(allocator := context.allocator) -> [Platonic_Solid]Mesh {
	solids := [Platonic_Solid]Mesh{}
	for &s, t in solids {
		s = mesh_generate_platonic_solid(t, allocator)
	}
	return solids
}

mesh_generate_truncated_tetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_tetrahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_octahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_rhombicuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_snub_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_classical_snub(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_icosahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_rhombicosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_truncated_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc = temp_alloc)
	mesh_convay_truncate(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_snub_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_classical_snub(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_archimedean_solid :: proc (solid : Archimedean_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	switch solid {
	case .Snub_Cube:					return mesh_generate_snub_cube(allocator, temp_alloc)
	case .Cuboctahedron:		 		return mesh_generate_cuboctahedron(allocator, temp_alloc)
	case .Truncated_Cube:		 		return mesh_generate_truncated_cube(allocator, temp_alloc)
	case .Snub_dodecahedron:			return mesh_generate_snub_dodecahedron(allocator, temp_alloc)
	case .Icosidodecahedron:			return mesh_generate_icosidodecahedron(allocator, temp_alloc)
	case .Rhombicuboctahedron:	 		return mesh_generate_rhombicuboctahedron(allocator, temp_alloc)
	case .Truncated_Octahedron:	 		return mesh_generate_truncated_octahedron(allocator, temp_alloc)
	case .Truncated_Tetrahedron: 		return mesh_generate_truncated_tetrahedron(allocator, temp_alloc)
	case .Truncated_icosahedron:		return mesh_generate_truncated_icosahedron(allocator, temp_alloc)
	case .Truncated_dodecahedron:		return mesh_generate_truncated_dodecahedron(allocator, temp_alloc)
	case .Rhombicosidodecahedron:		return mesh_generate_rhombicosidodecahedron(allocator, temp_alloc)
	case .Truncated_Cuboctahedron:		return mesh_generate_truncated_cuboctahedron(allocator, temp_alloc)
	case .Truncated_Icosidodecahedron:	return mesh_generate_truncated_icosidodecahedron(allocator, temp_alloc)
	}
	unreachable()
}

mesh_generate_all_archimedean_solids :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> [Archimedean_Solid]Mesh {
	solids := [Archimedean_Solid]Mesh{}
	for &s, t in solids {
		s = mesh_generate_archimedean_solid(t, allocator, temp_alloc)
	}
	return solids
}

mesh_generate_triakis_tetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_tetrahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_TRI_TETRAHEDRON_KIS_HEIGHT, temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_rhombic_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_join(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_triakis_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_octahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_TRI_OCTAHEDRON_KIS_HEIGHT, temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_tetrakis_hexahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_kis(&mesh, CATALAN_TETRA_HEXAHEDRON_KIS_HEIGHT, temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_deltoidal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ortho(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_disdyakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_meta(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_pentagonal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_classical_gyro(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_rhombic_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_join(&mesh)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_triakis_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_icosahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_TRI_ICOSAHEDRON_KIS_HEIGHT, temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_pentakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_PENTA_DODECAHEDRON_KIS_HEIGHT, temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_deltoidal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ortho(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_disdyakis_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_meta(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_pentagonal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_classical_gyro(&mesh, temp_alloc = temp_alloc)
	mesh_normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_catalan_solid :: proc(solid : Catalan_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	switch solid {
	case .Triakis_Octahedron:		 	return mesh_generate_triakis_octahedron(allocator, temp_alloc)
	case .Triakis_Tetrahedron:			return mesh_generate_triakis_tetrahedron(allocator, temp_alloc)
	case .Triakis_Icosahedron:			return mesh_generate_triakis_icosahedron(allocator, temp_alloc)
	case .Tetrakis_Hexahedron:			return mesh_generate_tetrakis_hexahedron(allocator, temp_alloc)
	case .Rhombic_Dodecahedron:			return mesh_generate_rhombic_dodecahedron(allocator, temp_alloc)
	case .Pentakis_Dodecahedron:		return mesh_generate_pentakis_dodecahedron(allocator, temp_alloc)
	case .Disdyakis_Dodecahedron:		return mesh_generate_disdyakis_dodecahedron(allocator, temp_alloc)
	case .Rhombic_Triacontahedron:		return mesh_generate_rhombic_triacontahedron(allocator, temp_alloc)
	case .Disdyakis_Triacontahedron:	return mesh_generate_disdyakis_triacontahedron(allocator, temp_alloc)
	case .Deltoidal_Hexecontahedron:	return mesh_generate_deltoidal_hexecontahedron(allocator, temp_alloc)
	case .Deltoidal_Icositetrahedron:	return mesh_generate_deltoidal_icositetrahedron(allocator, temp_alloc)
	case .Pentagonal_Hexecontahedron:	return mesh_generate_pentagonal_hexecontahedron(allocator, temp_alloc)
	case .Pentagonal_Icositetrahedron:	return mesh_generate_pentagonal_icositetrahedron(allocator, temp_alloc)
	}
	unreachable()
}

mesh_generate_all_catalan_solids :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> [Catalan_Solid]Mesh {
	solids := [Catalan_Solid]Mesh{}
	for &s, t in solids {
		s = mesh_generate_catalan_solid(t, allocator, temp_alloc)
	}
	return solids
}

mesh_generate_sierpinski_tetrahedron :: proc(depth: int, allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	p0 := [3]f32{1, 1, 1}
	p1 := [3]f32{-1, -1, 1}
	p2 := [3]f32{-1, 1, -1}
	p3 := [3]f32{1, -1, -1}

	subdivide :: proc(mesh: ^Mesh, a, b, c, d: [3]f32, depth: int) {
		if depth == 0 {
			ia := mesh_add_vertex(mesh, a)
			ib := mesh_add_vertex(mesh, b)
			ic := mesh_add_vertex(mesh, c)
			id := mesh_add_vertex(mesh, d)

			mesh_add_faces(mesh, {ia, ib, ic}, {ia, ic, id}, {ia, id, ib}, {ib, id, ic})
			return
		} else {
			mab := (a + b) * 0.5
			mac := (a + c) * 0.5
			mad := (a + d) * 0.5
			mbc := (b + c) * 0.5
			mbd := (b + d) * 0.5
			mcd := (c + d) * 0.5

			subdivide(mesh, a, mab, mac, mad, depth - 1)
			subdivide(mesh, mab, b, mbc, mbd, depth - 1)
			subdivide(mesh, mac, mbc, c, mcd, depth - 1)
			subdivide(mesh, mad, mbd, mcd, d, depth - 1)
		}
	}

	subdivide(&mesh, p0, p1, p2, p3, depth)

	mesh_add_boundaries(&mesh)
	return mesh
}

// Validations
mesh_validate :: proc(mesh: Mesh, caller_expression := #caller_expression) -> (invalid: bool) {
	log.info("--- Validation :", caller_expression)
	base_ok := mesh_validate_base_constraints(mesh)
	indices_ok := mesh_validate_indicies(mesh)
	opposites_ok := mesh_validate_opposites(mesh)
	link_ok := mesh_validate_links(mesh)
	face_loops_ok := mesh_validate_face_loops(mesh)
	vertex_references_ok := mesh_validate_vertex_references(mesh)
	edge_uniqueness_ok := mesh_validate_global_uniqueness(mesh)

	if !base_ok {log.info("Base constraint test passed")} else {log.errorf("Base constraint test failed")}
	if !link_ok {log.info("Edge link test passed")} else {log.error("Edge link test failed")}
	if !indices_ok {log.info("Indices test passed")} else {log.error("Indicies test failed")}
	if !opposites_ok {log.info("Opposites test passed")} else {log.error("Opposites test failed")}
	if !face_loops_ok {log.info("Face loops test passed")} else {log.error("Face loops test failed")}
	if !vertex_references_ok {log.info("Vertex reference test passed")} else {log.error("Vertex reference test failed")}
	if !edge_uniqueness_ok {log.info("Global half-edge uniqueness test passed")} else {log.error("Global half-edge uniqueness test failed")}

	return(
		base_ok ||
		indices_ok ||
		opposites_ok ||
		face_loops_ok ||
		vertex_references_ok ||
		edge_uniqueness_ok \
	)
}

mesh_validate_base_constraints :: proc (mesh: Mesh) -> (invalide: bool) {
	for i in mesh.active_edges {
		next, ok_next := mesh_get_edge_next_safe(mesh, i)
		prev, ok_prev := mesh_get_edge_prev_safe(mesh, i)
		opposite, ok_opposite := mesh_get_edge_opposite_safe(mesh, i)

		if !(ok_next || ok_prev || ok_opposite) {
			log.errorf("Invalid Next: %v, Invalid Prev: %v, Invalid Opposite: %v. Edge: %v", !ok_next, !ok_prev, !ok_opposite, mesh_get_edge(mesh, i))
			return true
		}
	}

	for i in mesh.active_faces {
		face := mesh_get_face(mesh, i)
		face_edge, ok := mesh_get_edge_safe(mesh, face.edge)

		if !ok {
			log.errorf("Face does not reference an existing edge. Face: %v", face)
			return true
		}
	}

	for i in mesh.active_verts {
		vert := mesh_get_vertex(mesh, i)
		vert_edge, ok := mesh_get_edge_safe(mesh, vert.edge)

		if !ok {
			log.errorf("Vertex does not reference an existing edge. Vertex Index: %i, Vertex: %v", i, vert)
			return true
		}
	}

	return false
}

mesh_validate_opposites :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.active_edges {
		e := mesh_get_edge(mesh, i)
		e_op := mesh_get_edge_opposite(mesh, i)

		if e_op.face == e.face && e.face != -1 && e_op.face != -1 { // Opposites must belong to differnet/unique faces. Except boundaries
			log.errorf("Half-edge opposite must belong to a different face. Edge's Face : %i, Opposite's Face : %i", e.face, e_op.face)
			return true
		}

		// Invarient : Opposite of the Opposite must point to current edge
		if (e_op.opposite != i) {
			log.errorf("opposite(opposite(h)) property violated. Edge : %v, Opposite : %v", e, e_op)
			return true
		}
	}
	return false
}

mesh_validate_indicies :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.active_edges {
		e := mesh_get_edge(mesh, i)
		if e.opposite < 0 || int(e.opposite) >= len(mesh.edges) {
			log.errorf("Invalide opposite half-edge index. Got %i, expected between [0-%i)", e.opposite, len(mesh.edges))
			return true
		}

		if e.next < 0 || int(e.next) >= len(mesh.edges) {
			log.errorf("Invalide next half-edge index. Got %i, expected between [0-%i)", e.next, len(mesh.edges))
			return true
		}

		if e.face < -1 || int(e.face) >= len(mesh.faces) {
			log.errorf("Invalide face index. Got %i, expected between [-1-%i)", e.face, len(mesh.faces))
			return true
		}

		if e.vertex < 0 || int(e.vertex) >= len(mesh.verts) {
			log.errorf("Invalide vertex index. Got %i, expected between [0-%i)", e.vertex, len(mesh.verts))
			return true
		}

		if e.prev < 0 || int(e.prev) >= len(mesh.edges) {
			log.errorf("Invalid prev half-edge index. Got %i, expected between [0-%i)", e.prev, len(mesh.edges))
		}
	}
	return false
}

mesh_validate_links :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.active_edges {
		e := mesh_get_edge(mesh, i)

		next := mesh_get_edge_next(mesh, i)
		prev := mesh_get_edge_prev(mesh, i)

		if next.prev != i {
			log.errorf("Linkage Invariance Violated: edges[%i].next.prev is %i, expected %i", i, next.prev, i)
			return true
		}

		if prev.next != i {
			log.errorf("Linkage Invariance Violated: edges[%i].prev.next is %i, expected %i", i, prev.next, i)
			return true
		}

		e_op := mesh_get_edge_opposite(mesh, i)
		if e_op.vertex != prev.vertex {
			log.errorf( "current.source == prev.target property violated. Current, Prev : %v, %v", e, prev)
			return true
		}
	}
	return false
}

mesh_validate_face_loops :: proc( mesh: Mesh, allocator := context.temp_allocator, ) -> ( invalid: bool, ) {
	lookup := make(map[Half_Edge_Index]struct{}, allocator)
	defer delete(lookup)

	for i in mesh.active_faces {
		f := mesh.faces[i]
		defer clear(&lookup)
		start := f.edge

		if mesh.edges[start].face != i { 	// The edge must point/reference current face
			log.errorf("Half-edge doesn't point to current face being walked. Face Index : %i, Face : %v, Edge : %v", i, f, mesh.edges[start])
			return true
		}

		curr := start
		walk_distance := int(0)

		lookup[start] = {}

		for {
			prev := curr
			curr = mesh.edges[curr].next

			if mesh.edges[curr].face != i {
				log.errorf( "Half-edge doesn't point to current face being walked. Face Index : %i, Face : %v, Edge : %v", i, f, mesh.edges[curr])
				return true
			}

			walk_distance += 1

			if curr == start {
				break
			}

			_, exists := lookup[curr]

			if exists { 	// Check locally that no repition happens in the face walk
				log.errorf("Half-edges repeated during face walk. Edge : %v", mesh.edges[curr])
				return true
			}

			lookup[curr] = {}
		}

		if walk_distance < 2 { 	// degenerate face
			log.errorf( "Degenerate Face consisting of only single half-edge. Face Index : %i, Face : %v", i, f, )
			return true
		}
	}
	return false
}

mesh_validate_vertex_references :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.active_verts {
		v := mesh.verts[i]
		if v.edge < 0 || int(v.edge) >= len(mesh.edges) {
			log.errorf("Vertex references invalid half-edge. Vertex Index : %i, Vertex : %v", i, v)
			return true
		}

		if mesh.edges[v.edge].vertex != i {
			log.errorf("Half-edge vertex references does not have the vertex as its target. Vertex Index : %i, Vertex : %v, Edge : %v", i, v, mesh.edges[v.edge])
			return true
		}
	}
	return false
}

mesh_validate_global_uniqueness :: proc(mesh: Mesh, allocator := context.temp_allocator) -> (failed: bool) {
	visited := make(map[Half_Edge_Index]struct{}, allocator)
	defer delete(visited)

	for i in mesh.active_edges {
		curr_idx := Half_Edge_Index(i)
		if curr_idx in visited do continue

		start_idx := curr_idx
		walk_idx := curr_idx

		for {
			if walk_idx in visited {
				log.errorf("Half-edge %i visited more than once. Mesh has non-manifold or corrupt links.", walk_idx)
				return true
			}
			visited[walk_idx] = {}

			walk_idx = mesh.edges[walk_idx].next

			if walk_idx == start_idx do break

			if len(visited) > len(mesh.edges) {
				log.errorf("Face walk exceeded edge count. Circularity broken.")
				return true
			}
		}
	}

	if len(visited) != len(mesh.active_edges) {
		log.errorf("Orphaned edges detected: %i / %i reached", len(visited), len(mesh.edges))
		return true
	}

	return false
}
