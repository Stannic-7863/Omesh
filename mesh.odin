package mesh

import "vendor:glfw/bindings"
import "core:encoding/csv"
import "core:container/queue"
import "core:slice"
import "core:image/netpbm"
import "core:math/linalg"
import "base:intrinsics"
import "base:runtime"
import "core:log"

/*

	Gonna paste some good stuff here
	-	https://levskaya.github.io/polyhedronisme
	-

*/

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

Free_List :: struct ($type: typeid, $type_index: typeid) where intrinsics.type_is_integer(type_index) {
	active: [dynamic]type_index,
	free:	[dynamic]type_index,
	all:	[dynamic]type,
}

Mesh :: struct {
	faces:      Free_List(Face, Face_Index),
	verts:      Free_List(Vertex, Vertex_Index),
	edges:      Free_List(Half_Edge, Half_Edge_Index),
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
	free_list_create(&m.faces, allocator)
	free_list_create(&m.edges, allocator)
	free_list_create(&m.verts, allocator)
	m.lookup = make(map[Lookup_Pair]Half_Edge_Index, allocator)
	return m
}

mesh_destroy :: proc(mesh: Mesh) {
	delete(mesh.lookup)

	free_list_destroy(mesh.faces)
	free_list_destroy(mesh.edges)
	free_list_destroy(mesh.verts)
}

meshes_destroy :: proc(meshes: ..Mesh) {
	for mesh in meshes {
		mesh_destroy(mesh)
	}
}

mesh_create_face_edge_iterator :: proc(mesh: ^Mesh, face: Face_Index) -> Face_Edge_Iterator {
	f := mesh_get_face_unsafe(mesh^, face)
	return {
		current = f.edge,
		start = f.edge,
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
	e := mesh_get_edge_ptr_unsafe(iter.mesh^, iter.current)
	iter.current = e.next

	return e, prev, true
}

mesh_face_edge_backward_iter :: proc(iter: ^Face_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	prev := iter.current
	e := mesh_get_edge_ptr_unsafe(iter.mesh^, iter.current)
	iter.current = e.prev

	return e, prev, true
}

mesh_create_vertex_edge_iterator :: proc(mesh: ^Mesh, vertex: Vertex_Index) -> Vertex_Edge_Iterator {
	v := mesh_get_vertex_unsafe(mesh^, vertex)
	edge := mesh_get_edge_unsafe(mesh^, v.edge)

	return {
		mesh = mesh,
		start = v.edge,
		current = v.edge,
		step = 0,
	}
}

mesh_vertex_incomming_edge_iter :: proc(iter: ^Vertex_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	e := mesh_get_edge_ptr_unsafe(iter.mesh^, iter.current)
	prev := iter.current
	iter.current = mesh_get_edge_unsafe(iter.mesh^, e.next).opposite

	return e, prev, true
}

mesh_vertex_outgoing_edge_iter :: proc(iter: ^Vertex_Edge_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1

	e := mesh_get_edge_ptr_unsafe(iter.mesh^, iter.current)
	iter.current = mesh_get_edge_next_unsafe(iter.mesh^, iter.current).opposite

	return mesh_get_edge_ptr_unsafe(iter.mesh^, e.opposite), e.opposite, true
}

mesh_create_triangle_emitter_iter :: proc(mesh: ^Mesh) -> Triangle_Emitter_Iterator {
	face_index := mesh.faces.active[0]
	face := mesh_get_face_unsafe(mesh^, face_index)
	normal := mesh_calculate_face_normal(mesh, face_index)
	return {
		mesh = mesh,
		face = face_index,
		start = face.edge,
		edge = mesh_get_edge_unsafe(mesh^, face.edge).next,
		normal = normal
	}
}

mesh_triangle_emitter_indexed_flat_iter :: proc(iter: ^Triangle_Emitter_Iterator) -> (count: i32, positions: [3]Vec3f32, normal: Vec3f32, indices: [3]i32, ok: bool) {
	if mesh_get_edge_unsafe(iter.mesh^, iter.edge).next == iter.start { // loop till < n-1
		if iter.face_step < i32(len(iter.mesh.faces.active) - 1) {
			iter.walk_step = 0
			iter.face_step += 1
			iter.face = iter.mesh.faces.active[iter.face_step]
			iter.normal = mesh_calculate_face_normal(iter.mesh, iter.face)
			iter.start =  mesh_get_face_unsafe(iter.mesh^, iter.face).edge
			iter.edge = mesh_get_edge_unsafe(iter.mesh^, iter.start).next
			iter.vertex_base = iter.vertex_step
		} else {
			return 0, 0, 0, 0, false
		}
	}

	edge := mesh_get_edge_unsafe(iter.mesh^, iter.edge)

	first := mesh_get_edge_target_unsafe(iter.mesh^, iter.start) // 0
	n := mesh_get_edge_target_unsafe(iter.mesh^, iter.edge) // n
	n_next := mesh_get_edge_target_unsafe(iter.mesh^, edge.next) // n + 1

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
		v := mesh_get_vertex_unsafe(mesh^, vertex)

		incomming_index := v.edge
		incomming := mesh_get_edge_ptr_unsafe(mesh^, incomming_index)
		outgoing_index := incomming.next
		outgoing := mesh_get_edge_ptr_unsafe(mesh^, outgoing_index)

		incomming_op := mesh_get_edge_ptr_unsafe(mesh^, incomming.opposite)
		outgoing_op := mesh_get_edge_ptr_unsafe(mesh^, outgoing.opposite)

		source, target := incomming_op.vertex, outgoing_op.vertex

		outgoing.prev = incomming.prev
		outgoing_op.next = incomming_op.next
		outgoing_op.vertex = incomming_op.vertex

		mesh_get_edge_ptr_unsafe(mesh^, incomming.prev).next = outgoing_index
		mesh_get_edge_ptr_unsafe(mesh^, outgoing_op.next).prev = outgoing.opposite

		incomming_op_vertex := mesh_get_vertex_ptr_unsafe(mesh^, incomming_op.vertex)
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
		vert := mesh_get_vertex_ptr_unsafe(iter.mesh^, e.vertex)
		if vert.edge == i {
			vert.edge = mesh_get_edge_next_unsafe(iter.mesh^, i).opposite
		}

		append(&outgoing, i)
	}

	for i := 0; i < len(outgoing); i += 1 {
		u := outgoing[i]
		v := outgoing[(i + 1) % len(outgoing)]

		u_e := mesh_get_edge_unsafe(mesh^, u)
		v_e := mesh_get_edge_unsafe(mesh^, v)

		if u_e.next == v || v_e.next == u { continue } // Skip adjacent vertices. This split function handles that but this is here to avoid logging the warnings

		mesh_split_face(mesh, v_e.face, u_e.vertex, v_e.vertex)
	}


	face_edge := mesh_get_edge_opposite_unsafe(mesh^, iter.start).next
	face := mesh_alloc_face(mesh, {-1})

	iter = mesh_create_vertex_edge_iterator(mesh, vertex)
	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
		v := mesh_get_vertex_ptr_unsafe(mesh^, e.vertex)
		if v.edge == i {
			v.edge = mesh_get_edge_next_unsafe(iter.mesh^, i).opposite
		}

		next := mesh_get_edge_next_ptr_unsafe(iter.mesh^, i)
		op := mesh_get_edge_ptr_unsafe(iter.mesh^, e.opposite)
		op_prev := mesh_get_edge_prev_ptr_unsafe(iter.mesh^, e.opposite)

		next.prev = op.prev
		op_prev.next = e.next
		delete_key(&mesh.lookup, Lookup_Pair{vertex, e.vertex})
		delete_key(&mesh.lookup, Lookup_Pair{e.vertex, vertex})
		mesh_free_half_edge(mesh, i)
		mesh_free_half_edge(mesh, e.opposite)
		mesh_free_face(mesh, e.face)
	}


	v := mesh_get_vertex_unsafe(mesh^, vertex)
	start := mesh_get_edge_opposite_unsafe(mesh^, v.edge).next
	curr := start

	for {
		edge := mesh_get_edge_ptr_unsafe(mesh^, curr)
		edge.face = face
		curr = edge.next
		if curr == start { break }
	}

	mesh_get_face_ptr_unsafe(mesh^, face).edge = face_edge
    mesh_free_vertex(mesh, vertex)

    return face
}

mesh_dissolve_half_edge :: proc(mesh: ^Mesh, edge: Half_Edge_Index) -> (kept_face: Face_Index) {
	if edge < 0 {return -1}

	// Wire edge.incomming-edge to the edge.opposite.outgoing edge and vice versa
	// Wire edge.outgoing-edge to the edge.opposite.incomming-edge and vice versa

	// This is blender's edge dissolve with the "dissolve vertex" option unselected

	e := mesh_get_edge_ptr_unsafe(mesh^, edge)
	e_op := mesh_get_edge_ptr_unsafe(mesh^, e.opposite)

	e_prev_index := e.prev
	e_op_prev_index := e_op.prev

	e_next_index := e.next
	e_op_next_index := e_op.next

	mesh_get_edge_ptr_unsafe(mesh^, e.next).prev = e_op_prev_index
	mesh_get_edge_ptr_unsafe(mesh^, e.prev).next = e_op_next_index
	mesh_get_edge_ptr_unsafe(mesh^, e_op.next).prev = e_prev_index
	mesh_get_edge_ptr_unsafe(mesh^, e_op.prev).next = e_next_index

	selected_edge_index := edge

	target_index := e.vertex
	source_index := e_op.vertex
	target, source := mesh_get_vertex_ptr_unsafe(mesh^, target_index), mesh_get_vertex_unsafe(mesh^, source_index)

	if target.edge == edge {
		target.edge = e_op.prev
	}

	if source.edge == e.opposite {
		source.edge = e.prev
	}

	if e.face == -1 {
		selected_edge_index = e.opposite
	}

	selected_edge := mesh_get_edge_unsafe(mesh^, selected_edge_index)
	to_be_deleted_index := selected_edge.opposite
	to_be_deleted := mesh_get_edge_unsafe(mesh^, selected_edge.opposite)

	if selected_edge.face != -1 {
		f := mesh_get_face_ptr_unsafe(mesh^, selected_edge.face)
		if f.edge == selected_edge_index {
			f.edge = selected_edge.next
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
		op := mesh_get_edge_unsafe(mesh^, e.opposite)
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

		curr_vert := mesh_get_vertex_ptr_unsafe(mesh^, curr_vert_index)
		next_vert := mesh_get_vertex_ptr_unsafe(mesh^, next_vert_index)

		lookup_pair := Lookup_Pair{curr_vert_index, next_vert_index}
		lookup_pair_op := Lookup_Pair{next_vert_index, curr_vert_index}

		op_index := mesh.lookup[lookup_pair_op] or_else -1
		mesh.lookup[lookup_pair] = curr_edge_index

		if op_index >= 0 {
			mesh_get_edge_ptr_unsafe(mesh^, op_index).opposite = curr_edge_index
		}

		if next_vert.edge < 0 {
			next_vert.edge = curr_edge_index
		}

		next_edge_index := Half_Edge_Index(-1)
		if i < n - 1 { next_edge_index = mesh_alloc_half_edge(mesh, {})
		} else { next_edge_index = first_edge_index }

		mesh_get_edge_ptr_unsafe(mesh^, curr_edge_index)^ = Half_Edge {
			face     = face_index,
			next     = next_edge_index,
			prev     = prev_edge_index,
			vertex   = next_vert_index,
			opposite = op_index,
		}

		prev_edge_index = curr_edge_index
		curr_edge_index = next_edge_index
	}

	mesh_get_edge_ptr_unsafe(mesh^, first_edge_index).prev = prev_edge_index
	return face_index
}

mesh_split_edges_all :: proc(mesh: ^Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.edges.active), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.edges.active), temp_alloc)
	copy(prev_edges[:], mesh.edges.active[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge(mesh, i, factor)
			lookup[mesh_get_edge_unsafe(mesh^, i).opposite] = {}
		}
	}
}

mesh_split_edges_twice_all :: proc(mesh: ^Mesh, factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	// Todo : Make a split varient for creating N-splits
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.edges.active), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.edges.active), temp_alloc)
	copy(prev_edges[:], mesh.edges.active[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge_twice(mesh, i, factor)
			lookup[mesh_get_edge_unsafe(mesh^, i).opposite] = {}
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

	e := mesh_get_edge_ptr_unsafe(mesh^, half_edge_index)
	e_op := mesh_get_edge_ptr_unsafe(mesh^, e.opposite)

	e_index := e_op.opposite
	e_op_index := e.opposite
	e_next_index := e.next
	e_op_next_index := e_op.next
	e_prev_index := e.prev
	e_op_prev_index := e_op.prev

	target_index, source_index := e.vertex, e_op.vertex
	target, source := mesh_get_vertex_ptr_unsafe(mesh^, target_index), mesh_get_vertex_ptr_unsafe(mesh^, source_index)

	new_vertex := mesh_get_vertex_ptr_unsafe(mesh^, new_vertex_index)
	new_vertex1 := mesh_get_vertex_ptr_unsafe(mesh^, new_vertex1_index)
	mid_point := (source.position + target.position) / 2
	new_vertex.position = target.position + (mid_point - target.position) * factor
	new_vertex1.position = source.position + (mid_point - source.position) * factor

	new_e := mesh_get_edge_ptr_unsafe(mesh^, new_e_index)
	new_e_op := mesh_get_edge_ptr_unsafe(mesh^, new_e_op_index)
	new_e1 := mesh_get_edge_ptr_unsafe(mesh^, new_e1_index)
	new_e1_op := mesh_get_edge_ptr_unsafe(mesh^, new_e1_op_index)

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

	mesh_get_edge_ptr_unsafe(mesh^, e_next_index).prev = new_e_index
	mesh_get_edge_ptr_unsafe(mesh^, e_op.prev).next = new_e_op_index

	mesh_get_edge_ptr_unsafe(mesh^, e_prev_index).next = new_e1_index
	mesh_get_edge_ptr_unsafe(mesh^, e_op_next_index).prev = new_e1_op_index

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

	e := mesh_get_edge_ptr_unsafe(mesh^, half_edge_index)
	e_op := mesh_get_edge_ptr_unsafe(mesh^, e.opposite)

	e_index := e_op.opposite
	e_op_index := e.opposite
	e_next_index := e.next
	e_op_next_index := e_op.next

	e_next := mesh_get_edge_ptr_unsafe(mesh^, e_next_index)
	e_op_next := mesh_get_edge_ptr_unsafe(mesh^, e_op_next_index)

	target_index, source_index := e.vertex, e_op.vertex
	target, source := mesh_get_vertex_ptr_unsafe(mesh^, target_index), mesh_get_vertex_ptr_unsafe(mesh^, source_index)
	new_vertex := mesh_get_vertex_ptr_unsafe(mesh^, new_vertex_index)
	new_vertex.position = source.position + (target.position - source.position) * factor

	new_e := mesh_get_edge_ptr_unsafe(mesh^, new_e_index)
	new_e_op := mesh_get_edge_ptr_unsafe(mesh^, new_e_op_index)

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

	e_next.prev = new_e_index
	mesh_get_edge_ptr_unsafe(mesh^, e_op.prev).next = new_e_op_index

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
	face := mesh_get_face_ptr_unsafe(mesh^, face_index)

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

	outgoing_a_index = mesh_get_edge_unsafe(mesh^, incomming_a_index).next
	outgoing_b_index = mesh_get_edge_unsafe(mesh^, incomming_b_index).next

	if mesh_get_edge_unsafe(mesh^, outgoing_a_index).next == outgoing_b_index || mesh_get_edge_unsafe(mesh^, outgoing_b_index).next == outgoing_a_index {
		log.warnf("Two adjacent vertices cannot be used to split a face. Vertex a : %i, Vertex b : %i", a_index, b_index)
		return true
	}

	a_to_b_index := mesh_alloc_half_edge(mesh, {})
	b_to_a_index := mesh_alloc_half_edge(mesh, {})
	a_to_b := mesh_get_edge_ptr_unsafe(mesh^, a_to_b_index)
	b_to_a := mesh_get_edge_ptr_unsafe(mesh^, b_to_a_index)

	incomming_a, incomming_b := mesh_get_edge_ptr_unsafe(mesh^, incomming_a_index), mesh_get_edge_ptr_unsafe(mesh^, incomming_b_index)
	outgoing_a, outgoing_b := mesh_get_edge_ptr_unsafe(mesh^, outgoing_a_index), mesh_get_edge_ptr_unsafe(mesh^, outgoing_b_index)

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
	total_non_boundary_edges := len(mesh.edges.active)
	for i in mesh.edges.active {
		e := mesh_get_edge_unsafe(mesh^, i)
		if e.opposite != -1 {continue}

		op := Half_Edge {
			face     = -1,
			opposite = Half_Edge_Index(i),
			vertex   = mesh_get_edge_unsafe(mesh^, e.prev).vertex,
			next     = -1,
			prev     = -1,
		}

		op_index := mesh_alloc_half_edge(mesh, op)
		mesh_get_edge_ptr_unsafe(mesh^, i).opposite = op_index
		mesh.lookup[Lookup_Pair{op.vertex, e.vertex}] = op_index
	}

	for i in mesh.edges.active[total_non_boundary_edges:] {
		b_index := Half_Edge_Index(i)
		b := mesh_get_edge_ptr_unsafe(mesh^, b_index)

		curr := b.opposite
		for {
			curr = mesh_get_edge_prev_unsafe(mesh^, curr).opposite
			c := mesh_get_edge_ptr_unsafe(mesh^, curr)
			if c.face == -1 {
				b.next = curr
				c.prev = b_index
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
        centroid += mesh_get_vertex_unsafe(mesh^, e.vertex).position
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

		edge_to_centroid_index := mesh_alloc_half_edge(mesh, {})
		edge_from_centroid_index := mesh_alloc_half_edge(mesh, {})
		new_face := mesh_alloc_face(mesh, {edge_index})

		mesh_get_vertex_ptr_unsafe(mesh^, centroid_vertex_index).edge = edge_to_centroid_index

		edge_ptr := mesh_get_edge_ptr_unsafe(mesh^, edge_index)

		target_index := edge_ptr.vertex
		source_index := mesh_get_edge_unsafe(mesh^, edge_ptr.opposite).vertex

		edge_to_centroid := Half_Edge{face = new_face, next = edge_from_centroid_index, prev = edge_index, vertex = centroid_vertex_index}
		edge_from_centroid := Half_Edge{face = new_face, prev = edge_to_centroid_index, next = edge_index, vertex = source_index}

		edge_to_centroid_op := mesh.lookup[Lookup_Pair{centroid_vertex_index, target_index}] or_else -1
		edge_from_centroid_op := mesh.lookup[Lookup_Pair{source_index, centroid_vertex_index}] or_else -1

		edge_to_centroid.opposite = edge_to_centroid_op
		edge_from_centroid.opposite = edge_from_centroid_op

		if edge_from_centroid_op != -1 {
			mesh_get_edge_ptr_unsafe(mesh^, edge_from_centroid_op).opposite = edge_from_centroid_index
		}

		if edge_to_centroid_op != -1 {
			mesh_get_edge_ptr_unsafe(mesh^, edge_to_centroid_op).opposite = edge_to_centroid_index
		}

		edge_ptr.face = new_face
		edge_ptr.next = edge_to_centroid_index
		edge_ptr.prev = edge_from_centroid_index

		mesh_get_face_ptr_unsafe(mesh^, new_face).edge = edge_to_centroid_index
		mesh_get_edge_ptr_unsafe(mesh^, edge_to_centroid_index)^ = edge_to_centroid
		mesh_get_edge_ptr_unsafe(mesh^, edge_from_centroid_index)^ = edge_from_centroid
		mesh.lookup[Lookup_Pair{target_index, centroid_vertex_index}] = edge_to_centroid_index
		mesh.lookup[Lookup_Pair{centroid_vertex_index, source_index}] = edge_from_centroid_index
	}

	mesh_free_face(mesh, face)
	return centroid_vertex_index
}

mesh_triangulate_face_from_vertex :: proc(mesh: ^Mesh, face: Face_Index, vertex: Vertex_Index) {}

mesh_calculate_face_normal :: proc(mesh: ^Mesh, face: Face_Index) -> Vec3f32 { // Newells algo to find normal
	normal := Vec3f32{}
	iter := mesh_create_face_edge_iterator(mesh, face)
	for e_c in mesh_face_edge_forward_iter(&iter) {
		v_c := mesh_get_vertex_unsafe(mesh^, e_c.vertex).position
		v_n := mesh_get_vertex_unsafe(mesh^, mesh_get_edge_unsafe(mesh^, e_c.next).vertex).position

		normal.x += (v_n.y - v_c.y) * (v_n.z + v_c.z)
		normal.y += (v_n.z - v_c.z) * (v_n.x + v_c.x)
		normal.z += (v_n.x - v_c.x) * (v_n.y + v_c.y)
	}
	return linalg.normalize0(normal)
}

mesh_normalize_onto_sphere :: proc(mesh: ^Mesh) {
	length := f32(0)
	for v in mesh.verts.active {
		length = max(length, linalg.length(mesh_get_vertex_unsafe(mesh^, v).position))
	}

	for v in mesh.verts.active {
		mesh_get_vertex_ptr_unsafe(mesh^, v).position /= length
	}
}
