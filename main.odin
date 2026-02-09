package main

import "core:slice"
import "core:math/linalg"
import "base:intrinsics"
import "base:runtime"
import "core:log"

import rl "vendor:raylib"

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
	lookup:     map[Lookup_Pair]Half_Edge_Index, // Source - Target
}

Face_Walk_Iterator :: struct {
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

Convey_Operation :: enum {
	Ambo,
	Bevel,
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
	Zip,
}

main :: proc() {
	context.logger = log.create_console_logger(.Debug, {.Procedure, .Level, .Line, .Terminal_Color})
	defer log.destroy_console_logger(context.logger)

	platonic_solids := mesh_generate_all_platonic_solids()
	defer meshes_destroy(..slice.enumerated_array(&platonic_solids))

	// archimedean_solids := [Archimedean_Solid]Mesh{}
	archimedean_solids := mesh_generate_all_archimedean_solids()
	defer meshes_destroy(..slice.enumerated_array(&archimedean_solids))

	// catalan_solids := [Catalan_Solid]Mesh{}
	catalan_solids := mesh_generate_all_catalan_solids()
	defer meshes_destroy(..slice.enumerated_array(&catalan_solids))

	selected_type : union #no_nil {Platonic_Solid, Archimedean_Solid, Catalan_Solid} = Platonic_Solid.Cube

	rl.InitWindow(1000, 1000, "Mesh")
	defer rl.CloseWindow()

	rl.SetTargetFPS(60)

	camera := rl.Camera3D{}
	camera.position = {1, 0, 0}
	camera.target = {0, 0, 0}
	camera.fovy = 45
	camera.projection = .PERSPECTIVE
	camera.up = {0, 1, 0}

	draw_debug := false

	mesh_split_edges_twice_all(&platonic_solids[.Cube])
	mesh_convay_kis(&platonic_solids[.Cube], 0)
	mesh_validate(platonic_solids[.Cube])

	for !rl.WindowShouldClose() {

		if rl.IsKeyPressed(.E) {draw_debug = ~draw_debug}
		if rl.IsKeyPressed(.U) { selected_type = .Cube	}
		if rl.IsKeyPressed(.I) { selected_type = .Truncated_Tetrahedron }
		if rl.IsKeyPressed(.P) { selected_type = .Triakis_Tetrahedron }

		if rl.IsKeyPressed(.LEFT) {
			switch &v in selected_type {
			case Platonic_Solid:	v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case Archimedean_Solid:	v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case Catalan_Solid:		v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		if rl.IsKeyPressed(.RIGHT) {
			switch &v in selected_type {
			case Platonic_Solid:	v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case Archimedean_Solid:	v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case Catalan_Solid:		v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		start := rl.Vector2{40, 40}
		size := rl.Vector2{100, 50}
		margin := rl.Vector2{5, 5}
		for o in Convey_Operation {
			if o == .Gyro || o == .Snub {
				rl.DrawRectangleV(start, size, {85, 50, 55, 255})
			} else {
				rl.DrawRectangleV(start, size, {45, 50, 55, 255})
			}
			if rl.CheckCollisionPointRec(rl.GetMousePosition(), {start.x, start.y, size.x, size.y}) {
				rl.DrawRectangleLinesEx({start.x, start.y, size.x, size.y}, 2, {80, 85, 100, 255})
				if rl.IsMouseButtonPressed(.LEFT) {
					switch v in selected_type {
					case Catalan_Solid:		mesh_convay_operation(&catalan_solids[v], o)
					case Platonic_Solid: 	mesh_convay_operation(&platonic_solids[v], o, kis_height = 0.5)
					case Archimedean_Solid: mesh_convay_operation(&archimedean_solids[v], o)
					}
				}
			}

			text := rl.TextFormat("%v", o)
			text_size := rl.MeasureTextEx(rl.GetFontDefault(), text, 20, 1)
			rl.DrawTextEx(rl.GetFontDefault(), text, start - text_size / 2 + size / 2, 20, 1, rl.WHITE)
			start.y += size.y + margin.y
		}

		{
			t := rl.TextFormat("%v", selected_type)
			rl.DrawText(t, rl.GetScreenWidth() / 2 - rl.MeasureText(t, 20) / 2, 10, 20, rl.WHITE)
		}

		{
			rl.DrawRectangleV(start, size, {45, 50, 55, 255})
			if rl.CheckCollisionPointRec(rl.GetMousePosition(), {start.x, start.y, size.x, size.y}) {
				rl.DrawRectangleLinesEx({start.x, start.y, size.x, size.y}, 2, {80, 85, 100, 255})
				if rl.IsMouseButtonPressed(.LEFT) {
					switch v in selected_type {
					case Catalan_Solid:		mesh_destroy(catalan_solids[v]);		catalan_solids[v] = mesh_generate_catalan_solid(v)
					case Platonic_Solid: 	mesh_destroy(platonic_solids[v]);		platonic_solids[v] = mesh_generate_platonic_solid(v)
					case Archimedean_Solid: mesh_destroy(archimedean_solids[v]);	archimedean_solids[v] = mesh_generate_archimedean_solid(v)
					}
				}
			}
			text : cstring = "Reset"
			text_size := rl.MeasureTextEx(rl.GetFontDefault(), text, 20, 1)
			rl.DrawTextEx(rl.GetFontDefault(), text, start - text_size / 2 + size / 2, 20, 1, rl.WHITE)
			start.y += size.y + margin.y
		}

		rl.UpdateCamera(&camera, .ORBITAL)
		rl.BeginDrawing()
		rl.ClearBackground(rl.BLACK)
		switch v in selected_type {
		case Catalan_Solid:		draw_mesh_edges(catalan_solids[v], camera, draw_debug)
		case Platonic_Solid:	draw_mesh_edges(platonic_solids[v], camera, draw_debug)
		case Archimedean_Solid: draw_mesh_edges(archimedean_solids[v], camera, draw_debug)
		}
		rl.EndDrawing()
	}
}

draw_mesh_edges :: proc(mesh: Mesh, camera: rl.Camera3D, draw_debug: bool) {
	OFFSET_MAGNITUDE: f32 = 0.025 if draw_debug else 0
	SHORTEN_AMOUNT: f32 = 0.025 if draw_debug else 0
	ARROW_SIZE: f32 = 0.01 if draw_debug else 0

	mesh := mesh

	rl.BeginMode3D(camera)

	for i in mesh.active_edges {
		e := mesh.edges[i]
		target_v := mesh.verts[e.vertex].position
		source_v := mesh.verts[mesh.edges[e.opposite].vertex].position

		dir := rl.Vector3Normalize(target_v - source_v)

		up := rl.Vector3{0, 1, 0}
		if abs(rl.Vector3DotProduct(dir, up)) > 0.9 {
			up = {1, 0, 0}
		}

		perp := rl.Vector3Normalize(rl.Vector3CrossProduct(dir, up))

		offset := perp * OFFSET_MAGNITUDE

		p1 := source_v + offset + (dir * SHORTEN_AMOUNT)
		p2 := target_v + offset - (dir * SHORTEN_AMOUNT)

		col_v := 255 * ((target_v + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}

		col = e.face == -1 ? rl.WHITE : col

		rl.DrawLine3D(p1, p2, col)
	}
	rl.EndMode3D()

	if !draw_debug { return }
	for i in mesh.active_faces {
		iter := mesh_create_face_walk_iterator(&mesh, i)
		centroid := Vec3f32{}
		for e, i in mesh_face_walk_iter(&iter) {
			centroid += mesh.verts[e.vertex].position
		}
		centroid /= f32(iter.step)

		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("F: %i, E: %i", i, mesh.faces[i].edge), rl.GetWorldToScreen(centroid, camera), 20, 1, rl.YELLOW)
	}

	for i in mesh.active_edges {
		e := mesh.edges[i]
		target, source := e.vertex, mesh.edges[e.opposite].vertex
		s, t := mesh.verts[source], mesh.verts[target]
		dir := linalg.normalize0(s.position - t.position)
		col_v := 255 * ((t.position + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}
		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("E: %i N: %i P: %i V: %i", i, e.next, e.prev, e.vertex), rl.GetWorldToScreen(((s.position + t.position + dir / 2) / 2), camera), 20, 2, col)
	}

	for i in mesh.active_verts {
		v := mesh_get_vertex(mesh, i)
		col_v := 255 * ((v.position + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}
		rl.BeginMode3D(camera)
		rl.DrawSphere(v.position, 0.025, col)
		rl.EndMode3D()
		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("V: %i E: %i", i, v.edge), rl.GetWorldToScreen(v.position, camera), 20, 2, rl.ORANGE)
	}

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

mesh_create_face_walk_iterator :: proc(mesh: ^Mesh, face: Face_Index) -> Face_Walk_Iterator {
	return {
		current = mesh.faces[face].edge,
		start = mesh.faces[face].edge,
		mesh = mesh,
		step = 0,
	}
}

mesh_face_walk_iter :: proc(iter: ^Face_Walk_Iterator) -> (^Half_Edge, Half_Edge_Index, bool) {
	if iter.step > 0 && iter.current == iter.start {
		return nil, -1, false
	}

	iter.step += 1
	prev := iter.current
	e := &iter.mesh.edges[iter.current]
	iter.current = e.next

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
			unordered_remove(&mesh.active_edges, j)
			return
		}
	}
}

mesh_free_vertex :: proc(mesh: ^Mesh, vertex: Vertex_Index) {
	if vertex < 0 { return }
	append(&mesh.free_verts, vertex)
	for i, j in mesh.active_verts {
		if i == vertex {
			unordered_remove(&mesh.active_verts, j)
			return
		}
	}
}

mesh_free_face :: proc(mesh: ^Mesh, face: Face_Index) {
	if face < 0 { return }
	append(&mesh.free_faces, face)
	for i, j in mesh.active_faces {
		if i == face {
			unordered_remove(&mesh.active_faces, j)
			return
		}
	}
}

mesh_dissolve_vertex_face_split :: proc(mesh: ^Mesh, vertex: Vertex_Index, temp_alloc := context.temp_allocator) -> (new_face: Face_Index) {
	if vertex < 0 { return -1 }

	iter := mesh_create_vertex_edge_iterator(mesh, vertex)

	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
        if e.face == -1 {
            log.error("Cannot dissolve boundary vertex")
            return -1
        }
	}

	if iter.step < 3 {
		v := mesh.verts[vertex]

		incomming_index := v.edge
		outgoing_index := mesh.edges[incomming_index].next

		incomming := &mesh.edges[incomming_index]
		outgoing := &mesh.edges[outgoing_index]
		incomming_op := &mesh.edges[incomming.opposite]
		outgoing_op := &mesh.edges[outgoing.opposite]

		s, t := incomming_op.vertex, outgoing_op.vertex

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
		delete_key(&mesh.lookup, Lookup_Pair{s, vertex})
		delete_key(&mesh.lookup, Lookup_Pair{vertex, s})
		delete_key(&mesh.lookup, Lookup_Pair{t, vertex})
		delete_key(&mesh.lookup, Lookup_Pair{vertex, t})

		mesh.lookup[Lookup_Pair{s, t}] = outgoing_index
		mesh.lookup[Lookup_Pair{t, s}] = outgoing.opposite
		return
	}

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

		mesh_split_face(mesh, mesh.edges[v].face, mesh.edges[u].vertex, mesh.edges[v].vertex)
	}

	face := mesh_dissolve_vertex_rim(mesh, vertex)

	return face
}

mesh_dissolve_vertex_rim :: proc(mesh: ^Mesh, vertex: Vertex_Index) -> (new_face: Face_Index) {
    if vertex < 0 { return -1 }

	iter := mesh_create_vertex_edge_iterator(mesh, vertex)

	for e, i in mesh_vertex_outgoing_edge_iter(&iter) {
        if e.face == -1 {
            log.error("Cannot dissolve boundary vertex")
            return -1
        }
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

mesh_dissolve_half_edge_pair :: proc(mesh: ^Mesh, edge: Half_Edge_Index) -> (kept_face: Face_Index) {
	if edge < 0 {return -1}

	e := &mesh.edges[edge]
	e_op := &mesh.edges[e.opposite]

	e_prev := e.prev
	e_op_prev := e_op.prev

	e_next := e.next
	e_op_next := e_op.next

	mesh.edges[e.next].prev = e_op_prev
	mesh.edges[e.prev].next = e_op_next

	mesh.edges[e_op.next].prev = e_prev
	mesh.edges[e_op.prev].next = e_next

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

		iter := mesh_create_face_walk_iterator(mesh, selected_edge.face)
		for face_e in mesh_face_walk_iter(&iter) {
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

mesh_remove_face :: proc(mesh: ^Mesh, face: Face_Index) {
	iter := mesh_create_face_walk_iterator(mesh, face)
	for e in mesh_face_walk_iter(&iter) {
		e.face = -1
	}

	mesh_free_face(mesh, face)
}

mesh_merge_face :: proc(mesh: ^Mesh, face_a: Face_Index, face_b: Face_Index) -> (kept: Face_Index) {
	if face_a < 0 || face_b < 0 {
		return -1
	}

	common_edge := Half_Edge_Index(-1)

	iter := mesh_create_face_walk_iterator(mesh, face_a)
	for e, i in mesh_face_walk_iter(&iter) {
		op := mesh.edges[e.opposite]
		if op.face == face_b {
			common_edge = i
			break
		}
	}

	return mesh_dissolve_half_edge_pair(mesh, common_edge)
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

mesh_split_edges_all :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.active_edges), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.active_edges), temp_alloc)
	copy(prev_edges[:], mesh.active_edges[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge(mesh, i)
			lookup[mesh.edges[i].opposite] = {}
		}
	}
}

mesh_split_edges_twice_all :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	prev_edges := make([dynamic]Half_Edge_Index, len(mesh.active_edges), temp_alloc)
	lookup := make(map[Half_Edge_Index]struct{}, len(mesh.active_edges), temp_alloc)
	copy(prev_edges[:], mesh.active_edges[:])

	for i in prev_edges {
		_, done := lookup[i]
		if !done {
			mesh_split_edge_twice(mesh, i)
			lookup[mesh.edges[i].opposite] = {}
		}
	}
}

mesh_split_edge_twice :: proc(mesh: ^Mesh, half_edge_index: Half_Edge_Index) {
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
	new_vertex.position = source.position + (target.position - source.position) * 0.666666666
	new_vertex1.position = source.position + (target.position - source.position) * 0.333333333

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

mesh_split_edge :: proc(mesh: ^Mesh, half_edge_index: Half_Edge_Index) -> (Half_Edge_Index) {
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
	new_vertex.position = (source.position + target.position) / 2

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
		log.errorf("Same vertex cannot be used to split a face")
		return true
	}

	{ 	// Validate inputs first
		found_a, found_b := false, false

		iter := mesh_create_face_walk_iterator(mesh, face_index)
		for e, i in mesh_face_walk_iter(&iter) {
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
			log.errorf("Vertices do not belong to provided face. Face Index : %i, Face : %v, Vertex a : %i, Vertex b : %i", face_index, face, a_index, b_index)
			return true
		}
	}

	outgoing_a_index = mesh.edges[incomming_a_index].next
	outgoing_b_index = mesh.edges[incomming_b_index].next

	if mesh.edges[outgoing_a_index].next == outgoing_b_index || mesh.edges[outgoing_b_index].next == outgoing_a_index {
		log.errorf("Two adjacent vertices cannot be used to split a face. Vertex a : %i, Vertex b : %i", a_index, b_index)
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
		iter := mesh_create_face_walk_iterator(mesh, face_index)
		for e in mesh_face_walk_iter(&iter) {
			e.face = face_index
		}

		iter = mesh_create_face_walk_iterator(mesh, new_face_index)
		for e in mesh_face_walk_iter(&iter) {
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

mesh_whirl_face :: proc(mesh: ^Mesh, face: Face_Index, twist_factor := f32(0.33), inset_factor := f32(0.33), temp_alloc := context.temp_allocator) {
    original_edges := make([dynamic]Half_Edge_Index, temp_alloc)
    original_verts := make([dynamic]Vertex_Index, temp_alloc)

    centroid := Vec3f32{}
    iter := mesh_create_face_walk_iterator(mesh, face)
    for e, e_idx in mesh_face_walk_iter(&iter) {
        v_idx := mesh.edges[e_idx].vertex
        centroid += mesh.verts[v_idx].position
        append(&original_edges, e_idx)
        append(&original_verts, v_idx)
    }

    count := len(original_verts)
    if count < 3 do return
    centroid /= f32(count)

    // 2. Generate Inner Vertices (The "Whirled" Ring)
    inner_verts := make([dynamic]Vertex_Index, count, temp_alloc)

    for i in 0..<count {
        // Current edge goes from source -> target (verts[i])
        // To twist, we usually blend between Source and Target.
        // Let's grab the actual Source Vertex for clarity:
        prev_idx := (i - 1 + count) % count
        source_v := original_verts[prev_idx] // Source of edge leading to verts[i]
        target_v := original_verts[i]

        p_source := mesh.verts[source_v].position
        p_target := mesh.verts[target_v].position

        // Twist: Move along the edge
        p_twist  := linalg.lerp(p_source, p_target, twist_factor)

        // Inset: Move towards centroid
        p_final  := linalg.lerp(p_twist, centroid, inset_factor)

        inner_verts[i] = mesh_add_vertex(mesh, p_final)
    }

    // 3. Topology Construction
    // We are replacing 1 face with:
    // - 1 Inner Face (connects all inner_verts)
    // - 2 * N Triangles (the skirt)

    // Pre-allocate inner edges for the center face
    inner_face_edges := make([dynamic]Half_Edge_Index, count, temp_alloc)
    inner_face_idx := mesh_alloc_face(mesh, {}) // We fill data later

    // Create the Inner Loop connectivity
    for i in 0..<count {
        curr_v := inner_verts[i]
        next_v := inner_verts[(i + 1) % count]

        he := mesh_alloc_half_edge(mesh, Half_Edge{
            vertex = next_v,
            face   = inner_face_idx,
        })
        inner_face_edges[i] = he

        // Set vertex edge pointer to an outgoing edge (this one is on the inner ring)
        mesh.verts[curr_v].edge = he
    }

    // Link Inner Loop (Next/Prev)
    for i in 0..<count {
        curr := inner_face_edges[i]
        prev := inner_face_edges[(i - 1 + count) % count]
        next := inner_face_edges[(i + 1) % count]

        mesh.edges[curr].prev = prev
        mesh.edges[curr].next = next

        // Setup lookup for the inner loop
        v_curr := inner_verts[i]
        v_next := inner_verts[(i + 1) % count]
        mesh.lookup[Lookup_Pair{v_curr, v_next}] = curr
    }
    mesh.faces[inner_face_idx].edge = inner_face_edges[0]

    // 4. Create the "Skirt" (Triangulate the gap)
    // The gap is between original_verts and inner_verts.
    // We treat the original edges as the boundary.
    // For each side i, we have a quad: Source_i -> Target_i -> Inner_i -> Inner_{i-1}
    // We split this quad into 2 triangles.
    // Choice of split determines chirality (CW/CCW visual).

    for i in 0..<count {
        // Indices relative to the loop
        prev_i   := (i - 1 + count) % count

        // Vertices
        v_source := original_verts[prev_i]
        v_target := original_verts[i]
        v_inner  := inner_verts[i]
        v_inner_prev := inner_verts[prev_i]

        // The original boundary edge (we reuse it, but change its face)
        edge_boundary := original_edges[i]

        // We will form 2 triangles:
        // T1: Source -> Target -> Inner ( v_source, v_target, v_inner )
        // T2: Source -> Inner -> Inner_Prev ( v_source, v_inner, v_inner_prev )

        // -- Triangle 1 --
        t1 := mesh_alloc_face(mesh, {})

        // Edges for T1
        // E1: Original Boundary (Source -> Target)
        // E2: Target -> Inner (New)
        // E3: Inner -> Source (New, Diagonal)

        e_target_inner := mesh_alloc_half_edge(mesh, Half_Edge{ vertex = v_inner, face = t1 })
        e_inner_source := mesh_alloc_half_edge(mesh, Half_Edge{ vertex = v_source, face = t1 })

        // Link T1
        mesh.edges[edge_boundary].face = t1
        mesh.edges[edge_boundary].next = e_target_inner
        mesh.edges[edge_boundary].prev = e_inner_source

        mesh.edges[e_target_inner].prev = edge_boundary
        mesh.edges[e_target_inner].next = e_inner_source

        mesh.edges[e_inner_source].prev = e_target_inner
        mesh.edges[e_inner_source].next = edge_boundary

        mesh.faces[t1].edge = edge_boundary

        // -- Triangle 2 --
        t2 := mesh_alloc_face(mesh, {})

        // Edges for T2
        // E1: Source -> Inner (Opposite of T1's E3)
        // E2: Inner -> Inner_Prev (Opposite of Inner Face Edge [prev_i])
        // E3: Inner_Prev -> Source (New)

        e_source_inner      := mesh_alloc_half_edge(mesh, Half_Edge{ vertex = v_inner, face = t2 })
        e_inner_inner_prev  := mesh_alloc_half_edge(mesh, Half_Edge{ vertex = v_inner_prev, face = t2 })
        e_inner_prev_source := mesh_alloc_half_edge(mesh, Half_Edge{ vertex = v_source, face = t2 })

        // Link T2
        mesh.edges[e_source_inner].next = e_inner_inner_prev
        mesh.edges[e_source_inner].prev = e_inner_prev_source

        mesh.edges[e_inner_inner_prev].next = e_inner_prev_source
        mesh.edges[e_inner_inner_prev].prev = e_source_inner

        mesh.edges[e_inner_prev_source].next = e_source_inner
        mesh.edges[e_inner_prev_source].prev = e_inner_inner_prev

        mesh.faces[t2].edge = e_source_inner

        // -- Update Lookups & Opposites --

        // 1. Diagonal (Source <-> Inner)
        mesh.edges[e_inner_source].opposite = e_source_inner
        mesh.edges[e_source_inner].opposite = e_inner_source

        // 2. Target -> Inner (connects to T2 of the NEXT edge iteration? No, connects to T2 of NEXT iteration's 'Inner_Prev -> Source')
        // Actually, let's just populate lookup map and let a repair pass handle or do it manually:
        mesh.lookup[Lookup_Pair{v_target, v_inner}] = e_target_inner
        mesh.lookup[Lookup_Pair{v_inner, v_source}] = e_inner_source

        mesh.lookup[Lookup_Pair{v_source, v_inner}] = e_source_inner
        mesh.lookup[Lookup_Pair{v_inner, v_inner_prev}] = e_inner_inner_prev
        mesh.lookup[Lookup_Pair{v_inner_prev, v_source}] = e_inner_prev_source

        // Link Inner_Inner_Prev to the Inner Face edge
        inner_edge_idx := inner_face_edges[prev_i]
        mesh.edges[e_inner_inner_prev].opposite = inner_edge_idx
        mesh.edges[inner_edge_idx].opposite = e_inner_inner_prev

        // Link e_target_inner to e_inner_prev_source of the NEXT iteration
        // (Target of current == Source of next)
        // (Inner of current == Inner_Prev of next)
        // This is tricky to do in one pass without lookups.
        // Reliance on `mesh_stitch_opposites` later is recommended,
        // OR:
        op_candidate := mesh.lookup[Lookup_Pair{v_inner, v_target}] or_else -1
        if op_candidate != -1 {
            mesh.edges[e_target_inner].opposite = op_candidate
            mesh.edges[op_candidate].opposite = e_target_inner
        }

        op_candidate_2 := mesh.lookup[Lookup_Pair{v_source, v_inner_prev}] or_else -1
        if op_candidate_2 != -1 {
            mesh.edges[e_inner_prev_source].opposite = op_candidate_2
            mesh.edges[op_candidate_2].opposite = e_inner_prev_source
        }
    }

    // Free the original face (it was replaced by the inner face and skirt)
    mesh_free_face(mesh, face)
}

mesh_triangulate_face_from_centroid :: proc (mesh: ^Mesh, face: Face_Index, height := f32(0), temp_alloc := context.temp_allocator) -> Vertex_Index {
    collected_edges := make([dynamic]Half_Edge_Index, temp_alloc)
	centroid := Vec3f32{}
    iter := mesh_create_face_walk_iterator(mesh, face)
    for e, i in mesh_face_walk_iter(&iter) {
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
	iter := mesh_create_face_walk_iterator(mesh, face)
	for e_c in mesh_face_walk_iter(&iter) {
		v_c := mesh.verts[e_c.vertex].position
		v_n := mesh.verts[mesh.edges[e_c.next].vertex].position

		normal.x += (v_n.y - v_c.y) * (v_n.z + v_c.z)
		normal.y += (v_n.z - v_c.z) * (v_n.x + v_c.x)
		normal.z += (v_n.x - v_c.x) * (v_n.y + v_c.y)
	}
	return linalg.normalize0(normal)
}

mesh_normalize :: proc(mesh: ^Mesh) {
	length := f32(0)
	for v in mesh.active_verts {
		length = max(0, linalg.length(mesh.verts[v].position))
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
		iter := mesh_create_face_walk_iterator(mesh, f)
		for e, i in mesh_face_walk_iter(&iter) {
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
        mesh_add_face(&dual, verts[:])
		clear(&verts)
	}

	mesh_destroy(mesh^)
	mesh^ = dual
}

mesh_convay_kis :: proc(mesh: ^Mesh, height := f32(1), temp_alloc := context.temp_allocator) {
	faces := make([dynamic]Face_Index, len(mesh.active_faces), temp_alloc)
	copy(faces[:], mesh.active_faces[:])

	for f in faces {
		mesh_triangulate_face_from_centroid(mesh, f, height, temp_alloc)
	}
}

mesh_convay_ambo :: proc (mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	verts := make([dynamic]Vertex_Index, len(mesh.active_verts), temp_alloc)
	copy(verts[:], mesh.active_verts[:])

	mesh_split_edges_all(mesh, temp_alloc)

	for i in verts {
		mesh_dissolve_vertex_face_split(mesh, i, temp_alloc)
	}
}

mesh_convay_truncate :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	verts :=  make([dynamic]Vertex_Index, len(mesh.active_verts), temp_alloc)
	copy(verts[:], mesh.active_verts[:])

	mesh_split_edges_twice_all(mesh, temp_alloc)

	for v in verts {
		mesh_dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

mesh_convay_snub :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
}

mesh_convay_gyro :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
}

mesh_convay_bevel :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_ambo(mesh, temp_alloc)
	mesh_convay_truncate(mesh, temp_alloc)
}

mesh_convay_expand :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_ambo(mesh, temp_alloc)
	mesh_convay_ambo(mesh, temp_alloc)
}

mesh_convay_join :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_ambo(mesh, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_ortho :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_expand(mesh, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_meta :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_bevel(mesh, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_needle :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_dual(mesh, temp_alloc)
	mesh_convay_kis(mesh, 0, temp_alloc)
}

mesh_convay_zip :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	mesh_convay_kis(mesh, 0, temp_alloc)
	mesh_convay_dual(mesh, temp_alloc)
}

mesh_convay_operation :: proc(mesh: ^Mesh, operation: Convey_Operation, temp_alloc := context.temp_allocator, kis_height := f32(0)) {
	switch operation {
		case .Kis:		mesh_convay_kis(mesh, kis_height, temp_alloc)
		case .Zip:		mesh_convay_zip(mesh, temp_alloc)
		case .Ambo:		mesh_convay_ambo(mesh, temp_alloc)
		case .Dual:		mesh_convay_dual(mesh, temp_alloc)
		case .Snub:		mesh_convay_snub(mesh, temp_alloc)
		case .Join:		mesh_convay_join(mesh, temp_alloc)
		case .Meta:		mesh_convay_meta(mesh, temp_alloc)
		case .Gyro:		mesh_convay_gyro(mesh, temp_alloc)
		case .Ortho:	mesh_convay_ortho(mesh, temp_alloc)
		case .Bevel:	mesh_convay_bevel(mesh, temp_alloc)
		case .Needle:	mesh_convay_needle(mesh, temp_alloc)
		case .Expand:	mesh_convay_expand(mesh, temp_alloc)
		case .Truncate:	mesh_convay_truncate(mesh, temp_alloc)
	}
}


// Polygon generation
mesh_generate_tetrahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, {1, 1, 1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1})
	mesh_add_faces(&mesh, {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2})
	mesh_normalize(&mesh)
	return mesh
}

mesh_generate_cube :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1})
	mesh_add_faces(&mesh, {3, 2, 1, 0}, {3, 0, 5, 4}, {4, 5, 6, 7}, {7, 6, 1, 2}, {6, 5, 0, 1}, {2, 3, 4, 7})
	mesh_normalize(&mesh)
	return mesh
}

mesh_generate_octahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	mesh_add_vertices(&mesh, { 1,  0,  0},{-1,  0,  0},{ 0,  1,  0},{ 0, -1,  0},{ 0,  0,  1},{ 0,  0, -1})
	mesh_add_faces(&mesh,{0, 4, 2},{2, 4, 1},{1, 4, 3},{3, 4, 0},{2, 5, 0},{1, 5, 2},{3, 5, 1},{0, 5, 3})
	mesh_normalize(&mesh)
	return mesh
}

mesh_generate_icosahedron :: proc(allocator := context.allocator) -> Mesh {
	mesh := mesh_create(allocator)
	phi :: 1.618033988749894
	mesh_add_vertices(&mesh,
		{-1,  phi,  0}, { 1,  phi,  0}, {-1, -phi,  0}, { 1, -phi,  0},
		{ 0, -1,  phi}, { 0,  1,  phi}, { 0, -1, -phi}, { 0,  1, -phi},
		{ phi,  0, -1}, { phi,  0,  1}, {-phi,  0, -1}, {-phi,  0,  1},
	)
	mesh_add_faces(&mesh,
		{0, 11, 5},  {0, 5, 1},   {0, 1, 7},   {0, 7, 10},  {0, 10, 11},
		{1, 5, 9},   {5, 11, 4},  {11, 10, 2}, {10, 7, 6},  {7, 1, 8},
		{3, 9, 4},   {3, 4, 2},   {3, 2, 6},   {3, 6, 8},   {3, 8, 9},
		{4, 9, 5},   {2, 4, 11},  {6, 2, 10},  {8, 6, 7},   {9, 8, 1},
	)
	mesh_normalize(&mesh)
	return mesh
}

mesh_generate_dodecahedron :: proc(allocator := context.allocator) -> Mesh {
    mesh := mesh_create(allocator)
    phi     :: 1.618033988749894
    inv_phi :: 0.618033988749894

    mesh_add_vertices(&mesh,
        { 1,  1,  1}, { 1,  1, -1}, { 1, -1,  1}, { 1, -1, -1}, // 0-3: Cube vertices
        {-1,  1,  1}, {-1,  1, -1}, {-1, -1,  1}, {-1, -1, -1}, // 4-7: Cube vertices
        { 0, inv_phi,  phi}, { 0, inv_phi, -phi}, { 0, -inv_phi,  phi}, { 0, -inv_phi, -phi}, // 8-11
        { inv_phi,  phi, 0}, { inv_phi, -phi, 0}, {-inv_phi,  phi, 0}, {-inv_phi, -phi, 0}, // 12-15
        { phi, 0,  inv_phi}, { phi, 0, -inv_phi}, {-phi, 0,  inv_phi}, {-phi, 0, -inv_phi}, // 16-19
    )

    mesh_add_faces(&mesh,
        {0, 16, 2, 10, 8},   {0, 8, 4, 14, 12},   {0, 12, 1, 17, 16},
        {3, 17, 1, 9, 11},   {3, 11, 7, 15, 13},  {3, 13, 2, 16, 17},
        {5, 9, 1, 12, 14},   {5, 14, 4, 18, 19},  {5, 19, 7, 11, 9},
        {6, 10, 2, 13, 15},  {6, 15, 7, 19, 18},  {6, 18, 4, 8, 10},
    )
    mesh_normalize(&mesh)
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
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_octahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_rhombicuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	mesh_convay_ambo(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_snub_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	return {}
}

mesh_generate_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_icosahedron(allocator)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_rhombicosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	mesh_convay_ambo(&mesh, temp_alloc)
	return mesh
}

mesh_generate_truncated_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ambo(&mesh, temp_alloc)
	mesh_convay_truncate(&mesh, temp_alloc)
	return mesh
}

mesh_generate_snub_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	return {}
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
	return mesh
}

mesh_generate_rhombic_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_join(&mesh, temp_alloc)
	return mesh
}

mesh_generate_triakis_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_octahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_TRI_OCTAHEDRON_KIS_HEIGHT, temp_alloc)
	return mesh
}

mesh_generate_tetrakis_hexahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_kis(&mesh, CATALAN_TETRA_HEXAHEDRON_KIS_HEIGHT, temp_alloc)
	return mesh
}

mesh_generate_deltoidal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_ortho(&mesh, temp_alloc)
	return mesh
}

mesh_generate_disdyakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_cube(allocator)
	mesh_convay_meta(&mesh, temp_alloc)
	return mesh
}

mesh_generate_pentagonal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	return {} // gyro of cube
}

mesh_generate_rhombic_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_join(&mesh)
	return mesh
}

mesh_generate_triakis_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_icosahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_TRI_ICOSAHEDRON_KIS_HEIGHT, temp_alloc)
	return mesh
}

mesh_generate_pentakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_kis(&mesh, CATALAN_PENTA_DODECAHEDRON_KIS_HEIGHT, temp_alloc)
	return mesh
}

mesh_generate_deltoidal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_ortho(&mesh, temp_alloc)
	return mesh
}

mesh_generate_disdyakis_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	mesh := mesh_generate_dodecahedron(allocator)
	mesh_convay_meta(&mesh, temp_alloc)
	return mesh
}

mesh_generate_pentagonal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> Mesh {
	return {} // gyro of dodecahedron
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
