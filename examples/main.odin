package example

import "core:math/linalg"
import "core:slice"
import "core:log"

import rl "vendor:raylib"
import m "../"

main :: proc() {
	context.logger = log.create_console_logger(.Debug, {.Procedure, .Level, .Line, .Terminal_Color})
	defer log.destroy_console_logger(context.logger)

	platonic_solids := m.mesh_generate_all_platonic_solids()
	defer m.meshes_destroy(..slice.enumerated_array(&platonic_solids))

	// archimedean_solids := [Archimedean_Solid]Mesh{}
	archimedean_solids := m.mesh_generate_all_archimedean_solids()
	defer m.meshes_destroy(..slice.enumerated_array(&archimedean_solids))

	// catalan_solids := [Catalan_Solid]Mesh{}
	catalan_solids := m.mesh_generate_all_catalan_solids()
	defer m.meshes_destroy(..slice.enumerated_array(&catalan_solids))

	selected_type : union #no_nil {m.Platonic_Solid, m.Archimedean_Solid, m.Catalan_Solid} = .Tetrahedron

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

	m.mesh_split_edges_twice_all(&platonic_solids[.Cube])

	m.mesh_validate(platonic_solids[.Cube])

	for !rl.WindowShouldClose() {
		if rl.IsKeyPressed(.E) { draw_debug = ~draw_debug }
		if rl.IsKeyPressed(.U) { selected_type = .Tetrahedron }
		if rl.IsKeyPressed(.I) { selected_type = .Truncated_Tetrahedron }
		if rl.IsKeyPressed(.P) { selected_type = .Triakis_Tetrahedron }

		if rl.IsKeyPressed(.LEFT) {
			switch &v in selected_type {
			case m.Platonic_Solid:	v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case m.Archimedean_Solid:	v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case m.Catalan_Solid:		v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		if rl.IsKeyPressed(.RIGHT) {
			switch &v in selected_type {
			case m.Platonic_Solid:	v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case m.Archimedean_Solid:	v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case m.Catalan_Solid:		v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		start := rl.Vector2{70, 40}
		size := rl.Vector2{100, 50}
		margin := rl.Vector2{5, 5}
		for o in m.Convay_Operation {
			text := rl.TextFormat("%v", o)
			text_size := rl.MeasureTextEx(rl.GetFontDefault(), text, 20, 1)
			d_pos := start + {-text_size.x / 2 + size.x / 2 - 8 if text_size.x > size.x else 0, 0}
			d_size := rl.Vector2{max(size.x, text_size.x + 16), size.y}
			rl.DrawRectangleV(d_pos, d_size, {45, 50, 55, 255})
			rl.DrawTextEx(rl.GetFontDefault(), text, start - text_size / 2 + size / 2, 20, 1, rl.WHITE)
			if rl.CheckCollisionPointRec(rl.GetMousePosition(), {d_pos.x, d_pos.y, d_size.x, d_size.y}) {
				rl.DrawRectangleLinesEx({d_pos.x, d_pos.y, d_size.x, d_size.y}, 2, {80, 85, 100, 255})
				if rl.IsMouseButtonPressed(.LEFT) {
					switch v in selected_type {
					case m.Catalan_Solid:		m.mesh_convay_operation(&catalan_solids[v], o)
					case m.Platonic_Solid: 		m.mesh_convay_operation(&platonic_solids[v], o)
					case m.Archimedean_Solid: 	m.mesh_convay_operation(&archimedean_solids[v], o)
					}
				}
			}
			start.y += d_size.y + margin.y
		}

		{
			t := rl.TextFormat("%v", selected_type)
			rl.DrawText(t, rl.GetScreenWidth() / 2 - rl.MeasureText(t, 20) / 2, 10, 20, rl.WHITE)
			t = cstring("")
			switch v in selected_type {
			case m.Catalan_Solid: 		t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(catalan_solids[v].active_verts), len(catalan_solids[v].active_faces), len(catalan_solids[v].active_edges))
			case m.Platonic_Solid: 		t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(platonic_solids[v].active_verts), len(platonic_solids[v].active_faces), len(platonic_solids[v].active_edges))
			case m.Archimedean_Solid: 	t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(archimedean_solids[v].active_verts), len(archimedean_solids[v].active_faces), len(archimedean_solids[v].active_edges))
			}
			rl.DrawText(t, rl.GetScreenWidth() / 2 - rl.MeasureText(t, 20) / 2, 40, 20, rl.WHITE)
		}

		{
			rl.DrawRectangleV(start, size, {45, 50, 55, 255})
			if rl.CheckCollisionPointRec(rl.GetMousePosition(), {start.x, start.y, size.x, size.y}) {
				rl.DrawRectangleLinesEx({start.x, start.y, size.x, size.y}, 2, {80, 85, 100, 255})
				if rl.IsMouseButtonPressed(.LEFT) {
					switch v in selected_type {
					case m.Catalan_Solid:		m.mesh_destroy(catalan_solids[v]);		catalan_solids[v] = m.mesh_generate_catalan_solid(v)
					case m.Platonic_Solid: 		m.mesh_destroy(platonic_solids[v]);		platonic_solids[v] = m.mesh_generate_platonic_solid(v)
					case m.Archimedean_Solid: 	m.mesh_destroy(archimedean_solids[v]);	archimedean_solids[v] = m.mesh_generate_archimedean_solid(v)
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
		case m.Catalan_Solid:		draw_mesh_edges(catalan_solids[v], camera, draw_debug)
		case m.Platonic_Solid:		draw_mesh_edges(platonic_solids[v], camera, draw_debug)
		case m.Archimedean_Solid: 	draw_mesh_edges(archimedean_solids[v], camera, draw_debug)
		}
		rl.EndDrawing()
	}
 }

draw_mesh_edges :: proc(mesh: m.Mesh, camera: rl.Camera3D, draw_debug: bool) {
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
		iter := m.mesh_create_face_edge_iterator(&mesh, i)
		centroid := m.Vec3f32{}
		for e, i in m.mesh_face_edge_forward_iter(&iter) {
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
		v := m.mesh_get_vertex(mesh, i)
		col_v := 255 * ((v.position + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}
		rl.BeginMode3D(camera)
		rl.DrawSphere(v.position, 0.025, col)
		rl.EndMode3D()
		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("V: %i E: %i", i, v.edge), rl.GetWorldToScreen(v.position, camera), 20, 2, rl.ORANGE)
	}
}
