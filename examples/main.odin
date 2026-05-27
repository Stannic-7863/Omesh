package example

import "core:math/linalg"
import "core:slice"
import "core:log"

import rl "vendor:raylib"
import m "../"
import c "../convay"
import p "../convay/polygon"

main :: proc() {
	context.logger = log.create_console_logger(.Debug, {.Procedure, .Level, .Line, .Terminal_Color})
	defer log.destroy_console_logger(context.logger)

	platonic_solids := p.generate_all_platonic_solids()
	defer m.destroy(..slice.enumerated_array(&platonic_solids))
	defer m.validate(platonic_solids[.Cube])

	// archimedean_solids := [c.Archimedean_Solid]m.Mesh{}
	archimedean_solids := p.generate_all_archimedean_solids()
	defer m.destroy(..slice.enumerated_array(&archimedean_solids))

	// catalan_solids := [m.Catalan_Solid]m.Mesh{}
	catalan_solids := p.generate_all_catalan_solids()
	defer m.destroy(..slice.enumerated_array(&catalan_solids))

	selected_type : union #no_nil {p.Platonic_Solid, p.Archimedean_Solid, p.Catalan_Solid} = .Tetrahedron

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

	for !rl.WindowShouldClose() {
		if rl.IsKeyPressed(.E) { draw_debug = ~draw_debug }
		if rl.IsKeyPressed(.U) { selected_type = .Tetrahedron }
		if rl.IsKeyPressed(.I) { selected_type = .Truncated_Tetrahedron }
		if rl.IsKeyPressed(.O) { selected_type = .Triakis_Tetrahedron }

		if rl.IsKeyPressed(.LEFT) {
			switch &v in selected_type {
			case p.Platonic_Solid:		v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case p.Archimedean_Solid:	v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			case p.Catalan_Solid:		v = type_of(v) ( (int(v) + 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		if rl.IsKeyPressed(.RIGHT) {
			switch &v in selected_type {
			case p.Platonic_Solid:		v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case p.Archimedean_Solid:	v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			case p.Catalan_Solid:		v = type_of(v) ( (int(v) - 1 + len(type_of(v))) % len(type_of(v)) )
			}
		}

		if rl.IsKeyPressed(.V) {
			switch &v in selected_type {
			case p.Platonic_Solid:		m.validate(platonic_solids[v])
			case p.Archimedean_Solid:  	m.validate(archimedean_solids[v])
			case p.Catalan_Solid:		m.validate(catalan_solids[v])
			}
		}

		start := rl.Vector2{70, 40}
		size := rl.Vector2{100, 50}
		margin := rl.Vector2{5, 5}
		for o in c.Operation {
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
					case p.Catalan_Solid:		c.operation(&catalan_solids[v], o)
					case p.Platonic_Solid: 		c.operation(&platonic_solids[v], o)
					case p.Archimedean_Solid: 	c.operation(&archimedean_solids[v], o)
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
			case p.Catalan_Solid: 		t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(catalan_solids[v].verts.active), len(catalan_solids[v].faces.active), len(catalan_solids[v].edges.active))
			case p.Platonic_Solid: 		t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(platonic_solids[v].verts.active), len(platonic_solids[v].faces.active), len(platonic_solids[v].edges.active))
			case p.Archimedean_Solid: 	t = rl.TextFormat("Verts %v\nFaces %v\nEdges %v", len(archimedean_solids[v].verts.active), len(archimedean_solids[v].faces.active), len(archimedean_solids[v].edges.active))
			}
			rl.DrawText(t, rl.GetScreenWidth() / 2 - rl.MeasureText(t, 20) / 2, 40, 20, rl.WHITE)
		}

		{
			rl.DrawRectangleV(start, size, {45, 50, 55, 255})
			if rl.CheckCollisionPointRec(rl.GetMousePosition(), {start.x, start.y, size.x, size.y}) {
				rl.DrawRectangleLinesEx({start.x, start.y, size.x, size.y}, 2, {80, 85, 100, 255})
				if rl.IsMouseButtonPressed(.LEFT) {
					switch v in selected_type {
					case p.Catalan_Solid:		m.destroy(catalan_solids[v]);		catalan_solids[v] = p.generate_catalan_solid(v)
					case p.Platonic_Solid: 		m.destroy(platonic_solids[v]);		platonic_solids[v] = p.generate_platonic_solid(v)
					case p.Archimedean_Solid: 	m.destroy(archimedean_solids[v]);	archimedean_solids[v] = p.generate_archimedean_solid(v)
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
		case p.Catalan_Solid:		draw_mesh_edges(catalan_solids[v], camera, draw_debug)
		case p.Platonic_Solid:		draw_mesh_edges(platonic_solids[v], camera, draw_debug)
		case p.Archimedean_Solid: 	draw_mesh_edges(archimedean_solids[v], camera, draw_debug)
		}
		rl.EndDrawing()
	}
 }

draw_mesh_edges :: proc(mesh: m.Mesh, camera: rl.Camera3D, draw_debug: bool) {
	OFFSET_MAGNITUDE: f32 = 0.025 if draw_debug else 0
	SHORTEN_AMOUNT: f32 = 0.025 if draw_debug else 0

	mesh := mesh

	rl.BeginMode3D(camera)

	for i in mesh.edges.active {
		e := m.get_edge_unsafe(mesh, i)
		target_v := m.get_vertex_unsafe(mesh, e.vertex).position
		source_v := m.get_vertex_unsafe(mesh, m.get_edge_opposite_unsafe(mesh, i).vertex).position

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
	for i in mesh.faces.active {
		iter := m.create_face_edge_iterator(&mesh, i)
		centroid := m.Vec3f32{}
		for e in m.face_edge_forward_iter(&iter) {
			centroid += m.get_vertex_unsafe(mesh, e.vertex).position
		}
		centroid /= f32(iter.step)
		face := m.get_face_unsafe(mesh, i)

		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("F: %i, E: %i", i, face.edge), rl.GetWorldToScreen(centroid, camera), 20, 1, rl.YELLOW)
	}

	for i in mesh.edges.active {
		e := m.get_edge_unsafe(mesh, i)
		target, source := e.vertex, m.get_edge_opposite_unsafe(mesh, i).vertex
		s, t := m.get_vertex_unsafe(mesh, target), m.get_vertex_unsafe(mesh, source)
		dir := linalg.normalize0(s.position - t.position)
		col_v := 255 * ((t.position + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}
		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("E: %i N: %i P: %i V: %i", i, e.next, e.prev, e.vertex), rl.GetWorldToScreen(((s.position + t.position + dir / 2) / 2), camera), 20, 2, col)
	}

	for i in mesh.verts.active {
		v := m.get_vertex_unsafe(mesh, i)
		col_v := 255 * ((v.position + 3) / 4)
		col := rl.Color{u8(col_v.r), u8(col_v.g), u8(col_v.b), 255}
		rl.BeginMode3D(camera)
		rl.DrawSphere(v.position, 0.025, col)
		rl.EndMode3D()
		rl.DrawTextEx(rl.GetFontDefault(), rl.TextFormat("V: %i E: %i", i, v.edge), rl.GetWorldToScreen(v.position, camera), 20, 2, rl.ORANGE)
	}
}
