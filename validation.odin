package mesh

import "core:log"

validate :: proc(mesh: Mesh, temp_alloc := context.temp_allocator, caller_expression := #caller_expression) -> (invalid: bool) {
	log.info("--- Validation : ", caller_expression)
	base_ok := validate_base_constraints(mesh)
	indices_ok := validate_indicies(mesh)
	opposites_ok := validate_opposites(mesh)
	link_ok := validate_links(mesh)
	face_loops_ok := validate_face_loops(mesh, temp_alloc)
	vertex_references_ok := validate_vertex_references(mesh)
	edge_uniqueness_ok := validate_global_uniqueness(mesh, temp_alloc)

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

validate_base_constraints :: proc (mesh: Mesh) -> (invalide: bool) {
	for i in mesh.edges.active {
		_, ok_next := get_edge_next(mesh, i)
		_, ok_prev := get_edge_prev(mesh, i)
		_, ok_opposite := get_edge_opposite(mesh, i)

		if !(ok_next || ok_prev || ok_opposite) {
			log.errorf("Invalid Next: %v, Invalid Prev: %v, Invalid Opposite: %v. Edge: %v", !ok_next, !ok_prev, !ok_opposite, get_edge_unsafe(mesh, i))
			return true
		}
	}

	for i in mesh.faces.active {
		face := get_face_unsafe(mesh, i)
		_, ok := get_edge(mesh, face.edge)

		if !ok {
			log.errorf("Face does not reference an existing edge. Face: %v", face)
			return true
		}
	}

	for i in mesh.verts.active {
		vert := get_vertex_unsafe(mesh, i)
		_, ok := get_edge(mesh, vert.edge)

		if !ok {
			log.errorf("Vertex does not reference an existing edge. Vertex Index: %i, Vertex: %v", i, vert)
			return true
		}
	}

	return false
}

validate_opposites :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.edges.active {
		e := get_edge_unsafe(mesh, i)
		e_op := get_edge_unsafe(mesh, e.opposite)

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

validate_indicies :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.edges.active {
		e := get_edge_unsafe(mesh, i)
		if e.opposite < 0 || int(e.opposite) >= len(mesh.edges.all) {
			log.errorf("Invalide opposite half-edge index. Got %i, expected between [0-%i)", e.opposite, len(mesh.edges.all))
			return true
		}

		if e.next < 0 || int(e.next) >= len(mesh.edges.all) {
			log.errorf("Invalide next half-edge index. Got %i, expected between [0-%i)", e.next, len(mesh.edges.all))
			return true
		}

		if e.face < -1 || int(e.face) >= len(mesh.faces.all) {
			log.errorf("Invalide face index. Got %i, expected between [-1-%i)", e.face, len(mesh.faces.all))
			return true
		}

		if e.vertex < 0 || int(e.vertex) >= len(mesh.verts.all) {
			log.errorf("Invalide vertex index. Got %i, expected between [0-%i)", e.vertex, len(mesh.verts.all))
			return true
		}

		if e.prev < 0 || int(e.prev) >= len(mesh.edges.all) {
			log.errorf("Invalid prev half-edge index. Got %i, expected between [0-%i)", e.prev, len(mesh.edges.all))
		}
	}
	return false
}

validate_links :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.edges.active {
		e := get_edge_unsafe(mesh, i)

		next := get_edge_next_unsafe(mesh, i)
		prev := get_edge_prev_unsafe(mesh, i)

		if next.prev != i {
			log.errorf("Linkage Invariance Violated: edges[%i].next.prev is %i, expected %i", i, next.prev, i)
			return true
		}

		if prev.next != i {
			log.errorf("Linkage Invariance Violated: edges[%i].prev.next is %i, expected %i", i, prev.next, i)
			return true
		}
		e_op := get_edge_unsafe(mesh, e.opposite)
		if e_op.vertex != prev.vertex {
			log.errorf( "current.source == prev.target property violated. Current, Prev : %v, %v", e, prev)
			return true
		}
	}
	return false
}

validate_face_loops :: proc( mesh: Mesh, allocator := context.temp_allocator) -> ( invalid: bool, ) {
	lookup := make(map[Half_Edge_Index]struct{}, allocator)
	defer delete(lookup)

	for i in mesh.faces.active {
		f := mesh.faces.all[i]
		defer clear(&lookup)
		start := f.edge


		if get_edge_unsafe(mesh, start).face != i { 	// The edge must point/reference current face
			log.errorf("Half-edge doesn't point to current face being walked. Face Index : %i, Face : %v, Edge : %v", i, f, get_edge_unsafe(mesh, start))
			return true
		}

		curr := start
		walk_distance := int(0)

		lookup[start] = {}

		for {
			c := get_edge_unsafe(mesh, curr)
			curr = c.next

			if c.face != i {
				log.errorf( "Half-edge doesn't point to current face being walked. Face Index : %i, Face : %v, Edge : %v", i, f, c)
				return true
			}

			walk_distance += 1

			if curr == start {
				break
			}

			_, exists := lookup[curr]

			if exists { 	// Check locally that no repition happens in the face walk
				log.errorf("Half-edges repeated during face walk. Edge : %v", c)
				return true
			}

			lookup[curr] = {}
		}

		if walk_distance < 2 { 	// degenerate face
			log.errorf( "Degenerate Face consisting of only single half-edge. Face Index : %i, Face : %v", i, f)
			return true
		}
	}
	return false
}

validate_vertex_references :: proc(mesh: Mesh) -> (invalid: bool) {
	for i in mesh.verts.active {
		v := mesh.verts.all[i]
		c := get_edge_unsafe(mesh, v.edge)
		if c.vertex != i {
			log.errorf("Half-edge vertex references does not have the vertex as its target. Vertex Index : %i, Vertex : %v, Edge : %v", i, v, c)
			return true
		}
	}
	return false
}

validate_global_uniqueness :: proc(mesh: Mesh, allocator := context.temp_allocator) -> (failed: bool) {
	visited := make(map[Half_Edge_Index]struct{}, allocator)
	defer delete(visited)

	for i in mesh.edges.active {
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

			c := get_edge_unsafe(mesh, walk_idx)
			walk_idx = c.next

			if walk_idx == start_idx do break

			if len(visited) > len(mesh.edges.active) {
				log.errorf("Face walk exceeded edge count. Circularity broken.")
				return true
			}
		}
	}

	if len(visited) != len(mesh.edges.active) {
		log.errorf("Orphaned edges detected: %i / %i reached", len(visited), len(mesh.edges.active))
		return true
	}

	return false
}
