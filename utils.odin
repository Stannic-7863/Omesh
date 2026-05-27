package mesh

// TODO: handle allocator errors?

free_list_create :: proc(free_list: ^Free_List($type, $type_index), allocator := context.allocator) {
	free_list.active = make([dynamic]type_index, allocator)
	free_list.free = make([dynamic]type_index, allocator)
	free_list.all = make([dynamic]type, allocator)
}

free_list_destroy :: proc(free_list: Free_List($type, $type_index)) {
	delete(free_list.all)
	delete(free_list.free)
	delete(free_list.active)
}

free_list_add :: proc(free_list: ^Free_List($type, $type_index), value: type) -> type_index {
	index, ok := pop_safe(&free_list.free)
	defer append(&free_list.active, index)
	if ok {
		free_list.all[index] = value
		return index
	}
	index = type_index(len(free_list.all))
	append(&free_list.all, value)
	return index
}

free_list_remove :: proc(free_list: ^Free_List($type, $type_index), index: type_index) {
	if index < 0 { return }
	found := false
	for i, j in free_list.active {
		if i == index {
			unordered_remove(&free_list.active, j)
			found = true
			break
		}
	}
	if found { append(&free_list.free, index) }
}

free_list_is_item_free :: proc(free_list: Free_List($type, $type_index), index: type_index) -> (ok: bool) {
	if index < 0 { return false }
	for i in free_list.free { if index == i { return false } }
	return true
}

free_list_get_item :: proc(free_list: Free_List($type, $type_index), index: type_index) -> (item: type, ok: bool) {
	free_list_is_item_free(free_list, index) or_return
	return free_list.all[index], true
}

free_list_get_item_unsafe :: proc(free_list: Free_List($type, $type_index), index: type_index) -> (item: type) {
	return free_list.all[index]
}

free_list_get_item_ptr :: proc(free_list: Free_List($type, $type_index), index: type_index) -> (item: ^type, ok: bool) {
	free_list_is_item_free(free_list, index) or_return
	return &free_list.all[index], true
}

free_list_get_item_ptr_unsafe :: proc(free_list: Free_List($type, $type_index), index: type_index) -> (item: ^type) {
	return &free_list.all[index]
}

get_face_unsafe :: proc (mesh: Mesh, index: Face_Index) -> Face {
	return free_list_get_item_unsafe(mesh.faces, index)
}

get_vertex_unsafe :: proc (mesh: Mesh, index: Vertex_Index) -> Vertex {
	return free_list_get_item_unsafe(mesh.verts, index)
}

get_edge_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return free_list_get_item_unsafe(mesh.edges, index)
}

get_face_ptr_unsafe :: proc (mesh: Mesh, index: Face_Index) -> ^Face {
	return free_list_get_item_ptr_unsafe(mesh.faces, index)
}

get_vertex_ptr_unsafe :: proc (mesh: Mesh, index: Vertex_Index) -> ^Vertex {
	return free_list_get_item_ptr_unsafe(mesh.verts, index)
}

get_edge_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return free_list_get_item_ptr_unsafe(mesh.edges, index)
}

get_edge_next_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return get_edge_unsafe(mesh, get_edge_unsafe(mesh, index).next)
}

get_edge_prev_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return get_edge_unsafe(mesh, get_edge_unsafe(mesh, index).prev)
}

get_edge_opposite_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Half_Edge {
	return get_edge_unsafe(mesh, get_edge_unsafe(mesh, index).opposite)
}

get_edge_next_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return get_edge_ptr_unsafe(mesh, get_edge_unsafe(mesh, index).next)
}

get_edge_prev_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return get_edge_ptr_unsafe(mesh, get_edge_unsafe(mesh, index).prev)
}

get_edge_opposite_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Half_Edge {
	return get_edge_ptr_unsafe(mesh, get_edge_unsafe(mesh, index).opposite)
}

get_edge_source_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Vertex {
	return get_vertex_unsafe(mesh, get_edge_prev_unsafe(mesh, index).vertex)
}

get_edge_target_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> Vertex {
	return get_vertex_unsafe(mesh, get_edge_unsafe(mesh, index).vertex)
}

get_edge_source_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Vertex {
	return get_vertex_ptr_unsafe(mesh, get_edge_prev_unsafe(mesh, index).vertex)
}

get_edge_target_ptr_unsafe :: proc (mesh: Mesh, index: Half_Edge_Index) -> ^Vertex {
	return get_vertex_ptr_unsafe(mesh, get_edge_unsafe(mesh, index).vertex)
}

get_face :: proc (mesh: Mesh, index: Face_Index) -> (face: Face, ok: bool) {
	return free_list_get_item(mesh.faces, index)
}

get_vertex :: proc (mesh: Mesh, index: Vertex_Index) -> (vertex: Vertex, ok: bool) {
	return free_list_get_item(mesh.verts, index)
}

get_edge :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	return free_list_get_item(mesh.edges, index)
}

get_face_ptr :: proc (mesh: Mesh, index: Face_Index) -> (face: ^Face, ok: bool) {
	return free_list_get_item_ptr(mesh.faces, index)
}

get_vertex_ptr :: proc (mesh: Mesh, index: Vertex_Index) -> (vertex: ^Vertex, ok: bool) {
	return free_list_get_item_ptr(mesh.verts, index)
}

get_edge_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	return free_list_get_item_ptr(mesh.edges, index)
}

get_edge_next :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := get_edge(mesh, index) or_return
	return get_edge(mesh, e.next)
}

get_edge_prev :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := get_edge(mesh, index) or_return
	return get_edge(mesh, e.prev)
}

get_edge_opposite :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: Half_Edge, ok: bool) {
	e := get_edge(mesh, index) or_return
	return get_edge(mesh, e.opposite)
}

get_edge_next_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := get_edge_ptr(mesh, index) or_return
	return get_edge_ptr(mesh, e.next)
}

get_edge_prev_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := get_edge_ptr(mesh, index) or_return
	return get_edge_ptr(mesh, e.prev)
}

get_edge_opposite_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (edge: ^Half_Edge, ok: bool) {
	e := get_edge_ptr(mesh, index) or_return
	return get_edge_ptr(mesh, e.opposite)
}

get_edge_source :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: Vertex, ok: bool) {
	edge := get_edge(mesh, index) or_return
	prev := get_edge(mesh, edge.prev) or_return
	return get_vertex(mesh, prev.vertex)
}

get_edge_target :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: Vertex, ok: bool) {
	edge := get_edge(mesh, index) or_return
	return get_vertex(mesh, edge.vertex)
}

get_edge_source_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: ^Vertex, ok: bool) {
	edge := get_edge(mesh, index) or_return
	prev := get_edge(mesh, edge.prev) or_return
	return get_vertex_ptr(mesh, prev.vertex)
}

get_edge_target_ptr :: proc (mesh: Mesh, index: Half_Edge_Index) -> (vertex: ^Vertex, ok: bool) {
	edge := get_edge(mesh, index) or_return
	return get_vertex_ptr(mesh, edge.vertex)
}

alloc_face :: proc(mesh: ^Mesh, face: Face) -> Face_Index {
	return free_list_add(&mesh.faces, face)
}

alloc_half_edge :: proc(mesh: ^Mesh, half_edge: Half_Edge) -> Half_Edge_Index {
	return free_list_add(&mesh.edges, half_edge)
}

alloc_vertex :: proc(mesh: ^Mesh, vertex: Vertex) -> Vertex_Index {
	return free_list_add(&mesh.verts, vertex)
}

free_half_edge :: proc(mesh: ^Mesh, edge: Half_Edge_Index) {
	free_list_remove(&mesh.edges, edge)
}

free_vertex :: proc(mesh: ^Mesh, vertex: Vertex_Index) {
	free_list_remove(&mesh.verts, vertex)
}

free_face :: proc(mesh: ^Mesh, face: Face_Index) {
	free_list_remove(&mesh.faces, face)
}
