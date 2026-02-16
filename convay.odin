package mesh

import "core:log"
import "core:slice"

mesh_convay_dual :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	dual := mesh_create(mesh.edges.all.allocator)

    dual_verts := make([dynamic]Vertex_Index, len(mesh.faces.all), temp_alloc)

    for f in mesh.faces.active {
		centroid := Vec3f32{}
		iter := mesh_create_face_edge_iterator(mesh, f)
		for e, i in mesh_face_edge_backward_iter(&iter) {
			centroid += mesh_get_vertex_unsafe(mesh^, e.vertex).position
		}
		centroid /= f32(iter.step)
		dual_verts[f] = mesh_add_vertex(&dual, centroid)
    }

	verts := make([dynamic]Vertex_Index, temp_alloc)
    for v in mesh.verts.active {
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
	faces := make([dynamic]Face_Index, len(mesh.faces.active), temp_alloc)
	copy(faces[:], mesh.faces.active[:])

	for f in faces {
		mesh_triangulate_face_from_centroid(mesh, f, kis_height, temp_alloc)
	}
}

mesh_convay_ambo :: proc (mesh: ^Mesh, ambo_factor := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

	mesh_split_edges_all(mesh, ambo_factor, temp_alloc)

	for v in verts {
		mesh_dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

mesh_convay_truncate :: proc(mesh: ^Mesh, truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	verts :=  make([dynamic]Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

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
	verts := make([dynamic]Vertex_Index, len(mesh.verts.active), temp_alloc)
	faces := make([dynamic]Face_Index, len(mesh.faces.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])
	copy(faces[:], mesh.faces.active[:])

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
			mesh_dissolve_half_edge(mesh, mesh_get_edge_unsafe(mesh^, e).next)
			mesh_dissolve_half_edge(mesh, e)
		}

		delete(to_dissolve)
	}
}

mesh_convay_classical_alternation :: proc(mesh: ^Mesh, temp_alloc := context.temp_allocator) {
	lookup := make(map[Vertex_Index]bool, temp_alloc)
	queue := make([dynamic]Vertex_Index, temp_alloc)

	append(&queue, mesh.verts.active[0])
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
