package convay

import "core:slice"
import "core:log"
import m "../"

Operation :: enum {
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

dual :: proc(mesh: ^m.Mesh, temp_alloc := context.temp_allocator) {
	dual := m.create(mesh.edges.all.allocator)

    dual_verts := make([dynamic]m.Vertex_Index, len(mesh.faces.all), temp_alloc)

    for f in mesh.faces.active {
		centroid := m.Vec3f32{}
		iter := m.create_face_edge_iterator(mesh, f)
		for e in m.face_edge_backward_iter(&iter) {
			centroid += m.get_vertex_unsafe(mesh^, e.vertex).position
		}
		centroid /= f32(iter.step)
		dual_verts[f] = m.add_vertex(&dual, centroid)
    }

	verts := make([dynamic]m.Vertex_Index, temp_alloc)
    for v in mesh.verts.active {
        iter := m.create_vertex_edge_iterator(mesh, v)
        for e in m.vertex_outgoing_edge_iter(&iter) {
            append(&verts, dual_verts[e.face])
        }
        slice.reverse(verts[:])
        m.add_face(&dual, verts[:])
		clear(&verts)
    }

	m.destroy(mesh^)
	mesh^ = dual
}

kis :: proc(mesh: ^m.Mesh, kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	faces := make([dynamic]m.Face_Index, len(mesh.faces.active), temp_alloc)
	copy(faces[:], mesh.faces.active[:])

	for f in faces {
		m.triangulate_face_from_centroid(mesh, f, kis_height, temp_alloc)
	}
}

ambo :: proc (mesh: ^m.Mesh, ambo_factor := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]m.Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

	m.split_edges_all(mesh, ambo_factor, temp_alloc)

	for v in verts {
		m.dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}


truncate :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	verts :=  make([dynamic]m.Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

	m.split_edges_twice_all(mesh, truncate_factor, temp_alloc)

	for v in verts {
		m.dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

snub :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	gyro(mesh, truncate_factor, gyro_height, temp_alloc)
	kis(mesh, kis_height, temp_alloc)
}

gyro :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), height := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]m.Vertex_Index, len(mesh.verts.active), temp_alloc)
	faces := make([dynamic]m.Face_Index, len(mesh.faces.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])
	copy(faces[:], mesh.faces.active[:])

	m.split_edges_twice_all(mesh, truncate_factor, temp_alloc)

	for f in faces {
		centroid_vert := m.triangulate_face_from_centroid(mesh, f, height, temp_alloc)

		to_dissolve := make([dynamic]m.Half_Edge_Index, temp_alloc)
		for v in verts {
			iter_v := m.create_vertex_edge_iterator(mesh, v)
			for e, i in m.vertex_outgoing_edge_iter(&iter_v) {
				if e.vertex == centroid_vert {
					append(&to_dissolve, i)
					break
				}
			}
		}

		for e in to_dissolve {
			m.dissolve_half_edge(mesh, m.get_edge_unsafe(mesh^, e).next)
			m.dissolve_half_edge(mesh, e)
		}
	}
}

classical_alternation :: proc(mesh: ^m.Mesh, temp_alloc := context.temp_allocator) {
	lookup := make(map[m.Vertex_Index]bool, temp_alloc)
	queue := make([dynamic]m.Vertex_Index, temp_alloc)

	append(&queue, mesh.verts.active[0])
	lookup[queue[0]] = false

	for len(queue) > 0 {
		v := pop_front(&queue)

		iter := m.create_vertex_edge_iterator(mesh, v)
		for e in m.vertex_outgoing_edge_iter(&iter) {
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
			m.dissolve_vertex_face_split(mesh, k, temp_alloc)
		}
	}
}

classical_snub :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	ambo(mesh, ambo_factor, temp_alloc)
	truncate(mesh, truncate_factor, temp_alloc)
	classical_alternation(mesh, temp_alloc)
}

classical_gyro :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
	dual(mesh, temp_alloc)
}

bevel :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	ambo(mesh, ambo_factor, temp_alloc)
	truncate(mesh, truncate_factor, temp_alloc)
}

expand :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	ambo(mesh, factor, temp_alloc)
	ambo(mesh, factor, temp_alloc)
}

join :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	dual(mesh, temp_alloc)
	ambo(mesh, factor, temp_alloc)
	dual(mesh, temp_alloc)
}

ortho :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	dual(mesh, temp_alloc)
	expand(mesh, factor, temp_alloc)
	dual(mesh, temp_alloc)
}

meta :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	dual(mesh, temp_alloc)
	bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
	dual(mesh, temp_alloc)
}

needle :: proc(mesh: ^m.Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	dual(mesh, temp_alloc)
	kis(mesh, height, temp_alloc)
}

zip :: proc(mesh: ^m.Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	kis(mesh, height, temp_alloc)
	dual(mesh, temp_alloc)
}

operation :: proc(mesh: ^m.Mesh, operation: Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	switch operation {
	case .Kis:						kis(mesh, kis_height, temp_alloc)
	case .Zip:						zip(mesh, kis_height, temp_alloc)
	case .Ambo:						ambo(mesh, ambo_factor, temp_alloc)
	case .Dual:						dual(mesh, temp_alloc)
	case .Snub:						snub(mesh, truncate_factor, gyro_height, kis_height, temp_alloc)
	case .Join:						join(mesh, ambo_factor, temp_alloc)
	case .Meta:						meta(mesh, ambo_factor, truncate_factor, temp_alloc)
	case .Gyro:						gyro(mesh, truncate_factor, gyro_height, temp_alloc)
	case .Ortho:					ortho(mesh, ambo_factor, temp_alloc)
	case .Bevel:					bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
	case .Needle:					needle(mesh, kis_height, temp_alloc)
	case .Expand:					expand(mesh, ambo_factor, temp_alloc)
	case .Truncate:					truncate(mesh, truncate_factor, temp_alloc)
	case .Classical_Snub:			classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
	case .Classical_Gyro:			classical_gyro(mesh, ambo_factor, truncate_factor, temp_alloc)
	case .Classical_Alternation:	classical_alternation(mesh, temp_alloc)
	}
}

operations :: proc(mesh: ^m.Mesh, ops: ..Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	for op in ops {
		operation(mesh, op, temp_alloc, ambo_factor, truncate_factor, gyro_height, kis_height)
	}
}
