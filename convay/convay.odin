package convay

import "core:slice"
import "core:log"
import m "../"

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

convay_dual :: proc(mesh: ^m.Mesh, temp_alloc := context.temp_allocator) {
	dual := m.create(mesh.edges.all.allocator)

    dual_verts := make([dynamic]m.Vertex_Index, len(mesh.faces.all), temp_alloc)

    for f in mesh.faces.active {
		centroid := m.Vec3f32{}
		iter := m.create_face_edge_iterator(mesh, f)
		for e, i in m.face_edge_backward_iter(&iter) {
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

convay_kis :: proc(mesh: ^m.Mesh, kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	faces := make([dynamic]m.Face_Index, len(mesh.faces.active), temp_alloc)
	copy(faces[:], mesh.faces.active[:])

	for f in faces {
		m.triangulate_face_from_centroid(mesh, f, kis_height, temp_alloc)
	}
}

convay_ambo :: proc (mesh: ^m.Mesh, ambo_factor := f32(0.5), temp_alloc := context.temp_allocator) {
	verts := make([dynamic]m.Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

	m.split_edges_all(mesh, ambo_factor, temp_alloc)

	for v in verts {
		m.dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}


convay_truncate :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	verts :=  make([dynamic]m.Vertex_Index, len(mesh.verts.active), temp_alloc)
	copy(verts[:], mesh.verts.active[:])

	m.split_edges_twice_all(mesh, truncate_factor, temp_alloc)

	for v in verts {
		m.dissolve_vertex_face_split(mesh, v, temp_alloc)
	}
}

convay_snub :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_gyro(mesh, truncate_factor, gyro_height, temp_alloc)
	convay_kis(mesh, kis_height, temp_alloc)
}

convay_gyro :: proc(mesh: ^m.Mesh, truncate_factor := f32(2.0/3.0), height := f32(0.5), temp_alloc := context.temp_allocator) {
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

		delete(to_dissolve)
	}
}

convay_classical_alternation :: proc(mesh: ^m.Mesh, temp_alloc := context.temp_allocator) {
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

convay_classical_snub :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	convay_ambo(mesh, ambo_factor, temp_alloc)
	convay_truncate(mesh, truncate_factor, temp_alloc)
	convay_classical_alternation(mesh, temp_alloc)
}

convay_classical_gyro :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	convay_classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
	convay_dual(mesh, temp_alloc)
}

convay_bevel :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	convay_ambo(mesh, ambo_factor, temp_alloc)
	convay_truncate(mesh, truncate_factor, temp_alloc)
}

convay_expand :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_ambo(mesh, factor, temp_alloc)
	convay_ambo(mesh, factor, temp_alloc)
}

convay_join :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_dual(mesh, temp_alloc)
	convay_ambo(mesh, factor, temp_alloc)
	convay_dual(mesh, temp_alloc)
}

convay_ortho :: proc(mesh: ^m.Mesh, factor := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_dual(mesh, temp_alloc)
	convay_expand(mesh, factor, temp_alloc)
	convay_dual(mesh, temp_alloc)
}

convay_meta :: proc(mesh: ^m.Mesh, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), temp_alloc := context.temp_allocator) {
	convay_dual(mesh, temp_alloc)
	convay_bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
	convay_dual(mesh, temp_alloc)
}

convay_needle :: proc(mesh: ^m.Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_dual(mesh, temp_alloc)
	convay_kis(mesh, height, temp_alloc)
}

convay_zip :: proc(mesh: ^m.Mesh, height := f32(0.5), temp_alloc := context.temp_allocator) {
	convay_kis(mesh, height, temp_alloc)
	convay_dual(mesh, temp_alloc)
}

mesh_convay_operation :: proc(mesh: ^m.Mesh, operation: Convay_Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	switch operation {
		case .Kis:						convay_kis(mesh, kis_height, temp_alloc)
		case .Zip:						convay_zip(mesh, kis_height, temp_alloc)
		case .Ambo:						convay_ambo(mesh, ambo_factor, temp_alloc)
		case .Dual:						convay_dual(mesh, temp_alloc)
		case .Snub:						convay_snub(mesh, truncate_factor, gyro_height, kis_height, temp_alloc)
		case .Join:						convay_join(mesh, ambo_factor, temp_alloc)
		case .Meta:						convay_meta(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Gyro:						convay_gyro(mesh, truncate_factor, gyro_height, temp_alloc)
		case .Ortho:					convay_ortho(mesh, ambo_factor, temp_alloc)
		case .Bevel:					convay_bevel(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Needle:					convay_needle(mesh, kis_height, temp_alloc)
		case .Expand:					convay_expand(mesh, ambo_factor, temp_alloc)
		case .Truncate:					convay_truncate(mesh, truncate_factor, temp_alloc)
		case .Classical_Snub:			convay_classical_snub(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Classical_Gyro:			convay_classical_gyro(mesh, ambo_factor, truncate_factor, temp_alloc)
		case .Classical_Alternation:	convay_classical_alternation(mesh, temp_alloc)
	}
}

mesh_convay_operations :: proc(mesh: ^m.Mesh, operations: ..Convay_Operation, temp_alloc := context.temp_allocator, ambo_factor := f32(0.5), truncate_factor := f32(2.0/3.0), gyro_height := f32(0.5), kis_height := f32(0.5)) {
	for operation in operations {
		mesh_convay_operation(mesh, operation, temp_alloc, ambo_factor, truncate_factor, gyro_height, kis_height)
	}
}

generate_tetrahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, {1, 1, 1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1})
	m.add_faces(&mesh, {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2})
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_cube :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1})
	m.add_faces(&mesh, {3, 2, 1, 0}, {3, 0, 5, 4}, {4, 5, 6, 7}, {7, 6, 1, 2}, {6, 5, 0, 1}, {2, 3, 4, 7})
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_octahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, { 1,  0,  0},{-1,  0,  0},{ 0,  1,  0},{ 0, -1,  0},{ 0,  0,  1},{ 0,  0, -1})
	m.add_faces(&mesh,{0, 4, 2},{2, 4, 1},{1, 4, 3},{3, 4, 0},{2, 5, 0},{1, 5, 2},{3, 5, 1},{0, 5, 3})
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_icosahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh,
		{-1,  PHI,  0}, { 1,  PHI,  0}, {-1, -PHI,  0}, { 1, -PHI,  0},
		{ 0, -1,  PHI}, { 0,  1,  PHI}, { 0, -1, -PHI}, { 0,  1, -PHI},
		{ PHI,  0, -1}, { PHI,  0,  1}, {-PHI,  0, -1}, {-PHI,  0,  1},
	)
	m.add_faces(&mesh,
		{5, 11, 0}, {1, 5, 0}, {7, 1, 0}, {10, 7, 0},
		{11, 10, 0}, {9, 5, 1}, {4, 11, 5}, {2, 10, 11},
	 	{6, 7, 10}, {8, 1, 7}, {4, 9, 3}, {2, 4, 3},
		{6, 2, 3}, {8, 6, 3}, {9, 8, 3}, {5, 9, 4},
		{11, 4, 2}, {10, 2, 6}, {7, 6, 8}, {1, 8, 9}
	)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_dodecahedron :: proc(allocator := context.allocator) -> m.Mesh {
    mesh := m.create(allocator)

    m.add_vertices(&mesh,
        { 1,  1,  1}, { 1,  1, -1}, { 1, -1,  1}, { 1, -1, -1}, // 0-3: Cube vertices
        {-1,  1,  1}, {-1,  1, -1}, {-1, -1,  1}, {-1, -1, -1}, // 4-7: Cube vertices
        { 0, INV_PHI,  PHI}, { 0, INV_PHI, -PHI}, { 0, -INV_PHI,  PHI}, { 0, -INV_PHI, -PHI}, // 8-11
        { INV_PHI,  PHI, 0}, { INV_PHI, -PHI, 0}, {-INV_PHI,  PHI, 0}, {-INV_PHI, -PHI, 0}, // 12-15
        { PHI, 0,  INV_PHI}, { PHI, 0, -INV_PHI}, {-PHI, 0,  INV_PHI}, {-PHI, 0, -INV_PHI}, // 16-19
    )

    m.add_faces(&mesh,
        {0, 16, 2, 10, 8},   {0, 8, 4, 14, 12},   {0, 12, 1, 17, 16},
        {3, 17, 1, 9, 11},   {3, 11, 7, 15, 13},  {3, 13, 2, 16, 17},
        {5, 9, 1, 12, 14},   {5, 14, 4, 18, 19},  {5, 19, 7, 11, 9},
        {6, 10, 2, 13, 15},  {6, 15, 7, 19, 18},  {6, 18, 4, 8, 10},
    )
    m.normalize_onto_sphere(&mesh)
    return mesh
}

generate_platonic_solid :: proc(solid : Platonic_Solid, allocator := context.allocator) -> m.Mesh {
	switch solid {
	case .Cube:			return generate_cube(allocator)
	case .Octahedron:	return generate_octahedron(allocator)
	case .Tetrahedron:	return generate_tetrahedron(allocator)
	case .Dodecaheron:	return generate_dodecahedron(allocator)
	case .Icosahedron:	return generate_icosahedron(allocator)
	}
	unreachable()
}

generate_all_platonic_solids :: proc(allocator := context.allocator) -> [Platonic_Solid]m.Mesh {
	solids := [Platonic_Solid]m.Mesh{}
	for &s, t in solids {
		s = generate_platonic_solid(t, allocator)
	}
	return solids
}

generate_truncated_tetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_tetrahedron(allocator)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_octahedron(allocator)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_rhombicuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_snub_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_classical_snub(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_icosahedron(allocator)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_rhombicosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_truncated_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_ambo(&mesh, temp_alloc = temp_alloc)
	convay_truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_snub_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_classical_snub(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_archimedean_solid :: proc (solid : Archimedean_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	#partial switch solid {
	case .Snub_Cube:					return generate_snub_cube(allocator, temp_alloc)
	case .Cuboctahedron:		 		return generate_cuboctahedron(allocator, temp_alloc)
	case .Truncated_Cube:		 		return generate_truncated_cube(allocator, temp_alloc)
	case .Snub_dodecahedron:			return generate_snub_dodecahedron(allocator, temp_alloc)
	case .Icosidodecahedron:			return generate_icosidodecahedron(allocator, temp_alloc)
	case .Rhombicuboctahedron:	 		return generate_rhombicuboctahedron(allocator, temp_alloc)
	case .Truncated_Octahedron:	 		return generate_truncated_octahedron(allocator, temp_alloc)
	case .Truncated_Tetrahedron: 		return generate_truncated_tetrahedron(allocator, temp_alloc)
	case .Truncated_icosahedron:		return generate_truncated_icosahedron(allocator, temp_alloc)
	case .Truncated_dodecahedron:		return generate_truncated_dodecahedron(allocator, temp_alloc)
	case .Rhombicosidodecahedron:		return generate_rhombicosidodecahedron(allocator, temp_alloc)
	case .Truncated_Cuboctahedron:		return generate_truncated_cuboctahedron(allocator, temp_alloc)
	case .Truncated_Icosidodecahedron:	return generate_truncated_icosidodecahedron(allocator, temp_alloc)
	}
	unreachable()
}

generate_all_archimedean_solids :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> [Archimedean_Solid]m.Mesh {
	solids := [Archimedean_Solid]m.Mesh{}
	for &s, t in solids {
		s = mesh_generate_archimedean_solid(t, allocator, temp_alloc)
	}
	return solids
}

generate_triakis_tetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_tetrahedron(allocator)
	convay_kis(&mesh, CATALAN_TRI_TETRAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_rhombic_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_join(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_triakis_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_octahedron(allocator)
	convay_kis(&mesh, CATALAN_TRI_OCTAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_tetrakis_hexahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_kis(&mesh, CATALAN_TETRA_HEXAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_deltoidal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_ortho(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_disdyakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_meta(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_pentagonal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	convay_classical_gyro(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_rhombic_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_join(&mesh)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_triakis_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_icosahedron(allocator)
	convay_kis(&mesh, CATALAN_TRI_ICOSAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_pentakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_kis(&mesh, CATALAN_PENTA_DODECAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_deltoidal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_ortho(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_disdyakis_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_meta(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

generate_pentagonal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	convay_classical_gyro(&mesh, temp_alloc = temp_alloc)
	m.normalize_onto_sphere(&mesh)
	return mesh
}

mesh_generate_catalan_solid :: proc(solid : Catalan_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	switch solid {
	case .Triakis_Octahedron:		 	return generate_triakis_octahedron(allocator, temp_alloc)
	case .Triakis_Tetrahedron:			return generate_triakis_tetrahedron(allocator, temp_alloc)
	case .Triakis_Icosahedron:			return generate_triakis_icosahedron(allocator, temp_alloc)
	case .Tetrakis_Hexahedron:			return generate_tetrakis_hexahedron(allocator, temp_alloc)
	case .Rhombic_Dodecahedron:			return generate_rhombic_dodecahedron(allocator, temp_alloc)
	case .Pentakis_Dodecahedron:		return generate_pentakis_dodecahedron(allocator, temp_alloc)
	case .Disdyakis_Dodecahedron:		return generate_disdyakis_dodecahedron(allocator, temp_alloc)
	case .Rhombic_Triacontahedron:		return generate_rhombic_triacontahedron(allocator, temp_alloc)
	case .Disdyakis_Triacontahedron:	return generate_disdyakis_triacontahedron(allocator, temp_alloc)
	case .Deltoidal_Hexecontahedron:	return generate_deltoidal_hexecontahedron(allocator, temp_alloc)
	case .Deltoidal_Icositetrahedron:	return generate_deltoidal_icositetrahedron(allocator, temp_alloc)
	case .Pentagonal_Hexecontahedron:	return generate_pentagonal_hexecontahedron(allocator, temp_alloc)
	case .Pentagonal_Icositetrahedron:	return generate_pentagonal_icositetrahedron(allocator, temp_alloc)
	}
	unreachable()
}

generate_all_catalan_solids :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> [Catalan_Solid]m.Mesh {
	solids := [Catalan_Solid]m.Mesh{}
	for &s, t in solids {
		s = mesh_generate_catalan_solid(t, allocator, temp_alloc)
	}
	return solids
}

generate_sierpinski_tetrahedron :: proc(depth: int, allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	p0 := [3]f32{1, 1, 1}
	p1 := [3]f32{-1, -1, 1}
	p2 := [3]f32{-1, 1, -1}
	p3 := [3]f32{1, -1, -1}

	subdivide :: proc(mesh: ^m.Mesh, a, b, c, d: [3]f32, depth: int) {
		if depth == 0 {
			ia := m.add_vertex(mesh, a)
			ib := m.add_vertex(mesh, b)
			ic := m.add_vertex(mesh, c)
			id := m.add_vertex(mesh, d)

			m.add_faces(mesh, {ia, ib, ic}, {ia, ic, id}, {ia, id, ib}, {ib, id, ic})
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

	m.add_boundaries(&mesh)
	return mesh
}
