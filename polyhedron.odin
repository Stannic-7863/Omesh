package mesh

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
	#partial switch solid {
	case .Snub_Cube:					return mesh_generate_snub_cube(allocator, temp_alloc)
	// case .Cuboctahedron:		 		return mesh_generate_cuboctahedron(allocator, temp_alloc)
	// case .Truncated_Cube:		 		return mesh_generate_truncated_cube(allocator, temp_alloc)
	// case .Snub_dodecahedron:			return mesh_generate_snub_dodecahedron(allocator, temp_alloc)
	// case .Icosidodecahedron:			return mesh_generate_icosidodecahedron(allocator, temp_alloc)
	// case .Rhombicuboctahedron:	 		return mesh_generate_rhombicuboctahedron(allocator, temp_alloc)
	// case .Truncated_Octahedron:	 		return mesh_generate_truncated_octahedron(allocator, temp_alloc)
	// case .Truncated_Tetrahedron: 		return mesh_generate_truncated_tetrahedron(allocator, temp_alloc)
	// case .Truncated_icosahedron:		return mesh_generate_truncated_icosahedron(allocator, temp_alloc)
	// case .Truncated_dodecahedron:		return mesh_generate_truncated_dodecahedron(allocator, temp_alloc)
	// case .Rhombicosidodecahedron:		return mesh_generate_rhombicosidodecahedron(allocator, temp_alloc)
	// case .Truncated_Cuboctahedron:		return mesh_generate_truncated_cuboctahedron(allocator, temp_alloc)
	// case .Truncated_Icosidodecahedron:	return mesh_generate_truncated_icosidodecahedron(allocator, temp_alloc)
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
