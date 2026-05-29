package polygon

import "core:math"
import "core:math/linalg"
import c "../" // convay
import m "../../" // mesh

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

generate_tetrahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, {1, 1, 1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1})
	m.add_faces(&mesh, {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2})
	m.normalize(&mesh)
	return mesh
}

generate_cube :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1})
	m.add_faces(&mesh, {3, 2, 1, 0}, {3, 0, 5, 4}, {4, 5, 6, 7}, {7, 6, 1, 2}, {6, 5, 0, 1}, {2, 3, 4, 7})
	m.normalize(&mesh)
	return mesh
}

generate_octahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, { 1,  0,  0},{-1,  0,  0},{ 0,  1,  0},{ 0, -1,  0},{ 0,  0,  1},{ 0,  0, -1})
	m.add_faces(&mesh,{0, 4, 2},{2, 4, 1},{1, 4, 3},{3, 4, 0},{2, 5, 0},{1, 5, 2},{3, 5, 1},{0, 5, 3})
	m.normalize(&mesh)
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
		{11, 4, 2}, {10, 2, 6}, {7, 6, 8}, {1, 8, 9},
	)
	m.normalize(&mesh)
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
    m.normalize(&mesh)
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
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_octahedron(allocator)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_rhombicuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_cuboctahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_snub_cube :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.classical_snub(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_icosahedron(allocator)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_rhombicosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_truncated_icosidodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.ambo(&mesh, temp_alloc = temp_alloc)
	c.truncate(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_snub_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.classical_snub(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_archimedean_solid :: proc (solid : Archimedean_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
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
		s = generate_archimedean_solid(t, allocator, temp_alloc)
	}
	return solids
}

generate_triakis_tetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_tetrahedron(allocator)
	c.kis(&mesh, CATALAN_TRI_TETRAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_rhombic_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.join(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_triakis_octahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_octahedron(allocator)
	c.kis(&mesh, CATALAN_TRI_OCTAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_tetrakis_hexahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.kis(&mesh, CATALAN_TETRA_HEXAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_deltoidal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.ortho(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_disdyakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.meta(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_pentagonal_icositetrahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_cube(allocator)
	c.classical_gyro(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_rhombic_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.join(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_triakis_icosahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_icosahedron(allocator)
	c.kis(&mesh, CATALAN_TRI_ICOSAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_pentakis_dodecahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.kis(&mesh, CATALAN_PENTA_DODECAHEDRON_KIS_HEIGHT, temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_deltoidal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.ortho(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_disdyakis_triacontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.meta(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_pentagonal_hexecontahedron :: proc(allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
	mesh := generate_dodecahedron(allocator)
	c.classical_gyro(&mesh, temp_alloc = temp_alloc)
	m.normalize(&mesh)
	return mesh
}

generate_catalan_solid :: proc(solid : Catalan_Solid, allocator := context.allocator, temp_alloc := context.temp_allocator) -> m.Mesh {
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
		s = generate_catalan_solid(t, allocator, temp_alloc)
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

generate_torus :: proc(major_radius, minor_radius: f32, segment_count_u, segment_count_v: int, generate_u, generate_v: int, allocator := context.allocator) -> m.Mesh {
    gen_u := min(segment_count_u, generate_u)
    gen_v := min(segment_count_v, generate_v)

    mesh := m.create(allocator)

    // Generate all vertices from the parametric torus equation.
    // x : major_radius * cos u + minor_radius * cos u * cos v
    // y : major_radius * cos v + minor_radius * cos u * sin v
    // z : minor_radius * sin v

    for u in 0..<gen_u {
        u_ratio := f32(u) / f32(segment_count_u)
        theta := u_ratio * 2.0 * f32(linalg.PI)

        cos_theta := math.cos(theta)
        sin_theta := math.sin(theta)

        for v in 0..<gen_v {
            v_ratio := f32(v) / f32(segment_count_v)
            phi := v_ratio * 2.0 * f32(linalg.PI)

            cos_phi := math.cos(phi)
            sin_phi := math.sin(phi)

            x := (major_radius + minor_radius * cos_phi) * cos_theta
            y := (major_radius + minor_radius * cos_phi) * sin_theta
            z := minor_radius * sin_phi

            m.add_vertices(&mesh, {x, y, z})
        }
    }

    index :: proc(u, v, seg_v: int) -> m.Vertex_Index {
        return m.Vertex_Index(u * seg_v + v)
    }

    wrap_u := gen_u == segment_count_u
    wrap_v := gen_v == segment_count_v

    max_u := gen_u if wrap_u else gen_u - 1
    max_v := gen_v if wrap_v else gen_v - 1

    for u in 0..<max_u {
        next_u := (u + 1) % gen_u

        for v in 0..<max_v {
            next_v := (v + 1) % gen_v

            a := index(u,      v,      gen_v)
            b := index(next_u, v,      gen_v)
            c := index(next_u, next_v, gen_v)
            d := index(u,      next_v, gen_v)

            m.add_faces(&mesh, {d, c, b, a})
        }
    }

    return mesh
}
