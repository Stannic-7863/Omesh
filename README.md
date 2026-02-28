A mesh library written in odin. The original motivation for write this library was to generate basic shapes/polyhedrons in a more formal way. 

# Progress
- All convay ops [from the wikipedia list](https://en.wikipedia.org/wiki/Conway_polyhedron_notation?useskin=vector#:~:text=Original%20operations)
- Functions to generate platonic, archimedean and catalan solids. 

# Usage 

```odin
import m "mesh"
import c "mesh/convay"
import p "mesh/convay/polygon"

// procedure from "mesh/convay/polygone" package
generate_tetrahedron :: proc(allocator := context.allocator) -> m.Mesh {
	mesh := m.create(allocator)
	m.add_vertices(&mesh, {1, 1, 1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1})
	m.add_faces(&mesh, {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2})
	m.normalize_onto_sphere(&mesh)
	// m.validate(mesh) validates some constraints for the mesh and the half-edge structure. 
	return mesh
}

do_convay_stuff_with_mesh :: proc(mesh: ^c.Mesh) {
	// All function that need some temporary storage will take in a temp_alloc paramter with context.temp_allocator as the default value.
	c.ambo(mesh, temp_alloc = context.temp_alloc) 
	// some convay operations are parameteried. In near future, all will be parameterized (I hope). 
	// These parameters have default values, so you don't have to specify em all the time
	c.kis(mesh, kis_height = 0.2) 
	
	// For convenience
	c.operations(mesh, .Kis, .Ambo, .Kis, .Truncate, kis_height = 0.7, ambo_factor = 0.2) 
}

main :: proc() {
	cube := p.generate_cube(context.allocator)
	tetrahedron := generate_tetrahedron()
	do_convay_stuff_with_mesh(&cube)	
	do_convay_stuff_with_mesh(&tetrahedron)

	// A mesh for rendering can be created like this: 
	render_indices := [dynamic]int{}
	render_vertex := [dynamic]Vertex{} // position + normal 
	iter := m.create_triangle_emitter_iter(&cube)
 	for count, position, normal, indices in m.triangle_emitter_indexed_flat_iter(&iter) { 
  	for i in indices { append(&render_indices, i) }
		for i := 0; i < count; i += 1 { append(&render_vertex, Vertex{position = position[i], normal = normal}) }
	}
	
	magically_render(render_indices[:], render_vertex[:])
	
	m.destroy(cube, tetrahedron)
}
```
