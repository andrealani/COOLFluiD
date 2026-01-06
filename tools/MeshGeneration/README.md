# Spherical Shell Q2 Mesh Generator

Generates 3D spherical shell meshes with 18-node prismatic Q2 elements using icosphere refinement or sphere triangulation, and radial extrusion.

The mesh is created by first generating a triangulated inner sphere (inlet boundary), then extruding radially outward to create prismatic layers until reaching the outer sphere (outlet boundary). 

Two methods are available for generating the sphere surface:
- **Icosphere**: Uses recursive refinement of an icosahedron. Triangle count is 20 × 4^level.
- **Fibonacci**: Uses Fibonacci spiral point distribution with Delaunay triangulation. Allows arbitrary point count (~2N-4 triangles for N points).

The radial refinement is controlled by the number of layers (`-n`) and the distribution method (`-d`), where layers can be equidistant, follow a geometric progression, or use a bump distribution for refinement at both boundaries.

Boundaries are annotated as "Inlet" (inner sphere) and "Outlet" (outer sphere). Output format is CFmesh by default.

## Usage

```bash
# Using icosphere method (default)
python generate_spherical_shell_Q2.py -r INNER_RADIUS -R OUTER_RADIUS -l REFINEMENT_LEVEL -n NUM_LAYERS [options]

# Using fibonacci method
python generate_spherical_shell_Q2.py -r INNER_RADIUS -R OUTER_RADIUS -s fibonacci -p NUM_POINTS -n NUM_LAYERS [options]
```

For detailed help and all available options:
```bash
python generate_spherical_shell_Q2.py --help
```

## Examples

```bash
# Icosphere with geometric (radial) distribution 
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 20 -k 1.2

# Fibonacci sphere with 500 points (~996 triangles)
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -s fibonacci -p 500 -n 20
 
# Equidistant layers
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 30 -d equidistant

# Bump distribution (refined at both inlet and outlet)
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 20 -d bump -b 0.2

# Manual layer heights
python generate_spherical_shell_Q2.py -r 1.0 -R 2.0 -l 2 -d manual -m 0.3 0.4 0.3

# Output both CFmesh and MSH formats
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 20 -k 1.2 -f both
```

### Required Parameters
- `-r`: Inner sphere radius
- `-R`: Outer sphere radius  
- `-l`: Refinement level for icosphere (0-5), OR
- `-p`: Number of points for fibonacci method
- `-n`: Number of radial layers (except for manual distribution)

### Optional Parameters
- `-s`: Sphere method: `icosphere` (default) or `fibonacci`
- `-d`: Distribution method: `equidistant`, `geometric` (default), `bump`, or `manual`
- `-k`: Growth ratio for geometric distribution (default: 1.0)
- `-b`: Bump ratio for bump distribution (0-1, default: 0.3). Controls boundary/middle layer ratio. Smaller values = more refinement at boundaries
- `-m`: Manual layer heights (space-separated)
- `-f`: Output format: `cfmesh` (default), `msh`, or `both`
- `-o`: Output filename (auto-generated if not specified)

## Sphere Methods

### Icosphere (default)
Uses recursive subdivision of an icosahedron. Triangle count increases by 4x per refinement level.

| Level | Triangles | Description      |
|-------|-----------|------------------|
| 0     | 20        | Base icosahedron |
| 1     | 80        | Coarse           |
| 2     | 320       | Medium           |
| 3     | 1,280     | Fine             |
| 4     | 5,120     | Very fine        |
| 5     | 20,480    | Extra fine       |

### Fibonacci
Uses Fibonacci spiral for uniform point distribution, then Delaunay triangulation. Provides fine-grained control over triangle count.

| Points | ~Triangles | 
|--------|------------|
| 50     | 96         |
| 100    | 196        |
| 200    | 396        |
| 500    | 996        |
| 1000   | 1996       |
| 2000   | 3996       |

Formula: triangles ≈ 2 × points - 4 

## Distribution Methods

- **Equidistant**: Uniform layer thickness
- **Geometric**: Exponential growth with ratio `k` (k>1 for boundary layers at inlet)
- **Bump**: Refined at both boundaries (inlet and outlet), coarser in the middle. Uses a sine-based distribution. Control with `-b` (bump ratio, 0-1). Smaller values give more refinement at boundaries.
- **Manual**: Custom layer heights

## Reference

For more information about the mesh generation procedure, please refer to the Mesh Generation section of:

**Dhib, R., Ben Ameur, F., Sharma, V., Lani, A., & Poedts, S. (2025).** *Toward High-order Solar Corona Simulations: A High-order Hyperbolized Poisson Approach for Magnetic Field Initialization.* The Astrophysical Journal, 980(2), 163. [DOI: 10.3847/1538-4357/adace5](https://doi.org/10.3847/1538-4357/adace5)

