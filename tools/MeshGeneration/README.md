# Spherical Shell Q2 Mesh Generator

Generates 3D spherical shell meshes with 18-node prismatic Q2 elements using icosphere refinement and radial extrusion.

The mesh is created by first generating a triangulated inner sphere (inlet boundary), then extruding radially outward to create prismatic layers until reaching the outer sphere (outlet boundary). The inner sphere surface resolution is controlled by the refinement level (`-l`). The radial refinement is controlled by the number of layers (`-n`) and the distribution method (`-d`), where layers can be equidistant or follow a geometric progression. 

Boundaries are annotated as "Inlet" (inner sphere) and "Outlet" (outer sphere). Output format is CFmesh by default.

## Usage

```bash
python generate_spherical_shell_Q2.py -r INNER_RADIUS -R OUTER_RADIUS -l REFINEMENT_LEVEL -n NUM_LAYERS [options]
```

For detailed help and all available options:
```bash
python generate_spherical_shell_Q2.py --help
```

## Examples

```bash
# Geometric (radial) distribution 
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 20 -k 1.2
 
# Equidistant layers
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 30 -d equidistant

# Manual layer heights
python generate_spherical_shell_Q2.py -r 1.0 -R 2.0 -l 2 -d manual -m 0.3 0.4 0.3

# Output both CFmesh and MSH formats
python generate_spherical_shell_Q2.py -r 21.5 -R 235 -l 3 -n 20 -k 1.2 -f both
```

### Required Parameters
- `-r`: Inner sphere radius
- `-R`: Outer sphere radius  
- `-l`: Refinement level (0-5)
- `-n`: Number of radial layers (except for manual distribution)

### Optional Parameters
- `-d`: Distribution method: `equidistant`, `geometric` (default), or `manual`
- `-k`: Growth ratio for geometric distribution (default: 1.0)
- `-m`: Manual layer heights (space-separated)
- `-f`: Output format: `cfmesh` (default), `msh`, or `both`
- `-o`: Output filename (auto-generated if not specified)

## Refinement Levels

| Level | Triangles | Description      |
|-------|-----------|------------------|
| 0     | 20        | Base icosahedron |
| 1     | 80        | Coarse           |
| 2     | 320       | Medium           |
| 3     | 1,280     | Fine             |
| 4     | 5,120     | Very fine        |
| 5     | 20,480    | Extra fine       |

## Distribution Methods

- **Equidistant**: Uniform layer thickness
- **Geometric**: Exponential growth with ratio `k` (k>1 for boundary layers)
- **Manual**: Custom layer heights

## Reference

For more information about the mesh generation procedure, please refer to the Mesh Generation section of:

**Dhib, R., Ben Ameur, F., Sharma, V., Lani, A., & Poedts, S. (2025).** *Toward High-order Solar Corona Simulations: A High-order Hyperbolized Poisson Approach for Magnetic Field Initialization.* The Astrophysical Journal, 980(2), 163. [DOI: 10.3847/1538-4357/adace5](https://doi.org/10.3847/1538-4357/adace5)

