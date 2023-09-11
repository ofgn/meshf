# meshf

> 🚧 **Warning**: This project is in its development stage. The current version only supports basic 2D Delaunay triangulation. Future plans include implementing Delaunay refinement, constrained Delaunay, 3D Delaunay, and structured mesh capabilities.

## Features:

- 2D Delaunay Triangulation (Basic)

## Supported File Formats:

### Current Support:

- Partial support for `.vtk`, `triangle/tetgen` file formats.

## Planned Features:

- 3D Delaunay triangulation
- 2D/3D Delaunay refinement
- 2D/3D Constrained Delaunay Triangulation
- Structured mesh generation
- Mesh file format conversion

## Prerequisites

1. **Fortran Compiler**: A Fortran compiler is required. It is recommended that the compiler used is compatible with Fortran 2018. 
    - Popular Fortran compilers include `gfortran`, `ifort`, and `ifx`.

2. **Make Utility**: GNU Make is recommended for building the project.

## Installation

1. **Clone the Repository**
   ```bash
    git clone https://github.com/ofgn/meshf.git
    cd meshf/
    ```
2. **Build the Project**

    To build the project, specify your Fortran compiler and run make:

    ```bash
    make FC=ifx
    ```
  
## Usage:

> 

## Development:

Contributions are welcome! To contribute, please follow these steps:

1. Fork this repository.
2. Create a new feature branch (`git checkout -b new-feature`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Create a Pull Request.


## License:

This project is licensed under the GNU Affero General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

