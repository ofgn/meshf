program test
    use vtk
    use tetrahedral_mesh
    use data_structures

    implicit none

    integer(int64) :: i
    character(len=256) :: path
    character(len=256) :: node_file, ele_file, sigma_file
    character(len=256) :: vtk_file
    type(VtkUnstructuredGrid) :: u_grid
    type(VtkUnstructuredGrid) :: u_grid2
    type(VtkData) :: point_data, cell_data
    logical, allocatable :: cell_mask(:), point_mask(:)

    type(TetrahedralMesh) :: mesh
    type(HashMap) :: map
    integer(int64), allocatable :: tetra(:)
    integer(int64), allocatable :: face(:)

    integer(int64) :: tetrahedron
    integer(int64) :: u
    integer(int64) :: v
    integer(int64) :: w
    integer(int64) :: x

    path = '/home/oliver/example_mesh/'

    node_file = trim(path) // 'e4d.1.node'
    ele_file = trim(path) // 'e4d.1.ele'
    vtk_file = trim(path) // 'mesh.vtk'

    call u_grid%read_tetgen_node(node_file, point_data)
    call u_grid%read_tetgen_ele(ele_file, cell_data)

    


    ! cell_mask = cell_data%scalar_real64(1, :) < 1.5d0
    ! call u_grid%mask_cells(cell_mask, cell_data=cell_data)
    
    call mesh%initialise(u_grid)
    print *, mesh%tetrahedron_neighbours(2)



    u_grid2 = mesh%export_u_grid()
    call u_grid2%write_legacy_vtk(vtk_file, cell_data=cell_data)



end program test