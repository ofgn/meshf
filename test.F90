program test
    use vtk
    use tetrahedral_mesh
    use data_structures

    implicit none

    integer(int64) :: i
    integer(int64), allocatable :: star(:)

    character(len=256) :: path
    character(len=256) :: node_file
    character(len=256) :: edge_file
    character(len=256) :: face_file
    character(len=256) :: ele_file
    character(len=256) :: neigh_file
    character(len=256) :: t2f_file

    character(len=256) :: sigma_file
    character(len=256) :: vtk_file

    type(VtkUnstructuredGrid) :: u_grid
    type(VtkUnstructuredGrid) :: u_grid2
    type(VtkData) :: point_data, cell_data
    logical, allocatable :: cell_mask(:), point_mask(:)

    integer(int64), allocatable :: value(:)

    type(MultiTesselation) :: mesh


    type(PriorityQueue) :: queue
    integer(int64), allocatable :: edge(:)
    real(real64) :: priority

    path = '/srv/wrk/temp/example_mesh/'

    node_file = trim(path) // 'e4d.1.node'
    edge_file = trim(path) // 'e4d.1.edge'
    face_file = trim(path) // 'e4d.1.face'
    ele_file = trim(path) // 'e4d.1.ele'
    neigh_file = trim(path) // 'e4d.1.neigh'
    t2f_file = trim(path) // 'e4d.1.t2f'
    vtk_file = trim(path) // 'mesh.vtk'

    ! call u_grid%read_tetgen_node(node_file, point_data)
    ! call u_grid%read_tetgen_ele(ele_file, cell_data)

    ! cell_mask = cell_data%scalar_real64(1, :) < 1.5d0
    ! call u_grid%mask_cells(cell_mask, cell_data=cell_data)

    ! call mesh%initialise(u_grid)
    call mesh%read_tetgen_node_file(node_file)
    call mesh%read_tetgen_edge_file(edge_file)
    call mesh%read_tetgen_face_file(face_file)
    call mesh%read_tetgen_ele_file(ele_file)
    call mesh%read_tetgen_neigh_file(neigh_file)
    call mesh%read_tetgen_t2f_file(t2f_file)

    
    call mesh%half_edge_collapse(47511,  28183, 0d0)
    call mesh%apply_updates(1)

    u_grid = mesh%export_vtk()
    call u_grid%write_legacy_vtk(vtk_file)


    ! call mesh%build_edge_queue()

    ! call u_grid2%write_legacy_vtk(vtk_file, cell_data=cell_data)



end program test