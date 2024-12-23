program test
    use vtk
    use multi_resolution_analysis
    use data_structures_algorithms
    use geometry
    use global

    implicit none

    integer(int64) :: i, j
    integer(int64), allocatable :: fine(:)
    integer(int64), allocatable :: coarse(:)
    real(real64) :: volume

    character(len=256) :: directory
    character(len=256) :: node_file
    character(len=256) :: edge_file
    character(len=256) :: face_file
    character(len=256) :: ele_file
    character(len=256) :: neigh_file
    character(len=256) :: t2f_file
    character(len=256) :: vtk_file

    type(TetrahedralMesh) :: coarse_mesh
    type(TetrahedralMesh) :: base_mesh
    type(VtkUnstructuredGrid) :: u_grid
    type(UpdateNode), pointer :: current_update  

    real(real64) :: tetra1(3, 4)
    real(real64) :: tetra2(3, 4)

    directory = '/srv/wrk/temp/example_mesh/'

    node_file = trim(directory) // 'e4d.1.node'
    edge_file = trim(directory) // 'e4d.1.edge'
    face_file = trim(directory) // 'e4d.1.face'
    ele_file = trim(directory) // 'e4d.1.ele'
    neigh_file = trim(directory) // 'e4d.1.neigh'
    t2f_file = trim(directory) // 'e4d.1.t2f'
    vtk_file = trim(directory) // 'mesh.vtk'

    !call configure_logging(enabled=.false.)

    call coarse_mesh%read_tetgen_node_file(node_file)
    call coarse_mesh%read_tetgen_ele_file(ele_file)
    call coarse_mesh%read_tetgen_neigh_file(neigh_file)

    base_mesh = coarse_mesh

    call coarse_mesh%simplify_mesh(0.1_real64)
    call coarse_mesh%apply_updates()

    current_update => coarse_mesh%first_update
    do while (associated(current_update))
        print *, current_update%u
        fine = [current_update%u_plus, current_update%u_minus]
        coarse = current_update%u_plus

        do i = 1, size(coarse)
            tetra1 = coarse_mesh%vertices(:, coarse_mesh%tetrahedra(:, coarse(i)))
            do j = 1, size(fine)
                tetra2 = base_mesh%vertices(:, base_mesh%tetrahedra(:, fine(j)))
                if (tetrahedra_intersect(tetra1, tetra2)) then
                    print *, "Coarse ", coarse(i), " intersects with fine ", fine(j)
                end if
            end do
        end do
        
        current_update => current_update%next
    end do


    ! u_grid = mesh%vtk_unstructured_grid()
    ! call u_grid%write_legacy_vtk(vtk_file, ascii=.false.)


    ! do i = 1, mesh%n_vertices
    !     print *, 'Vertex ', i, ': ', mesh%vertices(:, i)
    ! end do  

    




end program test