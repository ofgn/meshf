! ------------------------------------------------------------------------------
! @file tetrahedral.f90
! @brief Module for handling tetrahedral mesh operations
! @author ofgn
! @date 2024-10-03
! ------------------------------------------------------------------------------
module tetrahedral_mesh
    use iso_fortran_env, only: int32, int64, real32, real64
    use utility

    implicit none

    type :: MultiTesselation
        integer(int64) :: n_updates = 0                                         ! Number of updates              
        integer(int64) :: n_vertices = 0                                        ! Number of vertices in the grid
        integer(int64) :: n_edges = 0                                           ! Number of edges in the grid
        integer(int64) :: n_faces = 0                                           ! Number of faces in the grid
        integer(int64) :: n_tetrahedra = 0                                      ! Number of tetrahedra in the grid
        real(real64), allocatable :: vertices(:, :)                             ! Coordinates of each vertex (3 x n_vertices)
        integer(int64), allocatable :: edges(:, :)                              ! Edge connectivity (2 x n_edges)
        integer(int64), allocatable :: faces(:, :)                              ! Face connectivity (3 x n_faces)
        integer(int64), allocatable :: tetrahedra(:, :)                         ! Cell connectivity (4 x n_tetra)
        integer(int64), allocatable :: edge_boundary(:)                         ! Boundary markers for each edge
        integer(int64), allocatable :: face_boundary(:)                         ! Boundary markers for each face
        integer(int64), allocatable :: vertex_adjacency(:)                      ! A single tetrahedron adjacent to each vertex
        !integer(int64), allocatable :: edge_adjacency(:)                        !NOT USED! Tetrahedra adjacent to each edge
        !integer(int64), allocatable :: face_adjacency(:, :)                     !NOT USED ! Tetrahedra adjacent to each face
        integer(int64), allocatable :: tetrahedron_adjacency(:, :)              ! Tetrahedra adjacent to each tetrahedron
        integer(int64), allocatable :: tetrahedron_faces(:, :)                  ! Faces of each tetrahedron
        type(UpdateNode), pointer :: first_update => null()                     ! Head of updates list
        type(UpdateNode), pointer :: last_update => null()                      ! Tail of updates list
    contains
        procedure :: read_tetgen_node_file                                      ! Read node file
        procedure :: read_tetgen_edge_file                                      ! Read edge file
        procedure :: read_tetgen_face_file                                      ! Read face file    
        procedure :: read_tetgen_ele_file                                       ! Read ele file                
        procedure :: read_tetgen_neigh_file                                     ! Read neigh file
        procedure :: read_tetgen_t2f_file                                       ! Read t2f file
        procedure :: simplify_mesh                                              ! Simplify the mesh
        procedure :: export_vtk                                                 ! Export the mesh to a VtkUnstructuredGrid         
        procedure :: find_vertex_star                                           ! Find all tetrahedra containing a vertex
        procedure :: half_edge_collapse                                         ! Collapse an edge in the mesh
        procedure :: add_update
        procedure :: apply_updates
    end type MultiTesselation

    type :: UpdateNode
        integer(int64) :: id                                                    ! Unique node ID
        integer(int64) :: v                                                     ! Vertex being collapsed 
        integer(int64) :: w                                                     ! Target vertex
        integer(int64), allocatable :: u_minus(:)                               ! Tetrahedra removed
        integer(int64), allocatable :: u_plus(:)                                ! Tetrahedra modified
        real(real64) :: error                                                   ! Error estimate
        type(UpdateNode), pointer :: next => null()                             ! Next update in list
        type(UpdateNode), pointer :: prev => null()                             ! Previous update
    end type UpdateNode


contains

    ! --------------------------------------------------------------------------
    ! @brief Read a TetGen node file.
    ! @param[inout] self The tetrahedral mesh to populate.
    ! @param[in] file_path The path to the node file.
    ! @details Reads the node file and populates the vertices of the mesh.
    ! --------------------------------------------------------------------------
    subroutine read_tetgen_node_file(self, file_path)
        implicit none

        class(MultiTesselation), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to node file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index                  
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen node file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_vertices

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate(self%vertices(3, self%n_vertices))

        write(message, "(A, I0)") "Number of vertices: ", self%n_vertices
        call report(trim(message), 1)

        do j = 1, self%n_vertices
            read (unit, *, iostat=io_status) i, self%vertices(:, j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if
        end do 

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_node_file

    ! --------------------------------------------------------------------------
    ! @brief Read a TetGen edge file.
    ! @param[inout] self The tetrahedral mesh to populate.
    ! @param[in] file_path The path to the edge file.
    ! @details Reads the edge file and populates the edges of the mesh.
    ! --------------------------------------------------------------------------
    subroutine read_tetgen_edge_file(self, file_path)
        implicit none

        class(MultiTesselation), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to edge file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int32) :: j                                                     ! Loop index                
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen edge file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_edges

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate(self%edges(2, self%n_edges))
        allocate(self%edge_boundary(self%n_edges))
        !allocate(self%edge_adjacency(self%n_edges))

        write(message, "(A, I0)") "Number of edges: ", self%n_edges
        call report(trim(message), 1)

        do j = 1, self%n_edges
            read (unit, *, iostat=io_status) i, self%edges(:, j) &
                , self%edge_boundary(j)!, self%edge_adjacency(j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if
        end do 

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_edge_file

    subroutine read_tetgen_face_file(self, file_path)
        implicit none

        class(MultiTesselation), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to face file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index                  
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen face file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_faces

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate(self%faces(3, self%n_faces))
        allocate(self%face_boundary(self%n_faces))
        !allocate(self%face_adjacency(2, self%n_faces))


        write(message, "(A, I0)") "Number of faces: ", self%n_faces
        call report(trim(message), 1)

        do j = 1, self%n_faces
            read (unit, *, iostat=io_status) i, self%faces(:, j) &
                , self%face_boundary(j)!, self%face_adjacency(:, j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if
        end do 

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_face_file

    subroutine read_tetgen_ele_file(self, file_path)
        use geometry

        implicit none

        class(MultiTesselation), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to face file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index 
        integer(int64) :: u                                                     ! Vertex u
        integer(int64) :: v                                                     ! Vertex v
        integer(int64) :: w                                                     ! Vertex w
        integer(int64) :: x                                                     ! Vertex x                   
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        real(real128) :: orientation                                            ! Orientation of tetrahedron
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen ele file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_tetrahedra

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate(self%tetrahedra(4, self%n_tetrahedra))
        allocate(self%vertex_adjacency(self%n_vertices))

        do j = 1, self%n_tetrahedra
            read (unit, *, iostat=io_status) i, self%tetrahedra(:, j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if

            ! Get vertices of the tetrahedron
            u = self%tetrahedra(1, i)
            v = self%tetrahedra(2, i)
            w = self%tetrahedra(3, i)
            x = self%tetrahedra(4, i)

            ! Calculate orientation of tetrahedron
            orientation = orientation3d(self%vertices(:, u), &
                self%vertices(:, v), self%vertices(:, w), self%vertices(:, x))

            ! Ensure tetrahedron is positively oriented
            if (orientation < 0.0_real128) then
                ! Swap vertices to ensure positive orientation
                self%tetrahedra(:, i) = [u, v, x, w]
                u = self%tetrahedra(1, i)
                v = self%tetrahedra(2, i)
                w = self%tetrahedra(3, i)
                x = self%tetrahedra(4, i)
            end if
    
            ! Set vertices into adjacency list
            self%vertex_adjacency(u) = i
            self%vertex_adjacency(v) = i
            self%vertex_adjacency(w) = i
            self%vertex_adjacency(x) = i
        end do

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_ele_file

    subroutine read_tetgen_neigh_file(self, file_path)
        implicit none

        class(MultiTesselation), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to face file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index                  
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen neigh file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_tetrahedra

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate(self%tetrahedron_adjacency(4, self%n_tetrahedra))

        do j = 1, self%n_tetrahedra
            read (unit, *, iostat=io_status) i, self%tetrahedron_adjacency(:, j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if
        end do

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_neigh_file

    subroutine read_tetgen_t2f_file(self, file_path)
        implicit none

        class(MultiTesselation), intent(inout) :: self                          ! Tetrahedral mesh structure
        character(len=256), intent(in) :: file_path                             ! Path to face file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index                  
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: message                                          ! String buffer

        call cpu_time(start_time)

        call report("Reading TetGen t2f file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        allocate(self%tetrahedron_faces(4, self%n_tetrahedra))

        do j = 1, self%n_tetrahedra
            read (unit, *, iostat=io_status) i, self%tetrahedron_faces(:, j)
            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if
        end do

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Export the tetrahedral mesh to an unstructured grid.
    ! @param[in] self The tetrahedral mesh to export.
    ! @param[out] u_grid The resulting unstructured grid.
    ! @details Converts the internal representation of vertices and tetrahedra 
    !   to the VtkUnstructuredGrid.
    ! --------------------------------------------------------------------------
    function export_vtk(self) result(u_grid)
        use vtk
        implicit none

        class(MultiTesselation), intent(in) :: self                              ! The tetrahedral mesh

        type(VtkUnstructuredGrid) :: u_grid                                     ! Output unstructured grid
        integer(kind=int64) :: i                                                ! Loop variable
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        logical, allocatable :: mask(:), mask2(:, :)                            ! Masks for packing connectivity
        character(len=4096) :: message                                          ! String buffer

        write(message, '(A)') "Converting TetrahedralMesh to " &
            // "VtkUnstructuredGrid"
        call report(trim(message), 1)

        u_grid%n_points = self%n_vertices
        mask = .not. all(self%tetrahedra == 0, dim=1)
        mask2 = spread(mask, dim=1, ncopies=4)

        u_grid%n_cells = self%n_tetrahedra

        allocate (u_grid%points(3, u_grid%n_points))
        allocate (u_grid%cell_connectivity(4 * u_grid%n_cells))
        allocate (u_grid%cell_type(u_grid%n_cells))
        allocate (u_grid%cell_size(u_grid%n_cells))

        u_grid%points = self%vertices
        u_grid%cell_connectivity = pack(self%tetrahedra, mask2) - 1
        u_grid%cell_type = 10
        u_grid%cell_size = 4

        write(message, '(A, I0)') "Number of points: ", u_grid%n_points
        call report(trim(message), 1)
        write(message, '(A, I0)') "Number of cells: ", u_grid%n_cells
        call report(trim(message), 1)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1) 
    end function export_vtk

    ! --------------------------------------------------------------------------
    ! @brief Find all tetrahedra containing a vertex using depth-first search.
    ! @param[in] self The tetrahedral mesh.
    ! @param[in] vertex_index The vertex to search from.
    ! @return vertex_star Array of tetrahedra IDs containing the vertex.
    ! --------------------------------------------------------------------------
    function find_vertex_star(self, vertex_index) result(vertex_star)
        use algorithms

        implicit none
        
        class(MultiTesselation), intent(in) :: self
        integer(int64), intent(in) :: vertex_index
        integer(int64), allocatable :: vertex_star(:)
        
        logical :: visited(self%n_tetrahedra)  ! Stack allocated
        integer(int64) :: stack(self%n_tetrahedra)  ! Fixed size stack
        integer(int64) :: stack_top, result_top
        integer(int64) :: current_tet, next_tet, i
        integer(int64) :: start_tet
        
        ! Fast initialisation
        visited = .false.
        result_top = 0
        
        ! Get starting tetrahedron
        start_tet = self%vertex_adjacency(vertex_index)
        if (start_tet <= 0 .or. start_tet > self%n_tetrahedra) then
            allocate(vertex_star(0))
            return
        end if
        
        ! Initialise stack
        stack_top = 1
        stack(1) = start_tet
        visited(start_tet) = .true.
        
        ! Preallocate result (will trim later)
        allocate(vertex_star(self%n_tetrahedra))
        vertex_star(1) = start_tet
        result_top = 1
        
        ! Fast DFS traversal
        do while (stack_top > 0)
            current_tet = stack(stack_top)
            stack_top = stack_top - 1
            
            ! Check all 4 neighbors inline
            do i = 1, 4
                next_tet = self%tetrahedron_adjacency(i, current_tet)
                if (next_tet > 0) then
                    if (.not. visited(next_tet)) then
                        if (any(self%tetrahedra(:,next_tet) == vertex_index)) then
                            stack_top = stack_top + 1
                            stack(stack_top) = next_tet
                            visited(next_tet) = .true.
                            result_top = result_top + 1
                            vertex_star(result_top) = next_tet
                        end if
                    end if
                end if
            end do
        end do
        
        ! Trim result to actual size
        vertex_star = vertex_star(1:result_top)
        call quicksort_1d(vertex_star)
    end function find_vertex_star

    
    subroutine half_edge_collapse(self, v, w, error)
        implicit none

        class(MultiTesselation), intent(inout) :: self
        integer(int64), intent(in) :: v
        integer(int64), intent(in) :: w
        real(real64), intent(in) :: error
        
        type(UpdateNode) :: new_update
        integer(int64), allocatable :: vertex_star(:)
        logical, allocatable :: affected(:)
        character(len=4096) :: message
        integer(int64) :: i, tetra
        
        ! Initialise update
        new_update%id = self%n_updates + 1
        new_update%v = v
        new_update%w = w
        new_update%error = error
        
        allocate(new_update%u_minus(0))
        allocate(new_update%u_plus(0))
        
        ! Find affected tetrahedra
        vertex_star = self%find_vertex_star(new_update%v)
        allocate(affected(self%n_tetrahedra))
        affected = .false.
        
        ! Mark affected tetrahedra
        do i = 1, size(vertex_star)
            tetra = vertex_star(i)
            affected(tetra) = .true.
            if (any(self%tetrahedra(:,tetra) == new_update%w)) then
                new_update%u_minus = [new_update%u_minus, tetra]
            else
                new_update%u_plus = [new_update%u_plus, tetra]
            end if
        end do
        
        ! Store update
        call self%add_update(new_update)
        
        write(message, '(A,I0,A,I0,A,I0)') "Stored update ", new_update%id, &
            ": remove ", size(new_update%u_minus), " modify ", size(new_update%u_plus)
        call report(trim(message), 1)
    end subroutine

    subroutine add_update(self, update)
        class(MultiTesselation), intent(inout) :: self
        type(UpdateNode), intent(in) :: update
        type(UpdateNode), pointer :: new_update
        
        ! Allocate new update node
        allocate(new_update)
        new_update = update
        new_update%next => null()
        new_update%prev => null()
        
        if (.not. associated(self%first_update)) then
            ! First update in list
            self%first_update => new_update
            self%last_update => new_update
        else
            ! Add to end of list
            new_update%prev => self%last_update
            self%last_update%next => new_update
            self%last_update => new_update
        end if
        
        self%n_updates = self%n_updates + 1
    end subroutine

    subroutine apply_updates(self, n)
        class(MultiTesselation), intent(inout) :: self
        integer(int64), intent(in) :: n
        
        type(UpdateNode), pointer :: curr_update
        integer(int64) :: i, j, k, count
        logical, allocatable :: valid_tetra(:), valid_vertex(:), valid_edge(:)
        integer(int64), allocatable :: new_tetra_indices(:), new_vertex_indices(:), new_edge_indices(:)
        integer(int64), allocatable :: temp_tetra(:,:), temp_vertices(:,:), temp_edges(:,:)
        integer(int64), allocatable :: temp_adj(:,:), temp_faces(:,:)
        character(len=4096) :: message
        integer(int64) :: orig_vertices, orig_tetra, orig_edges

        ! Store original counts
        orig_vertices = self%n_vertices
        orig_tetra = self%n_tetrahedra
        orig_edges = self%n_edges

        if (n <= 0 .or. n > self%n_updates) then
            write(message, '(A,I0,A)') "Invalid number of updates: ", n
            call report(trim(message), 3)
            return
        endif

        ! Initialize validity masks
        allocate(valid_tetra(self%n_tetrahedra))
        allocate(valid_vertex(self%n_vertices))
        allocate(valid_edge(self%n_edges))
        valid_tetra = .true.
        valid_vertex = .true.
        valid_edge = .true.

        ! Apply updates sequentially
        curr_update => self%first_update
        count = 0

        do while (associated(curr_update) .and. count < n)
            ! Mark tetrahedra for removal
            do i = 1, size(curr_update%u_minus)
                valid_tetra(curr_update%u_minus(i)) = .false.
            end do

            ! Mark vertex being collapsed for removal
            valid_vertex(curr_update%v) = .false.

            ! Update tetrahedra vertex references
            do i = 1, size(curr_update%u_plus)
                k = curr_update%u_plus(i)
                where (self%tetrahedra(:,k) == curr_update%v)
                    self%tetrahedra(:,k) = curr_update%w
                end where
            end do

            count = count + 1
            curr_update => curr_update%next
        end do

        ! Mark edges connected to removed vertices for removal
        do i = 1, self%n_edges
            if (.not.(valid_vertex(self%edges(1,i)) .and. valid_vertex(self%edges(2,i)))) then
                valid_edge(i) = .false.
            end if
        end do

        ! Create new vertex indexing
        allocate(new_vertex_indices(self%n_vertices))
        j = 0
        do i = 1, self%n_vertices
            if (valid_vertex(i)) then
                j = j + 1
                new_vertex_indices(i) = j
            else
                new_vertex_indices(i) = 0
            end if
        end do

        ! Update edge connectivity with new vertex indices
        do i = 1, self%n_edges
            if (valid_edge(i)) then
                self%edges(:,i) = new_vertex_indices(self%edges(:,i))
            end if
        end do

        ! Create new edge indexing
        allocate(new_edge_indices(self%n_edges))
        k = 0
        do i = 1, self%n_edges
            if (valid_edge(i)) then
                k = k + 1
                new_edge_indices(i) = k
            else
                new_edge_indices(i) = 0
            end if
        end do

        ! Create new tetrahedra indexing
        allocate(new_tetra_indices(self%n_tetrahedra))
        k = 0
        do i = 1, self%n_tetrahedra
            if (valid_tetra(i)) then
                k = k + 1
                new_tetra_indices(i) = k
            else
                new_tetra_indices(i) = 0
            end if
        end do

        ! Allocate temporary arrays
        allocate(temp_vertices(3, j))
        allocate(temp_edges(2, maxval(new_edge_indices)))  ! Use maxval here
        allocate(temp_tetra(4, k))
        allocate(temp_adj(4, k))
        allocate(temp_faces(4, k))

        ! Pack vertices
        j = 0
        do i = 1, self%n_vertices
            if (valid_vertex(i)) then
                j = j + 1
                temp_vertices(:,j) = self%vertices(:,i)
            end if
        end do

        ! Pack edges
        k = 0
        do i = 1, self%n_edges
            if (valid_edge(i)) then
                k = k + 1
                temp_edges(:,k) = self%edges(:,i)
            end if
        end do

        ! Pack tetrahedra arrays
        k = 0
        do i = 1, self%n_tetrahedra
            if (valid_tetra(i)) then
                k = k + 1
                temp_tetra(:,k) = new_vertex_indices(self%tetrahedra(:,i))
                temp_adj(:,k) = self%tetrahedron_adjacency(:,i)
                temp_faces(:,k) = self%tetrahedron_faces(:,i)
            end if
        end do

        ! Deallocate and reallocate arrays
        deallocate(self%vertices, self%edges)
        deallocate(self%tetrahedra)
        deallocate(self%tetrahedron_adjacency)
        deallocate(self%tetrahedron_faces)

        allocate(self%vertices(3, size(temp_vertices,2)))
        allocate(self%edges(2, size(temp_edges,2)))
        allocate(self%tetrahedra(4, size(temp_tetra,2)))
        allocate(self%tetrahedron_adjacency(4, size(temp_adj,2)))
        allocate(self%tetrahedron_faces(4, size(temp_faces,2)))

        ! Copy back updated arrays
        self%vertices = temp_vertices
        self%edges = temp_edges
        self%tetrahedra = temp_tetra
        self%tetrahedron_adjacency = temp_adj
        self%tetrahedron_faces = temp_faces

        ! Update counts
        self%n_vertices = size(self%vertices, 2)
        self%n_edges = size(self%edges, 2)
        self%n_tetrahedra = size(self%tetrahedra, 2)

        ! Update adjacency references in tetrahedra
        do i = 1, self%n_tetrahedra
            where (self%tetrahedron_adjacency(:,i) > 0)
                self%tetrahedron_adjacency(:,i) = new_tetra_indices(self%tetrahedron_adjacency(:,i))
            end where
        end do

        ! ! Update edge_boundary if necessary
        ! if (allocated(self%edge_boundary)) then
        !     self%edge_boundary = self%edge_boundary(valid_edge)
        ! end if

        ! ! Update face_boundary if necessary
        ! if (allocated(self%face_boundary)) then
        !     ! Implement similar updates for faces if needed
        ! end if

        ! Report the changes
        write(message, '(A,I0,A,I0,A,I0,A)') "Applied ", count, " updates, removed ", &
            orig_vertices - self%n_vertices, " vertices and ", orig_tetra - self%n_tetrahedra, " tetrahedra"
        call report(trim(message), 1)

    end subroutine

end module tetrahedral_mesh
