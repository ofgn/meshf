! ------------------------------------------------------------------------------
! @file tetrahedral.f90
! @brief Module for handling tetrahedral mesh operations
! @author ofgn
! @date 2024-10-03
! ------------------------------------------------------------------------------
module tetrahedral_mesh
    use iso_fortran_env, only: int32, int64, real32, real64
    use utility
    use vtk
    use data_structures

    implicit none

    type :: TetrahedralMesh
        integer(int64) :: n_vertices                                            ! Number of vertices in the grid
        integer(int64) :: n_tetrahedra                                          ! Number of tetrahedra in the grid
        integer(int64) :: n_edges                                               ! Number of edges in the grid
        real(real64), allocatable :: vertices(:, :)                             ! Coordinates of each vertex (3 x n_vertices)
        integer(int64), allocatable :: tetrahedra(:, :)                         ! Cell connectivity (4 x n_tetra)
        logical, allocatable :: boundary(:)                                     ! Boundary vertices
        type(HashMap) :: adjacency_map                                          ! Hash map for O(1) cell adjacency lookups
        type(LinkedList), allocatable :: adjacency_list(:)                      ! Array of linked lists for O(1) vertex adjacency lookups
        type(PriorityQueue) :: edge_queue                                       ! Priority queue for edge collapses
    contains
        procedure :: initialise => initialise_tetrahedral_mesh                  ! Initialise the tetrahedral mesh
        procedure :: vertex_neighbours                                          ! Get neighbours of a vertex
        procedure :: tetrahedron_neighbours                                     ! Get neighbours of a tetrahedron
        procedure :: build_edge_queue                                           ! Build the edge queue
        procedure :: export_u_grid                                              ! Export to unstructured grid
    end type TetrahedralMesh

    type :: MultiTessellation
        integer(int64) :: count
        type(Update),  pointer :: head => null()
        type(Update),  pointer :: tail => null()
    contains
        procedure :: initialise => mt_initialise
        procedure :: add_update => mt_add_update
    end type MultiTessellation

    type :: Update
        integer(int64) :: id                                                         
        integer(int64) :: removed_vertex = 0
        integer(int64) :: n_created_tetrahedra = 0
        integer(int64) :: n_removed_tetrahedra = 0
        
        type(Update),  pointer :: previous => null()
        type(Update),  pointer :: next => null()                                            
    end type Update

contains

    ! --------------------------------------------------------------------------
    ! @brief Initialise the tetrahedral mesh using a VtkUnstructuredGrid.
    ! @param[inout] self The tetrahedral mesh to initialise.
    ! @param[in] u_grid The unstructured grid containing vertices and 
    !   connectivity.
    ! @details Copies vertices and cell connectivity from the unstructured grid 
    !   to set up the mesh.
    ! --------------------------------------------------------------------------
    subroutine initialise_tetrahedral_mesh(self, u_grid)
        use geometry
        
        implicit none

        class(TetrahedralMesh), intent(inout) :: self                           ! Tetrahedral mesh structure
        type(VtkUnstructuredGrid), intent(in) :: u_grid                         ! Input unstructured grid

        integer(int64) :: i                                                     ! Loop variable
        integer(int64) :: u                                                     ! Vertex u
        integer(int64) :: v                                                     ! Vertex v
        integer(int64) :: w                                                     ! Vertex w
        integer(int64) :: x                                                     ! Vertex x
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        real(real128) :: orientation                                            ! Orientation of tetrahedron
        character(len=4096) :: message                                          ! String buffer


        call cpu_time(start_time)
        call report("Initialising tetrahedral mesh", 1)

        self%n_vertices = u_grid%n_points
        self%n_tetrahedra = u_grid%n_cells

        write(message, '(A, I0)') "Vertices: ", self%n_vertices
        call report(trim(message), 1)
        write(message, '(A, I0)') "Tetrahedra: ", self%n_tetrahedra
        call report(trim(message), 1)

        allocate (self%vertices(3, self%n_vertices))
        allocate (self%tetrahedra(4, self%n_tetrahedra))

        self%vertices = u_grid%points
        self%tetrahedra = reshape(u_grid%cell_connectivity, &
            [4_int64, self%n_tetrahedra]) + 1

        call self%adjacency_map%initialise(4 * self%n_tetrahedra)
        allocate(self%adjacency_list(self%n_vertices))

        ! Populate adjacency map with each face of tetrahedra
        do i = 1, self%n_tetrahedra

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
    
            ! Insert faces into adjacency map
            call self%adjacency_map%insert([u, v, w], [i])
            call self%adjacency_map%insert([u, x, v], [i])
            call self%adjacency_map%insert([u, w, x], [i])
            call self%adjacency_map%insert([v, x, w], [i])
    
            ! Set vertices into adjacency list
            call self%adjacency_list(u)%append(i)
            call self%adjacency_list(v)%append(i)
            call self%adjacency_list(w)%append(i)
            call self%adjacency_list(x)%append(i)
        end do

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine initialise_tetrahedral_mesh

    ! --------------------------------------------------------------------------
    ! @brief Get the neighbours of a vertex.
    ! @param[in] self The tetrahedral mesh.
    ! @param[in] vertex The vertex to find neighbours for.
    ! @return The neighbours of the vertex.
    ! --------------------------------------------------------------------------
    function vertex_neighbours(self, vertex) result(neighbours)
        implicit none

        class(TetrahedralMesh), intent(in) :: self
        integer(int64), intent(in) :: vertex

        integer(int64), allocatable :: neighbours(:)

        neighbours = self%adjacency_list(vertex)%to_array()

    end function vertex_neighbours

    ! --------------------------------------------------------------------------
    ! @brief Get the neighbours of a tetrahedron.
    ! @param[in] self The tetrahedral mesh.
    ! @param[in] tetrahedron The tetrahedron to find neighbours for.
    ! @return The neighbours of the tetrahedron.
    ! --------------------------------------------------------------------------
    function tetrahedron_neighbours(self, tetrahedron) result(neighbours)
        implicit none
    
        class(TetrahedralMesh), intent(in) :: self
        integer(int64), intent(in) :: tetrahedron
    
        integer(int64) :: u, v, w, x                    ! Vertex indices
        integer(int64), allocatable :: value(:)         ! For map lookup results
        integer(int64), allocatable :: neighbours(:)    ! Result array
        integer(int64) :: face_idx                      ! Current face index
        integer(int64), dimension(3,6) :: face_order    ! All permutations per face
        integer(int64) :: i                            ! Loop counter
        logical :: found = .false.                     ! Found neighbour flag
    
        ! Get tetrahedron vertices
        u = self%tetrahedra(1, tetrahedron)
        v = self%tetrahedra(2, tetrahedron)
        w = self%tetrahedra(3, tetrahedron)
        x = self%tetrahedra(4, tetrahedron)
    
        allocate(neighbours(4))
        neighbours = 0
    
        ! Define all possible vertex permutations for each face
        ! Face 1 (u,v,w)
        face_order(:,1) = [u,v,w]
        face_order(:,2) = [u,w,v]
        face_order(:,3) = [v,w,u]
        face_order(:,4) = [v,u,w]
        face_order(:,5) = [w,u,v]
        face_order(:,6) = [w,v,u]
    
        ! Check each permutation for face 1
        do i = 1, 6
            if (self%adjacency_map%find(face_order(:,i), value)) then
                if (value(1) /= tetrahedron) then
                    neighbours(1) = value(1)
                    exit
                end if
            end if
        end do
    
        ! Face 2 (u,v,x)
        face_order(:,1) = [u,x,v]
        face_order(:,2) = [u,v,x]
        face_order(:,3) = [v,x,u]
        face_order(:,4) = [v,u,x]
        face_order(:,5) = [x,u,v]
        face_order(:,6) = [x,v,u]
    
        ! Check each permutation for face 2
        do i = 1, 6
            if (self%adjacency_map%find(face_order(:,i), value)) then
                if (value(1) /= tetrahedron) then
                    neighbours(2) = value(1)
                    exit
                end if
            end if
        end do
    
        ! Face 3 (u,w,x)
        face_order(:,1) = [u,w,x]
        face_order(:,2) = [u,x,w]
        face_order(:,3) = [w,x,u]
        face_order(:,4) = [w,u,x]
        face_order(:,5) = [x,u,w]
        face_order(:,6) = [x,w,u]
    
        ! Check each permutation for face 3
        do i = 1, 6
            if (self%adjacency_map%find(face_order(:,i), value)) then
                if (value(1) /= tetrahedron) then
                    neighbours(3) = value(1)
                    exit
                end if
            end if
        end do
    
        ! Face 4 (v,w,x)
        face_order(:,1) = [v,x,w]
        face_order(:,2) = [v,w,x]
        face_order(:,3) = [x,w,v]
        face_order(:,4) = [x,v,w]
        face_order(:,5) = [w,v,x]
        face_order(:,6) = [w,x,v]
    
        ! Check each permutation for face 4
        do i = 1, 6
            if (self%adjacency_map%find(face_order(:,i), value)) then
                if (value(1) /= tetrahedron) then
                    neighbours(4) = value(1)
                    exit
                end if
            end if
        end do
    end function tetrahedron_neighbours

    ! --------------------------------------------------------------------------
    ! @brief Build the edge queue for the tetrahedral mesh.
    ! @param[inout] self The tetrahedral mesh.
    ! @details Populates the edge queue with unique edges from each tetrahedron.
    ! --------------------------------------------------------------------------
    subroutine build_edge_queue(self)
        implicit none
    
        class(TetrahedralMesh), intent(inout) :: self                           ! The tetrahedral mesh

        integer(int64) :: i                                                     ! Loop variable
        integer(int64) :: u                                                     ! Vertex u
        integer(int64) :: v                                                     ! Vertex v
        integer(int64) :: w                                                     ! Vertex w
        integer(int64) :: x                                                     ! Vertex x               
        integer(int64), allocatable :: edge(:)                                  ! Edge
        integer(int64), allocatable :: edges(:, :)                              ! Edges
        integer(int64), allocatable :: value(:)                                 ! For map lookup results 
        real(real64) :: uv, uw, ux, vw, vx, wx                                  ! Edge priority
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Run time
        character(len=4096) :: message                                          ! String buffer
        type(HashMap) :: edge_map                                               ! Map of edges
    
        call cpu_time(start_time)
        call report("Building edge queue", 1)

        allocate(edge(2))
        call edge_map%initialise(6 * self%n_tetrahedra)

        ! Populate edge queue with unique edges from each tetrahedron
        do i = 1, self%n_tetrahedra
            u = self%tetrahedra(1, i)
            v = self%tetrahedra(2, i)
            w = self%tetrahedra(3, i)
            x = self%tetrahedra(4, i)

            edge = [min(u,v), max(u,v)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if

            edge = [min(u,w), max(u,w)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if

            edge = [min(u,x), max(u,x)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if

            edge = [min(v,w), max(v,w)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if

            edge = [min(v,x), max(v,x)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if

            edge = [min(w,x), max(w,x)]
            if (.not. edge_map%find(edge, value)) then
                call edge_map%insert(edge, [i])
            end if
        end do

        self%n_edges = edge_map%count

        edges = edge_map%keys()
        call edge_map%reset()
        call self%edge_queue%initialise(self%n_edges, reverse_priority=.true.)

        do i = 1, self%n_edges
            if ((mod(i, self%n_edges / 5)) == 0) then
                write(message, '(A, I0, A, I0)') "Inserting edge ", &
                    i, " of ", self%n_edges
                call report(trim(message), 1)
            end if
            call self%edge_queue%insert(edges(:, i), real(i, real64))
        end do

        deallocate(edge)
        
        write(message, '(A, I0)') "Number of unique edges: ", self%n_edges
        call report(trim(message), 1)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine build_edge_queue
    
    ! --------------------------------------------------------------------------
    ! @brief Export the tetrahedral mesh to an unstructured grid.
    ! @param[in] self The tetrahedral mesh to export.
    ! @param[out] u_grid The resulting unstructured grid.
    ! @details Converts the internal representation of vertices and tetrahedra 
    !   to the VtkUnstructuredGrid.
    ! --------------------------------------------------------------------------
    function export_u_grid(self) result(u_grid)
        implicit none

        class(TetrahedralMesh), intent(in) :: self                              ! The tetrahedral mesh

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
    end function export_u_grid

    ! --------------------------------------------------------------------------
    ! @brief Initializes a multi-tessellation structure.
    ! @param[inout] self The multi-tessellation instance.
    ! @param[in] n_nodes Number of nodes in tessellation.
    ! --------------------------------------------------------------------------
    subroutine mt_initialise(self, n_nodes)
        implicit none

        class(MultiTessellation), intent(inout) :: self
        integer(int64), intent(in) :: n_nodes

        integer(int64) :: i
        integer(int32) :: alloc_status
        
        ! Allocate nodes array
        if (allocated(self%nodes)) deallocate(self%nodes)
        allocate(self%nodes(n_nodes), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate multi-tessellation nodes", 3)
            return
        end if
        
        ! Initialize basic properties
        self%n_nodes = n_nodes
        
        ! Initialize each node
        do i = 1, n_nodes
            self%nodes(i)%id = i
            self%nodes(i)%remaining_vertex = 0
            self%nodes(i)%removed_vertex = 0
        end do
    end subroutine mt_initialise

    ! --------------------------------------------------------------------------
    ! @brief Adds or updates a node in the multi-tessellation.
    ! @param[inout] self The multi-tessellation instance.
    ! @param[in] node_id Node identifier.
    ! @param[in] remaining_vertex Vertex that remains after collapse.
    ! @param[in] removed_vertex Vertex removed in collapse.
    ! @param[in] ancestors Array of ancestor node IDs.
    ! --------------------------------------------------------------------------
    subroutine mt_add_update(self, remaining_vertex, removed_vertex, ancestors)
        implicit none

        class(MultiTessellation), intent(inout) :: self
        integer(int64), intent(in) :: node_id
        integer(int64), intent(in) :: remaining_vertex
        integer(int64), intent(in) :: removed_vertex
        integer(int64), intent(in) :: ancestors(:)

        integer(int64) :: i
        
        ! Validate input
        if (node_id < 1 .or. node_id > self%n_nodes) then
            call report("Invalid node ID in multi-tessellation update", 3)
            return
        end if
        
        ! Update node properties
        self%nodes(node_id)%remaining_vertex = remaining_vertex
        self%nodes(node_id)%removed_vertex = removed_vertex
        
        ! Create loop for dependencies
        call self%nodes(node_id)%loop%append(node_id)
        do i = 1, size(ancestors)
            call self%nodes(node_id)%loop%append(ancestors(i))
        end do
    end subroutine mt_add_update

end module tetrahedral_mesh
