! ------------------------------------------------------------------------------
! @file mra.f90
! @brief Module for handling tetrahedral mesh operations
! @author ofgn
! @date 2024-1003
! ------------------------------------------------------------------------------
module multi_resolution_analysis
    use global

    implicit none

    type :: TetrahedralMesh
        integer(int64) :: n_updates = 0                                         ! Number of updates              
        integer(int64) :: n_vertices = 0                                        ! Number of vertices in the grid
        integer(int64) :: n_tetrahedra = 0                                      ! Number of tetrahedra in the grid
        integer(int64), allocatable :: tetrahedra(:, :)                         ! Cell connectivity (4 x n_tetra)
        integer(int64), allocatable :: vertex_adjacency(:)                      ! A single tetrahedron adjacent to each vertex
        integer(int64), allocatable :: tetrahedron_adjacency(:, :)              ! Tetrahedra adjacent to each tetrahedron
        real(real64), allocatable :: vertices(:, :)                             ! Coordinates of each vertex (3 x n_vertices)
        type(UpdateNode), pointer :: first_update => null()                     ! Head of updates list
        type(UpdateNode), pointer :: last_update => null()                      ! Tail of updates list
    contains
        procedure :: read_tetgen_node_file                                      ! Read node file
        procedure :: read_tetgen_ele_file                                       ! Read ele file                
        procedure :: read_tetgen_neigh_file                                     ! Read neigh file
        procedure :: simplify_mesh                                              ! Apply half-edge collapse to compress the mesh
        procedure :: vtk_unstructured_grid                                      ! Export the mesh to a VtkUnstructuredGrid         
        procedure :: find_vertex_star                                           ! Find all tetrahedra containing a vertex
        procedure :: half_edge_collapse                                         ! Collapse an edge in the mesh
        procedure :: add_update
        procedure :: apply_updates
        procedure, private :: convex_hull
        procedure, private :: edge_error
        !procedure :: wavelet_compression
        procedure, private :: finalise
    end type TetrahedralMesh

    type :: MultiResolution
        type(TetrahedralMesh) :: base_mesh                                      ! Base mesh
    end type MultiResolution

    type :: UpdateNode
        integer(int64) :: id                                                    ! Unique node ID
        integer(int64) :: u                                                     ! Vertex being collapsed 
        integer(int64) :: v                                                     ! Target vertex
        integer(int64), allocatable :: u_minus(:)                               ! Tetrahedra removed
        integer(int64), allocatable :: u_plus(:)                                ! Tetrahedra modified
        integer(int64), allocatable :: definitions(:)                           ! Tetrahedra containing u
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

        class(TetrahedralMesh), intent(inout) :: self                           ! TetrahedralMesh type
        character(len=*), intent(in) :: file_path                               ! Path to node file
      
        integer(int64) :: i, j                                                  ! Loop index          
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        real(real64) :: start_time, end_time, run_time                          ! Timing variables                                             
        character(len=256) :: message                                           ! String buffer

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
    ! @brief Read a TetGen ele file.
    ! @param[inout] self The tetrahedral mesh to populate.
    ! @param[in] file_path The path to the ele file.
    ! @details Reads the ele file and populates the tetrahedra of the mesh.
    ! --------------------------------------------------------------------------
    subroutine read_tetgen_ele_file(self, file_path)
        use geometry, only: signed_volume
        use data_structures_algorithms, only: sort4, swap

        implicit none

        class(TetrahedralMesh), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=*), intent(in) :: file_path                               ! Path to face file
      
        integer(int64) :: i, j                                                  ! Loop index
        integer(int64) :: u, v, w, x                                            ! Vertex indices
        integer(int64) :: temp                                                  ! Temporary variable        
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        real(real64) :: orientation                                             ! Orientation of tetrahedron
        character(len=256) :: message                                           ! String buffer

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

        write(message, "(A, I0)") "Number of tetrahedra: ", self%n_tetrahedra
        call report(trim(message), 1)

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
        end do
        
        !$omp parallel do private(u, v, w, x, orientation)
        do i = 1, self%n_tetrahedra
            ! Sort vertices in ascending order
            call sort4(self%tetrahedra(:, i))

            u = self%tetrahedra(1, i)
            v = self%tetrahedra(2, i)
            w = self%tetrahedra(3, i)
            x = self%tetrahedra(4, i)

            ! Check orientation
            orientation = signed_volume(self%vertices(:, [u, v, w, x]))

            ! If orientation is negative, swap the last two vertices
            if (orientation < 0.0_real64) then
                call swap(w, x)    
                self%tetrahedra(:, i) = [u, v, w, x]
            end if

            self%vertex_adjacency(u) = i
            self%vertex_adjacency(v) = i
            self%vertex_adjacency(w) = i
            self%vertex_adjacency(x) = i
        end do
        !$omp end parallel do

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_ele_file

    subroutine read_tetgen_neigh_file(self, file_path)
        implicit none

        class(TetrahedralMesh), intent(inout) :: self                           ! Tetrahedral mesh structure
        character(len=*), intent(in) :: file_path                               ! Path to face file
      
        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index
        integer(int64) :: u, v, w, x                                            ! Vertex indices                 
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=256) :: message                                           ! String buffer

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

            where (self%tetrahedron_adjacency(:, j) == -1)
                self%tetrahedron_adjacency(:, j) = 0
            end where
        end do

        close (unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine read_tetgen_neigh_file

    function vtk_unstructured_grid(self)
        use vtk

        implicit none
    
        class(TetrahedralMesh), intent(in) :: self                              ! The tetrahedral mesh
    
        type(VtkUnstructuredGrid) :: vtk_unstructured_grid                      ! Output unstructured grid

        integer(kind=int64) :: i                                                ! Loop variable
        integer(kind=int64), allocatable :: flat_tetrahedra(:)                  ! Flattened tetrahedra array
        integer(kind=int64), allocatable :: packed_tetrahedra(:)                ! Packed tetrahedra array
        real(real64) :: start_time, end_time, run_time                          ! Timing variables                                              
        logical(bool8), allocatable :: mask(:)                                  ! Mask for tetrahedra
        logical(bool8), allocatable :: flat_mask(:)                             ! Flattened mask
        character(len=256) :: message                                           ! String buffer
    
        call cpu_time(start_time)
    
        write(message, '(A)') "Converting TetrahedralMesh to " &
            // "VtkUnstructuredGrid"
        call report(trim(message), 1)
    
        vtk_unstructured_grid%n_points = self%n_vertices
        mask = .not. all(self%tetrahedra == 0, dim=1)
    
        vtk_unstructured_grid%n_cells = count(mask)
    
        allocate &
            (vtk_unstructured_grid%points(3, vtk_unstructured_grid%n_points))
        allocate &
            (flat_tetrahedra(size(self%tetrahedra)))
        allocate &
            (flat_mask(size(flat_tetrahedra)))
        allocate &
            (vtk_unstructured_grid%cell_connectivity(4 * vtk_unstructured_grid%n_cells))
        allocate &
            (vtk_unstructured_grid%cell_type(vtk_unstructured_grid%n_cells))
        allocate &
            (vtk_unstructured_grid%cell_size(vtk_unstructured_grid%n_cells))
    
        ! Flatten the tetrahedra array
        flat_tetrahedra = reshape(self%tetrahedra, [size(self%tetrahedra)])
    
        ! Create a flattened mask
        flat_mask = reshape(spread(mask, dim=1, ncopies=4), [size(flat_tetrahedra)])
    
        ! Pack the flattened tetrahedra array
        packed_tetrahedra = pack(flat_tetrahedra, flat_mask)
    
        vtk_unstructured_grid%points = self%vertices
        vtk_unstructured_grid%cell_connectivity = packed_tetrahedra - 1
        vtk_unstructured_grid%cell_type = 10
        vtk_unstructured_grid%cell_size = 4
    
        write(message, '(A, I0)') "Number of points: ", vtk_unstructured_grid%n_points
        call report(trim(message), 1)
        write(message, '(A, I0)') "Number of cells: ", vtk_unstructured_grid%n_cells
        call report(trim(message), 1)
    
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A, F9.3, A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1) 
    end function vtk_unstructured_grid

    ! --------------------------------------------------------------------------
    ! @brief Find all tetrahedra containing a vertex using depth-first search.
    ! @param[in] self The tetrahedral mesh.
    ! @param[in] vertex_index The vertex to search from.
    ! @return vertex_star Array of tetrahedra IDs containing the vertex.
    ! --------------------------------------------------------------------------
    function find_vertex_star(self, vertex_index) result(vertex_star)
        use data_structures_algorithms

        implicit none
        
        class(TetrahedralMesh), intent(in) :: self
        integer(int64), intent(in) :: vertex_index
        integer(int64), allocatable :: vertex_star(:)
        
        logical(bool8) :: visited(self%n_tetrahedra)                            ! Stack allocated
        integer(int64) :: stack(self%n_tetrahedra)                              ! Fixed size stack
        integer(int64) :: stack_top                                             ! Stack pointer
        integer(int64) :: result_top                                            ! Result pointer
        integer(int64) :: start                                                 ! Starting tetrahedron
        integer(int64) :: current                                               ! Current tetrahedron
        integer(int64) :: next                                                  ! Next tetrahedron
        integer(int64) :: i                                                     ! Loop index
        
        ! Fast initialisation
        visited = .false.
        result_top = 0
        
        ! Get starting tetrahedron
        start = self%vertex_adjacency(vertex_index)
        if (start <= 0 .or. start > self%n_tetrahedra) then
            allocate(vertex_star(0))
            return
        end if
        
        ! Initialise stack
        stack_top = 1
        stack(1) = start
        visited(start) = .true.
        
        ! Preallocate result (will trim later)
        allocate(vertex_star(self%n_tetrahedra))
        vertex_star(1) = start
        result_top = 1
        
        ! Fast DFS traversal
        do while (stack_top > 0)
            current = stack(stack_top)
            stack_top = stack_top - 1
            
            ! Check all 4 neighbors inline
            do i = 1, 4
                next = self%tetrahedron_adjacency(i, current)
                if (next > 0) then
                    if (.not. visited(next)) then
                        if (any(self%tetrahedra(:,next) == vertex_index)) then
                            stack_top = stack_top + 1
                            stack(stack_top) = next
                            visited(next) = .true.
                            result_top = result_top + 1
                            vertex_star(result_top) = next
                        end if
                    end if
                end if
            end do
        end do
        
        ! Trim result to actual size
        vertex_star = vertex_star(1:result_top)
        call quicksort(vertex_star)
    end function find_vertex_star

    ! --------------------------------------------------------------------------
    ! @brief Apply half-edge collapse to simplify the mesh.
    ! @param[in] self The tetrahedral mesh.
    ! @param[in] v The vertex to collapse.
    ! @param[in] w The target vertex.
    ! @param[in] error The error threshold.
    ! --------------------------------------------------------------------------
    subroutine half_edge_collapse(self, u, v, error)
        implicit none
    
        class(TetrahedralMesh), intent(inout) :: self
        integer(int64), intent(in) :: u
        integer(int64), intent(in) :: v
        real(real64), intent(in) :: error
    
        type(UpdateNode) :: new_update
        integer(int64), allocatable :: vertex_star(:)
        integer(int64) :: i, tetra
    
        ! Initialise update node
        new_update%id = self%n_updates + 1
        new_update%u = u  ! Vertex to collapse
        new_update%v = v  ! Target vertex
        new_update%error = error
    
        ! Find affected tetrahedra
        vertex_star = self%find_vertex_star(u)
        
        ! Classify tetrahedra only - no modification
        allocate(new_update%u_minus(0))
        allocate(new_update%u_plus(0))
        
        do i = 1, size(vertex_star)
            tetra = vertex_star(i)
            if (any(self%tetrahedra(:,tetra) == v)) then
                ! Will be removed - contains both u and v
                new_update%u_minus = [new_update%u_minus, tetra]
            else
                ! Will be modified - contains only u
                new_update%u_plus = [new_update%u_plus, tetra]
            end if
        end do
    
        ! Store update
        call self%add_update(new_update)
    end subroutine half_edge_collapse
    

    ! --------------------------------------------------------------------------
    ! @brief Add an update node to the list.
    ! @param[inout] self The tetrahedral mesh.
    ! @param[in] update The update node to add.
    ! --------------------------------------------------------------------------
    subroutine add_update(self, update)
        class(TetrahedralMesh), intent(inout) :: self
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

    ! --------------------------------------------------------------------------
    ! @brief Apply a number of updates to the mesh.
    ! @param[inout] self The tetrahedral mesh.
    ! @param[in] n The number of updates to apply.
    ! --------------------------------------------------------------------------
    subroutine apply_updates(self, n)
        use data_structures_algorithms, only: sort4, swap
        use geometry, only: signed_volume

        implicit none
    
        class(TetrahedralMesh), intent(inout) :: self
        integer(int64), optional, intent(in) :: n
    
        integer(int64) :: n_iter
        integer(int64) :: i, j
        integer(int64) :: u, v, w, x
        integer(int64) :: updates_processed
        integer(int64) :: old_vertex, new_vertex
        integer(int64), allocatable :: vertex_map(:) 
        integer(int64), allocatable :: new_vertex_index(:)
        integer(int64) :: new_index
        integer(int64), allocatable :: flat_tetrahedra(:)
        real(real64), allocatable :: flat_vertices(:)
        logical(bool8), allocatable :: valid_vertex(:)
        logical(bool8), allocatable :: valid_tetra(:)
        logical(bool8), allocatable :: mask(:)
    
        real(real64) :: start_time, end_time, run_time
        character(len=256) :: message
        type(UpdateNode), pointer :: current_update
    
        call cpu_time(start_time)

        if (.not. associated(self%first_update)) then
            write(message, '(A)') "No updates to apply"
            call report(trim(message), 3)
            return
        end if
    
        if (present(n)) then
            n_iter = min(n, self%n_updates)
        else
            n_iter = self%n_updates
        end if
    
        if (n_iter <= 0) then
            write(message, '(A)') "No updates requested"
            call report(trim(message), 3)
            return
        end if
    
        write(message, '(A,I0,A)') "Applying ", n_iter, " updates"
        call report(trim(message), 1)
    
        allocate(valid_vertex(self%n_vertices))
        allocate(valid_tetra(self%n_tetrahedra))
    
        valid_vertex = .true.
        valid_tetra  = .true.
    
        allocate(vertex_map(self%n_vertices))
        do i = 1, self%n_vertices
            vertex_map(i) = i
        end do

        current_update => self%first_update
        updates_processed = 0
    
        do while (associated(current_update) .and. updates_processed < n_iter)
            ! Mark the collapsed vertex as mapping -> the final vertex
            vertex_map(current_update%u) = current_update%v
            valid_vertex(current_update%u) = .false.
    
            ! Mark all tetrahedra in u_minus as invalid
            valid_tetra(current_update%u_minus) = .false.
    
            updates_processed = updates_processed + 1
            current_update => current_update%next
        end do
    
        ! Count valid vertices and tetrahedra
        do i = 1, self%n_vertices
            if (.not. valid_vertex(i)) then
                do while (vertex_map(i) /= vertex_map(vertex_map(i)))
                    vertex_map(i) = vertex_map(vertex_map(i))
                end do
            end if
        end do

        do i = 1, self%n_tetrahedra
            if (valid_tetra(i)) then
                do j = 1, 4
                    ! Map each vertex of tetrahedron i
                    self%tetrahedra(j, i) = vertex_map(self%tetrahedra(j, i))
                end do
            end if
        end do
    
        allocate(flat_vertices(size(self%vertices)))
        flat_vertices = reshape(self%vertices, [size(self%vertices)])
    
        allocate(mask(size(flat_vertices)))
        mask = reshape(spread(valid_vertex, dim=1, ncopies=3), [size(flat_vertices)])
    
        ! Perform pack
        self%vertices = reshape(pack(flat_vertices, mask), [3, count(valid_vertex)])
        deallocate(flat_vertices, mask)
    
        allocate(new_vertex_index(self%n_vertices))
        new_vertex_index = 0
    
        new_index = 0
        do i = 1, self%n_vertices
            if (valid_vertex(i)) then
                new_index = new_index + 1
                new_vertex_index(i) = new_index
            else
                new_vertex_index(i) = 0
            end if
        end do
    
        ! Iterate through tetrahedra, reindexing vertices
        do i = 1, self%n_tetrahedra
            if (valid_tetra(i)) then
                do j = 1, 4
                    old_vertex = self%tetrahedra(j, i)
                    ! Convert from old_vertex ID -> new_vertex_index
                    new_vertex = new_vertex_index(old_vertex)
    
                    ! If it's 0, that means this tetra references a removed vertex
                    if (new_vertex == 0) then
                        valid_tetra(i) = .false.
                        exit
                    end if
    
                    self%tetrahedra(j, i) = new_vertex
                end do
            end if
        end do
    
        allocate(flat_tetrahedra(size(self%tetrahedra)))
        flat_tetrahedra = reshape(self%tetrahedra, [size(self%tetrahedra)])
    
        allocate(mask(size(flat_tetrahedra)))
        mask = reshape(spread(valid_tetra, dim=1, ncopies=4), [size(flat_tetrahedra)])
    
        self%tetrahedra = reshape(pack(flat_tetrahedra, mask), [4, count(valid_tetra)])
        deallocate(flat_tetrahedra, mask)
    
        self%n_vertices = size(self%vertices, 2)
        self%n_tetrahedra    = size(self%tetrahedra, 2)

        !$omp parallel do private(u, v, w, x)
        do i = 1, self%n_tetrahedra
            ! Sort vertices in ascending order
            call sort4(self%tetrahedra(:, i))

            u = self%tetrahedra(1, i)
            v = self%tetrahedra(2, i)
            w = self%tetrahedra(3, i)
            x = self%tetrahedra(4, i)

            if (signed_volume(self%vertices(:, [u, v, w, x])) < 0.0_real64) then
                call swap(w, x)    
                self%tetrahedra(:, i) = [u, v, w, x]
            end if

            self%vertex_adjacency(u) = i
            self%vertex_adjacency(v) = i
            self%vertex_adjacency(w) = i
            self%vertex_adjacency(x) = i
        end do
        !$omp end parallel do

    
        write(message, '(A, I0)') "Final vertex count: ", self%n_vertices
        call report(trim(message), 1)
        write(message, '(A, I0)') "Final tetra count : ", self%n_tetrahedra 
        call report(trim(message), 1)
    
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A,F9.3,A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(message), 1)
    
        deallocate(vertex_map, new_vertex_index, valid_vertex&
            , valid_tetra)
    end subroutine apply_updates
    
    
    ! --------------------------------------------------------------------------
    subroutine simplify_mesh(self, compression_ratio)
        use data_structures_algorithms
        
        implicit none
        
        class(TetrahedralMesh), intent(inout) :: self
        real(real64), intent(in) :: compression_ratio
        
        integer(int64), allocatable :: current_edge(:)
        integer(int64) :: target_tetrahedra
        integer(int64) :: count
        integer(int64) :: i, j
        integer(int64) :: u, v, w, x
        integer(int64) :: temp
        integer(int64) :: n_removals
        integer(int64) :: n_edges
        integer(int64) :: tetra
        integer(int64) :: new_edges(2, 3)
        integer(int64), allocatable :: edges(:, :)
        integer(int64), allocatable :: boundary(:)
        integer(int64) :: vertex_map(self%n_vertices)
        real(real64) :: error
        real(real64) :: start_time, end_time, run_time
        logical(bool8) :: valid_vertices(self%n_vertices)
        logical(bool8) :: valid_tetrahedra(self%n_tetrahedra)
        character(len=256) :: message
        type(MinHeap) :: edge_queue
        type(Set) :: edge_set

        call cpu_time(start_time)

        write(message, '(A)') "Mesh coarsening"
        call report(trim(message), 1)

        write(message, '(A)') "Finding unique edges"
        call report(trim(message), 1)

        call edge_set%initialise(12 * self%n_tetrahedra)
        
        ! Build edge set
        do i = 1, self%n_tetrahedra
            u = self%tetrahedra(1, i)
            v = self%tetrahedra(2, i)
            w = self%tetrahedra(3, i)
            x = self%tetrahedra(4, i)

            call edge_set%add([u, v])
            call edge_set%add([u, w])
            call edge_set%add([u, x])
            call edge_set%add([v, w])
            call edge_set%add([v, x])
            if (w > x) then
                call edge_set%add([x, w])
            else
                call edge_set%add([w, x])
            end if
        end do

        edges = edge_set%to_array()
        call edge_set%finalise()
        n_edges = size(edges, 2)

        boundary = self%convex_hull()
        valid_vertices = .true.
        valid_tetrahedra = .true.
        
        ! Initialise counts
        count = self%n_tetrahedra

        target_tetrahedra = ceiling(real(self%n_tetrahedra) &
            * compression_ratio, kind=int64)
        
        ! Build edge queue with lengths
        edge_queue = MinHeap()
        call edge_queue%initialise(n_edges)

        write(message, '(A)') "Queueing edges by error"
        call report(trim(message), 1)

        !$omp parallel do 
        do i = 1, n_edges
            if (boundary(u)) then
                if (boundary(v)) then
                    cycle
                else
                    temp = edges(1, i)
                    edges(1, i) = edges(2, i)
                    edges(2, i) = temp
                end if                
            end if
        end do
        !$omp end parallel do

        do i = 1, 3
            if (boundary(u)) then
                if (boundary(v)) then
                    cycle
                end if                
            end if

            call edge_queue%add(edges(:, i), &
                self%edge_error(edges(1, i), edges(2, i)))
        end do

        write(message, '(A)') "Collapsing edges"
        call report(trim(message), 1)

        do while ((count > target_tetrahedra) .and. &
            (edge_queue%count > 0))

            ! Get next edge
            current_edge = edge_queue%min(error)
            u = current_edge(1)
            v = current_edge(2)

            if (boundary(u)) cycle

            ! Check if edge is still valid
            if (.not. valid_vertices(u)) cycle
            if (.not. valid_vertices(v)) cycle

            call self%half_edge_collapse(u, v, error)

            valid_vertices(u) = .false.

            n_removals = 0
            ! Loop through u_minus (removed tetrahedra)
            do j = 1, size(self%last_update%u_minus)
                if (valid_tetrahedra(self%last_update%u_minus(j))) then
                    valid_tetrahedra(self%last_update%u_minus(j)) = .false.
                    n_removals = n_removals + 1
                end if
            end do
            count = count - n_removals
        end do

        write(message, '(A,I0)') "Half-edge collapses: ", self%n_updates
        call report(trim(message), 1)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(message, '(A,F9.3,A)') "Completed in ", &
            run_time, "s" // char(10)
        call report(trim(message), 1)
    end subroutine simplify_mesh

    
    function edge_error(self, u, v)
        implicit none

        class(TetrahedralMesh), intent(in) :: self                              ! TetrahedralMesh type
        integer(int64), intent(in) :: u                                         ! Vertex u
        integer(int64), intent(in) :: v                                         ! Vertex v

        real(real64) :: edge_error                                              ! Error value

        edge_error = norm2(self%vertices(:, u) - self%vertices(:, v))
    end function edge_error

    ! --------------------------------------------------------------------------
    ! @brief Calculate the boundary vertices of the mesh.
    ! @param[inout] self The tetrahedral mesh.
    ! --------------------------------------------------------------------------
    function convex_hull(self)
        implicit none

        class(TetrahedralMesh), intent(inout) :: self                           ! Tetrahedral mesh

        logical(bool8), allocatable :: convex_hull(:)                           ! Convex hull vertices

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index
        integer(int64) :: u                                                     ! Vertex u
        integer(int64) :: v                                                     ! Vertex v
        integer(int64) :: w                                                     ! Vertex w
        real(real64) :: run_time                                                ! Runtime
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        character(len=256) :: message                                           ! String buffer

        ! Start timer
        call cpu_time(start_time)
        write(message, '(A)') "Calculating boundary vertices"
        call report(trim(message), 1)

        ! Clean and allocate boundary array
        if (allocated(convex_hull)) deallocate(convex_hull)
        allocate(convex_hull(self%n_vertices))
        convex_hull = .false.
        
        ! Mark vertices that are part of boundary faces
        !$omp parallel do private(j, u, v, w)
        do i = 1, self%n_tetrahedra
            do j = 1, 4
                if (self%tetrahedron_adjacency(j,i) == 0) then
                    ! Calculate vertices of this face based on local connectivity
                    select case(j)
                        case(1)
                            u = self%tetrahedra(2,i)
                            v = self%tetrahedra(3,i) 
                            w = self%tetrahedra(4,i)
                        case(2)
                            u = self%tetrahedra(1,i)
                            v = self%tetrahedra(3,i)
                            w = self%tetrahedra(4,i)
                        case(3)
                            u = self%tetrahedra(1,i)
                            v = self%tetrahedra(2,i)
                            w = self%tetrahedra(4,i)
                        case(4)
                            u = self%tetrahedra(1,i)
                            v = self%tetrahedra(2,i)
                            w = self%tetrahedra(3,i)
                    end select
                    convex_hull(u) = .true.
                    convex_hull(v) = .true.
                    convex_hull(w) = .true.
                endif
            end do
        end do
        !$omp end parallel do

        write(message, '(A, I0)') "Boundary vertices: ", count(convex_hull)
        call report(trim(message), 1)

        call cpu_time(end_time)
        run_time = end_time - start_time
    end function convex_hull

    ! subroutine wavelet_compression(self, n_iterations)
    !     use geometry

    !     implicit none

    !     class(ExtendedTetrahedralMesh), intent(inout) :: self
    !     integer(int64), intent(in) :: n_iterations

    !     integer(int64) :: i
    !     type(TetrahedralMesh), allocatable :: levels(:)

    !     allocate(levels(n_iterations + 1))

    !     levels(1)%n_vertices = self%n_vertices
    !     levels(1)%n_tetrahedra = self%n_tetrahedra
    !     levels(1)%vertices = self%vertices
    !     levels(1)%tetrahedra = self%tetrahedra

    !     do i = 1, n_iterations
    !         call self%simplify_mesh(0.1d0)
    !         call self%apply_updates()
    !         levels(i+1)%n_vertices = self%n_vertices
    !         levels(i+1)%n_tetrahedra = self%n_tetrahedra
    !         levels(i+1)%vertices = self%vertices
    !         levels(i+1)%tetrahedra = self%tetrahedra
    !     end do

    !     !call self%finalise()

    ! end subroutine wavelet_compression

    subroutine finalise(self)
        implicit none
        
        class(TetrahedralMesh), intent(inout) :: self
        type(UpdateNode), pointer :: current, next
        
        ! Clean up linked list of updates
        current => self%first_update
        do while (associated(current))
            next => current%next
            if (allocated(current%u_minus)) deallocate(current%u_minus)
            if (allocated(current%u_plus)) deallocate(current%u_plus)
            deallocate(current)
            current => next
        end do
        
        ! Nullify update list pointers
        self%first_update => null()
        self%last_update => null()
        
        ! Deallocate arrays
        if (allocated(self%tetrahedra)) deallocate(self%tetrahedra)
        if (allocated(self%vertex_adjacency)) deallocate(self%vertex_adjacency)
        if (allocated(self%tetrahedron_adjacency)) deallocate(self%tetrahedron_adjacency)
        if (allocated(self%vertices)) deallocate(self%vertices)
        
        ! Reset counters
        self%n_updates = 0
        self%n_vertices = 0
        self%n_tetrahedra = 0
    end subroutine finalise
end module multi_resolution_analysis