! ------------------------------------------------------------------------------
! @file vtk.F90
! @brief Module contains routines to handle the VTK file format.
! @author ofgn
! @date 2024-07-01
! @todo - Complete u_grid_read_legacy_vtk.
!       - Complete support for all scalar data.
!       - Add support for vector data.
!       - Add support for other VTK dataset types.
!       - Add support for XML VTK files.
! ------------------------------------------------------------------------------
module vtk
    use iso_fortran_env, only: int16, int32, int64, real32, real64
    use global
    
    implicit none

    type :: VtkUnstructuredGrid
        integer(int64) :: n_points = 0                                          ! Number of points
        integer(int64) :: n_cells = 0                                           ! Number of cells
        integer(int32), allocatable :: cell_size(:)                             ! Number of points per cell
        integer(int32), allocatable :: cell_type(:)                             ! VTK cell type
        integer(int64), allocatable :: cell_connectivity(:)                     ! Cell connectivity array
        real(real64), allocatable :: points(:, :)                               ! Point coordinates
    contains
        procedure :: read_legacy_vtk => u_grid_read_legacy_vtk                  ! Read a VTK legacy file (.vtk) and extract point and cell definitions.
        procedure :: read_tetgen_node => u_grid_read_tetgen_node                ! Read a TetGen node file (.node) and extract point definitions.
        procedure :: read_tetgen_ele => u_grid_read_tetgen_ele                  ! Read a TetGen ele file (.ele) and extract cell definitions.
        procedure :: read_tetgen_face => u_grid_read_tetgen_face                ! Read a TetGen face file (.face) and extract face definitions.
        procedure :: write_legacy_vtk => u_grid_write_legacy_vtk                ! Write an unstructured grid to a VTK legacy file in either binary or ASCII format.
        procedure :: mask_points => u_grid_mask_points                          ! Mask the unstructured grid based on the given point mask.
        procedure :: mask_cells => u_grid_mask_cells                            ! Mask the unstructured grid based on the given cell mask.
        ! procedure :: mask_from_grd_file => u_grid_mask_from_grd_file          
        procedure :: clear => u_grid_clear                                      ! Clear the unstructured grid by deallocating all arrays.
    end type VtkUnstructuredGrid

    type :: VtkData
        integer(int64) :: n_data_points = 0                                     ! Number of points or cells
        integer(int32), allocatable :: scalar_int32(:, :)                       ! Scalar int32 data
        integer(int64), allocatable :: scalar_int64(:, :)                       ! Scalar int64 data
        real(real32), allocatable :: scalar_real32(:, :)                        ! Scalar real32 data
        real(real64), allocatable :: scalar_real64(:, :)                        ! Scalar real64 data
        character(len=32), allocatable :: scalar_int32_labels(:)                ! Labels for int32 scalar data
        character(len=32), allocatable :: scalar_int64_labels(:)                ! Labels for int64 scalar data
        character(len=32), allocatable :: scalar_real32_labels(:)               ! Labels for real32 scalar data
        character(len=32), allocatable :: scalar_real64_labels(:)               ! Labels for real64 scalar data
        integer(int32), allocatable :: scalar_int32_components(:)               ! Components of int32 scalar data
        integer(int32), allocatable :: scalar_int64_components(:)               ! Components of int64 scalar data
        integer(int32), allocatable :: scalar_real32_components(:)              ! Components of real32 scalar data
        integer(int32), allocatable :: scalar_real64_components(:)              ! Components of real64 scalar data
        integer(int32), allocatable :: vector_int32(:, :)                       ! Vector int32 data
        integer(int32), allocatable :: vector_int64(:, :)                       ! Vector int64 data
        real(real32), allocatable :: vector_real32(:, :)                        ! Vector real32 data
        real(real64), allocatable :: vector_real64(:, :)                        ! Vector real64 data
        character(len=32), allocatable :: vector_int32_labels(:)                ! Labels for int32 vector data
        character(len=32), allocatable :: vector_int64_labels(:)                ! Labels for int64 vector data
        character(len=32), allocatable :: vector_real32_labels(:)               ! Labels for real32 vector data
        character(len=32), allocatable :: vector_real64_labels(:)               ! Labels for real64 vector data
    contains
        procedure :: read_scalar_int32 => data_read_scalar_int32                ! Read int32 data from a file.
        procedure :: load_scalar_int32 => data_load_scalar_int32                ! Load scalar int32 data.
        procedure :: read_scalar_int64 => data_read_scalar_int64                ! Read int64 data from a file.
        procedure :: load_scalar_int64 => data_load_scalar_int64                ! Load scalar int64 data.
        procedure :: read_scalar_real32 => data_read_scalar_real32              ! Read real32 data from a file.
        procedure :: load_scalar_real32 => data_load_scalar_real32              ! Load scalar real32 data.
        procedure :: read_scalar_real64 => data_read_scalar_real64              ! Read real64 data from a file.
        procedure :: load_scalar_real64 => data_load_scalar_real64              ! Load scalar real64 data.
        procedure :: mask => data_mask                                          ! Mask data based on a logical mask.
    end type VtkData
contains

    ! -------------------------------------------------------------------------
    ! @brief Read a VtkUnstructuredGrid from a VTK legacy file (.vtk).
    ! @param[inout] self The unstructured grid type.
    ! @param[in] file_path Path to the VTK file.
    ! @param[out] point_data VtkData type to store point data (optional).
    ! @param[out] cell_data VtkData type to store cell data (optional).
    ! -------------------------------------------------------------------------
    subroutine u_grid_read_legacy_vtk(self, file_path, point_data, cell_data)
        use iso_fortran_env, only: int32, int64, real32, real64

        class(VtkUnstructuredGrid), intent(inout) :: self                       ! Unstructured grid
        character(len=*), intent(in) :: file_path                               ! Path to the VTK file
        type(VtkData), intent(inout), optional :: point_data                    ! Optional point data
        type(VtkData), intent(inout), optional :: cell_data                     ! Optional cell data

    end subroutine u_grid_read_legacy_vtk

    ! -------------------------------------------------------------------------
    ! @brief Read a TetGen node file (.node) and extract point definitions.
    ! @param[in] self The unstructured grid to store the point definitions.
    ! @param[in] file_path Path to the TetGen node file.
    ! @param[out] point_data VtkData type to store point data (optional).
    ! -------------------------------------------------------------------------
    subroutine u_grid_read_tetgen_node(self, file_path, point_data)
        implicit none

        class(VtkUnstructuredGrid), intent(inout) :: self                       ! The unstructured grid to store the node data
        character(len=*), intent(in) :: file_path                               ! Path to the TetGen node file
        type(VtkData), intent(inout), optional :: point_data                    ! Data structure for storing point data (optional)

        integer(int64) :: i                                                     ! Loop index                 
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: n_dimensions                                          ! Number of dimensions
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int32) :: id_flag                                               ! Boundary marker flag (1 = true, 0 = false)                    
        integer(int64) :: point_index                                           ! Point index
        real(real64) :: start_time                                              ! Start time                    
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
        call report("Reading TetGen node file", 1)
        call report("File: " // trim(file_path) , 1)

        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_points, n_dimensions, &
            n_fields, id_flag

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        if (n_dimensions .eq. 2) then
            call report("Format: 2D points (Triangle)", 1)
        else if (n_dimensions .eq. 3) then
            call report("Format: 3D points (TetGen)", 1)
        else
            write(string, "(A, I0)") "Unsupported number of dimensions: ", &
                n_dimensions
            call report(trim(string), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of points: ", self%n_points
        call report(trim(string), 1)

        allocate(self%points(3, self%n_points))

        if (present(point_data) .and. n_fields .gt. 0) then
            allocate(point_data%scalar_real64(n_fields, self%n_points))
            allocate(point_data%scalar_real64_labels(n_fields))
            allocate(point_data%scalar_real64_components(n_fields))

            write(string, "(A, I0)") "Number of scalar fields: ", n_fields
            call report(trim(string), 1)

            do i = 1, n_fields
                write(point_data%scalar_real64_labels(i), "(A, I0)") &
                    "scalar_real64_", i
            end do
            point_data%scalar_real64_components = 1
        end if

        if (present(point_data) .and. id_flag .gt. 0) then
            allocate(point_data%scalar_int32(1, self%n_points))
            allocate(point_data%scalar_int32_labels(1))
            allocate(point_data%scalar_int32_components(1))
            point_data%scalar_int32_labels(1) = "Boundary"
            point_data%scalar_int32_components = 1
        end if

        do i = 1, self%n_points
            if (present(point_data)) then
                if ((n_fields .gt. 0) .and. (id_flag .gt. 0)) then
                    read (unit, *, iostat=io_status) point_index, &
                        self%points(1:n_dimensions, i), &
                        point_data%scalar_real64(:, i), &
                        point_data%scalar_int32(1, i)
                    point_data%n_data_points = self%n_points
                else if (n_fields .gt. 0) then
                    read (unit, *, iostat=io_status) point_index, &
                        self%points(1:n_dimensions, i), &
                        point_data%scalar_real64(:, i)
                    point_data%n_data_points = self%n_points
                else if (id_flag .gt. 0) then
                    read (unit, *, iostat=io_status) point_index, &
                        self%points(1:n_dimensions, i), &
                        point_data%scalar_int32(1, i)
                    point_data%n_data_points = self%n_points
                else
                    read (unit, *, iostat=io_status) point_index, &
                    self%points(1:n_dimensions, i)
                end if
            else
                read (unit, *, iostat=io_status) point_index, &
                self%points(1:n_dimensions, i)
            end if

            if (io_status .ne. 0) then
                call report("File data read failed", 3)
                stop
            end if

            if (n_dimensions .eq. 2) then
                self%points(3, i) = 0.0d0
            end if
        end do

        close (unit)
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine u_grid_read_tetgen_node

    ! -------------------------------------------------------------------------
    ! @brief Read a TetGen ele file (.ele) and extract cell definitions.
    ! @param[inout] self The unstructured grid to store the cell definitions.
    ! @param[in] file_path Path to the TetGen element file.
    ! @param[out] cell_data VtkData type to store cell data (optional).
    ! -------------------------------------------------------------------------
    subroutine u_grid_read_tetgen_ele(self, file_path, cell_data)
        implicit none

        class(VtkUnstructuredGrid) :: self                                      ! The unstructured grid to store the cell data
        character(len=*), intent(in) :: file_path                               ! Path to the TetGen element file
        type(VtkData), optional :: cell_data                                    ! Data structure for storing cell data (optional)                 

        integer(int64) :: i                                                     ! Loop index                 
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)

        if (allocated(self%cell_connectivity)) then
            call report("Tetgen ele file read while existing cells exist", 2)
        end if

        if (.not. allocated(self%points)) then
            call report("Tetgen ele file read before node file", 2)
        end if

        call report("Reading TetGen ele file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_cells, cell_size, n_fields

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate (self%cell_connectivity(self%n_cells * cell_size))
        allocate (self%cell_size(self%n_cells))
        allocate (self%cell_type(self%n_cells))

        if (cell_size .eq. 3) then
            self%cell_type = 5
            call report("Format: 2D Unstructured grid (Triangle)", 1)
            call report("Cell type: Triangle ", 1)
        else if (cell_size .eq. 4) then
            self%cell_type = 10
            call report("Format: 3D Unstructured grid (TetGen)", 1)
            call report("Cell type: Tetrahedra (Linear)", 1)
        else if (cell_size .eq. 10) then
            self%cell_type = 24
            call report("Format: 3D Unstructured grid (TetGen)", 1)
            call report("Cell type: Tetrahedra (Quadratic)", 1)
        else
            write(string, "(A, I0)") &
                "Unsupported number of points per cell: ", cell_size
            call report(trim(string), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of cells: ", self%n_cells
        call report(trim(string), 1)

        self%cell_size = cell_size

        if (present(cell_data) .and. (n_fields > 0)) then
            cell_data%n_data_points = self%n_cells

            allocate (cell_data%scalar_real64(n_fields, self%n_cells))
            allocate (cell_data%scalar_real64_labels(n_fields))
            allocate (cell_data%scalar_real64_components(n_fields))
            cell_data%scalar_real64_components = 1
            write(string, "(A, I0)") "Number of scalar fields: ", n_fields
            call report(trim(string), 1)
        end if

        if (n_fields > 0) then
            if (present(cell_data)) then
                read (unit, *, iostat=io_status) &
                    (cell_index, self%cell_connectivity((i - 1) &
                    * cell_size + 1: i * cell_size), &
                    cell_data%scalar_real64(:, i), i = 1, self%n_cells)

                    do i = 1, n_fields
                        write(cell_data%scalar_real64_labels(i), "(A, I0)") &
                            "scalar_real64_", i
                    end do
            else 
                read (unit, *, iostat=io_status) &
                    (cell_index, self%cell_connectivity((i - 1) &
                    * cell_size + 1:i * cell_size), &
                    scalar_real64, i = 1, self%n_cells)
            end if
        else
            read (unit, *, iostat=io_status) &
                (cell_index, self%cell_connectivity((i - 1) * &
                cell_size + 1: i * cell_size), i=1, self%n_cells)
        end if

        if (io_status .ne. 0) then
            call report("File data read failed", 3)
            stop
        end if

        close (unit)
        self%cell_connectivity = self%cell_connectivity - 1
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine u_grid_read_tetgen_ele

    subroutine u_grid_read_tetgen_face(self, file_path, cell_data)
        implicit none

        class(VtkUnstructuredGrid) :: self                                      ! The unstructured grid to store the cell data
        character(len=*), intent(in) :: file_path                               ! Path to the TetGen element file
        type(VtkData), optional :: cell_data                                    ! Data structure for storing cell data (optional)                 

        integer(int64) :: i                                                     ! Loop index                 
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32) :: cell_size                                             ! Number of points per cell                   
        integer(int32) :: n_fields                                              ! Number of scalar fields
        integer(int64) :: cell_index                                            ! Cell index
        real(real64) :: scalar_real64                                           ! Real scalar value
        real(real64) :: start_time                                              ! Start time 
        real(real64) :: run_time                                                ! Runtime                   
        real(real64) :: end_time                                                ! End time
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)

        if (allocated(self%cell_connectivity)) then
            call report("Tetgen face file read while existing cells exist", 2)
        end if

        call cpu_time(start_time)

        if (.not. allocated(self%points)) then
            call report("Tetgen ele file read before node file", 2)
        end if

        call report("Reading TetGen ele file", 1)
        call report("File: " // trim(file_path) , 1)
        
        open (newunit=unit, file=file_path, status="old", action="read", &
            iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        read (unit, *, iostat=io_status) self%n_cells, n_fields
        cell_size = 3

        if (io_status .ne. 0) then
            call report("File header read failed", 3)
            stop
        end if

        allocate (self%cell_connectivity(self%n_cells * cell_size))
        allocate (self%cell_size(self%n_cells))
        allocate (self%cell_type(self%n_cells))

        if (cell_size .eq. 3) then
            self%cell_type = 5
            call report("Format: 3D Unstructured grid (Tetgen)", 1)
            call report("Cell type: Triangle ", 1)
        else
            write(string, "(A, I0)") &
                "Unsupported number of points per cell: ", cell_size
            call report(trim(string), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of cells: ", self%n_cells
        call report(trim(string), 1)

        self%cell_size = cell_size

        if (present(cell_data) .and. (n_fields > 0)) then
            cell_data%n_data_points = self%n_cells

            allocate (cell_data%scalar_real64(n_fields, self%n_cells))
            allocate (cell_data%scalar_real64_labels(n_fields))
            allocate (cell_data%scalar_real64_components(n_fields))
            cell_data%scalar_real64_components = 1
            write(string, "(A, I0)") "Number of scalar fields: ", n_fields
            call report(trim(string), 1)
        end if

        if (n_fields > 0) then
            if (present(cell_data)) then
                read (unit, *, iostat=io_status) &
                    (cell_index, self%cell_connectivity((i - 1) &
                    * cell_size + 1: i * cell_size), &
                    cell_data%scalar_real64(:, i), i = 1, self%n_cells)

                    do i = 1, n_fields
                        write(cell_data%scalar_real64_labels(i), "(A, I0)") &
                            "scalar_real64_", i
                    end do
            else 
                read (unit, *, iostat=io_status) &
                    (cell_index, self%cell_connectivity((i - 1) &
                    * cell_size + 1:i * cell_size), &
                    scalar_real64, i = 1, self%n_cells)
            end if
        else
            read (unit, *, iostat=io_status) &
                (cell_index, self%cell_connectivity((i - 1) * &
                cell_size + 1: i * cell_size), i=1, self%n_cells)
        end if

        if (io_status .ne. 0) then
            call report("File data read failed", 3)
            stop
        end if

        close (unit)
        self%cell_connectivity = self%cell_connectivity - 1
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)

    end subroutine u_grid_read_tetgen_face

    ! -------------------------------------------------------------------------
    ! @brief Write an unstructured grid to a VTK legacy file in 
    !   either binary or ASCII format.
    ! @param[inout] self The unstructured grid to store the cell data.
    ! @param[in] file_path Path to the output VTK file.
    ! @param[in] point_data VtkData type to store point data (optional).
    ! @param[in] cell_data VtkData type to store cell data (optional).
    ! @param[in] ascii Logical flag to write in ASCII format (default: false).
    ! --------------------------------------------------------------------------
    subroutine u_grid_write_legacy_vtk(self, file_path, point_data, cell_data, &
            ascii, big_endian)
        implicit none

        class(VtkUnstructuredGrid), intent(in) :: self                          ! Unstructured grid data structure
        character(len=*), intent(in) :: file_path                               ! Path to the output VTK file
        type(VtkData), optional, intent(in) :: point_data                       ! Data structure for point data (optional)
        type(VtkData), optional, intent(in) :: cell_data                        ! Data structure for cell data (optional)
        logical, intent(in), optional :: ascii                                  ! Write in ASCII format (optional, default: false)
        logical, intent(in), optional :: big_endian                             ! Write in big-endian format (optional, default: true)

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index                   
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int64), allocatable :: full_connectivity(:)                     ! Complete cell connectivity array
        real(real64) :: start_time                                              ! Start time                    
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        logical :: binary                                                       ! Flag to determine binary or ASCII format
        logical :: little_endian                                                ! Flag to determine little or big endian
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
        
        if (present(ascii)) then
            binary = .not. ascii
        else
            binary = .true.
        end if

        if (present(big_endian)) then
            little_endian = .not. big_endian
        else
            little_endian = .false.
        end if

        if (self%n_points > huge(0)) then
            binary = .false.
            call report("Legacy VTK does not support above " &
                // "2147483647 points in binary format.", 1)
            call report("Writing in ASCII format.", 1)
        end if

        if (self%n_points > huge(0)) then
            binary = .false.
            call report("Legacy VTK may not support above " &
                // "2147483647 cells in binary format.", 2)
            call report("Consider writing in ASCII format.", 2)
        end if


        if (binary) then
            if (little_endian) then
                open (newunit=unit, file=file_path, status="replace", &
                access="stream", action="write", form="unformatted", &
                iostat=io_status, convert="little_endian")
            else
                open (newunit=unit, file=file_path, status="replace", &
                access="stream", action="write", form="unformatted", &
                iostat=io_status, convert="big_endian")
            end if
        else
            open (newunit=unit, file=file_path, status="replace", &
            access="stream", action="write", form="formatted", &
            iostat=io_status)
        end if

        if (io_status .ne. 0) then
            call report("File open failed", 3)
            stop
        end if

        call report("Writing legacy VTK file", 1)
        call report("File: " // trim(file_path) , 1)

        if (binary) then
            call report("Format: Binary", 1)
        else
            call report("Format: ASCII", 1)
        end if

        call report("Type: Unstructured grid", 1)

        allocate (full_connectivity(size(self%cell_connectivity) &
            + self%n_cells))
        
        j = 1
        do i = 1, self%n_cells
            full_connectivity(j) = self%cell_size(i)
            full_connectivity(j + 1:j + self%cell_size(i)) = &
                self%cell_connectivity((i - 1) * self%cell_size(i) &
                + 1:i * self%cell_size(i))
            j = j + self%cell_size(i) + 1
        end do

        if (binary) then
            write(unit) "# vtk DataFile Version 3.0" // char(10)
            write(unit) "vtk" // char(10)
            write(unit) "BINARY" // char(10)
            write(unit) "DATASET UNSTRUCTURED_GRID" // char(10)
        else
            write(unit, "(a)") "# vtk DataFile Version 3.0"
            write(unit, "(a)") "vtk"
            write(unit, "(a)") "ASCII"
            write(unit, "(a)") "DATASET UNSTRUCTURED_GRID"
        end if

        if (binary) then
            write(string, "(A, I0, A)") "POINTS ", self%n_points, &
                " double" // char(10)
            write(unit) trim(string)
            write(unit) self%points
            write(unit) char(10)
        else
            write(unit, "(A, I0, A)") "POINTS ", self%n_points, " double"
            do i = 1, self%n_points
                write(unit, "(3E20.8)") self%points(:, i)
            end do
        end if

        if (binary) then
            write(string, "(A, I0, A, I0, A)") "CELLS ", &
                self%n_cells, " ", &
                size(full_connectivity), char(10)
            write(unit) trim(string)
            write(unit) int(full_connectivity, int32)
            write(unit) char(10)
        else
            write(unit, "(A, I0, A, I0)") "CELLS ", self%n_cells  , " ", &
                size(full_connectivity)
            do i = 1, size(full_connectivity)
                write(unit, "(I0)") full_connectivity(i)
            end do
        end if

        if (binary) then
            write(string, "(A, I0, A)") "CELL_TYPES ", &
                self%n_cells, char(10)
            write(unit) trim(string)
            write(unit) int(self%cell_type, int32)
            write(unit) char(10)
        else
            write(unit, "(A, I0)") "CELL_TYPES ", self%n_cells
            do i = 1, self%n_cells
                write(unit, "(I0)") self%cell_type(i)
            end do
        end if

        write(string, '(A, I0)') "Number of points: ", self%n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Number of cells: ", self%n_cells
        call report(trim(string), 1)

        if (present(point_data)) then
            if (self%n_points .ne. point_data%n_data_points) then
                write(string, "(A, I0, A, I0)") &
                "Mismatch in number of points and point data: ", &
                    self%n_points, " != ", point_data%n_data_points
                call report(trim(string), 3)
                write(string, "(A)") "Skipping point data"
                write(string, "(A)") "If using a mask, ensure it is applied" &
                    // " to both the unstructured grid and the data."
                call report(trim(string), 3)
                stop
            end if

            if (point_data%n_data_points .ne. 0) then
                if (binary) then
                    write(string, "(A, I0, A)") "POINT_DATA ", &
                    self%n_points, char(10)
                    write(unit) trim(string)
                else
                    write(unit, "(A, I0)") "POINT_DATA ", self%n_points
                end if

                if (allocated(point_data%scalar_int32)) then
                    do i = 1, size(point_data%scalar_int32, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(point_data%scalar_int32_labels(i)), &
                            " int ", &
                            point_data%scalar_int32_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) point_data%scalar_int32(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(point_data%scalar_int32_labels(i)), &
                            " int ", &
                            point_data%scalar_int32_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_points
                                write(unit, "(I0)") &
                                    point_data%scalar_int32(i, j)
                            end do
                        end if
                    end do
                end if

                if (allocated(point_data%scalar_int64)) then
                    do i = 1, size(point_data%scalar_int64, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(point_data%scalar_int64_labels(i)), &
                            " int ", &
                            point_data%scalar_int64_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) point_data%scalar_int64(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(point_data%scalar_int64_labels(i)), &
                            " int ", &
                            point_data%scalar_int64_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_points
                                write(unit, "(I0)") &
                                    point_data%scalar_int64(i, j)
                            end do
                        end if
                    end do
                end if

                if (allocated(point_data%scalar_real32)) then
                    do i = 1, size(point_data%scalar_real32, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(point_data%scalar_real32_labels(i)), &
                            " float ", &
                            point_data%scalar_real32_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) point_data%scalar_real32(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(point_data%scalar_real32_labels(i)), &
                            " float ", &
                            point_data%scalar_real32_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_points
                                write(unit, "(E20.8)") &
                                    point_data%scalar_real32(i, j)
                            end do
                        end if
                    end do
                end if

                if (allocated(point_data%scalar_real64)) then
                    do i = 1, size(point_data%scalar_real64, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(point_data%scalar_real64_labels(i)), &
                            " double ", &
                            point_data%scalar_real64_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) point_data%scalar_real64(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(point_data%scalar_real64_labels(i)), &
                            " double ", &
                            point_data%scalar_real64_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_points
                                write(unit, "(E20.8)") &
                                    point_data%scalar_real64(i, j)
                            end do
                        end if
                    end do
                end if
            end if
        end if

        if (present(cell_data)) then
            if (self%n_cells .ne. cell_data%n_data_points) then
                write(string, "(A, I0, A, I0)") &
                "Mismatch in number of cells and cell data: ", &
                    self%n_cells, " != ", cell_data%n_data_points
                call report(trim(string), 3)
                write(string, "(A)") "Skipping cell data"
                write(string, "(A)") "If using a mask, ensure it is applied" &
                    // " to both the unstructured grid and the data."
                call report(trim(string), 3)
                stop
            end if

            if (cell_data%n_data_points .ne. 0) then
                if (binary) then
                    write(string, "(A, I0, A)") "CELL_DATA ", &
                    self%n_cells, char(10)
                    write(unit) trim(string)
                else
                    write(unit, "(A, I0)") "CELL_DATA ", self%n_cells
                end if

                if (allocated(cell_data%scalar_int32)) then
                    do i = 1, size(cell_data%scalar_int32, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(cell_data%scalar_int32_labels(i)), &
                            " int ", &
                            cell_data%scalar_int32_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) cell_data%scalar_int32(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(cell_data%scalar_int32_labels(i)), &
                            " int ", &
                            cell_data%scalar_int32_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_cells
                                write(unit, "(I0)") &
                                    cell_data%scalar_int32(i, j)
                            end do
                        end if
                    end do
                end if

                if (allocated(cell_data%scalar_int64)) then
                    do i = 1, size(cell_data%scalar_int64, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(cell_data%scalar_int64_labels(i)), &
                            " int ", &
                            cell_data%scalar_int64_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) cell_data%scalar_int64(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(cell_data%scalar_int64_labels(i)), &
                            " int ", &
                            cell_data%scalar_int64_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_cells
                                write(unit, "(I0)") &
                                    cell_data%scalar_int64(i, j)
                            end do
                        end if
                    end do
                end if

                if (allocated(cell_data%scalar_real32)) then
                    do i = 1, size(cell_data%scalar_real32, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(cell_data%scalar_real32_labels(i)), &
                            " float ", &
                            cell_data%scalar_real32_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) cell_data%scalar_real32(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(cell_data%scalar_real32_labels(i)), &
                            " float ", &
                            cell_data%scalar_real32_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_cells
                                write(unit, "(E20.8)") &
                                    cell_data%scalar_real32(i, j)
                            end do
                        end if
                    end do
                end if
    
                if (allocated(cell_data%scalar_real64)) then
                    do i = 1, size(cell_data%scalar_real64, 1)
                        if (binary) then
                            write(string, "(A, A, A, I0, A)") "SCALARS ", &
                            trim(cell_data%scalar_real64_labels(i)), &
                            " double ", &
                            cell_data%scalar_real64_components(i), char(10)
                            write(unit) trim(string)
                            write(unit) "LOOKUP_TABLE default" // char(10)
                            write(unit) cell_data%scalar_real64(i, :)
                        else
                            write(unit, "(A, A, A, I0)") "SCALARS ", &
                            trim(cell_data%scalar_real64_labels(i)), &
                                " double ", &
                                cell_data%scalar_real64_components(i)
                            write(unit, "(A)") "LOOKUP_TABLE default"
                            do j = 1, self%n_cells
                                write(unit, "(E20.8)") &
                                    cell_data%scalar_real64(i, j)
                            end do
                        end if
                    end do
                end if
            end if
        end if

        close (unit)
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine u_grid_write_legacy_vtk

    ! --------------------------------------------------------------------------
    ! @brief Mask the unstructured grid based on the given point mask.
    ! @param[inout] self The unstructured grid to store the cell data.
    ! @param[in] point_mask Logical mask to determine which points to keep.
    ! --------------------------------------------------------------------------   
    subroutine u_grid_mask_points(self, point_mask, point_data, cell_data)
        implicit none

        class(VtkUnstructuredGrid), intent(inout) :: self                       ! The unstructured grid to modify
        logical, intent(in) :: point_mask(:)                                    ! Logical mask for filtering points
        type(VtkData), optional, intent(inout) :: point_data                    ! Data structure for storing point data (optional)
        type(VtkData), optional, intent(inout) :: cell_data                     ! Data structure for storing cell data (optional)

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index    
        integer(int64) :: k                                                     ! Loop index
        integer(int64) :: new_n_points                                          ! New number of points
        integer(int64) :: new_n_cells                                           ! New number of cells                  
        integer(int64), allocatable :: new_cell_connectivity(:)                 ! New cell connectivity array
        integer(int32), allocatable :: new_cell_size(:)                         ! New points-per-cell array
        integer(int32), allocatable :: new_cell_type(:)                         ! New cell type array
        integer(int64), allocatable :: point_map(:)                             ! Mapping from old points to new points
        real(real64), allocatable :: new_points(:, :)                           ! New point coordinates array
        real(real64) :: start_time                                              ! Start time                    
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        logical, allocatable :: cell_mask(:)                                    ! Masks for cells and connectivity
        logical, allocatable :: connectivity_mask(:)                            ! Masks for connectivity
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
        call report("Applying point mask", 1)

        if (size(point_mask) .ne. self%n_points) then
            call report("Mask size does not match number of points.", 3)
        end if

        allocate(cell_mask(self%n_cells), &
            connectivity_mask(size(self%cell_connectivity)))
        cell_mask = .true.
        connectivity_mask = .true.

        j = 1
        do i = 1, self%n_cells
            do k = 0, self%cell_size(i) - 1
                if (.not. point_mask(self%cell_connectivity(j + k) + 1)) then
                    cell_mask(i) = .false.
                    connectivity_mask(j:j + self%cell_size(i) - 1) &
                        = .false.
                    exit
                end if
            end do
            j = j + self%cell_size(i)
        end do

        allocate(point_map(self%n_points))
        point_map = -1
        j = 1
        do i = 1, self%n_points
            if (point_mask(i)) then
                point_map(i) = j
                j = j + 1
            end if
        end do

        new_n_points = count(point_mask)
        new_n_cells = count(cell_mask)

        allocate(new_points(3, new_n_points))
        new_points = self%points(:, pack([(i, i=1, self%n_points)], point_mask))

        new_cell_connectivity = pack(self%cell_connectivity, connectivity_mask)
        new_cell_size = pack(self%cell_size, cell_mask)
        new_cell_type = pack(self%cell_type, cell_mask)

        if (present(point_data)) then
            call point_data%mask(point_mask)
        end if

        if (present(cell_data)) then
            call cell_data%mask(cell_mask)
        end if

        new_cell_connectivity = point_map(new_cell_connectivity + 1) - 1
        deallocate(cell_mask, connectivity_mask)

        self%points = new_points
        self%cell_connectivity = new_cell_connectivity
        self%cell_size = new_cell_size
        self%cell_type = new_cell_type


        write(string, '(A, I0)') "Points prior to masking: ", self%n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Points after masking: ", new_n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Points removed: ", self%n_points - new_n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells prior to masking: ", self%n_cells
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells after masking: ", new_n_cells
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells removed: ", self%n_cells - new_n_cells
        call report(trim(string), 1)

        self%n_points = new_n_points
        self%n_cells = new_n_cells
        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine u_grid_mask_points

    ! --------------------------------------------------------------------------
    ! @brief Mask the unstructured grid based on the given cell mask.
    ! @param[inout] self The unstructured grid to modify.
    ! @param[in] cell_mask Logical mask to determine which cells to keep.
    ! --------------------------------------------------------------------------
    subroutine u_grid_mask_cells(self, cell_mask, point_data, cell_data)
        implicit none

        class(VtkUnstructuredGrid), intent(inout) :: self                       ! The unstructured grid to modify
        logical, intent(in) :: cell_mask(:)                                     ! Logical mask for filtering cells
        type(VtkData), optional, intent(inout) :: point_data                    ! Data structure for storing point data (optional)
        type(VtkData), optional, intent(inout) :: cell_data                     ! Data structure for storing cell data (optional)

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: j                                                     ! Loop index    
        integer(int64) :: k                                                     ! Loop index
        integer(int64) :: new_n_points                                          ! New number of points
        integer(int64) :: new_n_cells                                           ! New number of cells                  
        integer(int64), allocatable :: new_cell_connectivity(:)                 ! New cell connectivity array
        integer(int32), allocatable :: new_cell_size(:)                         ! New points-per-cell array
        integer(int32), allocatable :: new_cell_type(:)                         ! New cell type array
        integer(int64), allocatable :: point_map(:)                             ! Mapping from old points to new points
        real(real64), allocatable :: new_points(:, :)                           ! New point coordinates array
        real(real64) :: start_time                                              ! Start time                    
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        logical, allocatable :: point_mask(:)                                   ! Masks for cells and connectivity
        logical, allocatable :: connectivity_mask(:)                            ! Masks for connectivity
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
        call report("Applying cell mask", 1)

        if (size(cell_mask) .ne. self%n_cells) then
            call report("Mask size does not match number of cells.", 3)
            stop
        end if

        allocate(point_mask(self%n_points), &
            connectivity_mask(size(self%cell_connectivity)))

        point_mask = .false.
        connectivity_mask = .false.

        j = 1
        do i = 1, self%n_cells
            if (cell_mask(i)) then
                point_mask(self%cell_connectivity(j:j + self%cell_size(i) - 1) &
                    + 1) = .true.
                connectivity_mask(j:j + self%cell_size(i) - 1) = .true.
            end if
            j = j + self%cell_size(i)
        end do

        allocate(point_map(self%n_points))
        point_map = -1
        j = 1
        do i = 1, self%n_points
            if (point_mask(i)) then
                point_map(i) = j
                j = j + 1
            end if
        end do

        new_n_points = count(point_mask)
        new_n_cells = count(cell_mask)

        allocate(new_points(3, new_n_points))
        new_points = self%points(:, pack([(i, i=1, self%n_points)], point_mask))

        new_cell_connectivity = pack(self%cell_connectivity, connectivity_mask)
        new_cell_size = pack(self%cell_size, cell_mask)
        new_cell_type = pack(self%cell_type, cell_mask)

        if (present(point_data)) then
            call point_data%mask(point_mask)
        end if

        if (present(cell_data)) then
            call cell_data%mask(cell_mask)
        end if

        new_cell_connectivity = point_map(new_cell_connectivity + 1) - 1
        deallocate(point_mask, connectivity_mask)

        self%points = new_points
        self%cell_connectivity = new_cell_connectivity
        self%cell_size = new_cell_size
        self%cell_type = new_cell_type

        write(string, '(A, I0)') "Points prior to masking: ", self%n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Points after masking: ", new_n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Points removed: ", self%n_points - new_n_points
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells prior to masking: ", self%n_cells
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells after masking: ", new_n_cells
        call report(trim(string), 1)
        write(string, '(A, I0)') "Cells removed: ", self%n_cells - new_n_cells
        call report(trim(string), 1)

        self%n_points = new_n_points
        self%n_cells = new_n_cells

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine u_grid_mask_cells

    ! subroutine u_grid_mask_from_grd_file(self, file_path, point_mask, offset)
    !     use maths
    
    !     implicit none
    
    !     class(VtkUnstructuredGrid), intent(inout) :: self                       ! The unstructured grid to modify
    !     character(len=*), intent(in) :: file_path                               ! Path to the Surfer grid file
    !     logical, allocatable, intent(inout) :: point_mask(:)                    ! Logical mask for filtering points
    !     real(real64), optional, intent(in) :: offset                            ! Offset to apply to the z values
    
    !     integer(int64) :: i, j                                                  ! Loop indices
    !     integer(int32) :: unit                                                  ! File unit
    !     integer(int32) :: io_status                                             ! I/O status
    !     integer(int16) :: nx, ny                                                ! Grid dimensions
    !     real(real32), allocatable :: z(:, :)                                    ! Z values
    !     real(real64) :: x_min, x_max, y_min, y_max                              ! Grid bounds
    !     real(real64) :: c                                                       ! Offset value
    !     real(real64) :: start_time, end_time, run_time                          ! Timing variables
    !     real(real64), allocatable :: grid(:, :, :)
    !     real(real32), parameter :: no_data = 1.70141e38                         ! Surfer "no-data" value
    !     character(len=4) :: id                                                  ! File identifier
    !     character(len=4096) :: string                                           ! Log string
    !     real(real64) :: z_interp                                                ! Interpolated z-value
    
    !     call cpu_time(start_time)
    !     call report("Reading Surfer Grid file", 1)
    !     call report("File: " // trim(file_path), 1)
    
    !     ! Determine offset value
    !     if (present(offset)) then
    !         c = offset
    !     else
    !         c = 0.0d0
    !     end if
    
    !     ! Open the grid file
    !     open(newunit=unit, file=file_path, status="old", &
    !          access="stream", action="read", form="unformatted", &
    !          iostat=io_status)
    !     if (io_status /= 0) then
    !         call report("Error opening file: " // trim(file_path), 3)
    !         stop
    !     end if
    
    !     ! Read the file identifier
    !     read(unit) id
    !     if (id /= "DSBB") then
    !         call report("Unsupported file format: " // trim(id), 3)
    !         stop
    !     end if
    
    !     ! Read grid metadata
    !     read(unit) nx, ny, x_min, x_max, y_min, y_max
    !     allocate(z(ny, nx))
    !     read(unit) z
    !     close(unit)
    
    !     if (nx .eq. 1 .or. ny .eq. 1) then
    !         write(string, "(A)") "Grid dimensions must be greater than 1"
    !         call report(trim(string), 3)
    !         stop
    !     end if
    
    !     allocate(grid(nx, ny, 3))
    
    !     ! Populate the grid with x, y, z values
    !     do i = 1, nx
    !         do j = 1, ny
    !             if (z(j, i) == no_data) then
    !                 grid(i, j, 3) = no_data
    !             else
    !                 grid(i, j, 1) = x_min + (i - 1) * (x_max - x_min) / real(nx - 1, real64)
    !                 grid(i, j, 2) = y_min + (j - 1) * (y_max - y_min) / real(ny - 1, real64)
    !                 grid(i, j, 3) = z(j, i)
    !             end if
    !         end do
    !     end do
    
    !     ! Allocate point mask if not already allocated
    !     if (.not. allocated(point_mask)) then
    !         allocate(point_mask(self%n_points))
    !         point_mask = .true.
    !     end if
    
    !     ! Mask points based on interpolated z-values
    !     do i = 1, self%n_points
    !         if (self%points(1, i) < x_min .or. self%points(1, i) > x_max .or. &
    !             self%points(2, i) < y_min .or. self%points(2, i) > y_max) then
    !             point_mask(i) = .false.
    !         else
    !             z_interp = bicubic_interpolate(self%points(1, i), self%points(2, i), grid)
    !             if (z_interp == no_data) then
    !                 point_mask(i) = .false.
    !             else
    !                 point_mask(i) = self%points(3, i) < z_interp + c
    !             end if
    !         end if
    !     end do
    
    !     call cpu_time(end_time)
    !     run_time = end_time - start_time
    !     write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
    !     call report(trim(string), 1)
    ! end subroutine u_grid_mask_from_grd_file
    

    ! --------------------------------------------------------------------------
    ! @brief Reset the unstructured grid by deallocating all arrays.
    ! @param[inout] self The unstructured grid to clear.
    ! --------------------------------------------------------------------------
    subroutine u_grid_clear(self)
        implicit none

        class(VtkUnstructuredGrid), intent(inout) :: self                       ! The unstructured grid to clear

        deallocate(self%points)
        deallocate(self%cell_size)
        deallocate(self%cell_connectivity)
        deallocate(self%cell_type)

        self%n_points = 0
        self%n_cells = 0
    end subroutine u_grid_clear

    ! --------------------------------------------------------------------------
    ! @brief Load int32 scalar data into the VTK data structure.
    ! @param[inout] self The VtkData type to store the data
    ! @param[in] new_scalars_int32 New int32 scalar data to append.
    ! @param[in] new_scalars_int32_labels Labels for the new scalar data.
    ! @param[in] new_scalars_int32_components Number of components for the new
    !   scalar.
    ! --------------------------------------------------------------------------
    subroutine data_load_scalar_int32(self, new_scalars_int32, &
        new_scalars_int32_labels, new_scalars_int32_components)

        implicit none

        class(VtkData), intent(inout) :: self                                   ! Data structure containing int32 scalar data
        integer(int32), intent(in) :: new_scalars_int32(:, :)                  ! New int32 scalar data to append
        character(len=*), intent(in), optional :: new_scalars_int32_labels      ! Labels for the new scalar data
        integer(int64), intent(in), optional :: new_scalars_int32_components    ! Number of components for the new scalar

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: n_data                                                ! Number of elements per scalar
        integer(int64) :: n_new_scalars                                         ! Number of new scalar fields
        integer(int64), allocatable :: temp_components(:)                       ! Temporary array for components
        integer(int32), allocatable :: temp_scalars(:, :)                       ! Temporary array for scalar data
        logical :: unallocated                                                  ! Flag to check allocation state of arrays
        character(len=32), allocatable :: temp_labels(:)                        ! Temporary array for scalar labels

        n_new_scalars = size(new_scalars_int32, 1)
        n_data = size(new_scalars_int32, 2)
        self%n_data_points = n_data

        unallocated = .not. (allocated(self%scalar_int32) .and. &
            allocated(self%scalar_int32_labels) .and. &
            allocated(self%scalar_int32_components))

        if (unallocated) then
            self%scalar_int32 = new_scalars_int32

            if (present(new_scalars_int32_components)) then
                self%scalar_int32_components = new_scalars_int32_components
            else
                allocate (self%scalar_int32_components(n_new_scalars))
                self%scalar_int32_components = 1
            end if

            if (present(new_scalars_int32_labels)) then
                self%scalar_int32_labels = new_scalars_int32_labels
            else
                allocate (self%scalar_int32_labels&
                    (size(self%scalar_int32_components)))

                do i = 1, size(self%scalar_int32_components)
                    write(self%scalar_int32_labels(i), &
                        "(A, I0)") "scalar_int32_", i
                end do
            end if
        else
            
        end if
    end subroutine data_load_scalar_int32

    ! --------------------------------------------------------------------------
    ! @brief Load int64 scalar data into the VTK data structure.
    ! @param[inout] self The VtkData type to store the data
    ! @param[in] new_scalars_int64 New int64 scalar data to append.
    ! @param[in] new_scalars_int64_labels Labels for the new scalar data.
    ! @param[in] new_scalars_int64_components Number of components for the new
    !   scalar.
    ! --------------------------------------------------------------------------
    subroutine data_load_scalar_int64(self, new_scalars_int64, &
        new_scalars_int64_labels, new_scalars_int64_components)

        implicit none

        class(VtkData), intent(inout) :: self                                   ! Data structure containing int64 scalar data
        integer(int64), intent(in) :: new_scalars_int64(:, :)                   ! New int64 scalar data to append
        character(len=*), intent(in), optional :: new_scalars_int64_labels      ! Labels for the new scalar data
        integer(int64), intent(in), optional :: new_scalars_int64_components    ! Number of components for the new scalar

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: n_data                                                ! Number of elements per scalar
        integer(int64) :: n_new_scalars                                         ! Number of new scalar fields
        integer(int64), allocatable :: temp_components(:)                       ! Temporary array for components
        integer(int64), allocatable :: temp_scalars(:, :)                       ! Temporary array for scalar data
        logical :: unallocated                                                  ! Flag to check allocation state of arrays
        character(len=32), allocatable :: temp_labels(:)                        ! Temporary array for scalar labels

        n_new_scalars = size(new_scalars_int64, 1)
        n_data = size(new_scalars_int64, 2)
        self%n_data_points = n_data

        unallocated = .not. (allocated(self%scalar_int64) .and. &
            allocated(self%scalar_int64_labels) .and. &
            allocated(self%scalar_int64_components))

        if (unallocated) then
            self%scalar_int64 = new_scalars_int64

            if (present(new_scalars_int64_components)) then
                self%scalar_int64_components = new_scalars_int64_components
            else
                allocate (self%scalar_int64_components(n_new_scalars))
                self%scalar_int64_components = 1
            end if

            if (present(new_scalars_int64_labels)) then
                self%scalar_int64_labels = new_scalars_int64_labels
            else
                allocate (self%scalar_int64_labels&
                    (size(self%scalar_int64_components)))

                do i = 1, size(self%scalar_int64_components)
                    write(self%scalar_int64_labels(i), &
                        "(A, I0)") "scalar_int64_", i
                end do
            end if
        else
            allocate (temp_scalars(size(self%scalar_int64, 1) &
            + n_new_scalars, n_data))

            temp_scalars(1:size(self%scalar_int64, 1), :) = self%scalar_int64(:, :)

            temp_scalars(size(self%scalar_int64, 1) + 1: &
                size(self%scalar_int64, 1) &
                + n_new_scalars, :) = new_scalars_int64(1:n_new_scalars, :)

            deallocate (self%scalar_int64)
            self%scalar_int64 = temp_scalars

            allocate (temp_components(size(self%scalar_int64_components) &
                + n_new_scalars))

            temp_components(1:size(self%scalar_int64_components)) &
                = self%scalar_int64_components

            temp_components(size(self%scalar_int64_components) + &
                1:size(self%scalar_int64_components) + n_new_scalars) = 1

            deallocate (self%scalar_int64_components)

            self%scalar_int64_components = temp_components

            allocate (temp_labels(size(self%scalar_int64_labels) &
                + n_new_scalars))

            temp_labels(1:size(self%scalar_int64_labels)) &
                = self%scalar_int64_labels

            if (present(new_scalars_int64_labels)) then
                temp_labels(size(self%scalar_int64_labels) &
                    + 1:size(self%scalar_int64_labels) + n_new_scalars) &
                    = new_scalars_int64_labels
            else
                do i = 1, n_new_scalars
                    write(temp_labels(size(self%scalar_int64_labels) + i), &
                        "(A, I0)") "scalar_int64_", &
                        size(self%scalar_int64_labels) + i
                end do
            end if
        end if
    end subroutine data_load_scalar_int64

    ! --------------------------------------------------------------------------
    ! @brief Load real32 scalar data into the VTK data structure.
    ! @param[inout] self The VtkData type to store the data
    ! @param[in] new_scalars_real32 New real32 scalar data to append.
    ! @param[in] new_scalars_real32_labels Labels for the new scalar data.
    ! @param[in] new_scalars_real32_components Number of components for the new
    !   scalar.
    ! --------------------------------------------------------------------------
    subroutine data_load_scalar_real32(self, new_scalars_real32, &
        new_scalars_real32_labels, new_scalars_real32_components)

        implicit none

        class(VtkData), intent(inout) :: self                                   ! Data structure containing real(real32) scalar data
        real(real32), intent(in) :: new_scalars_real32(:, :)                    ! New real32 scalar data to append
        character(len=*), intent(in), optional :: new_scalars_real32_labels     ! Labels for the new scalar data
        integer(int64), intent(in), optional :: new_scalars_real32_components   ! Number of components for the new scalar

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: n_data                                                ! Number of elements per scalar
        integer(int64) :: n_new_scalars                                         ! Number of new scalar fields
        integer(int64), allocatable :: temp_components(:)                       ! Temporary array for components
        real(real32), allocatable :: temp_scalars(:, :)                         ! Temporary array for scalar data
        logical :: unallocated                                                  ! Flag to check allocation state of arrays
        character(len=32), allocatable :: temp_labels(:)                        ! Temporary array for scalar labels

        n_new_scalars = size(new_scalars_real32, 1)
        n_data = size(new_scalars_real32, 2)
        self%n_data_points = n_data

        unallocated = .not. (allocated(self%scalar_real32) .and. &
            allocated(self%scalar_real32_labels) .and. &
            allocated(self%scalar_real32_components))

        if (unallocated) then
            self%scalar_real32 = new_scalars_real32

            if (present(new_scalars_real32_components)) then
                self%scalar_real32_components = new_scalars_real32_components
            else
                allocate (self%scalar_real32_components(n_new_scalars))
                self%scalar_real32_components = 1
            end if

            if (present(new_scalars_real32_labels)) then
                self%scalar_real32_labels = new_scalars_real32_labels
            else
                allocate (self%scalar_real32_labels&
                    (size(self%scalar_real32_components)))

                do i = 1, size(self%scalar_real32_components)
                    write(self%scalar_real32_labels(i), &
                        "(A, I0)") "scalar_real32_", i
                end do
            end if
        else
            allocate (temp_scalars(size(self%scalar_real32, 1) &
            + n_new_scalars, n_data))

            temp_scalars(1:size(self%scalar_real32, 1), :) = self%scalar_real32(:, :)

            temp_scalars(size(self%scalar_real32, 1) + 1: &
                size(self%scalar_real32, 1) &
                + n_new_scalars, :) = new_scalars_real32(1:n_new_scalars, :)

            deallocate (self%scalar_real32)
            self%scalar_real32 = temp_scalars

            allocate (temp_components(size(self%scalar_real32_components) &
                + n_new_scalars))
                
            temp_components(1:size(self%scalar_real32_components)) &
                = self%scalar_real32_components

            temp_components(size(self%scalar_real32_components) + &
                1:size(self%scalar_real32_components) + n_new_scalars) = 1

            deallocate (self%scalar_real32_components)

            self%scalar_real32_components = temp_components

            allocate (temp_labels(size(self%scalar_real32_labels) &
                + n_new_scalars))

            temp_labels(1:size(self%scalar_real32_labels)) &
                = self%scalar_real32_labels

            if (present(new_scalars_real32_labels)) then
                temp_labels(size(self%scalar_real32_labels) &
                    + 1:size(self%scalar_real32_labels) + n_new_scalars) &
                    = new_scalars_real32_labels
            else
                do i = 1, n_new_scalars
                    write(temp_labels(size(self%scalar_real32_labels) + i), &
                        "(A, I0)") "scalar_real32_", &
                        size(self%scalar_real32_labels) + i
                end do
            end if
        end if
    end subroutine data_load_scalar_real32

    ! --------------------------------------------------------------------------
    ! @brief Load real64 scalar data into the VTK data structure.
    ! @param[inout] self The VtkData type to store the data
    ! @param[in] new_scalars_real64 New real64 scalar data to append.
    ! @param[in] new_scalars_real64_labels Labels for the new scalar data.
    ! @param[in] new_scalars_real64_components Number of components for the new 
    !   scalar.
    ! --------------------------------------------------------------------------
    subroutine data_load_scalar_real64(self, new_scalars_real64, &
        new_scalars_real64_labels, new_scalars_real64_components)

        implicit none

        class(VtkData), intent(inout) :: self                                   ! Data structure containing real(real64) scalar data
        real(real64), intent(in) :: new_scalars_real64(:, :)                    ! New real64 scalar data to append
        character(len=*), intent(in), optional :: new_scalars_real64_labels     ! Labels for the new scalar data
        integer(int64), intent(in), optional :: new_scalars_real64_components   ! Number of components for the new scalar

        integer(int64) :: i                                                     ! Loop index
        integer(int64) :: n_data                                                ! Number of elements per scalar
        integer(int64) :: n_new_scalars                                         ! Number of new scalar fields
        integer(int64), allocatable :: temp_components(:)                       ! Temporary array for components
        real(real64), allocatable :: temp_scalars(:, :)                         ! Temporary array for scalar data
        logical :: unallocated                                                  ! Flag to check allocation state of arrays
        character(len=32), allocatable :: temp_labels(:)                        ! Temporary array for scalar labels

        n_new_scalars = size(new_scalars_real64, 1)
        n_data = size(new_scalars_real64, 2)
        self%n_data_points = n_data

        unallocated = .not. (allocated(self%scalar_real64) .and. &
            allocated(self%scalar_real64_labels) .and. &
            allocated(self%scalar_real64_components))

        if (unallocated) then
            self%scalar_real64 = new_scalars_real64

            if (present(new_scalars_real64_components)) then
                self%scalar_real64_components = new_scalars_real64_components
            else
                allocate (self%scalar_real64_components(n_new_scalars))
                self%scalar_real64_components = 1
            end if

            if (present(new_scalars_real64_labels)) then
                self%scalar_real64_labels = new_scalars_real64_labels
            else
                allocate (self%scalar_real64_labels&
                    (size(self%scalar_real64_components)))

                do i = 1, size(self%scalar_real64_components)
                    write(self%scalar_real64_labels(i), &
                        "(A, I0)") "scalar_real64_", i
                end do
            end if
        else
            allocate (temp_scalars(size(self%scalar_real64, 1) &
            + n_new_scalars, n_data))

            temp_scalars(1:size(self%scalar_real64, 1), :) = self%scalar_real64(:, :)
            
            temp_scalars(size(self%scalar_real64, 1) + 1: &
                size(self%scalar_real64, 1) &
                + n_new_scalars, :) = new_scalars_real64(1:n_new_scalars, :)

            
            deallocate (self%scalar_real64)
            self%scalar_real64 = temp_scalars


            allocate (temp_components(size(self%scalar_real64_components) &
                + n_new_scalars))

            temp_components(1:size(self%scalar_real64_components)) &
                = self%scalar_real64_components

            temp_components(size(self%scalar_real64_components) + &
                1:size(self%scalar_real64_components) + n_new_scalars) = 1

            deallocate (self%scalar_real64_components)
            self%scalar_real64_components = temp_components

            allocate (temp_labels(size(self%scalar_real64_labels) &
                + n_new_scalars))

            temp_labels(1:size(self%scalar_real64_labels)) &
                = self%scalar_real64_labels
            
            if (present(new_scalars_real64_labels)) then
                temp_labels(size(self%scalar_real64_labels) &
                    + 1:size(self%scalar_real64_labels) + n_new_scalars) &
                    = new_scalars_real64_labels
            else
                do i = 1, n_new_scalars
                    write(temp_labels(size(self%scalar_real64_labels) + i), &
                        "(A, I0)") "scalar_real64_", &
                        size(self%scalar_real64_labels) + i
                end do
            end if
            deallocate (self%scalar_real64_labels)
            self%scalar_real64_labels = temp_labels
        end if
    end subroutine data_load_scalar_real64

    ! --------------------------------------------------------------------------
    ! @brief Read int32 data from a file.
    ! @param[in] self The VtkData type to store the data.
    ! @param[in] file_path Path to the input file.
    ! --------------------------------------------------------------------------
    subroutine data_read_scalar_int32(self, file_path)
        implicit none

        class(VtkData), intent(inout) :: self                                   ! The structure to store the scalar data
        character(len=*), intent(in) :: file_path                               ! Path to the input file

        integer(int64) :: n_data                                                ! Number of data points
        integer(int32) :: n_scalar_fields                                       ! Number of scalar fields
        integer(int64) :: i                                                     ! Loop index
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int32), allocatable :: scalar_int32(:, :)                       ! Int32 scalar data
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
    
        call report("Reading data file", 1)
        call report("File: " // trim(file_path) , 1)

        open(newunit=unit, file=file_path, status="old", &
            action="read", iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed" // char(10), 3)
            stop
        end if
    
        read(unit, *, iostat=io_status) n_data, n_scalar_fields

        if (io_status .ne. 0) then
            call report("File header read failed" // char(10), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of data points: ", n_data
        call report(trim(string), 1)
        write (string, "(A, I0)") "Number of fields: ", n_scalar_fields
        call report(trim(string), 1)
    
        allocate(scalar_int32(n_scalar_fields, n_data))
        read(unit, *, iostat=io_status) scalar_int32(:, :)
        
        call self%load_scalar_int32(scalar_int32)
    
        close(unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine data_read_scalar_int32

    ! --------------------------------------------------------------------------
    ! @brief Read int64 data from a file.
    ! @param[in] self The VtkData type to store the data.
    ! @param[in] file_path Path to the input file.
    ! --------------------------------------------------------------------------
    subroutine data_read_scalar_int64(self, file_path)
        implicit none

        class(VtkData), intent(inout) :: self                                   ! The structure to store the scalar data
        character(len=*), intent(in) :: file_path                               ! Path to the input file

        integer(int64) :: n_data                                                ! Number of data points
        integer(int32) :: n_scalar_fields                                       ! Number of scalar fields
        integer(int64) :: i                                                     ! Loop index
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        integer(int64), allocatable :: scalar_int64(:, :)                       ! Int64 scalar data
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
    
        call report("Reading data file", 1)
        call report("File: " // trim(file_path) , 1)

        open(newunit=unit, file=file_path, status="old", &
            action="read", iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed" // char(10), 3)
            stop
        end if
    
        read(unit, *, iostat=io_status) n_data, n_scalar_fields

        if (io_status .ne. 0) then
            call report("File header read failed" // char(10), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of data points: ", n_data
        call report(trim(string), 1)
        write (string, "(A, I0)") "Number of fields: ", n_scalar_fields
        call report(trim(string), 1)
    
        allocate(scalar_int64(n_scalar_fields, n_data))
        read(unit, *, iostat=io_status) scalar_int64(:, :)
        
        call self%load_scalar_int64(scalar_int64)
    
        close(unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine data_read_scalar_int64

    ! --------------------------------------------------------------------------
    ! @brief Read real32 data from a file.
    ! @param[in] self The VtkData type to store the data.
    ! @param[in] file_path Path to the input file.
    ! --------------------------------------------------------------------------
    subroutine data_read_scalar_real32(self, file_path)
        implicit none

        class(VtkData), intent(inout) :: self                                   ! The structure to store the scalar data
        character(len=*), intent(in) :: file_path                               ! Path to the input file

        integer(int64) :: n_data                                                ! Number of data points
        integer(int32) :: n_scalar_fields                                       ! Number of scalar fields
        integer(int64) :: i                                                     ! Loop index
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        real(real32), allocatable :: scalar_real32(:, :)                        ! Real32 scalar data
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
    
        call report("Reading data file", 1)
        call report("File: " // trim(file_path) , 1)

        open(newunit=unit, file=file_path, status="old", &
            action="read", iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed" // char(10), 3)
            stop
        end if
    
        read(unit, *, iostat=io_status) n_data, n_scalar_fields

        if (io_status .ne. 0) then
            call report("File header read failed" // char(10), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of data points: ", n_data
        call report(trim(string), 1)
        write (string, "(A, I0)") "Number of fields: ", n_scalar_fields
        call report(trim(string), 1)
    
        allocate(scalar_real32(n_scalar_fields, n_data))
        read(unit, *, iostat=io_status) scalar_real32(:, :)
        
        call self%load_scalar_real32(scalar_real32)
    
        close(unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine data_read_scalar_real32

    ! --------------------------------------------------------------------------
    ! @brief Read real64 data from a file.
    ! @param[in] self The VtkData type to store the data.
    ! @param[in] file_path Path to the input file.
    ! --------------------------------------------------------------------------
    subroutine data_read_scalar_real64(self, file_path)
        implicit none

        class(VtkData), intent(inout) :: self                                   ! The structure to store the scalar data
        character(len=*), intent(in) :: file_path                               ! Path to the input file

        integer(int64) :: n_data                                                ! Number of data points
        integer(int32) :: n_scalar_fields                                       ! Number of scalar fields
        integer(int64) :: i                                                     ! Loop index
        integer(int32) :: unit                                                  ! File unit
        integer(int32) :: io_status                                             ! I/O status
        real(real64), allocatable :: scalar_real64(:, :)                        ! Real64 scalar data
        real(real64) :: start_time                                              ! Start time
        real(real64) :: end_time                                                ! End time
        real(real64) :: run_time                                                ! Runtime
        character(len=4096) :: string                                           ! Character string

        call cpu_time(start_time)
    
        call report("Reading data file", 1)
        call report("File: " // trim(file_path) , 1)

        open(newunit=unit, file=file_path, status="old", &
            action="read", iostat=io_status)

        if (io_status .ne. 0) then
            call report("File open failed" // char(10), 3)
            stop
        end if
    
        read(unit, *, iostat=io_status) n_data, n_scalar_fields

        if (io_status .ne. 0) then
            call report("File header read failed" // char(10), 3)
            stop
        end if

        write(string, "(A, I0)") "Number of data points: ", n_data
        call report(trim(string), 1)
        write (string, "(A, I0)") "Number of fields: ", n_scalar_fields
        call report(trim(string), 1)
    
        allocate(scalar_real64(n_scalar_fields, n_data))
        read(unit, *, iostat=io_status) scalar_real64(:, :)
        
        call self%load_scalar_real64(scalar_real64)
    
        close(unit)

        call cpu_time(end_time)
        run_time = end_time - start_time
        write(string, '(A, F9.3, A)') "Completed in ", run_time, "s" // char(10)
        call report(trim(string), 1)
    end subroutine data_read_scalar_real64

    ! --------------------------------------------------------------------------
    ! @brief Mask data based on a logical mask.
    ! @param[inout] self The VTK data structure to mask.
    ! @param[in] mask Logical mask to apply to each data array.
    ! --------------------------------------------------------------------------
    subroutine data_mask(self, mask)
        implicit none

        class(VtkData), intent(inout) :: self                                   ! The VTK data structure to mask
        logical, intent(in) :: mask(:)                                          ! Logical mask for filtering data

        integer(int64) :: i                                                     ! Loop index
        integer(int32), allocatable :: temp_scalar_int32(:, :)                  ! Temporary array for int32 scalar data
        integer(int32), allocatable :: temp_vector_int32(:, :)                  ! Temporary array for int32 vector data
        integer(int64), allocatable :: temp_scalar_int64(:, :)                  ! Temporary array for int64 scalar data
        integer(int64), allocatable :: temp_vector_int64(:, :)                  ! Temporary array for int64 vector data
        real(real32), allocatable :: temp_scalar_real32(:, :)                   ! Temporary array for real32 scalar data
        real(real32), allocatable :: temp_vector_real32(:, :)                   ! Temporary array for real32 vector data
        real(real64), allocatable :: temp_scalar_real64(:, :)                   ! Temporary array for real64 scalar data
        real(real64), allocatable :: temp_vector_real64(:, :)                   ! Temporary array for real64 vector data

        self%n_data_points = count(mask)

        if (allocated(self%scalar_int32)) then
            allocate(temp_scalar_int32(size(self%scalar_int32, 1), &
                self%n_data_points))
            do i = 1, size(self%scalar_int32, 1)
                temp_scalar_int32(i, :) = pack(self%scalar_int32(i, :), &
                    mask)
            end do
            deallocate(self%scalar_int32)
            self%scalar_int32 = temp_scalar_int32
        end if

        if (allocated(self%scalar_int64)) then
            allocate(temp_scalar_int64(size(self%scalar_int64, 1), &
                self%n_data_points))
            do i = 1, size(self%scalar_int64, 1)
                temp_scalar_int64(i, :) = pack(self%scalar_int64(i, :), mask)
            end do
            deallocate(self%scalar_int64)
            self%scalar_int64 = temp_scalar_int64
        end if  

        if (allocated(self%scalar_real32)) then
            allocate(temp_scalar_real32(size(self%scalar_real32, 1), &
                self%n_data_points))
            do i = 1, size(self%scalar_real32, 1)
                temp_scalar_real32(i, :) = pack(self%scalar_real32(i, :), mask)
            end do
            deallocate(self%scalar_real32)
            self%scalar_real32 = temp_scalar_real32
        end if

        if (allocated(self%scalar_real64)) then
            allocate(temp_scalar_real64(size(self%scalar_real64, 1), &
                self%n_data_points))
            do i = 1, size(self%scalar_real64, 1)
                temp_scalar_real64(i, :) = pack(self%scalar_real64(i, :), mask)
            end do
            deallocate(self%scalar_real64)
            self%scalar_real64 = temp_scalar_real64
        end if

        if (allocated(self%vector_int32)) then
            allocate(temp_vector_int32(size(self%vector_int32, 1), &
                self%n_data_points))
            do i = 1, size(self%vector_int32, 1)
                temp_vector_int32(i, :) = pack(self%vector_int32(i, :), mask)
            end do
            deallocate(self%vector_int32)
            self%vector_int32 = temp_vector_int32
        end if

        if (allocated(self%vector_int64)) then
            allocate(temp_vector_int64(size(self%vector_int64, 1), &
                self%n_data_points))
            do i = 1, size(self%vector_int64, 1)
                temp_vector_int64(i, :) = pack(self%vector_int64(i, :), mask)
            end do
            deallocate(self%vector_int64)
            self%vector_int64 = temp_vector_int64
        end if

        if (allocated(self%vector_real32)) then
            allocate(temp_vector_real32(size(self%vector_real32, 1), &
                self%n_data_points))
            do i = 1, size(self%vector_real32, 1)
                temp_vector_real32(i, :) = pack(self%vector_real32(i, :), mask)
            end do
            deallocate(self%vector_real32)
            self%vector_real32 = temp_vector_real32
        end if

        if (allocated(self%vector_real64)) then
            allocate(temp_vector_real64(size(self%vector_real64, 1), &
                self%n_data_points))
            do i = 1, size(self%vector_real64, 1)
                temp_vector_real64(i, :) = pack(self%vector_real64(i, :), mask)
            end do
            deallocate(self%vector_real64)
            self%vector_real64 = temp_vector_real64
        end if
    end subroutine data_mask
    
end module vtk