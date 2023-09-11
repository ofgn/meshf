!> \file io.f90
!> \brief Read .poly file and populate mesh data.
!> \author ofgn
!> \version 0.1
!> \date 2023-09-09

module io
   use mesh_types
   use error_handling
   implicit none

   integer, parameter :: max_data_components = 8

contains

   ! ============================
   ! INPUT
   ! ============================

   !> \brief Read a .poly file to populate a mesh data structure.
!> \param[out] m The mesh data structure to be populated.
!> \param[in] file_path The path to the .poly file.
!>
!> This subroutine reads a .node file to populate a given mesh data structure.
!> It performs multiple checks to ensure the file and its contents are valid.
   subroutine read_poly(m, file_path)
      implicit none

      type(mesh) :: m
      integer :: i, j, io, stat
      character(len=*) :: file_path
      type(error_type) :: err
      integer :: boundary_mode

      ! Open the file
      open (newunit=io, file=file_path, status="old", action="read", iostat=stat)

      ! Error handling: unable to open file
      if (stat /= 0) then
         err%code = 1001
         err%procedure = "read_poly"
         write (err%message, "(3A)") "Unable to open file ", file_path, "."
         call handle_error(err)
         return
      end if

      ! Read header information
      read (io, *, iostat=stat) m%n_points, m%n_dims, m%n_pdata, boundary_mode

      ! Error handling for header
      if (stat /= 0) then
         err%code = 1002
         err%procedure = "read_poly"
         write (err%message, "(A)") "Failed to read .poly file header."
         call handle_error(err)
         close (io)
         return
      end if

      ! Validate header info
      if (m%n_dims == 2 .and. m%n_points < 3) then
         err%code = 1003
         err%procedure = "read_poly"
         err%message = "Error: Triangulation requires at least 3 points."
         call handle_error(err)
      else if (m%n_dims == 3 .and. m%n_points < 4) then
         err%code = 1004
         err%procedure = "read_poly"
         err%message = "Error: Tetrahedralisation requires at least 4 points."
         call handle_error(err)
      else if (m%n_dims < 2 .or. m%n_dims > 3) then
         err%code = 1005
         err%procedure = "read_poly"
         err%message = "Error: .poly file must be 2 or 3 dimensional."
         call handle_error(err)
      end if

      ! Memory allocation
      allocate (m%points(m%n_dims, m%n_points), m%pdata(m%n_pdata, m%n_points), m%pflags(m%n_points))

      ! Loop through points and read data
      do i = 1, m%n_points
         if (boundary_mode .and. m%n_pdata > 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pdata(1:m%n_pdata, i), m%pflags(i)
         else if (boundary_mode .and. m%n_pdata == 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pflags(i)
         else if (.not. boundary_mode .and. m%n_pdata > 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pdata(1:m%n_pdata, i)
            m%pflags(i) = 0.0d0
         else if (.not. boundary_mode .and. m%n_pdata == 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i)
            m%pflags(i) = 0.0d0
         else
            err%code = 1006
            err%procedure = "read_poly"
            err%message = "Error: Cannot read point."
            call handle_error(err)
            return
         end if

         if (i == 1) then
            m%x_min = m%points(1, i)
            m%x_max = m%points(1, i)
            m%y_min = m%points(2, i)
            m%y_max = m%points(2, i)
            if (m%n_dims .eq. 3) then
               m%z_min = m%points(3, i)
               m%z_max = m%points(3, i)
            else
               m%z_min = 0.0d0
               m%z_max = 0.0d0
            end if
         else
            m%x_min = min(m%x_min, m%points(1, i))
            m%x_max = max(m%x_max, m%points(1, i))
            m%y_min = min(m%y_min, m%points(2, i))
            m%y_max = max(m%y_max, m%points(2, i))
            if (m%n_dims .eq. 3) then
               m%z_min = min(m%z_min, m%points(3, i))
               m%z_max = max(m%z_max, m%points(3, i))
            end if
         end if
      end do
      ! Close file
      close (io)
   end subroutine read_poly

   !> \brief Read a .node file to populate a mesh data structure.
!> \param[out] m The mesh data structure to be populated.
!> \param[in] file_path The path to the .node file.
!>
!> This subroutine reads a .node file to populate a given mesh data structure.
!> It performs multiple checks to ensure the file and its contents are valid.
   subroutine read_node(m, file_path)
      implicit none

      type(mesh) :: m
      integer :: i, j, io, stat
      character(len=*) :: file_path
      type(error_type) :: err
      integer :: boundary_mode

      ! Open the file
      open (newunit=io, file=file_path, status="old", action="read", iostat=stat)

      ! Error handling: unable to open file
      if (stat /= 0) then
         err%code = 1001
         err%procedure = "read_node"
         write (err%message, "(3A)") "Unable to open file ", file_path, "."
         call handle_error(err)
         return
      end if

      ! Read header information
      read (io, *, iostat=stat) m%n_points, m%n_dims, m%n_pdata, boundary_mode

      ! Error handling for header
      if (stat /= 0) then
         err%code = 1002
         err%procedure = "read_node"
         write (err%message, "(A)") "Failed to read .node file header."
         call handle_error(err)
         close (io)
         return
      end if

      ! Validate header info
      if (m%n_dims == 2 .and. m%n_points < 3) then
         err%code = 1003
         err%procedure = "read_node"
         err%message = "Error: Triangulation requires at least 3 points."
         call handle_error(err)
      else if (m%n_dims == 3 .and. m%n_points < 4) then
         err%code = 1004
         err%procedure = "read_node"
         err%message = "Error: Tetrahedralisation requires at least 4 points."
         call handle_error(err)
      else if (m%n_dims < 2 .or. m%n_dims > 3) then
         err%code = 1005
         err%procedure = "read_node"
         err%message = "Error: .poly file must be 2 or 3 dimensional."
         call handle_error(err)
      end if

      if (m%n_pdata .lt. 0) then
         err%code = 1006
         err%procedure = "read_node"
         write (err%message, "(A)") "Invalid number of data components."
         call handle_error(err)
      else if (m%n_pdata .gt. max_data_components) then
         err%code = 1007
         err%procedure = "read_node"
         write (err%message, "(A, I0, A)") "Number of data components exceeds ", max_data_components, "."
         call handle_error(err)
      end if

      if ((boundary_mode .lt. 0) .or. (boundary_mode .gt. 1)) then
         err%code = 1008
         err%procedure = "read_node"
         write (err%message, "(A, I0, A)") &
            "Input file must denote boundary flags are present with 1, elsewise 0."
         call handle_error(err)
      end if

      ! Memory allocation
      allocate (m%points(m%n_dims, m%n_points), m%pdata(m%n_pdata, m%n_points), m%pflags(m%n_points))

      ! Loop through points and read data
      do i = 1, m%n_points
         if ((boundary_mode .eq. 1) .and. m%n_pdata > 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pdata(1:m%n_pdata, i), m%pflags(i)
         else if ((boundary_mode .eq. 1) .and. m%n_pdata == 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pflags(i)
         else if ((boundary_mode .eq. 0) .and. m%n_pdata > 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i), m%pdata(1:m%n_pdata, i)
            m%pflags(i) = 0.0d0
         else if ((boundary_mode .eq. 0) .and. m%n_pdata == 0) then
            read (io, *, iostat=stat) j, m%points(1:m%n_dims, i)
            m%pflags(i) = 0.0d0
         else
            err%code = 1006
            err%procedure = "read_node"
            err%message = "Error: Cannot read point."
            call handle_error(err)
            return
         end if

         if (i == 1) then
            m%x_min = m%points(1, i)
            m%x_max = m%points(1, i)
            m%y_min = m%points(2, i)
            m%y_max = m%points(2, i)
            if (m%n_dims .eq. 3) then
               m%z_min = m%points(3, i)
               m%z_max = m%points(3, i)
            else
               m%z_min = 0.0d0
               m%z_max = 0.0d0
            end if
         else
            m%x_min = min(m%x_min, m%points(1, i))
            m%x_max = max(m%x_max, m%points(1, i))
            m%y_min = min(m%y_min, m%points(2, i))
            m%y_max = max(m%y_max, m%points(2, i))
            if (m%n_dims .eq. 3) then
               m%z_min = min(m%z_min, m%points(3, i))
               m%z_max = max(m%z_max, m%points(3, i))
            end if
         end if
      end do

      ! Close file
      close (io)
   end subroutine read_node

   subroutine read_ele(m, file_path)
      implicit none

      type(mesh) :: m
      integer :: i, j, io, stat
      character(len=*) :: file_path
      type(error_type) :: err
      integer :: boundary_mode

      ! Open the file
      open (newunit=io, file=file_path, status="old", action="read", iostat=stat)

      ! Error handling: unable to open file
      if (stat /= 0) then
         err%code = 1001
         err%procedure = "read_ele"
         write (err%message, "(3A)") "Unable to open file ", file_path, "."
         call handle_error(err)
         return
      end if

      ! Read header information
      read (io, *, iostat=stat) m%n_cells, m%cell_type, m%n_cdata

      ! Error handling for header
      if (stat /= 0) then
         err%code = 1002
         err%procedure = "read_ele"
         write (err%message, "(A)") "Failed to read .ele file header."
         call handle_error(err)
         close (io)
         return
      end if

      ! Validate header info
      if (m%n_dims .ne. 2 .and. m%cell_type .eq. 3) then
         err%code = 1023
         err%procedure = "read_ele"
         err%message = "Error: Triangular mesh must be 2 dimensional."
         call handle_error(err)
      else if (m%n_dims .ne. 3 .and. m%cell_type .eq. 4) then
         err%code = 1024
         err%procedure = "read_ele"
         err%message = "Error: Tetrahedral mesh must be 3 dimensional."
         call handle_error(err)
      else if (m%cell_type .lt. 3 .or. m%cell_type .gt. 4) then
         err%code = 1025
         err%procedure = "read_ele"
         write (err%message, "(A, I0, A, I0, A)") "Error: Cells are incompatible with current mesh points. ", &
            m%n_dims, " dimensional mesh with type ", m%cell_type, " cells."
         call handle_error(err)
      end if

      if (m%n_cdata .lt. 0) then
         err%code = 1026
         err%procedure = "read_ele"
         write (err%message, "(A)") "Invalid number of data components."
         call handle_error(err)
      else if (m%n_pdata .gt. max_data_components) then
         err%code = 1027
         err%procedure = "read_ele"
         write (err%message, "(A, I0, A)") "Number of data components exceeds ", max_data_components, "."
         call handle_error(err)
      end if

      allocate (m%cells(m%cell_type, m%n_cells), m%cdata(1:m%n_cdata, m%n_cells))

      ! Loop through cells and read data
      do i = 1, m%n_cells
         if (m%n_pdata > 0) then
            read (io, *, iostat=stat) j, m%cells(1:m%cell_type, i), m%cdata(1:m%n_cdata, i)
         else if (m%n_cdata == 0) then
            read (io, *, iostat=stat) j, m%cells(1:m%cell_type, i)
         else
            err%code = 1036
            err%procedure = "read_ele"
            err%message = "Error: Cannot read cell."
            call handle_error(err)
         end if
      end do

      ! Close file
      close (io)

   end subroutine read_ele

   !> @brief Appends new cell data from an ASCII file to an existing mesh structure.
!>
!> @param m         The mesh structure to which new cell data will be appended.
!> @param filename  The name of the ASCII file containing the new cell data.
   subroutine read_cell_data_ascii(m, filename)
      type(mesh), intent(inout) :: m
      character(len=*), intent(in) :: filename

      ! Declare local variables
      integer :: io, i, n_cells, n_cdata_fields
      real, allocatable :: new_cdata(:, :)
      real, allocatable :: temp_cdata(:, :)
      type(error_type) :: err

      ! Open file
      open (newunit=io, file=filename, status="old", action="read")

      ! Read header to get the number of cells and data fields
      read (io, *) n_cells, n_cdata_fields

      ! Check for mismatch in number of cells
      if (n_cells /= m%n_cells) then
         err%code = 1040
         err%procedure = "read_cell_data_ascii"
         err%message = "Error: Number of cells in the ASCII file does not match the mesh."
         call handle_error(err)
      end if

      ! Allocate and read the new cell data
      allocate (new_cdata(n_cdata_fields, n_cells))
      do i = 1, n_cells
         read (io, *) new_cdata(:, i)
      end do
      close (io)

      ! Append new data to existing data
      if (allocated(m%cdata)) then
         ! Allocate space for the new concatenated array
         allocate (temp_cdata(m%n_cdata + n_cdata_fields, m%n_cells))

         ! Copy the existing data into the new array
         temp_cdata(1:m%n_cdata, :) = m%cdata

         ! Append the new data into the new array
         temp_cdata(m%n_cdata + 1:m%n_cdata + n_cdata_fields, :) = new_cdata

         ! De-allocate the old array and allocate space for the new concatenated array
         deallocate (m%cdata)
         allocate (m%cdata(m%n_cdata + n_cdata_fields, m%n_cells))

         ! Copy the concatenated data back into m%cdata
         m%cdata = temp_cdata

         ! De-allocate the temporary array
         deallocate (temp_cdata)
      else
         m%cdata = new_cdata
      end if

      ! Update the number of cell data fields in the mesh
      m%n_cdata = m%n_cdata + n_cdata_fields
   end subroutine read_cell_data_ascii

! ============================
! OUTPUT
! ============================

!> @brief Prints mesh vertices.
!>
!> This subroutine prints out the vertices of a given mesh. The output varies based
!> on the number of dimensions of the mesh.
!>
!> @param[in] m The mesh data.
   subroutine print_vertices(m)
      type(mesh), intent(in) :: m  !< Input mesh
      integer :: i                !< Loop variable for iterating through points

      ! Print general mesh information
      print *, "Mesh has ", m%n_points, " vertices."
      print *, "Dimensions: ", m%n_dims

      ! Loop through each vertex and print its coordinates
      do i = 1, m%n_points
         if (m%n_dims == 3) then
            ! For 3D meshes, print x, y, and z coordinates
            write (*, '(A, I4, A, 3F10.4)') "Vertex ", i, ": ", m%points(1:3, i)
         else
            ! For 2D meshes, print x and y coordinates
            write (*, '(A, I4, A, 2F10.4)') "Vertex ", i, ": ", m%points(1:2, i)
         end if
      end do

      ! Print range information
      print *, "X range: [", m%x_min, ",", m%x_max, "]"
      print *, "Y range: [", m%y_min, ",", m%y_max, "]"
      if (m%n_dims == 3) then
         print *, "Z range: [", m%z_min, ",", m%z_max, "]"
      end if
   end subroutine print_vertices

!> \brief Writes mesh to a legacy ASCII VTK file.
!> \param[in] m The mesh data structure.
!> \param[in] filename The output VTK filename.
!>
!> This subroutine writes a mesh to an ASCII VTK file, either as a
!> triangular unstructured grid or as a tetrahedral unstructured grid,
!> determined by m%cell_type.
   subroutine write_vtk_ascii(m, filename)
      type(mesh), intent(in) :: m
      character(len=*) :: filename
      integer :: io, i, j

      ! Open file for writing
      open (newunit=io, file=filename, status="replace", action="write")

      ! Write VTK header
      write (io, '(A)') "# vtk DataFile Version 3.0"
      write (io, '(A)') "Generated by MeshFort"
      write (io, '(A)') "ASCII"
      write (io, '(A)') "DATASET UNSTRUCTURED_GRID"

      ! Write vertices
      write (io, '("POINTS ", I0, " double")') m%n_points
      do i = 1, m%n_points
         if (m%n_dims == 2) then
            write (io, '(2F17.6, F17.6)') m%points(1, i), m%points(2, i), 0.d0
         else
            write (io, '(3F17.6)') m%points(:, i)
         end if
      end do

      ! Write cells
      if (m%cell_type == 3) then  ! Triangular grid
         write (io, '("CELLS ", I0, " ", I0)') m%n_cells, 4*m%n_cells
         do i = 1, m%n_cells
            write (io, '(I0, 3X, 3(I0, X))') 3, m%cells(:, i) - 1
         end do
         write (io, '("CELL_TYPES ", I0)') m%n_cells
         write (io, '(I0)') (5, i=1, m%n_cells)  ! VTK code 5 for triangles
      else if (m%cell_type == 4) then  ! Tetrahedral grid
         write (io, '("CELLS ", I0, " ", I0)') m%n_cells, 5*m%n_cells
         do i = 1, m%n_cells
            write (io, '(I0, 3X, 4(I0, X))') 4, m%cells(:, i) - 1
         end do
         write (io, '("CELL_TYPES ", I0)') m%n_cells
         write (io, '(I0)') (10, i=1, m%n_cells)  ! VTK code 10 for tetrahedra
      else
         ! Error handling for invalid cell_type
         print *, "Error: Invalid cell type."
         close (io)
         return
      end if

      ! Append point data if available
      if (allocated(m%pdata) .and. m%n_pdata > 0) then
         write (io, '("POINT_DATA ", I0)') m%n_points
         do j = 1, m%n_pdata
            write (io, '(A, I0)') "SCALARS point_data", j
            write (io, '(A)') "double 1"
            write (io, '(A)') "LOOKUP_TABLE default"
            do i = 1, m%n_points
               write (io, '(F17.6)') m%pdata(j, i)
            end do
         end do
      end if

      ! Append cell data if available
      if (allocated(m%cdata) .and. m%n_cdata > 0) then
         write (io, '("CELL_DATA ", I0)') m%n_cells
         do j = 1, m%n_cdata
            write (io, '(A, I0)') "SCALARS cell_data", j
            write (io, '(A)') "double 1"
            write (io, '(A)') "LOOKUP_TABLE default"
            do i = 1, m%n_cells
               write (io, '(F17.6)') m%cdata(j, i)
            end do
         end do
      end if

      ! Close the file
      close (io)

   end subroutine write_vtk_ascii

end module io
