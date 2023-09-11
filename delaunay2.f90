!> @file delaunay2.f90
!> @brief Delaunay triangulation in 2D.
!> This module provides types and operations for 2D Delaunay triangulation.

module delaunay2
   use mesh_types
   use data_structures
   use error_handling
   implicit none

contains

   !> Adds a triangle to the mesh.
   subroutine insert_triangle(m, u, v, w)
      implicit none

      !> Mesh instance.
      type(mesh), intent(inout) :: m
      !> Points of the new triangle.
      integer(4) :: u, v, w
      !> Keys for edge map.
      integer(4), allocatable :: uv(:), vw(:), wu(:)
      !> Variable for error handling
      type(error_type) :: err

      if (u < 0 .or. u > m%n_points .or. v < 0 .or. v > m%n_points .or. w < 0 .or. w > m%n_points) then
         err%code = 100
         err%procedure = "insert_triangle"
         write (err%message, "(A, I0, X, I0, X, I0, A)") &
            "Triangle ", u, v, w, " contains an undefined point."
         call handle_error(err)
      end if

      if (u .eq. v .or. v .eq. w .or. w .eq. u) then
         err%code = 101
         err%procedure = "insert_triangle"
         write (err%message, "(A, I0, X, I0, X, I0, A)") &
            "Triangle ", u, v, w, " is degenerate."
         call handle_error(err)
      end if

      m%n_cells = m%n_cells + 1

      !> Extend cells array if needed.
      if (m%n_cells > size(m%cells, dim=2)) then
         call array_resize(m%cells, size(m%cells, dim=1), size(m%cells, dim=2) + cell_block_size)
      end if

      !> Add new triangle points to cells array.
      m%cells(1:3, m%n_cells) = [u, v, w]

      !> Update edge map with new edges and associated triangle.
      uv = [u, v]
      call m%edge_map%map_insert(uv, m%n_cells)
      vw = [v, w]
      call m%edge_map%map_insert(vw, m%n_cells)
      wu = [w, u]
      call m%edge_map%map_insert(wu, m%n_cells)

   end subroutine insert_triangle

   !> Finds the point opposite to the edge defined by (u, v).
   !> Returns -2 if the opposite point is not found.
   function adjacent(m, u, v) result(w)
      implicit none

      !> Mesh instance.
      type(mesh), intent(in) :: m
      !> Input edge defined by u and v.
      integer(4), intent(in) :: u, v
      !> Triangle containing the edge.
      integer(4) :: a
      !> Point opposite to the edge.
      integer(4) :: w
      !> Store triangle points locally for fewer array accesses.
      integer(4), dimension(3) :: t
      !> Variable for error handling
      type(error_type) :: err

      ! Fetch triangle that contains the edge (u, v)
      a = fetch_triangle(m, u, v)

      ! Error check
      if (a .eq. -1) then
         w = -1
         return
      else if (a .lt. 1 .or. a .gt. m%n_cells) then
         err%code = 110
         err%procedure = "adjacent"
         write (err%message, "(A, I0, A)") &
            "Triangle (", a, ") is undefined."
         call handle_error(err)
      end if

      ! Fetch points of triangle 'a'
      t = m%cells(:, a)

      ! Find and return the opposite point.
      if (all(t(1:2) .eq. [u, v])) then
         w = t(3)
      else if (all(t(2:3) .eq. [u, v])) then
         w = t(1)
      else if (all(t([3, 1]) .eq. [u, v])) then
         w = t(2)
      else
         ! Opposite point not found.
         w = -1
      end if
   end function adjacent

   !> @brief Fetch the triangle that has the points u, v, and w in a given mesh.
   !> @param m The mesh instance where the triangle is to be looked for.
   !> @param u Point 1 of the triangle.
   !> @param v Point 2 of the triangle.
   !> @param w Point 3 of the triangle (optional).
   !> @return Returns the triangle identifier. Returns -1 if the triangle does not exist.
   !> @throws Prints error message to stdout and stops execution if triangle has been incorrectly stored in the hash table.
   function fetch_triangle(m, u, v, w) result(a)
      implicit none

      !> Mesh instance.
      type(mesh), intent(in) :: m
      !> Input points u, v, and optional w that define the triangle.
      integer(4), intent(in) :: u, v
      integer(4), intent(in), optional :: w
      !> Key for edge map.
      integer(4), allocatable :: uv(:), vw(:), wu(:)
      !> Triangle containing the points.
      integer(4) :: a, b, c
      !> Variable for error handling
      type(error_type) :: err

      if (u < 0 .or. u > m%n_points .or. v < 0 .or. v > m%n_points) then
         err%code = 120
         err%procedure = "fetch_triangle"
         write (err%message, "(A, I0, X, I0, X, I0, A)") &
            "Triangle ", u, v, w, " contains an undefined point."
         call handle_error(err)
      end if

      if (u .eq. v) then
         err%code = 121
         err%procedure = "fetch_triangle"
         write (err%message, "(A, I0, X, I0, X, I0, A)") &
            "Triangle ", u, v, w, " is degenerate."
         call handle_error(err)
      end if

      if (present(w)) then
         if (w < 0 .or. w > m%n_points) then
            err%code = 120
            err%procedure = "fetch_triangle"
            write (err%message, "(A, I0, X, I0, X, I0, A)") &
               "Triangle ", u, v, w, " contains an undefined point."
            call handle_error(err)
         end if
      end if

      if (present(w)) then
         if (v .eq. w .or. w .eq. u) then
            err%code = 121
            err%procedure = "fetch_triangle"
            write (err%message, "(A, I0, X, I0, X, I0, A)") &
               "Triangle ", u, v, w, " is degenerate."
            call handle_error(err)
         end if
      end if

      uv = [u, v]
      a = m%edge_map%map_find(uv)

      if (a .eq. -1) then
         return
      end if

      if (present(w)) then
         vw = [v, w]
         b = m%edge_map%map_find(vw)
         wu = [w, u]
         c = m%edge_map%map_find(wu)

         if (a .ne. b .or. b .ne. c .or. c .ne. a) then
            err%code = 122
            err%procedure = "fetch_triangle"
            write (err%message, "(A, 3(2A, I0, A, I0, A, I0))") &
               "Inconsistent key-value pairs in hash map.", &
               NEW_LINE('a'), "(", u, ",", v, "):", a, &
               NEW_LINE('a'), "(", v, ",", w, "):", b, &
               NEW_LINE('a'), "(", w, ",", u, "):", c
            call handle_error(err)
         end if
      end if

   end function fetch_triangle

   !> Deletes a triangle from the mesh.
   !>
   !> @param m The mesh to which the triangle belongs.
   !> @param a Identifier of the triangle to delete.
   subroutine delete_triangle(m, a)
      implicit none
      !> Mesh instance.
      type(mesh), intent(inout) :: m
      !> Triangle identifier.
      integer(4), intent(in) :: a
      !> Triangle points and neighbouring triangles.
      integer(4) :: u, v, w
      !> Keys for edge map.
      integer(4), allocatable :: uv(:), vw(:), wu(:)
      !> Variable for error handling
      type(error_type) :: err

      if (a < 1 .or. a > m%n_cells) then
         err%code = 130
         err%procedure = "delete_triangle"
         write (err%message, "(A, I0, A)") &
            "Triangle (", a, ") is undefined."
         call handle_error(err)
      end if

      !> Fetch points and neighbouring triangles of the triangle with identifier a.
      u = m%cells(1, a)
      v = m%cells(2, a)
      w = m%cells(3, a)

      !> Delete edges from the edge map.
      uv = [u, v]
      call m%edge_map%map_delete(uv)
      vw = [v, w]
      call m%edge_map%map_delete(vw)
      wu = [w, u]
      call m%edge_map%map_delete(wu)

      !> Replace this triangle with the last in the list and decrement the count.
      m%cells(:, a) = m%cells(:, m%n_cells)
      m%cells(:, m%n_cells) = -1
      m%n_cells = m%n_cells - 1

      u = m%cells(1, a)
      v = m%cells(2, a)
      w = m%cells(3, a)

      uv = [u, v]
      call m%edge_map%map_insert(uv, a)
      vw = [v, w]
      call m%edge_map%map_insert(vw, a)
      wu = [w, u]
      call m%edge_map%map_insert(wu, a)

   end subroutine delete_triangle

   !> @brief Determine if a point conflicts with a given triangle
!>
!> This function checks if a given point p lies within the circumcircle
!> of the triangle defined by points u, v, and w. If any point index is zero,
!> it checks if the point lies in the outer half-plane defined by the other two points.
!>
!> @param[in, out] m Mesh instance containing points and other mesh data
!> @param[in] p Index of the point to be checked
!> @param[in] u Index of the first point of the triangle
!> @param[in] v Index of the second point of the triangle
!> @param[in] w Index of the third point of the triangle
!>
!> @return True if there is a conflict, i.e., point p lies inside the circumcircle or outer half-plane, false otherwise
!>
   function is_conflict(m, p, u, v, w) result(conflict)
      use geometry  ! Assuming geometry module defines in_outer_halfplane and in_circle
      implicit none

      type(mesh), intent(inout) :: m  ! Mesh instance
      integer(4), intent(in) :: p, u, v, w  ! Point and triangle points
      real :: px, py, ux, uy, vx, vy, wx, wy  ! Coordinates for point and points
      logical :: conflict  ! Result: true if there is a conflict
      !> Variable for error handling
      type(error_type) :: err

      if (u < 0 .or. u > m%n_points .or. v < 0 .or. v > m%n_points .or. w < 0 .or. w > m%n_points) then
         err%code = 140
         err%procedure = "is_conflict"
         write (err%message, "(A, 3(I0, X), A)"), "Triangle ", u, v, w, "contains an undefined point."
         call handle_error(err)
      end if

      if (p < 0 .or. p > m%n_points) then
         err%code = 141
         err%procedure = "is_conflict"
         write (err%message, "(A, (I0, X), A)"), "Point ", p, " is undefined."
         call handle_error(err)
      end if

      if (p .eq. 0) then
         ux = m%points(1, u)
         uy = m%points(2, u)
         vx = m%points(1, v)
         vy = m%points(2, v)
         wx = m%points(1, w)
         wy = m%points(2, w)
         conflict = orient2([ux, uy], [wx, wy], [vx, vy]) .or. on_edge([ux, uy], [wx, wy], [vx, vy])
         return
      end if

      ! Fetch the coordinates of point p
      px = m%points(1, p)
      py = m%points(2, p)

      ! If any of the triangle points is zero, check for conflict using outer half-plane
      if (any([u, v, w] .eq. 0)) then
         ! Point u is zero; check using points v and w
         if (u .eq. 0) then
            vx = m%points(1, v)
            vy = m%points(2, v)
            wx = m%points(1, w)
            wy = m%points(2, w)
            conflict = orient2([px, py], [vx, vy], [wx, wy]) .or. on_edge([px, py], [vx, vy], [wx, wy])
            ! Point v is zero; check using points u and w
         else if (v .eq. 0) then
            ux = m%points(1, u)
            uy = m%points(2, u)
            wx = m%points(1, w)
            wy = m%points(2, w)
            conflict = orient2([px, py], [wx, wy], [ux, uy]) .or. on_edge([px, py], [wx, wy], [ux, uy])
            ! Point w is zero; check using points u and v
         else
            ux = m%points(1, u)
            uy = m%points(2, u)
            vx = m%points(1, v)
            vy = m%points(2, v)
            conflict = orient2([px, py], [ux, uy], [vx, vy]) .or. on_edge([px, py], [ux, uy], [vx, vy])
         end if
      else
         ! All points are non-zero; check for conflict using circumcircle
         ux = m%points(1, u)
         uy = m%points(2, u)
         vx = m%points(1, v)
         vy = m%points(2, v)
         wx = m%points(1, w)
         wy = m%points(2, w)
         conflict = in_circle([px, py], [ux, uy], [vx, vy], [wx, wy])
      end if

   end function is_conflict

   !> @brief Recursively dig a cavity in the mesh for the insertion of a new point.
   !>
   !> This subroutine is designed to prepare the mesh for the insertion of a new point
   !> by recursively digging a cavity. The initial triangle for the cavity is specified
   !> by its points u, v, and the point to be inserted is p.
   !>
   !> @param[in,out] m The mesh instance that will be modified by this subroutine.
   !> @param[in] p The index of the point to be inserted into the mesh.
   !> @param[in] u The index of the first point of the initial triangle.
   !> @param[in] v The index of the second point of the initial triangle.
   !>
   !> @note The subroutine uses recursion to continue the digging process if needed.
   !>
   recursive subroutine dig_cavity(m, p, u, v)
      implicit none

      !> Mesh instance to be modified.
      type(mesh), intent(inout) :: m
      !> Index of the point to be inserted.
      integer(4), intent(in) :: p
      !> Indices of the points of the initial triangle.
      integer(4), intent(in) :: u, v
      !> Additional variables
      integer(4) :: w, a
      !> Variable for error handling
      type(error_type) :: err

      !> Find the third point adjacent to u and v.
      w = adjacent(m, v, u)

      !> If the adjacent function returns -1, terminate the recursion.
      if (w .eq. -1) then
         return
      end if

      !> If there is a conflict with the third point, continue digging.

      if (is_conflict(m, w, p, u, v)) then
         !> Fetch the triangle for deletion.
         a = fetch_triangle(m, v, u, w)
         !> Delete the conflicting triangle.
         call delete_triangle(m, a)
         !> Continue digging recursively.
         call dig_cavity(m, p, u, w)
         call dig_cavity(m, p, w, v)
         return
      else

         !> If there is no conflict, add a new triangle.
         call insert_triangle(m, p, u, v)
         return
      end if

   end subroutine dig_cavity

   !> @brief Add a point to the mesh and retriangulate.
   !> @param m Mesh instance to modify.
   !> @param p Point to add.
   !> @param u, v, w The points of the triangle containing p.
   subroutine insert_point(m, p, u, v, w)
      implicit none

      !> Mesh instance.
      type(mesh), intent(inout) :: m
      !> Point to add.
      integer(4), intent(in) :: p
      !> Points of the triangle containing p.
      integer(4), intent(in) :: u, v, w
      integer(4) :: a
      !> Variable for error handling
      type(error_type) :: err

      if (u < 0 .or. u > m%n_points .or. v < 0 .or. v > m%n_points .or. w < 0 .or. w > m%n_points) then
         err%code = 160
         err%procedure = "insert_point"
         write (err%message, "(A, I0, X, I0, X, I0, A)") &
            "Triangle ", u, v, w, " contains an undefined point."
         call handle_error(err)
      end if

      a = fetch_triangle(m, u, v, w)

      if (a .lt. 1) then
         err%code = 161
         err%procedure = "insert_point"
         write (err%message, "(A, 3(I0, X), A)"), "Triangle ", u, v, w, "is undefined."
         call handle_error(err)
      end if

      call delete_triangle(m, a)
      call dig_cavity(m, p, u, v)
      call dig_cavity(m, p, v, w)
      call dig_cavity(m, p, w, u)
   end subroutine insert_point

   !> @brief Initialize the mesh with initial points and triangles.
!>
!> This subroutine initializes a mesh with given points u, v, and w.
!> It allocates memory for cell data and initializes the edge map.
!> It also checks the orientation of the points and adds initial triangles.
!>
!> @param m Mesh instance to initialize.
   subroutine initialise_mesh(m)
      use geometry
      use data_structures
      implicit none

      !> Mesh instance
      type(mesh), intent(inout) :: m
      integer :: u, v, w, i
      real, dimension(2) :: temp
      real :: ux, uy, vx, vy, wx, wy
      integer :: status  ! Variable to hold allocation status
      !> Variable for error handling
      type(error_type) :: err

      ! call permute_array2d(m%points, 2)

      ! Initialize initial points
      u = 1
      v = 2
      w = 3

      ux = m%points(1, u)
      uy = m%points(2, u)
      vx = m%points(1, v)
      vy = m%points(2, v)
      wx = m%points(1, w)
      wy = m%points(2, w)

      do while (.not. orient2([ux, uy], [vx, vy], [wx, wy]))

         ! call permute_array2d(m%points, 2)

         ux = m%points(1, u)
         uy = m%points(2, u)
         vx = m%points(1, v)
         vy = m%points(2, v)
         wx = m%points(1, w)
         wy = m%points(2, w)

      end do

      ! Allocate space for cells and check for successful allocation
      allocate (m%cells(3, cell_block_size), stat=status)
      if (status /= 0) then
         err%code = 170
         err%procedure = "initialise_mesh"
         write (err%message, "(A)") &
            "Failed to allocate cells."
         call handle_error(err)
      end if

      ! Initialize edge map
      call m%edge_map%map_initialize(map_block_size)

      ! Initialize cells to a negative value
      m%cells = -1

      ! Add initial triangles to the mesh
      call insert_triangle(m, u, v, w)
      call insert_triangle(m, v, u, 0)
      call insert_triangle(m, w, v, 0)
      call insert_triangle(m, u, w, 0)

   end subroutine initialise_mesh

   subroutine delete_dummy_triangles(m)
      use geometry
      implicit none
      type(mesh), intent(inout) :: m
      integer :: i

      i = 1
      do while (i <= m%n_cells)
         if (any(m%cells(1:3, i) .eq. 0)) then
            call delete_triangle(m, i)
            i = i - 1  ! Decrement the index to stay at the same place for the next iteration
         end if
         i = i + 1  ! Increment the index for the next iteration
      end do

   end subroutine delete_dummy_triangles

   !> @brief Implements the Bowyer-Watson algorithm for Delaunay triangulation.
!>
!> This subroutine iteratively adds each point to the existing triangulated mesh
!> and retriangulates as needed to maintain the Delaunay condition.
!>
!> @param m Mesh instance to modify.
   subroutine bowyer_watson2(m)
      use geometry  ! Use geometry module for additional functionalities like is_conflict
      implicit none

      !> Mesh instance
      type(mesh), intent(inout) :: m
      integer :: u, v, w, i, j

      call initialise_mesh(m)

      ! Iterate over all points starting from the 4th point
      do i = 4, m%n_points

         ! Check conflict for each existing triangle
         do j = 1, m%n_cells

            ! Get points of the j-th triangle
            u = m%cells(1, j)
            v = m%cells(2, j)
            w = m%cells(3, j)

            ! Check if the point lies inside the circumcircle of triangle uvw
            if (is_conflict(m, i, u, v, w)) then
               ! If so, modify mesh to insert new point
               call insert_point(m, i, u, v, w)
               exit  ! Exit the loop once point is added
            end if
         end do
      end do

      call delete_dummy_triangles(m)

   end subroutine bowyer_watson2

end module delaunay2
