module mesh_types
   use data_structures
   implicit none

   !> Size of cell blocks for array resizing.
   integer(4), parameter :: cell_block_size = 4096
   integer(4), parameter :: map_block_size = 12288

   !> @brief Mesh data type to represent the triangulation.
   type :: mesh
      integer(4) :: n_dims !> Number of dimensions.
      integer(4) :: n_points !> Number of points.
      integer(4) :: n_cells !> Number of cells.
      integer(4) :: cell_type !> Number of points per cell.
      integer(4) :: n_pdata !> Number of point data fields.
      integer(4) :: n_cdata !> Number of cell data fields.
      integer(4) :: n_constrs !> Number of constraints (for CDT).
    !   logical(1) :: b_flags !> Flag points on boundaries.
      integer(4) :: n_bounds !> Number of boundaries.
      integer(4) :: n_voids !> Number of voids (holes).
      real, allocatable, dimension(:, :) :: points !> Coordinates of points.
      real, allocatable, dimension (:, :) :: pdata
      character(len=16), allocatable, dimension (:) :: pdata_labels
      integer, allocatable, dimension (:) :: pflags
      integer(4), allocatable, dimension(:, :) :: cells !> Cells (triangles) defined by points.
      real, allocatable, dimension (:, :) :: cdata
      character(len=16), allocatable, dimension (:) :: cdata_labels
      type(hash_map) :: edge_map !> Hash map to store edges.
      real :: x_min, x_max !> Minimum and Maximum X coordinate.
      real :: y_min, y_max !> Minimum and Maximum Y coordinate.
      real :: z_min, z_max !> Minimum and Maximum Z coordinate.
   end type mesh
end module mesh_types
