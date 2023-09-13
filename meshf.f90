program meshf
   use data_structures
   use mesh_types
   use io
   use geometry
   use delaunay2
   implicit none

   type(mesh) :: m
   character(len=128) :: base_name, file_name
   integer :: i, n_args
   character(len=32) :: arg

   ! Get the number of command line arguments
   n_args = command_argument_count()
   call get_command_argument(1, arg)

   select case (arg)
   case ("-tetgen2vtk")
      call get_command_argument(2, base_name)
      call read_node(m, trim(base_name)//".node")
      call read_ele(m, trim(base_name)//".ele")
      call get_command_argument(3, file_name)
      call write_vtk_ascii(m, trim(file_name))
   end select


end program meshf
