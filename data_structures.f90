!> @file data_structures.f90
!> @brief Implementation of simple data structures in Fortran.
!> @author ofgn
!> @date 2023-09-01

module data_structures
   implicit none

   !> @brief Node for singly linked list.
   !>
   !> Each node contains a dynamic integer array and a pointer to the next node.
   type :: list_node
      integer, allocatable :: value(:)    !< Value held by this node.
      type(list_node), pointer :: next => null()  !< Pointer to the next node.
   end type list_node

   !> @brief Singly linked list of integer arrays.
   !>
   !> Contains a pointer to the head node and the list length.
   type :: list
      type(list_node), pointer :: head => null()  !< Pointer to the head node.
      integer :: length = 0  !< Number of nodes in the list.
   contains
      procedure :: list_append  !< Appends a value to the end.
      procedure :: list_pop     !< Removes and returns the last value.
   end type list

   !> @brief Defines the hash map node for storing key-value pairs.
   type :: hash_map_node
      integer, allocatable :: key(:) !< The key, represented as a 1D integer array of arbitrary length.
      integer :: value !< The value associated with the key.
      type(hash_map_node), pointer :: next => null()  !< Pointer to the next node in the bucket.
   end type hash_map_node

   !> @brief Defines a hash map bucket, essentially a linked list of hash_map_nodes.
   type :: hash_map_bucket
      type(hash_map_node), pointer :: head => null() !< Pointer to the head node of the bucket.
   end type hash_map_bucket

   !> @brief Defines the hash map data structure.
   type :: hash_map
      type(hash_map_bucket), allocatable :: buckets(:)  !< Array of hash_map_bucket types to contain the key-value pairs.
      integer :: count !< Number of inserted key-value pairs.
      integer :: size !< Current maximum number of key-value pairs
   contains
      procedure :: map_initialize, map_clear, map_insert, map_find, map_delete, map_resize  !< Procedures for hash_map.
   end type hash_map

contains

   ! ============================
   ! LINKED LIST
   ! ============================

   !> @brief Appends a value to the list end.
   !> @param self The list.
   !> @param value The integer array to add.
   subroutine list_append(self, value)
      class(list), intent(inout) :: self
      integer, intent(in) :: value(:)
      type(list_node), pointer :: new_node, curr_node

      allocate (new_node)
      if (.not. associated(new_node)) then
         print *, "Linked List Error: Memory allocation failed."
         return
      end if

      new_node%value = value
      new_node%next => null()

      if (associated(self%head)) then
         curr_node => self%head
         do while (associated(curr_node%next))
            curr_node => curr_node%next
         end do
         curr_node%next => new_node
      else
         self%head => new_node
      end if

      self%length = self%length + 1
   end subroutine list_append

   !> @brief Removes the last node and returns its value.
   !> @param self The list.
   !> @param value Value of the last node.
   subroutine list_pop(self, value)
      class(list), intent(inout) :: self
      integer, allocatable, intent(out) :: value(:)
      type(list_node), pointer :: curr_node, prev_node

      if (.not. associated(self%head)) then
         print *, "Linked List Error: The list is empty."
         return
      end if

      curr_node => self%head
      prev_node => null()
      do while (associated(curr_node%next))
         prev_node => curr_node
         curr_node => curr_node%next
      end do

      value = curr_node%value

      if (associated(prev_node)) then
         prev_node%next => null()
      else
         self%head => null()
      end if

      deallocate (curr_node)
      self%length = self%length - 1
   end subroutine list_pop

   ! ============================
   ! HASH MAP
   ! ============================

   !> @brief Initialize the hash map with a given size.
   subroutine map_initialize(self, size)
      class(hash_map), intent(inout) :: self
      integer, intent(in) :: size
      integer :: stat

      if (size <= 0) then
         write (*, *) "Error (hash_map): Size must be a positive integer."
         return
      end if

      allocate (self%buckets(size), stat=stat)
      self%size = size
      if (stat /= 0) then
         write (*, *) "Error (hash_map): Failed to initialize hash map."
         return
      end if
   end subroutine map_initialize

!> @brief Insert a key-value pair into the hash map.
   subroutine map_insert(self, key, value)
      class(hash_map), intent(inout) :: self
      integer, intent(in), allocatable :: key(:)
      integer, intent(in) :: value
      integer :: index, alloc_stat
      real(8) :: load_factor  ! New variable for load factor
      type(hash_map_node), pointer :: new_node

      index = map_hash(key, size(self%buckets))
      if (index == -1) return

      allocate (new_node, stat=alloc_stat)
      if (alloc_stat /= 0) then
         write (*, *) "Error (hash_map): Failed to insert new node."
         return
      end if

      allocate (new_node%key(size(key)), stat=alloc_stat)
      if (alloc_stat /= 0) then
         write (*, *) "Error (hash_map): Failed to allocate memory for key."
         return
      end if

      new_node%key = key
      new_node%value = value
      new_node%next => self%buckets(index)%head
      self%buckets(index)%head => new_node

      ! Increase the count of total items
      self%count = self%count + 1

      ! Calculate the load factor
      load_factor = real(self%count)/real(size(self%buckets))

      ! Check if resizing is needed
      if (load_factor > 0.7) then
         call map_resize(self, self%count + self%count/4)
      end if
   end subroutine map_insert
!> @brief Delete a key-value pair from the hash map.
   subroutine map_delete(self, key)
      class(hash_map), intent(inout) :: self
      integer, intent(in), allocatable :: key(:)
      integer :: index
      type(hash_map_node), pointer :: current, prev

      index = map_hash(key, size(self%buckets))
      if (index == -1) return

      prev => null()
      current => self%buckets(index)%head

      do while (associated(current))
         if (all(current%key == key)) then
            if (associated(prev)) then
               prev%next => current%next
            else
               self%buckets(index)%head => current%next
            end if
            deallocate (current%key)
            deallocate (current)
            exit
         end if
         prev => current
         current => current%next
      end do
   end subroutine map_delete

!> @brief Helper function for robust hash calculation.
   integer function hash(key, bucket_size)
      integer, intent(in), allocatable :: key(:)
      integer, intent(in) :: bucket_size
      integer :: i

      if (bucket_size <= 0) then
         write (*, *) "Error (hash_map): Bucket size must be a positive integer."
         hash = -1
         return
      end if

      hash = 0
      do i = 1, size(key)
         hash = mod(hash + key(i), bucket_size)
      end do
      hash = hash + 1
   end function hash

   !> @brief Find a value in the hash map by its key.
   function map_find(self, key) result(value)
      class(hash_map), intent(in) :: self
      integer, intent(in), allocatable :: key(:)
      integer :: value, index
      type(hash_map_node), pointer :: current

      index = map_hash(key, size(self%buckets))
      value = -1  ! Default value for not found

      current => self%buckets(index)%head
      do while (associated(current))
         if (all(current%key == key)) then
            value = current%value
            exit
         end if
         current => current%next
      end do
   end function map_find

   !> @brief Finalize the hash map, releasing all allocated memory.
   subroutine map_clear(self)
      class(hash_map), intent(inout) :: self
      integer :: i
      type(hash_map_node), pointer :: current, tmp

      do i = 1, size(self%buckets)
         current => self%buckets(i)%head
         do while (associated(current))
            tmp => current
            current => current%next
            deallocate (tmp)
         end do
         self%buckets(i)%head => null()
      end do
   end subroutine map_clear

   !> @brief Resize the hash map.
   subroutine map_resize(self, new_size)
      class(hash_map), intent(inout) :: self
      integer, intent(in) :: new_size
      type(hash_map_bucket), allocatable :: new_buckets(:)
      type(hash_map_node), pointer :: current, new_node, tmp
      integer :: i, new_index, alloc_stat

      ! Safety checks: Do not resize if the new size is not positive.
      if (new_size <= 0) then
         write (*, *) "Error: new_size must be a positive integer."
         return
      end if

      ! Step 1: Allocate a new array of buckets.
      allocate (new_buckets(new_size), stat=i)
      if (i /= 0) then
         write (*, *) "Error: Could not allocate memory for new_buckets."
         return
      end if

      ! Step 2 and 3: Rehash existing items and move them to new buckets.
      do i = 1, size(self%buckets)
         current => self%buckets(i)%head
         do while (associated(current))
            ! Rehash the key to find its place in the new array of buckets.
            new_index = map_hash(current%key, new_size)

            ! Allocate and initialize a new node.
            allocate (new_node)
            allocate (new_node%key(size(current%key)), stat=alloc_stat)  !< Use alloc_stat instead of i
            if (alloc_stat /= 0) then
               write (*, *) "Error: Could not allocate memory for new_node's key."
               return
            end if

            new_node%key = current%key
            new_node%value = current%value

            ! Insert the new node into the appropriate bucket in new_buckets.
            new_node%next => new_buckets(new_index)%head
            new_buckets(new_index)%head => new_node

            ! Move to the next node in the current bucket.
            tmp => current
            current => current%next

            ! Deallocate the current node from the old bucket.
            deallocate (tmp)
         end do
      end do

      ! Step 4: Replace old buckets with new_buckets.
      self%buckets = new_buckets
      self%size = size(self%buckets)
   end subroutine map_resize

   !> @brief Helper function for robust hash calculation.
   integer function map_hash(key, bucket_size)
      integer, intent(in), allocatable :: key(:)
      integer, intent(in) :: bucket_size
      integer :: i

      map_hash = 0
      do i = 1, size(key)
         map_hash = mod(map_hash + key(i), bucket_size)
      end do
      map_hash = map_hash + 1
   end function map_hash

   ! ============================
   ! ARRAY
   ! ============================

   !> @brief Resize an array.
   subroutine array_resize(arr, new_rows, new_cols)
      implicit none
      integer, allocatable, dimension(:, :), intent(inout) :: arr
      integer, intent(in) :: new_rows, new_cols
      integer, allocatable, dimension(:, :) :: temp_arr

      ! Store old data if array is already allocated
      if (allocated(arr)) then
         allocate (temp_arr(size(arr, dim=1), size(arr, dim=2)))
         temp_arr = arr
         deallocate (arr)
      end if

      ! Allocate new array with new dimensions
      allocate (arr(new_rows, new_cols))

      ! Copy old data back into the resized array
      if (allocated(temp_arr)) then
         arr(:, 1:min(size(temp_arr, dim=2), new_cols)) = temp_arr(:, 1:min(size(temp_arr, dim=2), new_cols))
         deallocate (temp_arr)
      end if
   end subroutine array_resize

!> @brief Randomly permute the rows or columns of a 2D array along a specified dimension.
   !>
   !> The subroutine uses the Fisher-Yates shuffle algorithm to randomly permute
   !> the rows or columns of the input 2D array `arr` in-place along the specified dimension.
   !>
   !> @param arr The 2D array to be permuted.
   !> @param dim The dimension along which to perform the permutation (1 for rows, 2 for columns).
   subroutine permute_array2d(arr, dim)
      implicit none
      real, intent(inout), dimension(:, :) :: arr  ! 2D array
      integer, intent(in) :: dim  ! Dimension to shuffle
      integer :: n, i, j
      real, allocatable :: temp_row(:), temp_col(:)

      ! Determine the size of the array
      n = size(arr, dim=dim)

      ! Initialize random seed
      call random_seed()

      if (dim == 1) then
         ! Shuffle along the first dimension (rows)
         allocate (temp_row(size(arr, dim=2)))
         do i = n, 2, -1
            call random_number(temp_row(1))
            j = int(temp_row(1)*i) + 1

            ! Swap rows i and j
            temp_row = arr(i, :)
            arr(i, :) = arr(j, :)
            arr(j, :) = temp_row
         end do
         deallocate (temp_row)

      elseif (dim == 2) then
         ! Shuffle along the second dimension (columns)
         allocate (temp_col(size(arr, dim=1)))
         do i = n, 2, -1
            call random_number(temp_col(1))
            j = int(temp_col(1)*i) + 1

            ! Swap columns i and j
            temp_col = arr(:, i)
            arr(:, i) = arr(:, j)
            arr(:, j) = temp_col
         end do
         deallocate (temp_col)

      else
         print *, "Error: Invalid dimension specified. Please use 1 or 2."
         return
      end if

   end subroutine permute_array2d

end module data_structures
