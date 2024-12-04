! ------------------------------------------------------------------------------
! @file data_structures.F90
! @brief Implementation of simple data structures in Fortran.
! @date 2024-04-27
! ------------------------------------------------------------------------------
module data_structures
    use iso_fortran_env, only: int32, int64, real64
    use utility

    implicit none

    ! --------------------------------------------------------------------------
    ! @brief Defines a hash map with dynamic key-value storage.
    ! --------------------------------------------------------------------------
    type :: HashMap
        integer(int64) :: count = 0                                             ! Number of key-value pairs
        integer(int64) :: capacity = 0                                          ! Size of the hash map (number of buckets)
        type(HashMapBucket), allocatable :: buckets(:)                          ! Array of buckets
    contains
        procedure :: initialise => hash_map_initialise                          ! Initialise the hash map
        procedure :: insert => hash_map_insert                                  ! Insert a key-value pair
        procedure :: delete => hash_map_delete                                  ! Delete a key-value pair
        procedure :: find => hash_map_find                                      ! Retrieve a value by key
        procedure :: get => hash_map_get                                        ! Retrieve a value by key
        procedure :: resize => hash_map_resize                                  ! Resize the hash map
        procedure :: clear => hash_map_clear                                    ! Clear the hash map
    end type HashMap

    ! --------------------------------------------------------------------------
    ! @brief Represents a single node in the hash map.
    ! --------------------------------------------------------------------------
    type :: HashMapNode
        integer(int64), allocatable :: key(:)                                   ! Key array
        integer(int64), allocatable :: value(:)                                 ! Value array
        type(HashMapNode), pointer :: next => null()                            ! Pointer to the next node in the bucket
    end type HashMapNode

    ! --------------------------------------------------------------------------
    ! @brief Represents a bucket in the hash map, containing a list of nodes.
    ! --------------------------------------------------------------------------
    type :: HashMapBucket
        type(HashMapNode), pointer :: head => null()                            ! Pointer to the head node
    end type HashMapBucket

    ! --------------------------------------------------------------------------
    ! @brief Defines a node in the linked list.
    ! --------------------------------------------------------------------------
    type :: ListNode
        integer(int64) :: value                                                 ! Value of the node
        type(ListNode), pointer :: next => null()                               ! Pointer to the next node
    end type ListNode

    ! --------------------------------------------------------------------------
    ! @brief Defines the linked list.
    ! --------------------------------------------------------------------------
    type :: LinkedList
        type(ListNode), pointer :: head => null()                               ! Pointer to the head of the list
        integer(kind=int64) :: count = 0                                        ! Number of elements in the list
    contains
        procedure :: append => linked_list_append                               ! Insert an element at the end
        procedure :: to_array => linked_list_to_array                           ! Retrieve the elements as an integer array
    end type LinkedList

contains

    ! --------------------------------------------------------------------------
    ! @brief Initialise the hash map with a specified size.
    ! @param[inout] self The hash map to initialise.
    ! @param[in] map_size Initial number of buckets.
    ! --------------------------------------------------------------------------
    subroutine hash_map_initialise(self, map_size)
        implicit none
        class(HashMap), intent(inout) :: self                                   ! Hash map instance
        integer(int64), intent(in) :: map_size                                  ! Initial number of buckets
        integer(int32) :: alloc_status                                          ! Allocation status

        if (map_size <= 0) then
            call report("HashMap size must be positive.", 3)
            return
        end if

        allocate(self%buckets(map_size), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for HashMap.", 3)
            return
        end if

        self%capacity = map_size
        self%count = 0
    end subroutine hash_map_initialise

    ! --------------------------------------------------------------------------
    ! @brief Insert a key-value pair into the hash map.
    ! @param[inout] self The hash map instance.
    ! @param[in] key The key to insert.
    ! @param[in] value The value associated with the key.
    ! --------------------------------------------------------------------------
    subroutine hash_map_insert(self, key, value)
        implicit none
        class(HashMap), intent(inout) :: self                                   ! Hash map instance
        integer(int64), intent(in) :: key(:)                                    ! Key array
        integer(int64), intent(in) :: value(:)                                  ! Value array
        integer(int64) :: index                                                 ! Bucket index
        type(HashMapNode), pointer :: new_node                                  ! New node to insert
        real(real64) :: load_factor                                             ! Current load factor
        integer(int32) :: alloc_status                                          ! Allocation status

        ! Compute the hash index
        index = compute_hash(key, self%capacity)

        ! Allocate new node
        allocate(new_node, stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for HashMapNode.", 3)
            return
        end if

        ! Assign key
        allocate(new_node%key(size(key)), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for key.", 3)
            deallocate(new_node)
            return
        end if
        new_node%key = key

        ! Assign value
        allocate(new_node%value(size(value)), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for value.", 3)
            deallocate(new_node%key)
            deallocate(new_node)
            return
        end if
        new_node%value = value

        ! Insert the node at the head of the bucket's linked list
        new_node%next => self%buckets(index)%head
        self%buckets(index)%head => new_node
        self%count = self%count + 1

        ! Check load factor and resize if necessary
        load_factor = real(self%count, real64) / real(self%capacity, real64)
        if (load_factor > 0.7) then
            call self%resize(self%capacity * 2)
        end if
    end subroutine hash_map_insert

    ! --------------------------------------------------------------------------
    ! @brief Retrieve the value associated with a key.
    ! @param[in] self The hash map instance.
    ! @param[in] key The key to search for.
    ! @param[out] value The value associated with the key.
    ! @return found Logical indicating if the key was found.
    ! --------------------------------------------------------------------------
    function hash_map_find(self, key, value) result(found)
        implicit none

        class(HashMap), intent(in) :: self                                       ! Hash map instance
        integer(int64), intent(in) :: key(:)                                     ! Key array
        integer(int64), allocatable, intent(out) :: value(:)                     ! Value array

        logical :: found                                                         ! Key found flag
        integer(int64) :: index                                                  ! Bucket index
        type(HashMapNode), pointer :: current_node                               ! Current node

        found = .false.
        index = compute_hash(key, self%capacity)
        current_node => self%buckets(index)%head

        do while (associated(current_node))
            if (all(current_node%key == key)) then
                found = .true.
                value = current_node%value
                return
            end if
            current_node => current_node%next
        end do
    end function hash_map_find

    ! --------------------------------------------------------------------------
    ! @brief Retrieve the value associated with a key.
    ! @param[in] self The hash map instance.
    ! @param[in] key The key to search for.
    ! @param[out] value The value associated with the key.
    ! --------------------------------------------------------------------------
    function hash_map_get(self, key) result(value)
        implicit none
        class(HashMap), intent(in) :: self                                       ! Hash map instance
        integer(int64), intent(in) :: key(:)                                     ! Key array
        integer(int64), allocatable :: value(:)                                  ! Value array

        if (.not. hash_map_find(self, key, value)) then
            call report("Key not found in the hash map.", 2)
        end if
    end function hash_map_get

    ! --------------------------------------------------------------------------
    ! @brief Delete a key-value pair from the hash map.
    ! @param[inout] self The hash map instance.
    ! @param[in] key The key to delete.
    ! --------------------------------------------------------------------------
    subroutine hash_map_delete(self, key)
        implicit none
        class(HashMap), intent(inout) :: self                                   ! Hash map instance
        integer(int64), intent(in) :: key(:)                                    ! Key array
        integer(int64) :: index                                                 ! Bucket index
        type(HashMapNode), pointer :: current_node, prev_node                   ! Current and previous nodes

        index = compute_hash(key, self%capacity)
        if (index < 0 .or. index >= self%capacity) then
            call report("Invalid hash index during delete operation.", 3)
            return
        end if

        prev_node => null()
        current_node => self%buckets(index)%head

        ! Traverse the linked list to find the key
        do while (associated(current_node))
            if (all(current_node%key == key)) then
                ! Key found: remove the node
                if (associated(prev_node)) then
                    prev_node%next => current_node%next
                else
                    self%buckets(index)%head => current_node%next
                end if
                ! Deallocate key, value, and node
                deallocate(current_node%key)
                deallocate(current_node%value)
                deallocate(current_node)
                self%count = self%count - 1
                return
            end if
            prev_node => current_node
            current_node => current_node%next
        end do

        call report("Key not found in the hash map.", 2)
    end subroutine hash_map_delete

    ! --------------------------------------------------------------------------
    ! @brief Clear the hash map by deallocating all nodes.
    ! @param[inout] self The hash map instance.
    ! --------------------------------------------------------------------------
    subroutine hash_map_clear(self)
        implicit none
        class(HashMap), intent(inout) :: self                                   ! Hash map instance
        integer(int64) :: i                                                     ! Loop index
        type(HashMapNode), pointer :: current_node, next_node                   ! Nodes to deallocate

        do i = 1, self%capacity
            current_node => self%buckets(i)%head
            do while (associated(current_node))
                next_node => current_node%next
                deallocate(current_node%key)
                deallocate(current_node%value)
                deallocate(current_node)
                current_node => next_node
            end do
        end do
        self%count = 0
    end subroutine hash_map_clear

    ! --------------------------------------------------------------------------
    ! @brief Resize the hash map.
    ! @param[inout] self The hash map instance.
    ! @param[in] new_size The new size for the hash map.
    ! --------------------------------------------------------------------------
    subroutine hash_map_resize(self, new_size)
        implicit none
        class(HashMap), intent(inout) :: self              ! Hash map instance
        integer(int64), intent(in) :: new_size             ! New size
        type(HashMapBucket), allocatable :: new_buckets(:) ! New buckets
        type(HashMapNode), pointer :: current, next        ! For traversal
        integer(int64) :: i, new_index                     ! Loop variables
        integer(int32) :: alloc_status                     ! Allocation status

        ! Validate size
        if (new_size <= 0) then
            call report("Resize size must be positive.", 3)
            return
        end if

        ! Allocate new buckets
        allocate(new_buckets(new_size), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate new buckets.", 3)
            return
        end if

        ! Rehash all existing elements
        do i = 1, self%capacity
            current => self%buckets(i)%head
            do while (associated(current))
                ! Save next pointer before moving node
                next => current%next
                
                ! Calculate new index
                new_index = compute_hash(current%key, new_size)
                
                ! Move node to new bucket
                current%next => new_buckets(new_index)%head
                new_buckets(new_index)%head => current
                
                ! Move to next node
                current => next
            end do
        end do

        ! Replace old buckets with new ones
        deallocate(self%buckets)
        self%buckets = new_buckets
        self%capacity = new_size
    end subroutine hash_map_resize

    ! --------------------------------------------------------------------------
    ! @brief Compute the hash for a given key.
    ! @param[in] key The key array.
    ! @param[in] bucket_size The size of the bucket array.
    ! @return index The hash index.
    ! --------------------------------------------------------------------------
    integer(int64) function compute_hash(key, bucket_size)
        implicit none
        integer(int64), intent(in) :: key(:), bucket_size
        integer(int64) :: i, hash
        
        ! Use smaller prime numbers to avoid overflow
        hash = 17_int64  ! Initial seed
        do i = 1, size(key)
            ! Combine using prime multiplication and XOR
            hash = ieor(hash * 31_int64, key(i))
        end do
        
        ! Ensure positive value and proper modulo
        compute_hash = modulo(abs(hash), bucket_size)
    end function compute_hash

    ! --------------------------------------------------------------------------
    ! @brief Insert an element at the end of the linked list.
    ! @param[inout] self The linked list.
    ! @param[in] value The value to insert.
    ! --------------------------------------------------------------------------
    subroutine linked_list_append(self, value)
        implicit none

        class(LinkedList), intent(inout) :: self
        integer(int64), intent(in) :: value
        
        type(ListNode), pointer :: new_node
        type(ListNode), pointer :: current
        integer :: status
        character(len=4096) :: message

        ! Allocate a new node for the value
        allocate (new_node, stat=status)
        if (status .ne. 0) then
            write(message, '(A)') "Failed to allocate memory for ListNode."
            call report(trim(message), 3)
            return
        end if

        new_node%value = value
        new_node%next => null()

        ! If the list is empty, set the new node as the head
        if (.not. associated(self%head)) then
            self%head => new_node
        else
            ! Traverse the list to find the last node
            current => self%head
            do while (associated(current%next))
                current => current%next
            end do
            current%next => new_node
        end if

        self%count = self%count + 1
    end subroutine linked_list_append

    ! --------------------------------------------------------------------------
    ! @brief Retrieve the elements of the linked list as an integer array.
    ! @param[in] self The linked list.
    ! @return array The array of elements.
    ! --------------------------------------------------------------------------
    function linked_list_to_array(self) result(array)
        implicit none

        class(LinkedList), intent(in) :: self                                   ! Linked list instance
        integer(int64), allocatable :: array(:)                                 ! Array of elements

        type(ListNode), pointer :: current                                      ! Current node
        integer(int64) :: i                                                     ! Loop index

        ! Check if the list is empty
        if (self%count .eq. 0) then
            allocate(array(0))  ! Return an empty array
            return
        end if

        ! Allocate the array with the size of the linked list
        allocate(array(self%count))

        ! Traverse the linked list and copy values into the array
        current => self%head
        do i = 1, self%count
            array(i) = current%value
            current => current%next
        end do
    end function linked_list_to_array

end module data_structures
