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
    ! @brief Defines a hash map with key-value storage.
    ! --------------------------------------------------------------------------
    type :: HashMap
        integer(int64) :: count = 0                                             ! Number of key-value pairs
        integer(int64) :: capacity = 0                                          ! Capacity of the hash map (number of buckets)
        type(HashMapBucket), allocatable :: buckets(:)                          ! Array of buckets
    contains
        procedure :: initialise => hash_map_initialise                          ! Initialise the hash map
        procedure :: keys => hash_map_keys                                      ! Retrieve all keys as a 2D matrix
        procedure :: insert => hash_map_insert                                  ! Insert a key-value pair
        procedure :: delete => hash_map_delete                                  ! Delete a key-value pair
        procedure :: find => hash_map_find                                      ! Retrieve a value by key
        procedure :: get => hash_map_get                                        ! Retrieve a value by key
        procedure :: resize => hash_map_resize                                  ! Resize the hash map
        procedure :: reset => hash_map_reset                                    ! Clear the hash map
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
    ! @brief Defines the linked list.
    ! --------------------------------------------------------------------------
    type :: LinkedList
        type(ListNode), pointer :: head => null()                               ! Pointer to the head of the list
        type(ListNode), pointer :: tail => null()                               ! Pointer to the tail of the list
        integer(kind=int64) :: count = 0                                        ! Number of elements in the list
        logical :: is_cyclic = .false.                                          ! True if the list is cyclic
    contains
        procedure :: append => linked_list_append                               ! Insert an element at the end
        procedure :: cyclic => linked_list_cyclic                               ! Make the list cyclic
        procedure :: to_array => linked_list_to_array                           ! Retrieve the elements as an integer array
    end type LinkedList

    ! --------------------------------------------------------------------------
    ! @brief Defines a node in the linked list.
    ! --------------------------------------------------------------------------
    type :: ListNode
        integer(int64) :: value                                                 ! Value of the node
        type(ListNode), pointer :: next => null()                               ! Pointer to the next node
    end type ListNode

    ! --------------------------------------------------------------------------
    ! @brief Defines a priority queue with configurable priority ordering.
    ! --------------------------------------------------------------------------
    type :: PriorityQueue
        type(PriorityNode), allocatable :: heap(:)    ! Array of PriorityNodes
        integer(int64) :: count = 0                   ! Current number of elements
        integer(int64) :: capacity = 0                ! Maximum capacity of the heap
        logical :: reverse_priority = .false.         ! True for low-priority-first
    contains
        procedure :: initialise => priority_queue_initialise
        procedure :: insert => priority_queue_insert
        procedure :: extract => priority_queue_extract
        procedure, private :: compare => priority_queue_compare
        procedure, private :: swap_nodes => priority_queue_swap_nodes
        procedure, private :: cleanup_heap => priority_queue_cleanup_heap
    end type PriorityQueue

    ! --------------------------------------------------------------------------
    ! @brief Node structure for the priority queue.
    ! --------------------------------------------------------------------------
    type :: PriorityNode
        integer(int64), allocatable :: value(:)    ! Allocatable array for the value
        real(real64) :: priority                   ! Priority of the node
    end type PriorityNode

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

        allocate (self%buckets(map_size), stat=alloc_status)
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
        type(HashMapNode), pointer :: current_node, new_node                    ! Nodes for traversal and insertion
        integer(int32) :: alloc_status                                          ! Allocation status
    
        if (.not. allocated(self%buckets) .or. self%capacity <= 0) then
            call report("HashMap is not initialised or has invalid capacity.", 3)
            return
        end if
    
        ! Compute the hash index
        index = compute_hash(key, self%capacity)
        if (index < 1 .or. index > self%capacity) then
            call report("Invalid hash index in insert operation.", 3)
            return
        end if
    
        ! Traverse the bucket's linked list to find if the key already exists
        current_node => self%buckets(index)%head
        do while (associated(current_node))
            if (all(current_node%key == key)) then
                ! Key already exists; update its value
                deallocate(current_node%value)
                allocate(current_node%value(size(value)), stat=alloc_status)
                if (alloc_status /= 0) then
                    call report("Failed to allocate memory for value.", 3)
                    return
                end if
                current_node%value = value
                return
            end if
            current_node => current_node%next
        end do
    
        ! Key not found; add a new node to the bucket
        allocate(new_node, stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for HashMapNode.", 3)
            return
        end if
    
        allocate(new_node%key(size(key)), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for key.", 3)
            deallocate(new_node)
            return
        end if
        new_node%key = key
    
        allocate(new_node%value(size(value)), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for value.", 3)
            deallocate(new_node%key)
            deallocate(new_node)
            return
        end if
        new_node%value = value
    
        ! Insert the new node at the head of the bucket
        new_node%next => self%buckets(index)%head
        self%buckets(index)%head => new_node
    
        ! Increment the count of unique key-value pairs
        self%count = self%count + 1
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
                deallocate (current_node%key)
                deallocate (current_node%value)
                deallocate (current_node)
                self%count = self%count - 1
                return
            end if
            prev_node => current_node
            current_node => current_node%next
        end do

        call report("Key not found in the hash map.", 2)
    end subroutine hash_map_delete

    ! --------------------------------------------------------------------------
    ! @brief Retrieve all keys from the hash map as a 2D matrix.
    ! @param[in] self The hash map instance.
    ! @return keys_matrix A 2D array where each row is a key.
    ! --------------------------------------------------------------------------
    function hash_map_keys(self) result(keys_matrix)
        class(HashMap), intent(in) :: self                   ! Hash map instance
        integer(int64), allocatable :: keys_matrix(:,:)      ! 2D matrix of keys

        type(HashMapNode), pointer :: current_node           ! Pointer to traverse nodes
        integer(int64), allocatable :: key_sizes(:)          ! Array to store sizes of keys
        integer(int64) :: i, bucket_index, max_key_length    ! Loop and indexing variables
        integer(int64) :: num_keys, current_key_index        ! Number of keys and current index

        ! Initialise variables
        num_keys = self%count
        allocate(key_sizes(num_keys))
        current_key_index = 0
        max_key_length = 0

        ! Determine the maximum key length and collect key sizes
        do bucket_index = 1, self%capacity
            current_node => self%buckets(bucket_index)%head
            do while (associated(current_node))
                current_key_index = current_key_index + 1
                key_sizes(current_key_index) = size(current_node%key)
                max_key_length = max(max_key_length, size(current_node%key))
                current_node => current_node%next
            end do
        end do

        ! Allocate the keys matrix and initialise to zero
        allocate(keys_matrix(max_key_length, num_keys))
        keys_matrix = 0_int64

        ! Populate the keys matrix
        current_key_index = 0
        do bucket_index = 1, self%capacity
            current_node => self%buckets(bucket_index)%head
            do while (associated(current_node))
                current_key_index = current_key_index + 1
                keys_matrix(1:key_sizes(current_key_index), current_key_index) = current_node%key
                current_node => current_node%next
            end do
        end do
    end function hash_map_keys


    ! --------------------------------------------------------------------------
    ! @brief Clear the hash map by deallocating all nodes.
    ! @param[inout] self The hash map instance.
    ! --------------------------------------------------------------------------
    subroutine hash_map_reset(self)
        implicit none
        class(HashMap), intent(inout) :: self                                   ! Hash map instance
        integer(int64) :: i                                                     ! Loop index
        type(HashMapNode), pointer :: current_node, next_node                   ! Nodes to deallocate

        do i = 1, self%capacity
            current_node => self%buckets(i)%head
            do while (associated(current_node))
                next_node => current_node%next
                deallocate (current_node%key)
                deallocate (current_node%value)
                deallocate (current_node)
                current_node => next_node
            end do
        end do
        self%count = 0
    end subroutine hash_map_reset

    ! --------------------------------------------------------------------------
    ! @brief Resize the hash map.
    ! @param[inout] self The hash map instance.
    ! @param[in] new_size The new size for the hash map.
    ! --------------------------------------------------------------------------
    subroutine hash_map_resize(self, new_size)
        implicit none

        class(HashMap), intent(inout) :: self
        integer(int64), intent(in) :: new_size
        type(HashMapBucket), allocatable :: new_buckets(:)
        type(HashMapNode), pointer :: current, next
        integer(int64) :: i, new_index
        integer(int32) :: alloc_status

        if (new_size <= 0) then
            call report("Resize size must be positive.", 3)
            return
        end if

        allocate (new_buckets(new_size), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate new buckets.", 3)
            return
        end if

        do i = 1, self%capacity
            current => self%buckets(i)%head
            do while (associated(current))
                next => current%next
                new_index = compute_hash(current%key, new_size)
                if (new_index < 1 .or. new_index > new_size) then
                    call report("Invalid rehash index during resize.", 3)
                    return
                end if
                current%next => new_buckets(new_index)%head
                new_buckets(new_index)%head => current
                current => next
            end do
        end do

        deallocate (self%buckets)
        self%buckets = new_buckets
        self%capacity = new_size
    end subroutine hash_map_resize

    ! --------------------------------------------------------------------------
    ! @brief Compute the hash for a given key.
    ! @param[in] key The key array.
    ! @param[in] bucket_size The size of the bucket array.
    ! @return index The hash index.
    ! --------------------------------------------------------------------------
    function compute_hash(key, bucket_size) result(index)
        implicit none

        integer(int64), intent(in) :: key(:), bucket_size
        integer(int64) :: i, hash, index

        if (bucket_size <= 0) then
            call report("Invalid bucket size in compute_hash.", 3)
            index = 1
            return
        end if

        hash = 17_int64
        do i = 1, size(key)
            hash = ieor(hash*31_int64, key(i))
        end do

        index = modulo(abs(hash), bucket_size) + 1
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
        integer :: status

        ! Allocate new node
        allocate(new_node, stat=status)
        if (status /= 0) then
            call report("Failed to allocate memory for ListNode.", 3)
            return
        end if

        new_node%value = value
        new_node%next => null()

        ! Empty list case
        if (.not. associated(self%head)) then
            self%head => new_node
            self%tail => new_node
        else
            ! Append to tail
            self%tail%next => new_node
            self%tail => new_node
        end if

        self%count = self%count + 1
    end subroutine linked_list_append

    ! --------------------------------------------------------------------------
    ! @brief Makes the linked list cyclic by connecting last node to head.
    ! @param[inout] self The linked list instance.
    ! --------------------------------------------------------------------------
    subroutine linked_list_cyclic(self)
        implicit none
        class(LinkedList), intent(inout) :: self
        type(ListNode), pointer :: current

        ! Check for empty list or single element
        if (.not. associated(self%head) .or. self%count <= 1) then
            return
        end if

        ! Already cyclic
        if (self%is_cyclic) then
            return
        end if

        ! Find last node
        current => self%head
        do while (associated(current%next))
            current => current%next
        end do

        ! Connect last node to head
        current%next => self%head
        self%is_cyclic = .true.

    end subroutine linked_list_cyclic

    ! --------------------------------------------------------------------------
    ! @brief Retrieve the elements of the linked list as an integer array.
    ! @param[in] self The linked list.
    ! @return array The array of elements.
    ! --------------------------------------------------------------------------
    function linked_list_to_array(self) result(array)
        implicit none

        class(LinkedList), intent(in) :: self                                   ! Linked list instance

        integer(int64), allocatable :: array(:)                                 ! Output array

        type(ListNode), pointer :: current                                      ! Current node pointer
        integer(int64) :: i                                                     ! Loop counter

        ! Handle empty list case
        if (self%count == 0) then
            allocate(array(0))
            return
        end if

        ! Allocate array and copy values
        allocate(array(self%count))
        current => self%head
        do i = 1, self%count
            if (.not. associated(current)) then
                call report("Corrupted linked list structure", 3)
                return
            end if
            array(i) = current%value
            current => current%next
        end do
    end function linked_list_to_array

    ! --------------------------------------------------------------------------
    ! @brief Initialise a priority queue with given capacity.
    ! @param[inout] self The priority queue instance.
    ! @param[in] capacity Maximum number of elements.
    ! @param[in] reverse_priority Optional - if true, creates min-heap instead of max-heap.
    ! --------------------------------------------------------------------------
    subroutine priority_queue_initialise(self, capacity, reverse_priority)
        implicit none

        class(PriorityQueue), intent(inout) :: self                             ! Priority queue instance
        integer(int64), intent(in) :: capacity                                  ! Maximum capacity
        logical, intent(in), optional :: reverse_priority                       ! Priority direction

        integer(int32) :: alloc_status                                          ! Allocation status

        ! Validate input
        if (capacity <= 0) then
            call report("Priority queue capacity must be positive", 3)
            return
        end if

        ! Clean up existing heap if any
        if (allocated(self%heap)) then
            call priority_queue_cleanup_heap(self)
            deallocate(self%heap)
        end if

        ! Allocate new heap
        allocate(self%heap(capacity), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate priority queue memory", 3)
            return
        end if

        ! Initialise queue state
        self%capacity = capacity
        self%count = 0
        self%reverse_priority = .false.
        if (present(reverse_priority)) then
            self%reverse_priority = reverse_priority
        end if
    end subroutine priority_queue_initialise

    ! --------------------------------------------------------------------------
    ! @brief Insert a value with associated priority into the queue.
    ! @param[inout] self The priority queue instance.
    ! @param[in] value The value to insert.
    ! @param[in] priority The priority associated with the value.
    ! --------------------------------------------------------------------------
    subroutine priority_queue_insert(self, value, priority)
        implicit none

        class(PriorityQueue), intent(inout) :: self                             ! Priority queue instance
        integer(int64), intent(in) :: value(:)                                  ! Value to insert
        real(real64), intent(in) :: priority                                    ! Priority value

        integer(int64) :: idx, parent                                           ! Node indices
        integer(int32) :: alloc_status                                          ! Allocation status

        ! Check capacity
        if (self%count >= self%capacity) then
            call report("Priority queue is full", 3)
            return
        end if

        ! Insert at bottom
        self%count = self%count + 1
        idx = self%count

        ! Allocate and copy value
        if (allocated(self%heap(idx)%value)) deallocate(self%heap(idx)%value)
        allocate(self%heap(idx)%value(size(value)), stat=alloc_status)
        if (alloc_status /= 0) then
            call report("Failed to allocate memory for new value", 3)
            self%count = self%count - 1
            return
        end if

        self%heap(idx)%value = value
        self%heap(idx)%priority = priority

        ! Bubble up to maintain heap property
        do while (idx > 1)
            parent = idx / 2
            if (self%compare(self%heap(idx)%priority, self%heap(parent)%priority)) then
                call self%swap_nodes(idx, parent)
                idx = parent
            else
                exit
            end if
        end do
    end subroutine priority_queue_insert

    ! Helper subroutine to clean up heap memory
    subroutine priority_queue_cleanup_heap(self)
        implicit none
        class(PriorityQueue), intent(inout) :: self
        integer(int64) :: i

        do i = 1, self%capacity
            if (allocated(self%heap(i)%value)) deallocate(self%heap(i)%value)
        end do
    end subroutine priority_queue_cleanup_heap

    ! --------------------------------------------------------------------------
    ! @brief Extract the value with the highest priority from the queue.
    ! @param[inout] self The priority queue instance.
    ! @param[out] value The extracted value.
    ! @param[out] priority The priority of the extracted value.
    ! @return success True if extraction was successful.
    ! --------------------------------------------------------------------------
    function priority_queue_extract(self, value, priority) result(success)
        implicit none

        class(PriorityQueue), intent(inout) :: self                             ! Priority queue instance
        integer(int64), allocatable, intent(out) :: value(:)                    ! Extracted value
        real(real64), intent(out) :: priority                                   ! Priority of extracted value

        logical :: success                                                      ! Extraction success flag

        integer(int64) :: idx, child, left, right                               ! Node indices
        
        ! Initialise result
        success = .false.
        
        ! Check empty queue
        if (self%count <= 0) then
            call report("PriorityQueue is empty", 2)
            return
        end if
        
        ! Check array bounds
        if (.not. allocated(self%heap) .or. self%count > self%capacity) then
            call report("PriorityQueue heap corruption detected", 3)
            return
        end if
        
        ! Save root value
        if (.not. allocated(self%heap(1)%value)) then
            call report("Root node value not allocated", 3)
            return
        end if
        
        allocate(value(size(self%heap(1)%value)))
        value = self%heap(1)%value
        priority = self%heap(1)%priority
        
        ! Handle last element special case
        if (self%count == 1) then
            deallocate(self%heap(1)%value)
            self%count = 0
            success = .true.
            return
        end if
        
        ! Move last element to root
        if (allocated(self%heap(1)%value)) deallocate(self%heap(1)%value)
        allocate(self%heap(1)%value(size(self%heap(self%count)%value)))
        self%heap(1)%value = self%heap(self%count)%value
        self%heap(1)%priority = self%heap(self%count)%priority
        
        ! Clean up last element
        deallocate(self%heap(self%count)%value)
        self%count = self%count - 1
        
        ! Heapify down
        idx = 1
        do while (idx <= self%count/2)  ! Only need to check until parent nodes
            child = idx
            left = 2 * idx
            right = 2 * idx + 1
            
            ! Find smallest/largest child based on priority
            if (left <= self%count .and. &
                self%compare(self%heap(left)%priority, self%heap(child)%priority)) then
                child = left
            end if
            
            if (right <= self%count .and. &
                self%compare(self%heap(right)%priority, self%heap(child)%priority)) then
                child = right
            end if
            
            ! If no swap needed, heap property restored
            if (child == idx) exit
            
            ! Swap with appropriate child
            call self%swap_nodes(idx, child)
            idx = child
        end do
        
        success = .true.
    end function priority_queue_extract

    ! --------------------------------------------------------------------------
    ! @brief Compare two priorities based on the queue's ordering.
    ! @param[in] self The priority queue instance.
    ! @param[in] a The first priority.
    ! @param[in] b The second priority.
    ! @return True if a has higher priority than b.
    ! --------------------------------------------------------------------------
    logical function priority_queue_compare(self, a, b)
        class(PriorityQueue), intent(in) :: self
        real(real64), intent(in) :: a, b

        if (self%reverse_priority) then
            ! Reverse priority: lower priority comes first
            priority_queue_compare = a > b
        else
            ! Default: higher priority comes first
            priority_queue_compare = a < b 
        end if
    end function priority_queue_compare

    ! --------------------------------------------------------------------------
    ! @brief Swap two nodes in the priority queue.
    ! @param[inout] self The priority queue instance.
    ! @param[in] i The first node index.
    ! @param[in] j The second node index.
    ! --------------------------------------------------------------------------
    subroutine priority_queue_swap_nodes(self, i, j)
        implicit none

        class(PriorityQueue), intent(inout) :: self
        integer(int64), intent(in) :: i, j
        
        type(PriorityNode) :: temp
        
        if (allocated(temp%value)) then
            deallocate(temp%value)
        end if

        allocate(temp%value(size(self%heap(i)%value)))
        
        temp%value = self%heap(i)%value
        temp%priority = self%heap(i)%priority
        
        if (allocated(self%heap(i)%value)) then
            deallocate(self%heap(i)%value)
        end if

        allocate(self%heap(i)%value(size(self%heap(j)%value)))
        self%heap(i)%value = self%heap(j)%value
        self%heap(i)%priority = self%heap(j)%priority
        
        if (allocated(self%heap(j)%value)) then
            deallocate(self%heap(j)%value)
        end if

        allocate(self%heap(j)%value(size(temp%value)))
        self%heap(j)%value = temp%value
        self%heap(j)%priority = temp%priority
    end subroutine priority_queue_swap_nodes
end module data_structures
