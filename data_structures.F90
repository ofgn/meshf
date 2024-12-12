! ------------------------------------------------------------------------------
! @file data_structures.F90
! @brief Implementation of simple data structures in Fortran.
! @date 2024-04-27
! ------------------------------------------------------------------------------
module data_structures
    use iso_fortran_env, only: int32, int64, real64
    use utility

    implicit none

    integer, parameter :: HASH_TABLE_SIZE = 100003  ! A large prime number

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

    type :: Set
        integer(int64), allocatable :: table(:)     ! Hash table
        logical, allocatable :: occupied(:)         ! Occupancy flags
        integer(int64) :: size = 0                  ! Number of elements
    contains
        procedure :: initialise => hashset_initialize
        procedure :: insert => hashset_insert
        procedure :: remove => hashset_remove
        procedure :: exists => hashset_exists
    end type Set

contains

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

    subroutine hashset_initialize(self)
        type(Set), intent(inout) :: self  ! Declare self as type(Set)
        allocate(self%table(HASH_TABLE_SIZE))
        allocate(self%occupied(HASH_TABLE_SIZE))
        self%occupied = .false.
        self%size = 0
    end subroutine hashset_initialize

    integer function hash_function(key) result(index)
        integer(int64), intent(in) :: key
        index = mod(key, HASH_TABLE_SIZE) + 1
    end function hash_function

    subroutine hashset_insert(self, key)
        type(Set), intent(inout) :: self  ! Declare self as type(Set)
        integer(int64), intent(in) :: key
        integer :: index, i

        index = hash_function(key)
        do i = 0, HASH_TABLE_SIZE - 1
            if (.not. self%occupied(index)) then
                self%table(index) = key
                self%occupied(index) = .true.
                self%size = self%size + 1
                return
            else if (self%table(index) == key) then
                ! Key already present
                return
            else
                index = mod(index, HASH_TABLE_SIZE) + 1
            end if
        end do
        print *, 'Set is full, cannot insert key:', key
    end subroutine hashset_insert

    subroutine hashset_remove(self, key)
        type(Set), intent(inout) :: self  ! Declare self as type(Set)
        integer(int64), intent(in) :: key
        integer :: index, i

        index = hash_function(key)
        do i = 0, HASH_TABLE_SIZE - 1
            if (.not. self%occupied(index)) then
                ! Key not found
                return
            else if (self%table(index) == key) then
                self%occupied(index) = .false.
                self%size = self%size - 1
                return
            else
                index = mod(index, HASH_TABLE_SIZE) + 1
            end if
        end do
    end subroutine hashset_remove

    logical function hashset_exists(self, key)
        type(Set), intent(in) :: self  ! Declare self as type(Set)
        integer(int64), intent(in) :: key
        integer :: index, i

        index = hash_function(key)
        do i = 0, HASH_TABLE_SIZE - 1
            if (.not. self%occupied(index)) then
                hashset_exists = .false.
                return
            else if (self%table(index) == key) then
                hashset_exists = .true.
                return
            else
                index = mod(index, HASH_TABLE_SIZE) + 1
            end if
        end do
        hashset_exists = .false.
    end function hashset_exists

end module data_structures
