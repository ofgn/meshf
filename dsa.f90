! ------------------------------------------------------------------------------
! @file dsa.F90
! @brief Fortran implementation of common data structures and algorithms.
! @details Optimized Set implementation for arbitrary key sizes.
! @date 2024-07-27
! ------------------------------------------------------------------------------
module data_structures_algorithms
    use global

    implicit none

    public :: Set, MinHeap
    public :: quicksort, sort4

    ! --------------------------------------------------------------------------
    ! @brief A min-heap data structure for storing integer values with priority.
    ! @details The heap is implemented as an array of nodes with integer values
    ! and real priorities. The heap is a complete binary tree where the parent
    ! node has a lower priority than its children.
    ! --------------------------------------------------------------------------
    type :: MinHeap
        type(MinHeapNode), allocatable :: values(:)                             ! The array of heap nodes
        integer(int64) :: count = 0                                             ! The number of elements
        integer(int64) :: capacity = 0                                          ! The capacity of the heap
    contains
        procedure :: initialise => min_heap_initialise
        procedure :: add => min_heap_add
        procedure :: min => min_heap_min
        procedure, private :: heapify => min_heap_heapify
        procedure :: finalise => min_heap_finalise
        final :: min_heap_finalise_auto
    end type MinHeap

    ! --------------------------------------------------------------------------
    ! @brief A min-heap node with an integer value and real priority.
    ! --------------------------------------------------------------------------
    type :: MinHeapNode
        integer(int64), allocatable :: value(:)                                 ! The integer value  
        real(real64) :: priority                                                ! The priority number
    end type MinHeapNode

    ! --------------------------------------------------------------------------
    ! @brief Interface for MinHeap constructor
    ! --------------------------------------------------------------------------
    interface MinHeap
        module procedure min_heap_constructor                                   
    end interface

    ! --------------------------------------------------------------------------
    ! @brief A set data structure for storing unique integer arrays.
    ! @details The set is implemented using a hash table with quadratic probing.       
    ! --------------------------------------------------------------------------
    type :: Set
        integer(int64) :: count = 0                                             ! The number of elements.
        integer(int64) :: capacity = 0                                          ! The capacity of the set.
        logical(bool8), allocatable :: occupied(:)                              ! The slot occupancy.
        type(SetNode), allocatable :: keys(:)                                   ! The array of keys.
    contains
        procedure :: initialise => set_initialise
        procedure :: add => set_add
        procedure :: contains => set_contains
        procedure :: remove => set_remove
        procedure :: reset => set_reset
        procedure :: size => set_size
        procedure :: to_array => set_to_array
        procedure :: finalise => set_finalise
        procedure, private :: resize => set_resize
        procedure, private :: hash => murmur3_hash
        procedure, private :: quadratic_probe => set_quadratic_probe
        final :: set_finalise_auto
    end type Set

    ! --------------------------------------------------------------------------
    ! @brief Interface for Set constructor
    ! --------------------------------------------------------------------------
    interface Set
        module procedure set_constructor
    end interface

    ! --------------------------------------------------------------------------
    ! @brief A derived type to store individual keys of arbitrary size.
    ! --------------------------------------------------------------------------
    type :: SetNode
        integer(int64), allocatable :: key(:)                                   ! The key array
    end type SetNode

contains

    ! --------------------------------------------------------------------------
    ! @brief Constructor for the MinHeap type.
    ! @return The new MinHeap instance.
    ! --------------------------------------------------------------------------
    function min_heap_constructor()
        implicit none

        type(MinHeap) :: min_heap_constructor
    end function

    ! --------------------------------------------------------------------------
    ! @brief Initialise a new min-heap with an optional initial capacity.
    ! @param[in] initial_capacity The initial capacity of the heap.
    ! --------------------------------------------------------------------------
    subroutine min_heap_initialise(self, initial_capacity)
        implicit none

        class(MinHeap) :: self                                                  ! The min-heap instance
        integer(int64), optional, intent(in) :: initial_capacity                ! The initial capacity of the heap

        if (present(initial_capacity)) then
            self%capacity = initial_capacity
        else 
            self%capacity = 1024
        endif

        allocate(self%values(self%capacity))
        self%count = 0
    end subroutine min_heap_initialise

    ! --------------------------------------------------------------------------
    ! @brief Add a new element to the min-heap with safety checks
    ! @param[in] key The integer array to add
    ! @param[in] value The priority value
    ! --------------------------------------------------------------------------
    subroutine min_heap_add(self, key, value)
        implicit none
        class(MinHeap), intent(inout) :: self
        integer(int64), intent(in) :: key(:)
        real(real64), intent(in) :: value

        integer(int64) :: i, parent
        
        ! Resize if needed
        if (self%count >= self%capacity) then
            block
                type(MinHeapNode), allocatable :: temp(:)
                allocate(temp(self%capacity * 2))
                
                ! Copy existing values
                do i = 1, self%count
                    if (allocated(self%values(i)%value)) then
                        allocate(temp(i)%value(size(self%values(i)%value)))
                        temp(i)%value = self%values(i)%value
                        temp(i)%priority = self%values(i)%priority
                        deallocate(self%values(i)%value)
                    end if
                end do
                
                call move_alloc(temp, self%values)
                self%capacity = self%capacity * 2
            end block
        endif
        
        ! Add new element
        self%count = self%count + 1
        if (.not. allocated(self%values(self%count)%value)) then
            allocate(self%values(self%count)%value(size(key)))
        else if (size(self%values(self%count)%value) /= size(key)) then
            deallocate(self%values(self%count)%value)
            allocate(self%values(self%count)%value(size(key)))
        end if
        
        self%values(self%count)%value = key
        self%values(self%count)%priority = value
        
        ! Bubble up with thread safety
        !$omp critical
        i = self%count
        do while (i > 1)
            parent = i / 2
            if (self%values(i)%priority < self%values(parent)%priority) then
                call swap_elements(self%values(i), self%values(parent))
                i = parent
            else
                exit
            end if
        end do
        !$omp end critical
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Remove and return the minimum element from the min-heap.
    ! @param[out] priority Optional priority value of the minimum element
    ! @return The value array of the minimum element
    ! --------------------------------------------------------------------------
    function min_heap_min(self, priority) result(min_val)
        implicit none
        class(MinHeap), intent(inout) :: self                                   ! The min-heap instance
        real(real64), optional, intent(out) :: priority                         ! The priority of the min element

        integer(int64), allocatable :: min_val(:)                               ! The min element value
        
        integer(int64) :: i                                                     ! The loop index
        integer(int64) :: smallest
        integer(int64) :: left
        integer(int64) :: right
        
        if (self%count < 1) return
        
        ! Save min element value
        allocate(min_val(size(self%values(1)%value)))
        min_val = self%values(1)%value
        
        ! Set priority if requested
        if (present(priority)) then
            priority = self%values(1)%priority
        endif

        ! Before removing the min element
        if (allocated(self%values(1)%value)) then
            deallocate(self%values(1)%value)
        end if
        
        ! Move last element to root
        if (self%count > 1) then
            self%values(1) = self%values(self%count)
        endif
        self%count = self%count - 1
        
        ! Heapify from root
        i = 1
        do while (.true.)
            smallest = i
            left = 2 * i
            right = 2 * i + 1
            
            if (left <= self%count) then
                if (self%values(left)%priority &
                    < self%values(smallest)%priority) then

                    smallest = left
                endif
            endif
            
            if (right <= self%count) then
                if (self%values(right)%priority &
                    < self%values(smallest)%priority) then
                        
                    smallest = right
                endif
            endif
            
            if (smallest /= i) then
                call swap_elements(self%values(i), self%values(smallest))
                i = smallest
            else
                exit
            endif
        end do
    end function

    ! --------------------------------------------------------------------------
    ! @brief Heapify the min-heap.
    ! --------------------------------------------------------------------------
    subroutine min_heap_heapify(self)
        implicit none

        class(MinHeap), intent(inout) :: self

        integer(int64) :: i
        
        !$omp parallel do if(self%count > PARALLEL_THRESHOLD) schedule(dynamic)
        do i = self%count/2, 1, -1
            call sift_down(self, i)
        end do
        !$omp end parallel do
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Sift down the min-heap from a given index.
    ! @param[inout] heap The min-heap.
    ! @param[in] start_idx The starting index.
    ! --------------------------------------------------------------------------
    subroutine sift_down(heap, start_idx)
        implicit none
        class(MinHeap), intent(inout) :: heap
        integer(int64), intent(in) :: start_idx

        integer(int64) :: current, smallest, left, right, n
        logical :: continue_sift

        n = heap%count
        current = start_idx
        continue_sift = .true.

        do while (continue_sift)
            smallest = current
            left = 2 * current
            right = left + 1

            ! Check left child
            if (left <= n) then
                if (heap%values(left)%priority &
                    < heap%values(smallest)%priority) then

                    smallest = left
                end if
            end if

            ! Check right child
            if (right <= n) then
                if (heap%values(right)%priority &
                    < heap%values(smallest)%priority) then

                    smallest = right
                end if
            end if

            ! Swap if needed and continue
            if (smallest /= current) then
                call swap_elements(heap%values(current), heap%values(smallest))
                current = smallest
            else
                continue_sift = .false.
            end if
        end do
    end subroutine sift_down

    ! --------------------------------------------------------------------------
    ! @brief Swap two MinHeapNode elements.
    ! @param[inout] a First element.
    ! @param[inout] b Second element.
    ! --------------------------------------------------------------------------
    subroutine swap_elements(a, b)
        implicit none
        type(MinHeapNode), intent(inout) :: a, b
        type(MinHeapNode) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_elements

    ! --------------------------------------------------------------------------
    ! @brief Quicksort algorithm for sorting a 1D integer array.
    ! @param[inout] arr The integer array to sort.
    ! --------------------------------------------------------------------------
    subroutine quicksort(array)
        implicit none

        integer(int64), intent(inout) :: array(:)                               ! The integer array to sort

        integer(int64) :: n                                                     ! The array size                    
        
        n = size(array)
        
        !$omp parallel
        !$omp single
        call quicksort_internal(array, 1, n)
        !$omp end single
        !$omp end parallel

    contains
        ! ----------------------------------------------------------------------
        ! @brief Internal quicksort subroutine.
        ! @param[inout] arr The integer array to sort.
        ! @param[in] low The lower bound of the array.
        ! @param[in] high The upper bound of the array.
        ! ----------------------------------------------------------------------
        recursive subroutine quicksort_internal(arr, low, high)
            implicit none

            integer(int64), intent(inout) :: arr(:)                             ! The integer array to sort
            integer(int64), intent(in) :: low                                   ! The lower bound of the array    
            integer(int64), intent(in) :: high                                  ! The upper bound of the array
            
            integer(int64) :: pivot_index                                       ! The pivot index    
            integer(int64) :: i                                                 ! The loop index
            integer(int64) :: j                                                 ! The loop index
            
            if (low >= high) return

            i = low
                do j = low + 1, high
                    if (arr(j) <= arr(low)) then
                        i = i + 1
                        block
                            integer :: temp
                            temp = arr(i)
                            arr(i) = arr(j)
                            arr(j) = temp
                        end block
                    end if
                end do
                block
                    integer :: temp
                    temp = arr(i)
                    arr(i) = arr(low)
                    arr(low) = temp
                end block

            pivot_index = i
            
            !$omp task shared(arr) if(high - low > PARALLEL_THRESHOLD)
            call quicksort_internal(arr, low, pivot_index - 1)
            !$omp end task
            
            !$omp task shared(arr) if(high - low > PARALLEL_THRESHOLD)
            call quicksort_internal(arr, pivot_index + 1, high)
            !$omp end task
            
            !$omp taskwait
        end subroutine
    end subroutine

    pure subroutine swap(a, b)
        implicit none

        integer(int64), intent(inout) :: a, b

        integer(int64) :: temp

        !$omp declare simd
        temp = a
        a = b
        b = temp
    end subroutine

    pure subroutine sort4(array)
        implicit none

        integer(int64), intent(inout) :: array(4)

        integer(int64) :: low1, low2, lowest
        integer(int64) :: middle1, middle2, middle
        integer(int64) :: high1, high2, highest

        !$omp declare simd
        if (array(1) < array(2)) then
            low1 = array(1)
            high1 = array(2)
        else
            low1 = array(2)
            high1 = array(1)
        end if

        if (array(3) < array(4)) then
            low2 = array(3)
            high2 = array(4)
        else
            low2 = array(4)
            high2 = array(3)
        end if

        if (low1 < low2) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        end if

        if (high1 > high2) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        end if

        if (middle1 < middle2) then
            middle = middle1
            high1 = middle2
        else
            middle = middle2
            high1 = middle1
        end if

        array = (/lowest, middle, high1, highest/)
    end subroutine sort4

    ! --------------------------------------------------------------------------
    ! @brief Constructor for the Set type.
    ! @return The new Set instance.
    ! --------------------------------------------------------------------------
    function set_constructor()
        implicit none

        type(Set) :: set_constructor
    end function

    ! --------------------------------------------------------------------------
    ! @brief Initialise a new set with an optional initial capacity.
    ! @param[in] initial_capacity The initial capacity of the set.
    ! --------------------------------------------------------------------------
    subroutine set_initialise(self, initial_capacity)
        implicit none
        
        class(Set) :: self
        integer(int64), optional, intent(in) :: initial_capacity


        if (present(initial_capacity)) then
            self%capacity = initial_capacity
        else
            self%capacity = 1024
        endif

        allocate(self%keys(self%capacity))
        allocate(self%occupied(self%capacity))
        self%occupied = .false.
        self%count = 0
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief MurmurHash3 64-bit hash function for integer arrays.
    ! @param[in] key The integer array to hash.
    ! @return The hash value.
    ! --------------------------------------------------------------------------
    function murmur3_hash(self, key) result(hash)
        implicit none

        class(Set), intent(in) :: self
        integer(int64), intent(in) :: key(:)
        integer(int64) :: hash
        integer(int64) :: c1, c2, k
        integer :: i
        integer(int64) :: len

        c1 = int(Z'87C37B91114253D5', int64)
        c2 = int(Z'4CF5AD432745937F', int64)

        hash = 104729_int64
        len = int(size(key) * 8_int64, int64)  ! Length in bytes

        do i = 1, size(key)
            k = key(i)
            k = k * c1
            k = rotl64(k, 31)
            k = k * c2

            hash = ieor(hash, k)
            hash = rotl64(hash, 27)
            hash = hash * 5_int64 + int(Z'52DCE729', int64)
        end do

        ! Finalization
        hash = ieor(hash, len)

        ! Avalanche operations
        hash = ieor(hash, hash / int(Z'100000000', int64))
        hash = hash * int(Z'FF51AFD7ED558CCD', int64)
        hash = ieor(hash, hash / int(Z'100000000', int64))
        hash = hash * int(Z'C4CEB9FE1A85EC53', int64)
        hash = ieor(hash, hash / int(Z'100000000', int64))

        ! Ensure positive hash value within table bounds
        hash = modulo(abs(hash), self%capacity)
    end function

    ! --------------------------------------------------------------------------
    ! @brief Rotate left function for 64-bit integers.
    ! @param[in] x The integer to rotate.
    ! @param[in] n The number of bits to rotate.
    ! @return The rotated integer.
    ! --------------------------------------------------------------------------
    integer(int64) function rotl64(x, n)
        implicit none
        integer(int64), intent(in) :: x
        integer, intent(in) :: n
        integer :: bits, s
        integer(int64) :: left_shifted, right_shifted

        bits = 64
        s = modulo(n, bits)

        if (s == 0) then
            rotl64 = x
        else
            left_shifted = ishft(x, s)
            right_shifted = lshift(x, bits - s)
            rotl64 = ior(left_shifted, right_shifted)
        end if
    end function

    ! --------------------------------------------------------------------------
    ! @brief Logical shift right function for 64-bit integers.
    ! @param[in] x The integer to shift.
    ! @param[in] n The number of bits to shift.
    ! @return The shifted integer.
    ! --------------------------------------------------------------------------
    integer(int64) function lshift(x, n)
        implicit none
        integer(int64), intent(in) :: x
        integer, intent(in) :: n
        integer(int64) :: mask

        if (n <= 0) then
            lshift = x
        else
            mask = ishft(not(0_int64), n)
            lshift = iand(ishft(x, -n), not(mask))
        end if
    end function

    ! --------------------------------------------------------------------------
    ! @brief Find an empty slot in the set using quadratic probing.
    ! @param[in] key The integer array to find.
    ! @return The slot index.
    ! --------------------------------------------------------------------------
    function set_quadratic_probe(self, key) result(slot)
        implicit none
        class(Set), intent(in) :: self
        integer(int64), intent(in) :: key(:)
        integer(int64) :: slot
        integer(int64) :: base_slot
        integer(int64) :: i

        base_slot = self%hash(key) + 1

        do i = 0_int64, self%capacity - 1
            slot = modulo(base_slot + i * i - 1, self%capacity) + 1

            if (.not. self%occupied(slot)) exit
            if (allocated(self%keys(slot)%key)) then
                if (size(self%keys(slot)%key) == size(key)) then
                    if (all(self%keys(slot)%key == key)) exit
                end if
            end if
        end do

        if (i >= self%capacity) then
            slot = -1_int64
        end if
    end function

    ! --------------------------------------------------------------------------
    ! @brief Add an integer array to the set.
    ! @param[in] key The integer array to add.
    ! --------------------------------------------------------------------------
    subroutine set_add(self, key)
        implicit none

        class(Set), intent(inout) :: self                                       ! The set type 
        integer(int64), intent(in) :: key(:)                                    ! The integer array to add  

        integer(int64) :: slot                                                  ! The slot index 
        integer(int64) :: key_len                                               ! The key length 
        
        key_len = size(key)

        ! Resize at 75% capacity
        if (real(self%count) / self%capacity > 0.75) then
            call self%resize(self%capacity * 2)
        endif

        slot = self%quadratic_probe(key)
        if (slot < 0 .or. self%occupied(slot)) return

        if (.not. self%occupied(slot)) then
            allocate(self%keys(slot)%key(size(key)))
            self%keys(slot)%key = key
            self%occupied(slot) = .true.
            self%count = self%count + 1
        endif
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Check if the set contains an integer array.
    ! @param[in] key The integer array to check.
    ! @return True if the set contains the key, false otherwise.
    ! --------------------------------------------------------------------------
    function set_contains(self, key) result(found)
        implicit none

        class(Set), intent(in) :: self                                          ! The set type        
        integer(int64), intent(in) :: key(:)                                    ! The integer array to check

        logical(bool8) :: found                                      ! The found flag

        integer(int64) :: slot                                                  ! The slot index
        
        slot = self%quadratic_probe(key)
        found = slot >= 0 .and. self%occupied(slot)
    end function

    ! --------------------------------------------------------------------------
    ! @brief Remove an integer array from the set.
    ! @param[in] key The integer array to remove.
    ! --------------------------------------------------------------------------
    subroutine set_remove(self, key)
        implicit none

        class(Set), intent(inout) :: self                                       ! The set type
        integer(int64), intent(in) :: key(:)                                    ! The integer array to remove

        integer(int64) :: slot                                                  ! The slot index
        
        slot = self%quadratic_probe(key)
        if (slot >= 0 .and. self%occupied(slot)) then
            deallocate(self%keys(slot)%key)
            self%occupied(slot) = .false.
            self%count = self%count - 1
        end if
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Resize the set to a new capacity.
    ! @param[in] new_capacity The new capacity of the set.
    ! --------------------------------------------------------------------------
    subroutine set_resize(self, new_capacity)
        implicit none

        class(Set), intent(inout) :: self                                       ! The set type
        integer(int64), intent(in) :: new_capacity                              ! The new capacity

        integer(int64) :: i, slot                                               ! Loop indices and slot index
        type(SetNode), allocatable :: temp_keys(:)                              ! Temporary keys array 
        logical(bool8), allocatable :: temp_occupied(:)              ! Temporary occupied array
        integer(int64) :: hash_val, base_slot, j
        integer(int64), allocatable :: key(:)

        allocate(temp_keys(new_capacity))
        allocate(temp_occupied(new_capacity))
        temp_occupied = .false.

        ! Rehash existing keys into new arrays
        do i = 1, self%capacity
            if (self%occupied(i)) then
                key = self%keys(i)%key
                hash_val = self%hash(key)
                base_slot = hash_val + 1
                ! Quadratic probing to find new slot
                do j = 0, new_capacity - 1
                    slot = modulo(base_slot + j * j - 1, new_capacity) + 1
                    if (.not. temp_occupied(slot)) then
                        allocate(temp_keys(slot)%key(size(key)))
                        temp_keys(slot)%key = key
                        temp_occupied(slot) = .true.
                        exit
                    end if
                end do
                deallocate(self%keys(i)%key)
            end if
        end do

        ! Deallocate old arrays and move new ones
        deallocate(self%keys)
        call move_alloc(temp_keys, self%keys)
        deallocate(self%occupied)
        call move_alloc(temp_occupied, self%occupied)
        self%capacity = new_capacity
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Get the number of elements in the set.
    ! @return The number of elements in the set.
    ! --------------------------------------------------------------------------
    function set_size(self) result(count)
        implicit none

        class(Set), intent(in) :: self                                          ! The set type

        integer(int64) :: count                                                 ! The number of elements in the set

        count = self%count
    end function

    ! --------------------------------------------------------------------------
    ! @brief Reset the set to an empty state.
    ! --------------------------------------------------------------------------
    subroutine set_reset(self)
        implicit none
        
        class(Set), intent(inout) :: self                                       ! The set type

        self%occupied = .false.
        self%count = 0
    end subroutine

    ! --------------------------------------------------------------------------
    ! @brief Convert the set to an array of keys.
    ! @param[out] array The array of keys.
    ! --------------------------------------------------------------------------
    function set_to_array(self) result(array)
        implicit none

        class(Set), intent(in) :: self

        integer(int64), allocatable :: array(:, :)                              ! The array of keys

        integer(int64) :: max_key_size                                          ! The maximum key size
        integer(int64) :: i                                                     ! The loop index
        integer(int64) :: j                                                     ! The loop index
        
        ! Find largest key size
        max_key_size = 0
        do i = 1, self%capacity
            if (self%occupied(i)) then
                max_key_size = max(max_key_size, size(self%keys(i)%key))
            end if
        end do
        
        ! Allocate output array with correct dimensions
        allocate(array(max_key_size, self%count))
        array = 0  ! Initialize to zero
        
        ! Fill array with keys
        j = 1
        do i = 1, self%capacity
            if (self%occupied(i)) then
                array(1:size(self%keys(i)%key), j) = self%keys(i)%key
                j = j + 1
            end if
        end do
    end function set_to_array

    ! --------------------------------------------------------------------------
    ! @brief Finalise the set by deallocating all allocated memory.
    ! --------------------------------------------------------------------------
    subroutine set_finalise(self)
        implicit none

        class(Set), intent(inout) :: self
        integer :: i

        do i = 1, self%capacity
            if (self%occupied(i)) then
                if (allocated(self%keys(i)%key)) then
                    deallocate(self%keys(i)%key)
                end if
            end if
        end do

        if (allocated(self%keys)) then
            deallocate(self%keys)
        end if

        if (allocated(self%occupied)) then
            deallocate(self%occupied)
        end if
    end subroutine set_finalise

    ! --------------------------------------------------------------------------
    ! @brief Finalise the set by deallocating all allocated memory.
    ! --------------------------------------------------------------------------
    subroutine set_finalise_auto(self)
        implicit none

        type(Set), intent(inout) :: self
        integer :: i

        if (allocated(self%occupied) .and. allocated(self%keys)) then
            do i = 1, self%capacity
                if (self%occupied(i)) then
                    if (allocated(self%keys(i)%key)) then
                        deallocate(self%keys(i)%key)
                    end if
                end if
            end do
        end if

        if (allocated(self%keys)) then
            deallocate(self%keys)
        end if

        if (allocated(self%occupied)) then
            deallocate(self%occupied)
        end if
    end subroutine set_finalise_auto

    ! --------------------------------------------------------------------------
    ! @brief Finalise the min-heap by deallocating all allocated memory.
    ! --------------------------------------------------------------------------
    subroutine min_heap_finalise(self)
        implicit none

        class(MinHeap), intent(inout) :: self
        integer :: i

        do i = 1, self%count
            if (allocated(self%values(i)%value)) then
                deallocate(self%values(i)%value)
            end if
        end do

        if (allocated(self%values)) then
            deallocate(self%values)
        end if
    end subroutine min_heap_finalise

    ! --------------------------------------------------------------------------
    ! @brief Finalise the min-heap by deallocating all allocated memory.
    ! --------------------------------------------------------------------------
    subroutine min_heap_finalise_auto(self)
        implicit none

        type(MinHeap), intent(inout) :: self
        integer :: i

        do i = 1, self%count
            if (allocated(self%values(i)%value)) then
                deallocate(self%values(i)%value)
            end if
        end do

        if (allocated(self%values)) then
            deallocate(self%values)
        end if
    end subroutine min_heap_finalise_auto
end module data_structures_algorithms