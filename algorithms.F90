module algorithms
    use iso_fortran_env, only: int64
    use omp_lib

    implicit none

    private
    public :: quicksort_1d
    
    integer(int64), parameter :: SERIAL_THRESHOLD = 2**16
    
contains
    subroutine quicksort_1d(arr)
        integer(int64), allocatable, intent(inout) :: arr(:)
        integer(int64) :: n
        
        if (.not. allocated(arr)) return
        if (size(arr) <= 1) return
        
        n = size(arr, kind=int64)
        
        !$omp parallel
        !$omp single
        call quicksort_recursive(arr, 1_int64, n)
        !$omp end single
        !$omp end parallel
    end subroutine quicksort_1d
    
    recursive subroutine quicksort_recursive(arr, low, high)
        integer(int64), intent(inout) :: arr(:)
        integer(int64), intent(in) :: low, high
        integer(int64) :: pivot_idx
        
        if (low < high) then
            if (high - low <= SERIAL_THRESHOLD) then
                ! Do serial sort inline
                pivot_idx = partition(arr, low, high)
                call quicksort_recursive(arr, low, pivot_idx - 1)
                call quicksort_recursive(arr, pivot_idx + 1, high)
                return
            end if
            
            pivot_idx = partition(arr, low, high)
            
            !$omp task shared(arr) if(high - low > SERIAL_THRESHOLD)
            call quicksort_recursive(arr, low, pivot_idx - 1)
            !$omp end task
            
            !$omp task shared(arr) if(high - low > SERIAL_THRESHOLD)
            call quicksort_recursive(arr, pivot_idx + 1, high)
            !$omp end task
            
            !$omp taskwait
        end if
    end subroutine quicksort_recursive
    
    function partition(arr, low, high) result(pivot_idx)
        integer(int64), intent(inout) :: arr(:)
        integer(int64), intent(in) :: low, high
        integer(int64) :: pivot_idx
        
        integer(int64) :: pivot, i, j, temp
        
        pivot = arr(low)
        i = low
        
        do j = low + 1, high
            if (arr(j) <= pivot) then
                i = i + 1
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
            end if
        end do
        
        temp = arr(i)
        arr(i) = arr(low)
        arr(low) = temp
        
        pivot_idx = i
    end function partition

end module algorithms