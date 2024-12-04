! ------------------------------------------------------------------------------
! @file maths.F90
! @brief Contains mathematical routines for cubic interpolation.
! @details This module provides cubic and bicubic interpolation routines
!   and supporting utility functions, with OpenMP parallelisation.
! @date 2024-11-28
! ------------------------------------------------------------------------------
module maths
    use iso_fortran_env, only: real64
    use omp_lib
    use utility

    implicit none

contains
    ! ------------------------------------------------------------------------------
    ! @brief Perform bicubic interpolation for multiple points using OpenMP.
    ! @param[in] x_array Array of x-coordinates to interpolate.
    ! @param[in] y_array Array of y-coordinates to interpolate.
    ! @param[in] grid 3D array containing x, y, and z values.
    ! @param[out] z_array Interpolated z-values.
    ! ------------------------------------------------------------------------------
    subroutine bicubic_interpolate_parallel(x_array, y_array, grid, z_array)
        implicit none

        real(real64), intent(in) :: x_array(:), y_array(:)         ! Interpolation coordinates
        real(real64), allocatable, intent(in) :: grid(:, :, :)     ! Input grid (x, y, z)
        real(real64), intent(out) :: z_array(size(x_array))        ! Interpolated z-values

        integer :: i                                               ! Loop index
        integer :: n_points                                        ! Number of points to interpolate

        n_points = size(x_array)

        ! Parallelise the interpolation loop
        !$omp parallel do private(i) shared(x_array, y_array, grid, z_array, n_points)
        do i = 1, n_points
            z_array(i) = bicubic_interpolate(x_array(i), y_array(i), grid)
        end do
        !$omp end parallel do
    end subroutine bicubic_interpolate_parallel

    ! ------------------------------------------------------------------------------
    ! @brief Perform bicubic interpolation on a 2D grid.
    ! @param[in] x X-coordinate to interpolate.
    ! @param[in] y Y-coordinate to interpolate.
    ! @param[in] grid 3D array containing x, y, and z values.
    ! @return Interpolated z-value at (x, y).
    ! ------------------------------------------------------------------------------
    function bicubic_interpolate(x, y, grid) result(z)
        implicit none

        real(real64), intent(in) :: x, y                                        ! Interpolation coordinates
        real(real64), allocatable, intent(in) :: grid(:, :, :)                  ! Input grid (x, y, z)
        real(real64) :: z                                                       ! Interpolated z-value

        integer :: nx, ny                                                       ! Grid dimensions
        integer :: i, j                                                         ! Indices

        real(real64) :: tx, ty                                                  ! Normalised distances
        real(real64), dimension(4, 4) :: p                                      ! Neighbouring z-values
        real(real64), dimension(4) :: cx, cy                                    ! Cubic interpolation weights

        ! Get grid dimensions
        nx = size(grid, 1)
        ny = size(grid, 2)

        ! Locate the bottom-left corner of the cell containing (x, y)
        do i = 2, nx - 2
            if (x >= grid(i, 1, 1) .and. x < grid(i + 1, 1, 1)) exit
        end do
        do j = 2, ny - 2
            if (y >= grid(1, j, 2) .and. y < grid(1, j + 1, 2)) exit
        end do

        ! Handle edge cases
        if (i < 2 .or. i > nx - 2 .or. j < 2 .or. j > ny - 2) then
            z = grid(i, j, 3)
            return
        end if

        ! Compute normalised distances within the cell
        tx = (x - grid(i, 1, 1))/(grid(i + 1, 1, 1) - grid(i, 1, 1))
        ty = (y - grid(1, j, 2))/(grid(1, j + 1, 2) - grid(1, j, 2))

        ! Collect the 4x4 z-values for the surrounding grid points
        p(1, :) = grid(i - 1:i + 2, j - 1, 3)
        p(2, :) = grid(i - 1:i + 2, j, 3)
        p(3, :) = grid(i - 1:i + 2, j + 1, 3)
        p(4, :) = grid(i - 1:i + 2, j + 2, 3)

        ! Compute cubic interpolation weights
        cx = cubic_weights(tx)
        cy = cubic_weights(ty)

        ! Perform bicubic interpolation
        z = sum(cx*matmul(p, cy))

    contains
        ! ----------------------------------------------------------------------
        ! @brief Compute cubic interpolation weights for a given parameter t.
        ! @param[in] t Interpolation parameter (0 <= t <= 1).
        ! @return Weights for cubic interpolation.
        ! ----------------------------------------------------------------------
        function cubic_weights(t) result(weights)
            implicit none

            real(real64), intent(in) :: t                                       ! Interpolation parameter
            real(real64), dimension(4) :: weights                               ! Interpolation weights

            weights(1) = -0.5d0*t**3 + t**2 - 0.5d0*t
            weights(2) = 1.5d0*t**3 - 2.5d0*t**2 + 1.0d0
            weights(3) = -1.5d0*t**3 + 2.0d0*t**2 + 0.5d0*t
            weights(4) = 0.5d0*t**3 - 0.5d0*t**2
        end function cubic_weights
    end function bicubic_interpolate

end module maths
