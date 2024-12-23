! -----------------------------------------------------------------------------
! @file geometry.f90
! @brief Module containing geometry functions.
! @author ofgn
! @date 2024-08-10
! -----------------------------------------------------------------------------
module geometry
    use global

    implicit none

contains

    ! --------------------------------------------------------------------------
    ! @brief Calculates the cross product of two 3D vectors.
    ! @param[in] u The first vector (3D coordinates).
    ! @param[in] v The second vector (3D coordinates).
    ! @return The cross product of the two vectors.
    ! --------------------------------------------------------------------------
    pure function cross_product(u, v)
        implicit none

        real(real64), intent(in) :: u(3), v(3)

        real(real64) :: cross_product(3)

        cross_product(1) = u(2) * v(3) - u(3) * v(2)
        cross_product(2) = u(3) * v(1) - u(1) * v(3)
        cross_product(3) = u(1) * v(2) - u(2) * v(1)
    end function cross_product

    ! !--------------------------------------------------------------------------
    ! ! @brief Calculates the signed volume of a tetrahedron.
    ! ! @param[in] a The vertices of the tetrahedron (3x4).
    ! ! @return The signed volume of the tetrahedron.
    ! !--------------------------------------------------------------------------
    pure function signed_volume(a)
        implicit none

        real(real64), intent(in) :: a(3, 4)

        real(real64) :: signed_volume
        real(real64) :: u(3), v(3), w(3), cross(3)
    
        !$omp declare simd

        ! Compute relative vectors
        u = a(:, 1) - a(:, 4)
        v = a(:, 2) - a(:, 4)
        w = a(:, 3) - a(:, 4)
    
        ! Compute cross-product
        cross = cross_product(v, w)

        ! Compute signed volume
        signed_volume = dot_product(u, cross) / 6.0_real64
    end function signed_volume

    !--------------------------------------------------------------------------
    ! @brief Tests if two tetrahedra intersect
    ! @param[in] a Vertices of first tetrahedron (3x4)
    ! @param[in] b Vertices of second tetrahedron (3x4)
    ! @return .true. if tetrahedra intersect, .false. otherwise
    !--------------------------------------------------------------------------
    pure function tetrahedra_intersect(a, b) result(intersect)
        implicit none

        real(real64), intent(in) :: a(3, 4), b(3, 4)

        logical :: intersect

        real(real64) :: volumes(8)

        !$omp declare simd

        ! Test vertices of tet2 against tet1
        volumes(1) = signed_volume([a(:,1), a(:,2), a(:,3), b(:,1)])
        volumes(2) = signed_volume([a(:,1), a(:,2), a(:,3), b(:,2)])
        volumes(3) = signed_volume([a(:,1), a(:,2), a(:,3), b(:,3)])
        volumes(4) = signed_volume([a(:,1), a(:,2), a(:,3), b(:,4)])

        ! Early exit if all vertices are on same side
        if (all(volumes(1:4) > 0) .or. all(volumes(1:4) < 0)) then
            intersect = .false.
            return
        end if

        ! Test vertices of tet1 against tet2
        volumes(5) = signed_volume([b(:,1), b(:,2), b(:,3), a(:,1)])
        volumes(6) = signed_volume([b(:,1), b(:,2), b(:,3), a(:,2)])
        volumes(7) = signed_volume([b(:,1), b(:,2), b(:,3), a(:,3)])
        volumes(8) = signed_volume([b(:,1), b(:,2), b(:,3), a(:,4)])

        ! Check if all vertices are on same side
        intersect = .not. (all(volumes(5:8) > 0) .or. all(volumes(5:8) < 0))
    end function tetrahedra_intersect
end module geometry
