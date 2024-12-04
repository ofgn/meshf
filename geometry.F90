module geometry
    use iso_fortran_env, only: int32, int64, real64, real128

    implicit none
    
    contains
    ! --------------------------------------------------------------------------
    ! @brief Calculates the orientation of a tetrahedron.
    ! @param[in] a The first vertex of the tetrahedron (3D coordinates).
    ! @param[in] b The second vertex of the tetrahedron (3D coordinates).
    ! @param[in] c The third vertex of the tetrahedron (3D coordinates).
    ! @param[in] d The fourth vertex of the tetrahedron (3D coordinates).
    ! @return The orientation of the tetrahedron. Positive if vertices are
    !   ordered counter-clockwise, negative if clockwise, and zero if coplanar.
    ! --------------------------------------------------------------------------
    function orientation3d(a, b, c, d) result(orientation)
        implicit none

        real(real64), intent(in) :: a(3), b(3), c(3), d(3)
        real(real128) :: orientation
        real(real128) :: ax, ay, az, bx, by, bz, cx, cy, cz

        ! Ensure input vertices are valid
        if (size(a) /= 3 .or. size(b) /= 3 .or. size(c) /= 3 .or. size(d) /= 3) then
            error stop "All input vertices must be 3D coordinates."
        end if

        ! Compute vector differences
        ax = real(a(1), kind=real128) - real(d(1), kind=real128)
        ay = real(a(2), kind=real128) - real(d(2), kind=real128)
        az = real(a(3), kind=real128) - real(d(3), kind=real128)
        bx = real(b(1), kind=real128) - real(d(1), kind=real128)
        by = real(b(2), kind=real128) - real(d(2), kind=real128)
        bz = real(b(3), kind=real128) - real(d(3), kind=real128)
        cx = real(c(1), kind=real128) - real(d(1), kind=real128)
        cy = real(c(2), kind=real128) - real(d(2), kind=real128)
        cz = real(c(3), kind=real128) - real(d(3), kind=real128)

        ! Compute determinant
        orientation = ax * (by * cz - bz * cy) &
            - ay * (bx * cz - bz * cx) &
            + az * (bx * cy - by * cx)
    end function orientation3d

end module geometry
