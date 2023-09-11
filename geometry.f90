!> @file geometry.f90
!> @brief Provides geometry-related utility functions.
module geometry
   implicit none

   !> @var q_prec
   !> Kind parameter for quad-precision real numbers.
   integer, parameter :: q_prec = selected_real_kind(33, 4931)

   !> @var eps
   !> Small constant for numerical comparisons.
   real(q_prec) :: eps = 1.0E-32_q_prec

contains

   !> @brief Initialize the value of eps.
   !> @param val New value for eps.
   subroutine init_eps(val)
      real(q_prec), intent(in) :: val  !< New eps value.
      eps = val
   end subroutine init_eps

   !> @brief Check if a point lies inside the circumcircle of a triangle.
   !> @param p Test point.
   !> @param u First vertex of the triangle.
   !> @param v Second vertex of the triangle.
   !> @param w Third vertex of the triangle.
   !> @return True if the point is inside, false otherwise.
   function in_circle(p, u, v, w) result(is_inside)
      real, intent(in) :: p(2), u(2), v(2), w(2)  !< Input points.
      logical :: is_inside                        !< Output result.
      real(q_prec) :: det, px, py, ux, uy, vx, vy, wx, wy

      ! Convert to quad precision
      px = real(p(1), q_prec)
      py = real(p(2), q_prec)
      ux = real(u(1), q_prec)
      uy = real(u(2), q_prec)
      vx = real(v(1), q_prec)
      vy = real(v(2), q_prec)
      wx = real(w(1), q_prec)
      wy = real(w(2), q_prec)

      ! Calculate determinant in quad-precision
      det = (ux - px)*((vy - py)*(wx*wx + wy*wy - px*px - py*py) - (vy*vy + vx*vx - px*px - py*py)*(wy - py)) &
            - (uy - py)*((vx - px)*(wx*wx + wy*wy - px*px - py*py) - (vy*vy + vx*vx - px*px - py*py)*(wx - px)) &
            + (ux*ux + uy*uy - px*px - py*py)*((vx - px)*(wy - py) - (vy - py)*(wx - px))

      ! Compare against eps
      is_inside = det > eps
   end function in_circle

   function on_edge(p, u, v) result(is_on_edge)
      real, intent(in) :: p(2), u(2), v(2)  !< Input points.
      logical :: is_on_edge                  !< Output result.
      real(q_prec) :: det, px, py, ux, uy, vx, vy
      real(q_prec) :: pu(2), pv(2)

      ! Convert to quad precision
      px = real(p(1), q_prec)
      py = real(p(2), q_prec)
      ux = real(u(1), q_prec)
      uy = real(u(2), q_prec)
      vx = real(v(1), q_prec)
      vy = real(v(2), q_prec)

      ! Calculate vectors pu and pv
      pu = [px, py] - [ux, uy]
      pv = [px, py] - [vx, vy]

      ! Calculate determinant in quad-precision
      det = pu(1)*pv(2) - pu(2)*pv(1)

      ! Compare against eps
      is_on_edge = abs(det) < eps
   end function on_edge

   !> @brief Check orientation of a triangle uvw.
   !> @param u First vertex of the triangle.
   !> @param v Second vertex of the triangle.
   !> @param w Third vertex of the triangle.
   !> @return True if positively oriented, false otherwise.
   function orient2(u, v, w) result(is_pos_oriented)
      real, intent(in) :: u(2), v(2), w(2)  !< Input points.
      logical :: is_pos_oriented             !< Output result.
      real(q_prec) :: det, ux, uy, vx, vy, wx, wy

      ! Convert to quad precision
      ux = real(u(1), q_prec)
      uy = real(u(2), q_prec)
      vx = real(v(1), q_prec)
      vy = real(v(2), q_prec)
      wx = real(w(1), q_prec)
      wy = real(w(2), q_prec)

      ! Calculate determinant in quad-precision
      det = ux*(vy - wy) + vx*(wy - uy) + wx*(uy - vy)

      ! Compare against eps
      is_pos_oriented = det > eps
   end function orient2

end module geometry
