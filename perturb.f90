subroutine wave_prep(N, M, lambda_k, A, phi)

 ! =====================================================
 !  Accepts number of harmonics, N, and max wavelength, lambda_k,
 !  then returns randomized amplitude, A, and phase shift, phi
 ! =====================================================

implicit none

integer, intent(in)                                 :: lambda_k, N, M
real(8), intent(out), dimension(-N:N, -M:M)         :: phi
real(8), intent(out), dimension(-N:N, -M:M)         :: A
real(8)                                             :: random_normal, randn
integer                                             :: i, j

do i =  -N, N
    do j = -M, M
        if (((i**2 + j**2) .gt. 0) .and. ((i**2 + j**2) .lt. N**2 + 1)) then
            randn = random_normal()
            A(i, j) = (lambda_k - 5*randn) * sqrt((real(i)**2 + real(j)**2)/(2*N**2))
        else
            A(i, j) = 0
        end if
        phi(i, j) = 10 * rand()
    end do
end do

end subroutine



subroutine fourier_series(pts, phi, A, lambda_k, sample_points, N, M, pts_perturbed, double_sum, B)

 ! =====================================================
 !  Accepts x y z points in the pts array, with other
 !  parameters for the fourier series perturbation process
 ! =====================================================

implicit none

integer, intent(in)                                       :: sample_points, N, M
real(8), intent(in)                                       :: lambda_k
real(8), intent(in), dimension(-N:N, -M:M)                :: phi, A
real(8), intent(in), dimension(0:sample_points-1,3)       :: pts
real(8), intent(out), dimension(0:sample_points-1,3)      :: pts_perturbed
real(8), intent(out), dimension(0:sample_points-1)        :: double_sum
real(8), intent(out), dimension(-N:N, -N:N)               :: B
real(8), dimension(-N:N, -N:N)                            :: N_mat, M_mat
real(8), dimension(1:2)                                   :: bounds
real(8), dimension(0:sample_points-1)                     :: x, y, z, r, th
real(8)                                                   :: pi = 4.0 * ATAN(1.0)
real(8)                                                   :: max_height, half_height, inv_maxr, inv_lambda
integer                                                   :: i, j

do i = 0, sample_points-1
  x(i) = pts(i,1) * 1000     ! axial coordinates in mils
  y(i) = pts(i,2) * 1000     ! radial coordinates in mils
  z(i) = pts(i,3) * 1000     ! radial coordinates in mils
end do

! Set axial boundaries (stagnation + 10 mils to [end - X] mils)
bounds = (/minval(x(:))+50, maxval(x(:))-250/)

! Map radial coordinates into cylindrical coordinate system
do i = 0, sample_points-1
  th(i) = atan2(z(i),y(i))
  r(i) = sqrt(y(i)**2 + z(i)**2)
end do

! Create N_mat and M_mat to avoid a triple do loop
do i = -N, N
  N_mat(i,:) = i
end do

do j = -M, M
  M_mat(:,j) = j
end do

max_height = 100 * lambda_k ! perturbations greater than this quantity are not acceptable
half_height = (0.5 * lambda_k)
inv_maxr = 1 / maxval(r(:))
inv_lambda = 1 / lambda_k


do i = 0, sample_points-1
    if ((x(i).gt.bounds(1)).and.(x(i).le.(bounds(2)))) then
        B = cos((2 * pi * N_mat * th(i) * inv_lambda) + (2 * pi * M_mat * x(i) * inv_lambda) + phi)
        double_sum(i) = sum(A*B)
        if (abs(double_sum(i)) .gt. max_height) then
            double_sum(i) = 0
        end if
        r(i) = r(i) + abs(double_sum(i)) * r(i) * inv_maxr
    end if
end do

! Convert coordinates back to cartesian
x(:) = x(:) / 1000              ! axial coordinates, mils to inches
y(:) = r(:) * cos(th(:)) / 1000 ! radial coordinates, mils to inches
z(:) = r(:) * sin(th(:)) / 1000 ! radial coordinates, mils to inches

! Gather perturbed x y z coordinates
do i = 0, sample_points-1
  pts_perturbed(i,1) = x(i)
  pts_perturbed(i,2) = y(i)
  pts_perturbed(i,3) = z(i)
end do

end subroutine



real(8) function random_normal()

 ! =====================================================
 !  Calculates pseudo-random normal (gaussian) numbers
 ! =====================================================

implicit none

  real(8)                   :: r1
  real(8)                   :: r2
  real(8)                   :: x
  real(8)                   :: pi = 4.0 * ATAN(1.0)

  r1 = rand()
  r2 = rand()
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * pi * r2 )
  random_normal = x

return
end function random_normal
