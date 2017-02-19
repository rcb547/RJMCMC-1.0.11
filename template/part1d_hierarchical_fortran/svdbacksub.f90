
subroutine svdbacksub(U, S, V, B, W, X, n)

  double precision, dimension(n, n), intent(in) :: U
  double precision, dimension(n), intent(in) :: S
  double precision, dimension(n, n), intent(in) :: V
  double precision, dimension(n), intent(in) :: B
  double precision, dimension(n), intent(out) :: W
  double precision, dimension(n), intent(out) :: X

  integer, intent(in) :: n

  integer :: i
  integer :: j

  do i = 1, n

     W(i) = 0.0

     do j = 1, n
        W(i) = W(i) + U(j, i)*B(j)
     end do
     
     W(i) = W(i)/S(i)

  end do

  do i = 1, n
     
     X(i) = 0.0

     do j = 1, n
        X(i) = X(i) + W(j)*V(j, i)
     end do

  end do

end subroutine svdbacksub

