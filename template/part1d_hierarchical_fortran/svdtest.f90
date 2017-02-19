program svdtest

character*1 :: jobz
double precision, dimension(5,5) :: A
double precision, dimension(5) :: S
double precision, dimension(5,5) :: U, VT
double precision, dimension(200) :: WORK

double precision, dimension(5) :: B
double precision, dimension(5) :: X
double precision, dimension(5) :: W

integer, dimension(100) :: IWORK

integer :: i

jobz = "A"

A(1,1) = 0.658771  
A(1,2) = 0.992937   
A(1,3) = 0.215037   
A(1,4) = 0.764384   
A(1,5) = 0.828987

A(2,1) = 0.359896   
A(2,2) = 0.361115   
A(2,3) = 0.948443   
A(2,4) = 0.939531   
A(2,5) = 0.159497

A(3,1) = 0.055671   
A(3,2) = 0.112755   
A(3,3) = 0.792010   
A(3,4) = 0.784137   
A(3,5) = 0.812457

A(4,1) = 0.882290   
A(4,2) = 0.936349   
A(4,3) = 0.843152   
A(4,4) = 0.454622   
A(4,5) = 0.196823

A(5,1) = 0.916061   
A(5,2) = 0.744779   
A(5,3) = 0.866891   
A(5,4) = 0.834291   
A(5,5) = 0.931443

call dgesdd(jobz, 5, 5, A, 5, S, U, 5, VT, 5, WORK, -1, IWORK, INFO)

write (*,*) WORK(1)

call dgesdd(jobz, 5, 5, A, 5, S, U, 5, VT, 5, WORK, 200, IWORK, INFO)

write (*,*) S

write (*,*) "U"

do i = 1, 5
   write (*,*) U(i, :)
end do

write (*,*) "V"

do i = 1, 5
   write (*,*) VT(i, :)
end do

B = 0.0
B(1) = 1.0

call svdbacksub(U, S, VT, B, W, X, 5)

write (*,*) "X (B = [1 0 0 0 0 ])"

write (*,*) X

end program svdtest
