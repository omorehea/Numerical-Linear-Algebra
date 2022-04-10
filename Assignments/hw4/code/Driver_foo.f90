!Main driver program. To compile and run, use command make LinAl, and then ./LinAl which will print all 
!the desired results to the screen for:

!Reduction to Hessenberg (tridiagonal) form for a Normal & Singular (symmetric) matrix
!QR algorithm (i) without and (ii) with shift
!Inverse iteration algorithm to calculate corresponding eigenvectors given the eigenvalues

Program Driver_foo

  use utility, only: dp

  use LinAl, only: mTrace, vec2Norm, matPrint, readMat, Vandermonde, houseTri, QR_houseH, backSub,&
                   inv_Iter, QR_ews_noShift, QR_ews_Shift

  implicit none

  real(dp), allocatable :: A(:,:),A_s(:,:)
  real(dp), allocatable :: Q(:,:), eye(:,:)
  real(dp), allocatable :: ews(:)
  real(dp), allocatable :: x(:,:)

  integer(dp) :: iters

  integer, dimension (1:2) :: order1 = (/1, 2/) !variable to use when creating matrix 
  integer, dimension (1:2) :: order2 = (/2, 1/)  


  integer :: i,j
  Logical :: sTest


  sTest = .false.


  allocate(A(4,4))
  A = 0.0

  
  !-------- Hessenberg matrix reduction to tridiagonal form --------!
  print*, "---- Reduce matrix A to tridiagonal form ----"

  !shape A with row major using reshape and order2, standard is using order1, column major
  A(:,:) = reshape((/5,4,1,1, 4,5,1,1, 1,1,4,2, 1,1,2,4/), (/4,4/), order = order2) 
  
  print*, ""
  print*, "Matrix A"
  call matPrint(A,shape(A))

  call houseTri(A)
 
  print*, ""
  print*, "Reduced matrix A"
  call matPrint(A,shape(A))

  deallocate(A)
  
  !-------- QR algorithm to find eigenvalues -- With and without shift --------
  print*, ""
  print*, "---- QR algorithm (i) without shift and (ii) with shift to calculate eigenvalues of matrix A ----"

  allocate(A(3,3))
  allocate(A_s(size(A,1),size(A,2)))
  A = 0.0
  A_s = 0.0
  A(:,:) = reshape((/3,1,0, 1,2,1, 0,1,1/), (/3,3/), order = order2)
  A_s = A

  print*, ""
  print*, "Matrix A"
  call matPrint(A_s,shape(A_s))

  allocate(Q(size(A,1),size(A,2)))
  allocate(eye(size(A,1),size(A,2)))
!  allocate(R(size(A,1),size(A,2)))
  allocate(ews(size(A,1)))

  Q = 0.0
 ! R = 0.0
  ews = 0.0

  eye = 0.0 !fill in identity matrix with ones in diagonal
  do i = 1, size(eye,2)
     eye(i,i) = 1.0
  end do

  !-------- QR algorithm to find eigenvalues -- WITHOUT SHIFT --------
  
  print*, ""
  print*, "---- QR algorithm WITHOUT SHIFT to find eigenvalues ----"

  iters = 0
  call QR_ews_noShift(A,Q,eye,iters)

  print*, ""
  print*, "Updated diagonal matrix holding eigenvalues"

  call matPrint(A,shape(A))

  do i = 1,size(A,1)
     ews(i) = A(i,i)
  end do
  print*, ""
  print*, "Eigenvalues of A:", ews
  
  print*, "number of iterations:",iters



  !-------- QR algorithm to find eigenvalues -- WITH SHIFT --------
  print*, ""
  print*, "------ QR algorithm WITH SHIFT to find eigenvalues ------"

  A = A_s !restart matrix A back to original matrix
  Q = 0

  ews = 0

  iters = 0

  call QR_ews_Shift(A,Q,eye,iters)

  print*, ""
  print*, "Updated diagonal matrix holding eigenvalues"
  call matPrint(A,shape(A))

  do i = 1,size(A,1)
     ews(i) = A(i,i)
  end do
  print*, "number of iterations:",iters
  print*, ""
  print*, "Eigenvalues of A:", ews

!----------------------------------------------

  deallocate(A)
  deallocate(Q)
  deallocate(ews)
  deallocate(eye)

  !-------- Inverse iteration method to calculate eigenvectors(appx) given the eigenvalue(appx)
  print*, ""
  print*, "------ Inverse iteration method to calculate eigenvectors ------"

  allocate(A(4,4))
  allocate(ews(4))

  A = 0
  A(:,:) = reshape((/2,1,3,4, 1,-3,1,5, 3,1,6,-2, 4,5,-2,-1/), (/4,4/), order = order2)

  ews = 0
  ews = (/-8.0286,7.9329,5.6689,-1.5732/) !vector of provdied eigenvalues

  allocate(eye(size(A,1),size(A,2))) !fill in identity matrix
  eye = 0
  do i = 1, size(eye,2)
     eye(i,i) = 1.0
  end do

  allocate(x(size(A,1),1))
 
  do i = 1, size(ews)
     x(:,:) = reshape((/0,1,0,0/), (/size(A,1),1/), order = order1) !create column vector x, with ||x|| = 1
     call inv_Iter(A,ews(i),x,eye)

     do j = 1, (size(x,1))  !normalizing eigenvectors such that bottom-most eigenvector is 1
        x(j,1) = x(j,1)/x(size(x,1),1)
     end do

     print*, ""
     print*, "eigenvalue:",ews(i)
     print*, ""
     print*, "corresponding eigenvectors:"
     call matPrint(x,shape(x))

  end do


deallocate(A)
deallocate(ews)
deallocate(x)

end program Driver_foo
