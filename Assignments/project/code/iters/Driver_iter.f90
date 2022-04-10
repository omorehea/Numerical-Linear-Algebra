Program Driver_iter

  use utility, only: dp

  use LinAl, only:  mTrace, vec2Norm, matPrint, readMat, Vandermonde, &
                      mat_normalEqs, Cholesky, Cholesky_backSub,&
                      g_jacobi, g_seidel, CG_smart, CG_precond

  implicit none


  real(dp), allocatable :: A(:,:), M_precond(:,:)
  real(dp), allocatable :: b(:), x(:)
  real(dp) :: accuracy, D
  real(dp) :: upper_tol

  integer :: i, j, m, l_gj, l_gs, l_cg

  character (len = 5) :: user_input

  m = 10  !size of symmetric matrix
  upper_tol = 1.0d15 !upper limit value to determine divergence


  allocate(A(m,m))
  allocate(M_precond(m,m))
  allocate(b(m))
  allocate(x(m))
 

  A = 1.d0  !matrix A full of ones except on diagonal where A_ii = D

  !---- D is value entered by user ----
  print*, ""
  print*, "Please input a value for D: "
  read*, D
  print*, ""

  !D = 1000

  do i = 1,size(A,2)  !fill diagonals of A with value D
     A(i,i) = D
  end do

  
  do i = 1,size(b)    !fill vector b_i = i
     b(i) = i
  end do


  !---- Ask user whether to use Gauss-Jacobi or Gauss-Seidel ----


  


  !---- Calling the specfied routine to solve Ax = b      ----
  !---- with accuracy 10^-5, and prints results to screen ----

  l_gj = 1
  l_gs = 1
  l_cg = 1

  accuracy = 1.0d-5


  M = 0
  do i = 1, size(A,1)
     M_precond(i,i) = A(i,i)
  end do
  



  x = 0  !initial solution guess

!  call matPrint(A,shape(A))

  print*, ""
  print*, "Performing Gauss-Jacobi algorithm..."
  call g_jacobi(A,b,x,accuracy,l_gj,upper_tol)
  print*,""
  print*,"Solution vector x"
  do i = 1,size(x)
     print*,x(i)
  end do
  print*, ""
  print*, "Total number of iterations:", l_gj


  x = 0
  
  print*, ""
  print*, "Performing Gauss-Seidel algorithm..."
  call g_seidel(A,b,x,accuracy,l_gs,upper_tol)
  print*,""
  print*,"Solution vector x"
  do i = 1,size(x)
     print*,x(i)
  end do
  print*, ""
  print*, "Total number of iterations:", l_gs

  x = 0  

  print*, ""
  print*, "Performing Conjugate Gradient algorithm..."
  call CG_smart(A,b,x,accuracy,l_cg)
  print*,""
  print*,"Solution vector x"
  do i = 1,size(x)
     print*,x(i)
  end do
  print*, ""
  print*, "Total number of iterations:", l_cg


  x = 0
  l_cg = 0

  print*, ""
  print*, "Performing Conjugate Gradient with Diagonal Preconditioner..."
  call CG_precond(A,b,x,M_precond,accuracy,l_cg)
  print*,""
  print*,"Solution vector x"
  do i = 1,size(x)
     print*,x(i)
  end do
  print*, ""
  print*, "Total number of iterations:", l_cg



end program Driver_iter
