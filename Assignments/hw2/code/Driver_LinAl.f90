!Main program. To compile and run, use command make LinAl, and then ./LinAl which will print
!all results to the screen for problems 2,3,4,5.
!all subroutine inputs come from LinAl.f90

!output includes solving the system Ax=B by means of: 
!-Gaussian Elimination with Partial Pivoting
!-LU Decomposition with Partial Pivoting

!lastly, gaussian elim is used to solve for the equation of a plane given 3 cartesian points in 3d space.


Program Driver_LinAl

  use utility, only: dp
  use LinAl, only: msize, nsize, mTrace, vec2Norm, matPrint, readMat, gausElim_pp, backSub, LU_pp, LU_backSub

  implicit none

!  integer, save :: msize, nsize
!  real, dimension(:,:), allocatable, save :: mat


  real(kind=8), allocatable :: A(:,:),B(:,:)
  real(kind=8), allocatable :: A_s(:,:),B_s(:,:)
  real(kind=8), allocatable :: X_s(:,:),X(:,:)
  integer ::  dimsA(2), dimsB(2)
  real(kind=8), allocatable :: E(:,:)
  integer(kind=8), allocatable :: s(:)

  real(kind=8),allocatable :: points(:,:) !block of variables here for problem 5, plane connecting 3 points
  real(kind=8),allocatable :: b_vec(:,:)
  real(kind=8),allocatable :: sol_vec(:,:)
  real(kind=8) :: pi_ref = acos(-1.d0)
  real(kind=8) :: e1 = exp(1.d0)

  character(len=100) :: myFileName
  integer :: i,j
  real(kind=8) :: trace
  real(kind=8) :: norm2
  Logical :: sTest 
  sTest = .false.
  
  myFileName = 'Amat.dat'

  open(10,file=myFileName) !read dimensions of A
  read(10,*) msize,nsize
  close(10)


  allocate(A_s(msize,nsize))
  allocate(A(msize,nsize))
  A = 0.0
  A_s = 0.0
  

  call readMat(myFileName,A_s) !read and allocate matrix A into A_s and A
  A = A_s  !copy A as saved matrix A_s
   
  myFileName = 'Bmat.dat'

  open(10,file=myFileName) !read dimensions of B
  read(10,*) msize,nsize
  close(10)

  allocate(B_s(msize,nsize)) !read and allocate matrix  B into B_s and B
  allocate(B(msize,nsize))
  B = 0.0
  B_s = 0.0

  allocate(X_s(msize,nsize))  !allocate solution vector X_s and its duplicate X to all ones
  allocate(X(msize,nsize))
  X_s = 1.0
  X = X_s

  allocate(E(size(X_s,1),size(X_s,2))) !allocate dimensions of E matrix to same size as B_s and X_s


  call readMat(myFileName,B_s) !read and allocate matrix B into B   
  B = B_s  !copy B as saved matrix B_s


!  do i = 1, msize
!     write(*,*) (mat(i,j) , j = 1, nsize )
!  end do
  
  !filling matrix dimensions into arrays
  dimsA(1) = size(A_s,1)
  dimsA(2) = size(A_s,2)
  dimsB(1) = size(B_s,1)
  dimsB(2) = size(B_s,2)

!-------Problem 2--------
!Print matrix, trace, norm

  print *,"Matrix A"
  call matPrint(A_s,dimsA)
 
  trace = 0.0
  call mTrace(A_s,trace)
  print *, 'Trace =', trace

  do j = 1, size(A_s,2)
     norm2 = 0.0d0
     call vec2Norm(A_s(:,j),norm2)
     print *, 'Euclidean norm of column',j,'=',norm2
  end do


!  do j = 1, size(B_s,2)
!     norm2 = 0.0d0
!     call vec2Norm(B_s(:,j),size(B_s,1),norm2)
!     print *, 'Euclidean norm of column',j,'=',norm2
!  end do


!--------Problem 3--------
!Gaussian Elimination

  print *,"Matrix A"
  call matPrint(A_s,dimsA)

  print *, ""
  print *,"Matrix B"
  call matPrint(B_s,dimsB)

  print*,""
  print*,"...Performing Gaussian Elimination with Partial Pivoting..."


  call gausElim_pp(A,B,dimsA,dimsB,sTest)
  print *,""
  print *,"Matrix A"
  call matPrint(A,dimsA)

  print *,""
  print *,"Matrix B"
  call matPrint(B,dimsB)

  print*,""
  print*,"...Performing Backsubstitution..."

  call backSub(A,B,X)

  print *,""
  print *,"Solution matrix X"
  call matPrint(X,shape(X))

  !calculate error matrix
  E = matmul(A_s,X)
  E = E - B_s

  print *,""
  print *,"Error Matrix, E = A_sX - B_s"
  call matPrint(E,shape(E))

  !calculate norm of each of column vectors of Error matrix
  do j = 1, size(E,2)
     norm2 = 0.0d0
     call vec2Norm(E(:,j),norm2)
     print *,""
     print *, "Euclidean norm of Error matrix column",j,"=",norm2
  end do

  
!-------Problem 4-------
!LU factorization

  print *,""
  print *,"...Performing LU Decomposition with Partial Pivoting..."
  
!Reinitialize variables:
  E = 0.0
  X = X_s
  A = A_s
  B = B_s
  sTest = .false.
  allocate(s(size(A,1)))
  s = 0.0

  print *,""
  print *,"Matrix A"
  call matPrint(A,dimsA)

!  print *,""
!  print *,"Matrix B"
!  call matPrint(B,size(B,1),size(B,2))

  call LU_pp(A,size(A,1),sTest,s)
  print *,""
  print *,"Decomposted matrix A = LU"
  call matPrint(A,dimsA)
  
  print*, ""
  print*, "...Performing LU Backsubstitution..."

  call LU_backSub(A,size(A,1),B,s,X)
  print *,""
  print *,"Solution matrix X"
  call matPrint(X,shape(X))

  !calculate error matrix                                                                                                                                       
  E = matmul(A_s,X)
  E = E - B_s

  print *,""
  print *,"Error Matrix, E = A_sX - B_s"
  call matPrint(E,shape(E))

  !calculate norm of each of column vectors of Error matrix                                                                                                     
  do j = 1, size(E,2)
     norm2 = 0.0d0
     call vec2Norm(E(:,j),norm2)
     print *,""
     print *, "Euclidean norm of Error matrix column",j,"=",norm2
  end do

!-------Problem 5-------
!solve for constant used to define equation of plane passing through three points

  allocate(points(3,3))
  allocate(b_vec(size(points,1),1))
  allocate(sol_vec(size(points,1),1))
  sol_vec = 1.d0

  sTest = .false.

  points(:,:) = reshape((/1.d0,-3.d0,pi_ref,2.d0,2.d0,e1,3.d0,5.d0,-1.d0*sqrt(2.d0)/),(/3,3/))
  b_vec(:,:) = reshape((/1.d0,1.d0,1.d0/),(/3,1/))

  print*, ""
  print*, "...Finding equation of plane passing through three points..."
  print*, "Matrix of points and rhs vector, augmented matrix:"
  ! solving: a(x-x_0) + b(y-y_0) + c(z-z_0) = d
  ! for each point as a system of linear equations
  
  do i = 1, size(points,1)
     print*, points(i,:), "|", b_vec(i,:)
  end do

  print*, "...Solving using Gaussian Elimination..."
  
  call gausElim_pp(points,b_vec,shape(points),shape(b_vec),sTest)


  call backSub(points,b_vec,sol_vec)

  print *,""
  print *,"Solution matrix X"
  call matPrint(sol_vec,shape(sol_vec))

  print *,""
  print *,"Solution matrix X normalized w.r.t a"
  call matPrint(sol_vec/sol_vec(1,1),shape(sol_vec))
 
  E = 0.0
  E = matmul(points,sol_vec)
  E = E - b_vec

  print *,""
  print *,"Error Matrix, E"
  call matPrint(E,shape(E))

!deallocating 

  deallocate(points)
  deallocate(b_vec)
  deallocate(sol_vec)

  deallocate(A)
  deallocate(B)
  deallocate(X_s)
  deallocate(X)
  deallocate(E)
  deallocate(s)

End Program Driver_LinAl
