!Main program. To compile and run, use command make LinAl, and then ./ LinAl which will print                  
!all results to the screen for Cholesky Decomposition                                                               
!all subroutine inputs come from LinAl.f90                                                                     

!output includes calculating solution of the least squares problem                                             
!which can be viewed as a minimization problem on the norm of the residual, r = b-Ax over all x                

! - Cholesky Decomposition

Program Driver_Cholesky

  use utility, only: dp

  use LinAl, only:  mTrace, vec2Norm, matPrint, readMat, Vandermonde, &
                      mat_normalEqs, Cholesky, Cholesky_backSub

  implicit none

!  integer, save :: msize, nsize                                                                               
!  real, dimension(:,:), allocatable, save :: mat                                                              


  real(dp), allocatable :: A_col(:,:),B_col(:,:),A(:,:),A_n(:,:),B_n(:,:)
  real(dp), allocatable :: X_s(:,:),X(:,:),A_s(:,:),B_s(:,:)
  real(dp), allocatable :: resids(:)
  real(dp), allocatable :: E(:,:)

  character(len=100) :: myFileName
  integer :: i,j
  real(dp) :: trace
  real(dp) :: norm2
  Logical :: sTest 
  integer :: N 
  sTest = .false.
  norm2 = 0.0
  trace = 0.0
  
  allocate(A_col(21,1))
  allocate(B_col(21,1))

  A_col = 0.0
  B_col = 0.0

  !-------- Read in data fine into two seperate column vectors --------
  myFileName = 'least_squares_data.dat'
  open(10,file=myFileName)

  do i = 1, 21
     read(10,*) A_col(i,1),B_col(i,1)
  end do

  close(10)

  print*, ""
  print*, "Matrix column of x_i's from data file"
  call matPrint(A_col,shape(A_col))
  
  print*, ""
  print*, "Matrix column of y_i's from data file"
  call matPrint(B_col,shape(B_col))

  !-------- Fit with a 3rd degree polynomial -> N = 3 --------
  N = 3

  print*, ""
  print*, "----Fitting with polynomial of degree",N



  !-------- Allocating rest of variables for the problem --------
  allocate(A(size(A_col,1),N+1))
  allocate(A_n(N+1,N+1))
  allocate(B_n(N+1,1))
  allocate(X(size(B_n,1),size(B_n,2)))
  allocate(E(size(B_n,1),size(B_n,2)))
  allocate(A_s(size(A_n,1),size(A_n,2)))
  allocate(B_s(size(B_n,1),size(B_n,2)))
  allocate(resids(size(A_col,1)))

  A = 0.0
  A_s = 0.0
  B_s = 0.0
  A_n = 0.0
  B_n = 0.0
  X = 1.0
  E = 0.0


  !-------- Form Vandermonde style matrix A --------
  call Vandermonde(A_col,A,N)

!  print*, ""
!  print*, "Vandermonde style matrix A"
!  call matPrint(A,shape(A))


  ! -------- Form matrices A_n and B_n associated with normal equaitons --------
  call mat_normalEqs(A,B_col,A_n,B_n) 
  A_s = A_n !save matries associated with normal equations A_n and B_n as A_s and B_s
  B_s = B_n 

  print*, ""
  print*, "Matrix A_n = (A^T)A  associated with normal equations:"
  call matPrint(A_s,shape(A_s))
  
  print*, ""
  print*, "Matrix B_n = (A^T)B formed associated with normal equations:"    
  call matPrint(B_s,shape(B_s))

  ! -------- Cholesky Decomposition and Back-substitution --------
  print*, ""
  print*, "-------- Performing Cholesky Decomposition and back-substitution --------"

  call Cholesky(A_n,sTest)
  print*, ""
  print*, "Cholesky Decomposed matrix A "
  call matPrint(A_n,shape(A_n))

  call Cholesky_backSub(A_n,B_n,X)

  print*, ""
  print*, "Solution matrix X "
  call matPrint(X,shape(X))

  !---- calculate error matrix ---- 
  !and 2-norm of each column (only one column if B rhs only has one column vector)
  E = matmul(A_s,X)
  E = E - B_s

  do j = 1, size(E,2)
     norm2 = 0.0
     call vec2Norm(E(:,j),norm2)
     print*, ""
     print*, "Euclidean norm error of solution vector =",norm2
  end do

  !---- output the fitted curve to a file ----
  open(1, file = 'fittedCurve_3degree_Ch.dat')
!  write(1,*) X(1,1),'+',X(2,1),'x','+',X(3,1),'x^2','+',X(4,1),'x^3'
  do i = 1, size(X,1)
     write(1,*) X(i,:)
  end do
  close(1)

  !---- calculate 2norm error between fitted curve and data using vec2Norm subroutine
  !first form vector of difference between fit and data for each point ----
  resids = 0.0
  do i = 1,size(A_col(:,1))
     resids(i) = B_col(i,1) - (X(1,1) + X(2,1)*A_col(i,1) &
     + X(3,1)*A_col(i,1)**2.0 + X(4,1)*A_col(i,1)**3.0)
  end do

  call vec2Norm(resids,norm2)
  print*, ""
  print*, "2 norm error bewteen fitted polynomial curve and data:",norm2




  !-----------------!!!------------------

  !-------- Next, fit with a 5th degree polynomial --------
  !The only variables we want to keep from previous allocation are A_col and B_col

  deallocate(A)
  deallocate(A_s)
  deallocate(B_s) 
  deallocate(A_n)
  deallocate(B_n)
  deallocate(X)
  deallocate(E)
  deallocate(resids)

  !-------- Fit with a 5th degree polynomial -> N = 5 --------                                                                                 
  N = 5

  print*, ""
  print*, "----Fitting with polynomial of degree",N


  !-------- Allocating rest of variables for the problem --------                                                                              

  allocate(A(size(A_col,1),N+1))
  allocate(A_n(N+1,N+1))
  allocate(B_n(N+1,1))
  allocate(X(size(B_n,1),size(B_n,2)))
  allocate(E(size(B_n,1),size(B_n,2)))
  allocate(A_s(size(A_n,1),size(A_n,2)))
  allocate(B_s(size(B_n,1),size(B_n,2)))
  allocate(resids(size(A_col,1)))

  A = 0.0
  A_s = 0.0
  B_s = 0.0
  A_n = 0.0
  B_n = 0.0
  X = 1.0
  E = 0.0

  !-------- Form Vandermonde style matrix A --------                                                                                           

  call Vandermonde(A_col,A,N)
 ! print*, ""
 ! print*, "Vandermonde style matrix A"
 ! call matPrint(A,shape(A))

  ! -------- Form matrices A_n and B_n associated with normal equaitons --------                                                              

  call mat_normalEqs(A,B_col,A_n,B_n)
  A_s = A_n !save matries associated with normal equations A_n and B_n as A_s and B_s                                                          
  B_s = B_n

  print*, ""
  print*, "Matrix A_n = (A^T)A  associated with normal equations:"
  call matPrint(A_s,shape(A_s))

  print*, ""
  print*, "Matrix B_n = (A^T)B formed associated with normal equations:"
  call matPrint(B_s,shape(B_s))


  ! -------- Cholesky Decomposition and Back-substitution --------                                                                             
  print*, ""
  print*, "-------- Performing Cholesky Decomposition and back-substitution --------"

  call Cholesky(A_n,sTest)
  print*, ""
  print*, "Cholesky Decomposed matrix A "
  call matPrint(A_n,shape(A_n))

  call Cholesky_backSub(A_n,B_n,X)

  print*, ""
  print*, "Solution matrix X "
  call matPrint(X,shape(X))

  !calculate error matrix and 2 norm of each column (only one column if B rhs only has one column vector)                                      
  E = matmul(A_s,X)
  E = E - B_s

  do j = 1, size(E,2)
     norm2 = 0.0
     call vec2Norm(E(:,j),norm2)
     print*, ""
     print*, "Euclidean norm error of solution vector =",norm2  
  end do

  !output the fitted curve to a file                                                                                                          

  open(1, file = 'fittedCurve_5degree_Ch.dat')
!  write(1,*) X(1,1),'+',X(2,1),'x','+',X(3,1),'x^2','+',X(4,1),'x^3','+',X(5,1),'x^4','+',X(6,1),'x^5'                                         
  do i = 1, size(X,1)
     write(1,*) X(i,:)
  end do
  close(1)

  !calculate 2norm error between fitted curve and data using vec2Norm subroutine                                                              
  !first form vector of difference between fit and data for each point                                                                         

  resids = 0.0
  do i = 1,size(A_col(:,1))
     resids(i) = B_col(i,1) - (X(1,1) + X(2,1)*A_col(i,1) &
     + X(3,1)*A_col(i,1)**2.0 + X(4,1)*A_col(i,1)**3.0 &
     + X(5,1)*A_col(i,1)**4.0 + X(6,1)*A_col(i,1)**5.0)
  end do

  call vec2Norm(resids,norm2)
  print*, ""
  print*, "2 norm error bewteen fitted polynomial curve and data:",norm2


  
  deallocate(A)
  deallocate(A_s)
  deallocate(B_s)
  deallocate(A_n)
  deallocate(B_n)
  deallocate(X)
  deallocate(E)
  deallocate(resids)
  
  deallocate(A_col)
  deallocate(B_col)

  
  

End Program Driver_Cholesky
