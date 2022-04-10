!Main program. To compile and run, use command make LinAl, and then ./ LinAl which will print                  
!all results to the screen for Householder QR Factorization                                                               
!all subroutine inputs come from LinAl.f90                                                                     

!output includes calculating solution of the least squares problem                                             
!which can be viewed as a minimization problem on the norm of the residual, r = b-Ax over all x                

! - QR solution of the least-squares problem using Householder method                                                                     

Program Driver_QR

  use utility, only: dp
  use LinAl, only:  mTrace, vec2Norm, matPrint, readMat, Vandermonde, &
                    mat_normalEqs, QR_houseH, backSub 

  implicit none

  real(dp), allocatable :: A_col(:,:),B_col(:,:),A(:,:)
  real(dp), allocatable :: X_s(:,:),X(:,:),A_s(:,:),B_s(:,:)
  real(dp), allocatable :: resids(:)
  real(dp), allocatable :: Q(:,:), eye(:,:), AAt(:,:), Qterm(:,:), QtermQtermT(:,:),Qtb(:,:)
  real(dp), allocatable :: Aterm(:,:)
  real(dp), allocatable :: R_nlength(:,:), Qtb_nlength(:,:)
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


  allocate(A(size(A_col,1),N+1))
  A = 0.0

  !-------- Form Vandermonde style matrix A --------                                                                                 
  call Vandermonde(A_col,A,N)

!  print*, ""
!  print*, "Vandermonde style matrix A"
!  call matPrint(A,shape(A))


  ! -------- Householder QR Factorization --------
  ! Utilize Vandermonde matrix A, and y column vector in data file, B_col
  ! updated m x n matrix R from QR decomp (still labeled A after goes through QR_houseH)
  ! and m x m matrix Q also formed in the QR_houseH subroutine 
 
  allocate(A_s(size(A,1),size(A,2)))
  allocate(AAt(size(A,1),size(A,1)))
  allocate(Q(size(A,1),size(A,1)))
  allocate(eye(size(A,1),size(A,1)))
  allocate(Qterm(size(Q,1),size(Q,2)))
  allocate(QtermQtermT(size(Q,1),size(Q,2)))
  allocate(Aterm(size(A,1),size(A,2)))
  allocate(Qtb(size(B_col,1),size(B_col,2)))
  allocate(X(size(A,2),1))
  allocate(R_nlength(size(A,2),size(A,2)))
  allocate(Qtb_nlength(size(A,2),1))
  allocate(E(size(A,2),1))
  allocate(resids(size(A_col,1)))

  A_s = 0.0
  A_s = A
  Aterm = 0.0
  AAt = 0.0
  Q = 0.0
  Qterm = 0.0
  QtermQtermT = 0.0
  R_nlength = 0.0
  Qtb_nlength = 0.0
  X = 0.0
  E = 0.0

  eye = 0.0  !fill in ideneity matrix with ones in diagonal
  do i = 1, size(eye,2)
     eye(i,i) = 1.0
  end do
!  print*, "Identity matrix formed has shape: ", shape(eye)
  

  !-------- perform QR decomposition on A --------
  call QR_houseH(A,Q,eye)  !run this algorithm on the Vandermonde style matrix A
  print*, ""
  print*, "A decomposed into R:"
  call matPrint(A,shape(A))


  !-------- compute and print A - QR --------
  Aterm = A_s - matmul(Q,A)

  print*, ""
  print*, "A - QR : "
  call matPrint(Aterm,shape(Aterm))

  !---- compute Frobenius-norm of A-QR (matrix frob norm = sqrt(Trace(AA^*)) ----
  AAt = matmul(Aterm,transpose(Aterm))  !A^T(A) used in frobenius norm
  call mTrace(AAt,trace)

  print*, ""
  print*, "Frobenius norm of A - QR:",sqrt(trace)

  !-------- compute and print Q^T(Q) - I and its Frobenius norm ---------
  !NOTE --- I dont truncate to the 6 row matrices until calculating the solution vector
  !     --- Therefore Q^TQ - I is 21x21 matrix so I wont print it to screen. I just print its error norm

  Qterm = matmul(transpose(Q),Q) - eye  !21x21 matrix unless I truncate it to first n rows beforehand
 ! call matPrint(Qterm,shape(Qterm))

  QtermQtermT = matmul(Qterm,transpose(Qterm)) !to use for frobenius norm
  call mTrace(QtermQtermT,trace)

  print*, ""
  print*, "Frobenius norm of Q^T(Q) - I:",sqrt(trace)

  ! ------!!!------
  !solve least square equation projected onto span of A, namely Rx = Q^T(b)
  !print solution vector x
  !verify solution correct by calculating its 2-norm error

  !use backsub subroutine to solve Rx = Q^T(b) = Qtb
  !only need to use the first n rows of the matrices to solve for x!

  Qtb = matmul(transpose(Q),B_col)

!  print*, ""
!  print*, "Full size Q^Tb:"
!  call matPrint(Qtb,shape(Qtb))

!  R_nlength = A(1:size(A,2),:)
!  Qtb_nlength = Qtb(1:size(A,2),:)
!  Dont need these 2 variables above. Can just input the truncated matrices into backSub 
  call backSub(A(1:size(A,2),:),Qtb(1:size(A,2),:), X)

  print*, ""
  print*, "Solution matrix X "
  call matPrint(X,shape(X))

  !calculate error matrix and 2 norm of each column (only one column if B rhs only has one column vector)                                                                          
  E = matmul(A(1:size(A,2),:),X)
  E = E - Qtb(1:size(A,2),:)

  do j = 1, size(E,2)
     norm2 = 0.0
     call vec2Norm(E(:,j),norm2)
     print*, ""
     print*, "Euclidean norm error of solution vector =",norm2
  end do

  !---- output the fitted curve to a file ----                                                                                                                                    

  open(1, file = 'fittedCurve_3degree_QR.dat')
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
  deallocate(AAt)
  deallocate(Q)
  deallocate(eye)
  deallocate(Qterm)
  deallocate(QtermQtermT)
  deallocate(Aterm)
  deallocate(Qtb)
  deallocate(X)
  deallocate(R_nlength)
  deallocate(Qtb_nlength)
  deallocate(E)
  deallocate(resids)

  !-------- Fit with a 5th degree polynomial -> N = 5 --------                                                          
                                                                                                                          
  N = 5

  print*, ""
  print*, "----Fitting with polynomial of degree",N

  allocate(A(size(A_col,1),N+1))
  A = 0.0

  !-------- Form Vandermonde style matrix A --------                                                                                              
  
  call Vandermonde(A_col,A,N)
 ! print*, ""
 ! print*, "Vandermonde style matrix A"
 ! call matPrint(A,shape(A))

  !-------- Allocating rest of variables for the problem --------      

  allocate(A_s(size(A,1),size(A,2)))
  allocate(AAt(size(A,1),size(A,1)))
  allocate(Q(size(A,1),size(A,1)))
  allocate(eye(size(A,1),size(A,1)))
  allocate(Qterm(size(Q,1),size(Q,2)))
  allocate(QtermQtermT(size(Q,1),size(Q,2)))
  allocate(Aterm(size(A,1),size(A,2)))
  allocate(Qtb(size(B_col,1),size(B_col,2)))
  allocate(X(size(A,2),1))
  allocate(R_nlength(size(A,2),size(A,2)))
  allocate(Qtb_nlength(size(A,2),1))
  allocate(E(size(A,2),1))
  allocate(resids(size(A_col,1)))

  A_s = 0.0
  A_s = A

  Aterm = 0.0
  AAt = 0.0
  Q = 0.0
  Qterm = 0.0
  Qtb = 0.0
  QtermQtermT = 0.0
  R_nlength = 0.0
  Qtb_nlength = 0.0
  X = 0.0
  E = 0.0

  eye = 0.0  !fill in ideneity matrix with ones in diagonal                                                                                       
  do i = 1, size(eye,2)
     eye(i,i) = 1.0
  end do
!  print*, "Identity matrix formed has shape: ", shape(eye)                                                                                                                          

  !-------- perform QR decomposition on A --------                                                                                           

  call QR_houseH(A,Q,eye)  !run this algorithm on the Vandermonde style matrix A                                                                  

  print*, ""
  print*, "A decomposed into R:"
  call matPrint(A,shape(A))

  !-------- compute and print A - QR --------                                                                                                      
  Aterm = A_s - matmul(Q,A)

  print*, ""
  print*, "A - QR : "
  call matPrint(Aterm,shape(Aterm))

  !---- compute Frobenius-norm of A-QR (matrix frob norm = sqrt(Trace(AA^*)) ----                                                                  
  AAt = matmul(Aterm,transpose(Aterm))  !A^T(A) used in frobenius norm                                                                           
  call mTrace(AAt,trace)

  print*, ""
  print*, "Frobenius norm of A - QR:",sqrt(trace)

  !-------- compute and print Q^T(Q) - I and its Frobenius norm ---------                                                                          
  Qterm = matmul(transpose(Q),Q) - eye  !print this 21 x 21 matrix?                                                                              
 ! call matPrint(Qterm,shape(Qterm))                                                                                                             
  QtermQtermT = matmul(Qterm,transpose(Qterm)) !to use for frobenius norm                                                                          
  call mTrace(QtermQtermT,trace)

  print*, ""
  print*, "Frobenius norm of Q^T(Q) - I:",sqrt(trace)
  
  ! ------!!!------                                                                                                                                 
  !solve least square equation projected onto span of A, namely Rx = Q^T(b)                                                                        
  !print solution vector x                                                                                                                          
  !verify solution correct by calculating its 2-norm error                                                                                          
  !use backsub subroutine to solve Rx = Q^T(b) = Qtb                                                                                               
  !only need to use the first n rows of the matrices to solve for x                                                                                
  
  Qtb = matmul(transpose(Q),B_col)
!  R_nlength = A(1:size(A,2),:)                                                                                                                
!  Qtb_nlength = Qtb(1:size(A,2),:)                                                                                                                

  call backSub(A(1:size(A,2),:),Qtb(1:size(A,2),:), X)

  print*, ""
  print*, "Solution matrix X "
  call matPrint(X,shape(X))

  !calculate error matrix and 2-norm of each column (only one column if B rhs has one column vector)
  E = matmul(A(1:size(A,2),:),X)
  E = E - Qtb(1:size(A,2),:)

  do j = 1, size(E,2)
     norm2 = 0.0
     call vec2Norm(E(:,j),norm2)
     print*, ""
     print*, "Euclidean norm error of solution vector =",norm2  
  end do

  !---- output the fitted curve to a file ----                                                                                                     
  open(1, file = 'fittedCurve_5degree_QR.dat')
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


  deallocate(A_col)
  deallocate(B_col)
  deallocate(A)
  deallocate(A_s)
  deallocate(AAt)
  deallocate(Q)
  deallocate(eye)
  deallocate(Qterm)
  deallocate(QtermQtermT)
  deallocate(Aterm)
  deallocate(Qtb)
  deallocate(X)
  deallocate(R_nlength)
  deallocate(Qtb_nlength)
  deallocate(E)
  deallocate(resids)

End Program Driver_QR
