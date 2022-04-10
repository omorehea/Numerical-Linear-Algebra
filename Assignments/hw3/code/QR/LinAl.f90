module LinAl
 
  use  utility, only: dp

  implicit none
  
!  integer, parameter :: dp = kind(0.d0)  !same declaration as in utility.f90
  integer, save :: msize, nsize

  
contains

  !********************************************************

  subroutine mTrace(mat,tr)
    ! subroutine: mtrace
    ! purpose: Takes matrix mat and calcualtes its trace.
    ! inputs: mat --  mxm square matrix.
    ! outputs: tr -- trace of A.

    implicit none
    real(dp),intent(in) :: mat(:,:)
    real(dp),intent(out) :: tr
    integer :: i

    if (.not. size(mat,1) == size(mat,2)) then
       print *, "Matrix is not symmetric!"
       stop
    end if

    tr = 0.0
    do i = 1, size(mat,1)
       tr = tr + mat(i,i)
    end do
  end subroutine mTrace

  subroutine vec2Norm(V,normTwo)
    ! subroutine: vec2Norm
    ! purpose: Takes vector V and its dimension dimV and calculates its Euclidean norm (2-norm)
    ! inputs: V -- Rank 1 array and its dimension, size(V).
    ! outputs: normTwo -- Euclidean norm of V
    implicit none
    real(dp),intent(in) :: V(:)
    real(dp),intent(out) :: normTwo
    integer :: i

    normTwo = 0.0
    do i = 1, size(V)
       normTwo = normTwo + (abs(V(i)**2))
    end do
    normTwo = sqrt(normTwo)
  end subroutine vec2Norm

  subroutine matPrint(mat,mdims)
    ! subroutine: matPrint
    ! purpose: Prints matrix and its dimensions
    ! inputs: mat -- mxn matrix
    !         mdims -- m and n dimensions of matrix mat respectively
   
    implicit none
    integer :: i
    integer,intent(in) :: mdims(2)
    real(dp),intent(in) :: mat(:,:)
    print "('Matrix has dimensions:',i2,' x',i2)",mdims(1),mdims(2)
    do i = 1,mdims(1)
       print*, mat(i,:)
      ! 600 format(F30.3)
    end do
  end subroutine matPrint
    
  subroutine gausElim_pp(A,B,dimsA,dimsB,sTest)
    !subroutine: gausElim_pp
    !purpose: Perform Gaussian Elimination with partial pivoting on matrix A with 
    !         corresponding operations on matrix B containing n rhs vectors.
    !inputs/outputs: A -- m x m real square matrix.
    !        B -- non-square real matrix containing n rhs vectors.
    !        dimsA / dimsB -- respective dimensions of A and B.
    !        sTest - logical flag that indicates whether the problem is singular or not.

    implicit none
    real(dp),intent(inout) :: A(:,:)
    real(dp),intent(inout) :: B(:,:)
    integer,intent(in) :: dimsA(2),dimsB(2)
    real(dp) :: tempA(size(A,2)), tempB(size(B,2))
    Logical, intent(inout) :: sTest
    real(dp) :: factor, tol
    integer :: i,j,K
!    real(kind=8) :: p

    tol = 1.0d-15 !tolerance for singular matrix check

    !----- gaussian elimination with  partial pivoting -----
    do j = 1, dimsA(2)-1
       !find index K  and pivot p s.t p = |a_{Kj}| = max_{k = j,...,m}|a_{Kj}|
       K = j
       do i = K+1, dimsA(1) !scan  over all rows i in column j and store row index of largest row element
          if (abs(A(i,j)) > abs(A(K,j))) then
             K = i
          end if
       end do
!       print*, "K = ",K
       if (.not. K == j) then  !swap rows of A and B if necessary
          tempA = A(j,:)
          A(j,:) = A(K,:)
          A(K,:) = tempA
          
          tempB = B(j,:)
          B(j,:) = B(K,:)
          B(K,:) = tempB
       end if

       if (abs(A(j,j)) < tol) then
          sTest = .true.
          print *, "Matrix is singular!"
          stop
       end if

       do i = j+1, dimsA(1) !loop over rows below j, transformation of remaining submatrix and RHS vector
          factor = A(i,j)/A(j,j)  
          A(i,:) = A(i,:) - factor*A(j,:)
          B(i,:) = B(i,:) - factor*B(j,:)

       end do

    end do

  end subroutine gausElim_pp


  subroutine backSub(U, B2, X)
    !subroutine: backSub
    !purpose: Performs a backsubstitution to return the solution to the problem AX = B
    !         where x is  amatrix containing the n solution vectors.
    !inputs/outputs: U -- m x m upper-diagonal matrix, for instance, the output, A, of the gausElim_pp subroutine.
    !        B2 -- a rhs matrix of n column vectors, i.e., the output, B, of the gausElim_pp subroutine.
    !        X -- matrix that  will contain the n solution vectors upon output.
    implicit none
    real(dp),intent(inout) :: U(:,:)
    real(dp),intent(inout) :: B2(:,:)  
    real(dp),intent(inout) :: X(:,:)
    real(dp) :: sum(size(B2,2))
    real(dp) :: tol
    integer :: i,j,K

    tol = 1.0d-15 !tolerance to check for singular matrix

    if (abs(U(size(U,1),size(U,2))) < tol) then !check for singular matrix
       print*,"Matrix U is singular!"
       stop
    end if

    do j = 1,size(B2,2)     ! X(m,:) = B2(m,:)/U(m,m) where m = size(U,1)
       X(size(B2,1),j) = B2(size(B2,1),j)/U(size(U,1),size(U,2))
    end do

    do i = size(U,1)-1, 1, -1  !loop over lines, bottom to top
       if (abs(U(i,i)) < tol) then
          print*,"Matrix U is singular!"
          stop
       end if
       
       sum = 0.0

       do j = 1, size(B2,2) !loop over columns          
          do K = i+1, size(U,1) 
             sum(j) = sum(j) + U(i,K)*X(K,j) !add to sum variable
          end do
       end do
       
       do j = 1, size(B2,2)
          X(i,j) = (B2(i,j)-sum(j))/U(i,i)  !update solution matrix
       end do
    end do

  end subroutine backSub


  subroutine LU_pp(A,dimA,sTest,s)
    !subroutine: LU_pp
    !purpose: Performs an LU decomposition routine with partial pivoting
    !inputs/outputs: A -- real square m x m matrix
    !        dimA -- the first dimension  of  A (i.e., m)
    !        sTest -- logical flag that indicates whether problem is singular or not
    !        s -- permutation vector of length m

    implicit none
    real(dp), intent(inout) :: A(:,:)
    integer, intent(in) :: dimA
    integer(dp), intent(inout) :: s(dimA)
    Logical, intent(inout) :: sTest
    real(dp) :: tempA(size(A,2))
    real(dp) :: tol
    integer(dp) :: tempS
    integer :: i,j,K
    

    tol= 1.0d-15 !tolerance for zero pivot check
    
    !---- initialize permutation vector ----
    do j = 1, dimA
       s(j) = j
    end do

    do j = 1, dimA    !loop over columns
       K = j
       !---- first perform partial pivoting ----
       do i = K, dimA !scan  over all rows i in column j and store row index of largest row element                                                             
          if (abs(A(i,j)) > abs(A(K,j))) then
             K = i
          end if
       end do
      ! print*, "K = ",K
       if (.not. K == j) then  !swap rows if necessary                                                                                                               
          tempA = A(j,:)
          A(j,:) = A(K,:)
          A(K,:) = tempA

          tempS = s(j)  !record permutation if rows switch   
          s(j) = s(K)
          s(K) = tempS
       end if

       if (abs(A(j,j)) < tol) then
          sTest = .true.
          print*, "pivot is zero, stopping program"
          stop
       end if

       do i = j+1,dimA
          A(i,j) = A(i,j)/A(j,j)  !create the l(i,j) and stores them in A(i,j)
          do K = j+1,dimA
             A(i,K) = A(i,K) - A(i,j)*A(j,K) !updates A into LU matrix
          end do
       end do
    end do

  end subroutine LU_pp


  subroutine LU_backSub(A,dimA,B,s,X)
    !subroutine: LU_backSub
    !purpose: LU backsubstitution routine that solves Ax = B
    !inputs/outputs: A -- real square matrix A decomposed into LU using the subroutine LU_pp
    !        dimA -- first dimension of A
    !        B -- real matrix B containing n rhs vectors
    !        s -- permutation vector
    !        X -- matrix of n solution vectors upon output

    implicit none
    real(dp), intent(inout) :: A(:,:)
    integer, intent(in) :: dimA
    real(dp), intent(inout) :: B(:,:)
    real(dp), intent(inout) :: X(:,:)
    integer(dp), intent(inout) :: s(dimA)
!    real(dp) :: tempB(1,size(B,2))
    real(dp) :: sum(size(B,2))
    real(dp) :: tol
    integer :: i,j,K
    
    tol = 1.0d-15 !tolerance for zero diagonal check

   ! print*,"s =", s
    
    !to initialize y (new B variable) with P*B, y(i,:) = B(s(i),:) variable notation here
    !instead of doing this directly I just replace all B(i,:) with B(s(i),:) within this routine 

 !   do i = 1, dimA
 !      tempB = B(s(i),:)
 !      B(s(i),:) = B(i,:)
 !      B(i,:) = tempB
       
 !   end do

    !---- forward substitution, y = L^{-1}*P*B ----
    do j = 1, dimA-1
       !code to do y = M_j*y
       do i = j+1, dimA
          B(s(i),:) = B(s(i),:) - B(s(j),:)*A(i,j)
       end do
    end do

    !---- backwards substitution, Ux = y ----
    do i = dimA, 1, -1
       !loop over rows, bottom to top
       if (abs(A(i,i)) < tol) then
          print*, "A = (LU)_{i,i} = 0"
          stop
       end if
       sum = 0.0

       do j = 1, size(B,2) !run through each column vector in B
          
          do K = i+1, dimA
!             print*, "i = ", i, "k = ",K
!             print*, "A(i,K) =", A(i,K)
             sum(j) = sum(j) + A(i,K)*X(K,j)
          end do
       end do

       do j = 1, size(B,2)
          X(i,j) = (B(s(i),j)-sum(j))/A(i,i) !update solution matrix
       end do      

    end do
  end subroutine LU_backSub


  subroutine Cholesky(A, sTest)
    !subroutine: Cholesky                                                                                                             
    !purpose: Cholesky decomposition routine that replaces square matrix A by its Cholesky decomposition upon output
    !inputs/outputs: A -- real square matrix A that is replaced with its Cholesky decomposition in  the lower triangular
    !      sTest -- logical flag that indicates whether problem is singular or not
   

    implicit none
    real(dp), intent(inout) :: A(:,:)
!    real(dp), intent(inout) :: A_diag(:,:)
    Logical, intent(inout) :: sTest
    real(dp) :: factor, tol
    integer :: i,j,K
    factor = 0.0
    tol= 1.0d-15

    if (size(A,1) .NE. size(A,2)) then
       sTest = .true.
       print*, "Matrix not symmetric!"
       stop
    end if

    do j = 1, size(A,2)  !loop over columns
       
       if (abs(A(j,j)) .LE. tol) then  !Matrix A must be SPD, symmetric positive definite (all diagonal elements > 0)
          sTest = .true.
          print *, "Matrix diagonal less than or equal to 0"
          stop
       end if

       !---- calculate new diagonal element ----
       do k = 1, j-1
          factor = A(j,k)*A(j,k)  !A(j,k)*A(j,k)^T = A(j,k)*A(k,j), no need to use Transpose function since only one element
          A(j,j) = A(j,j) - factor
         
       end do
       A(j,j) = sqrt(A(j,j))
!       print*, "A(j,j) for j =",j," is ", A(j,j)
 
       factor = 0.0
       !---- calculate elements below diagonal ----
       do i = j+1, size(A,1)
          do k = 1, j-1
             factor = A(i,k)*A(j,k)
             A(i,j) = A(i,j) - factor
          end do

          A(i,j) = A(i,j)/A(j,j)
       end do
    end do

  end subroutine Cholesky

  subroutine Cholesky_backSub(A, B, X)
    !subroutine: Cholesky_backSub
    !purpose: Perform Cholesky backsubstitution that solves LL^Tx = b
    !inputs/outputs: A -- real square matrix already Cholesky-decomposed from the Cholesky subroutine
    !      B -- matrix containing n rhs vectors
    !      X -- solution matrix containing n rhs vectors 
    
    implicit none
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: B(:,:)
    real(dp), intent(inout) :: X(:,:)
    real(dp) :: sum(size(B,2))
    real(dp) :: B_factor(size(B,2))
    real(dp) :: Y(size(B,1),size(B,2))
    real(dp) :: tol
    integer :: i,j,k

    tol= 1.0d-15

    !---- Forward substitution, solving Ly = b ----

    Y = B
    sum = 0.0
    do i = 1, size(A,1)
       
       sum(:) = B(i,:)
       do j = 1, i-1
          sum(:) = sum(:) - Y(j,:)*A(i,j)
       end do
       Y(i,:) = sum(:)/A(i,i)
    end do
    
    !---- Backward substitution, solving (L^*)x = y ----
    do i = size(A,1),1,-1
       if (abs(A(i,i)) < tol) then
          print*, "Diagonal element", i, "of A is zero"
          stop
       end if
       do k = i+1, size(A,1)
          Y(i,:) = Y(i,:) - A(k,i)*X(k,:)
       end do
       X(i,:) = Y(i,:)/A(i,i)
    end do

  end subroutine Cholesky_backSub

  subroutine QR_houseH(A,Q,eye)
   !subroutine: QR_houseH
   !purpose: Householder QR factorization algorithm
   !inputs/outputs: A -- input is n x n  matrix A, output is A transformed into R
   !                Q -- output is the n x n calculated matrix Q
   !                eye -- identity matrix (m x m)
    implicit none
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: Q(:,:),eye(:,:)
    real(dp) :: v(size(A,1),size(A,2))
    real(dp) :: s
    integer :: i,j
    integer :: sgn
    real(dp) :: Href(size(eye,1),size(eye,2))
    real(dp) :: vtA(1,size(A,2))
!    real(kind=8) :: v_transp(size(A,1),1)
    
    v = 0.0
    do j = 1, size(A,2) !loop over columns
       
       if (A(j,j) .GE. 0) then
          sgn = 1.0
       else
          sgn = -1.0
       end if
       s = sgn*norm2(A(j:size(A,1),j))  !---- compute signed norm

       v(j,j) = A(j,j) + s    !---- compute v(j,j), different from following row values
       do i = j+1, size(A,1)
          v(i,j) = A(i,j)     !---- row values below v(j,j) equal a(j+1,j)...a(m,j)
       end do
       v(:,j) = v(:,j)/norm2(v(:,j))  !---- normalized Householder vector

       vtA = matmul(transpose(v(:,j:j)),A) !---- calculate (v^T)A
       A = A - 2.0*matmul(v(:,j:j),vtA)  !---- update A -> A = A -2*v*((v^T)*A)
    end do

    !---- calculating matrix Q = (I-2v_1(v_1^T))...(I - 2v_n(v_n^T)) ----
    Q = (eye - 2.0*matmul(v(:,1:1),transpose(v(:,1:1))))
    do i = 2, size(v,2)
       Href = (eye - 2.0*matmul(v(:,i:i),transpose(v(:,i:i)))) 
       Q = matmul(Q,Href)
    end do

  end subroutine QR_houseH

  subroutine Vandermonde(A_col,A_v,N)
    !subroutine: Vandermonde
    !purpose: Form the vandermonde style matrix from the m x_i's in a column vector
    !inputs/outputs: A_col -- m x 1 column extracted from given data file (extracted in Driver program)
    !                A_v -- m x n Vandermonde matrix formed for a n-1 degree polynomial 
    !                N -- the degree of the polynomial we wish to fit the data to

    implicit none
    real(dp), intent(inout) :: A_col(:,:)
    real(dp), intent(inout) :: A_v(:,:)
    integer, intent(inout) :: N
    integer :: i
    
    A_v(:,1) = 1.0
    do i = 2,N+1   !N degree polynomial fitting                                                                                   
       A_v(:,i) = A_col(:,1)**(i-1)
    end do

  end subroutine Vandermonde

  subroutine mat_normalEqs(A,B,A_n,B_n)
    !subroutine: mat_normalEqs
    !purpose: First forms the Vandermonde style matrix A,
    !         then forms the matrices associated with the normal equation, (A^T)Ax = (A^T)B
    !inputs/outputs: A -- m x n Vandermonde style matrix A
    !                A_n -- A_n is an n x n matrix defined by the normal equatin, (A^T)A
    !                B -- B is a matrix of n column vectors with m rows
    !                B_n -- B_n is the new matrix accociated with the normal equation, (A^T)B                   

    implicit none
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: B(:,:)
    real(dp), intent(out) :: A_n(:,:)
    real(dp), intent(out) :: B_n(:,:)
    

    A_n = matmul(transpose(A),A)
    B_n = matmul(transpose(A),B)

!    print*, "shape A_n",shape(A_n)
!    print*, "shape B_n",shape(B_n)
  end subroutine mat_normalEqs

  subroutine houseTri(A)
    !subroutine: houseTri
    !purpose: Reduce a symmetric matrix to tridiagonal form using Householder matrices for the similarity transformations
    !         Otherwise known as, Reduction to Hessenberg form of square matrix
    !inputs/outputs: A - m x m symmetric matrix

    implicit none
    real(dp), intent(inout) :: A(:,:)
    real(dp) :: v(size(A,1),size(A,2))
    real(dp) :: s
    integer :: sgn
    real(dp) :: vtA(1,size(A,2))
    real(dp) :: Av(size(A,1),1)
    integer :: i,j

    
    v = 0.0
    do j = 1, (size(A,1)-1) !loop over columns                                                                                      
       if (A(j+1,j) .GE. 0) then
          sgn = 1.0
       else
          sgn = -1.0
       end if
       s = sgn*norm2(A(j+1:size(A,1),j))  !---- compute signed norm                                                                    
       
       v(j+1,j) = A(j+1,j) + s    !---- compute v(j+1,j), different from following row values                                         

       do i = j+2, size(A,1)
          v(i,j) = A(i,j)     !---- row values below v(j+1,j) equal a(j+1,j)...a(m,j)                                                  
       end do
       v(:,j) = v(:,j)/norm2(v(:,j))  !---- normalized Householder vector                                                           
       
       vtA = matmul(transpose(v(:,j:j)),A) !---- calculate (v^T)A                                                                   
       A = A - 2.0*matmul(v(:,j:j),vtA)  !---- update A from left -> A = A -2*v*((v^T)*A)                                            
       
       Av = matmul(A,v(:,j:j))  !---- calculate A(v)
       A = A - 2.0*matmul(Av,transpose(v(:,j:j)))  !---- update A from right -> A = A - 2*(A*v)*v^T
       
    end do

  end subroutine houseTri
    
  subroutine inv_Iter(A,ews,x,eye)
    !subroutine: inv_Iter
    !purpose: Inverse iteration algorithm to calculate eigenvector(approx) given eigenvalues(approx)
    !inputs/outputs: A - m x m symmetric matrix
               !ews - a given approx eigenvalue
               !x - m x 1 column vector of guesses for eigenvectors, on output returns the approx values
               !eye - m x m identity matrix
    implicit none
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: ews
    real(dp), intent(inout) :: x(:,:)
    real(dp), intent(inout) :: eye(:,:)
    logical :: sTest
    integer(dp) :: s(size(A,1))
    real(dp) :: B(size(A,1),size(A,2))
    real(dp) :: y(size(x),1)
    real(dp) :: mu, tol
    real(dp) :: r(size(x,1),1)

    integer :: l
    
    tol = 1.0d-12
    sTest = .false.
    y = 1
    r = 20
    l = 1

    mu = ews
    s = 0

    B = (A - mu*eye)

    call gausElim_pp(B,x,shape(B),shape(x),sTest) !---- factorize (A-mu*eye) outside of while loop
    do while (norm2(r(:,1)) > tol)
       call backSub(B,x,y)  !---- Perform backsubstitution to solve for y = (A-mu*eye)^{-1}(x) 
       y = y / norm2(y)
       r = y - x
       x = y
       l = l + 1
    end do
   ! print*, "norm2 of R",norm2(r(:,1))       
! end do

  end subroutine inv_Iter

  subroutine QR_ews_noShift(A,Q,eye,iters)
    !subroutine: QR_ews_noShift
    !purpose: Calculate eigenvalues of matrix A using QR algorithm WITHOUT SHIFTS
    !inputs/outputs: A -- real symmetric m x m matrix. On output, A holds the eigenvalues in diagonal
           !Q -- real symmetric m x m matrix Q used to compute QR factorization
           !eye -- m x m identity matrix
           !iters -- integer to hold total number of iterations on output

    implicit none
    real(dp), intent(inout) :: A(:,:),Q(:,:),eye(:,:)
    integer(dp), intent(inout) :: iters
    real(dp) :: R(size(A,1),size(A,2))
    real(dp) :: ew_n, ew_diff
    real(dp) :: tol
    integer :: i

    ew_n = 0
    ew_diff = 1
    
    !---- run one iteration outside of loop to establish variables
    call QR_houseH(A,Q,eye)
    R = A

    A = matmul(R,Q)
    ew_n = A(size(A,1)-1,size(A,2)-1)

    tol = 1.0d-12
    iters = 1
    do while (abs(ew_diff) > tol) !continue iteration steps until error small enough
       call QR_houseH(A,Q,eye)    !compute QR factorization of A
       R = A   !call this QR factorization of A new variable R

       A = matmul(R,Q)  !perform matrix multiplication, A = RQ
       ew_diff = abs(ew_n - A(size(A,1)-1,size(A,2)-1)) !difference between eigenvalues at each iteration
                                    !determines when while loop terminates
       ew_n = A(size(A,1)-1,size(A,2)-1)     !update ew_n
       iters = iters + 1 !update iteration count
    end do

  end subroutine QR_ews_noShift

  subroutine QR_ews_Shift(A,Q,eye,iters)
    !subroutine: QR_ews_Shift                                                                             
    !purpose: Calculate eigenvalues of matrix A using QR algorithm WITH SHIFTS
    !inputs/outputs: A -- real symmetric m x m matrix. On output, A holds the eigenvalues in diagonal      
           !Q -- real symmetric m x m matrix Q used to compute QR factorization                            
           !eye -- m x m identity matrix                                                                    
           !iters -- integer to hold total number of iterations on output                                    
    
    implicit none
    real(dp), intent(inout) :: A(:,:),Q(:,:),eye(:,:)
    integer(dp), intent(inout) :: iters
    real(dp) :: R(size(A,1),size(A,2)),Amu(size(A,1),size(A,2))
    real(dp) :: A_sub(size(A,1)-1,size(A,2)-1), Amu_sub(size(A,1)-1,size(A,2)-1)
    real(dp) :: eye_sub(size(eye,1)-1,size(eye,2)-1), Q_sub(size(Q,1)-1,size(Q,2)-1)
    real(dp) :: R_sub(size(eye,1)-1,size(eye,2)-1)
    real(dp) :: ew_n, ew_diff
    real(dp) :: tol, mu
    integer :: i,j

    mu = 0
    ew_n = 0
    ew_diff = 1

    !---- run one iteration outside of loop to establish variables
    
    mu = A(size(A,1),size(A,2))
    Amu = A - mu*eye
    call QR_houseH(Amu,Q,eye)
    R = Amu

    A = matmul(R,Q)
    A = A + mu*eye
    ew_n = A(size(A,1),size(A,2))
    
    tol = 1.0d-12
    iters = 1
    do while (abs(ew_diff) > tol)
       mu = A(size(A,1),size(A,2))
       Amu = A - mu*eye
       call QR_houseH(Amu,Q,eye)
       R = Amu
       
       A = matmul(R,Q)
       A = A + mu*eye
      
       ew_diff = abs(ew_n - A(size(A,1),size(A,2)))
       ew_n = A(size(A,1),size(A,2))
       iters = iters + 1
    end do
    
    !smallest eigenvalue is calculated to full precision
    !now reduce to m-1 x m-1 submatrix to finish calculating rest of eigenvalues
    !can also calculate the final subdiagonal entry

    mu = 0
    ew_n = 0
    ew_diff = 1

    !form submatrix
    A_sub = 0
    do i = 1, size(A_sub,1)
       do j = 1, size(A_sub,2)
          A_sub(i,j) = A(i,j)
       end do
    end do
    do i = 1, size(eye_sub,1)
       do j = 1, size(eye_sub,2)
          eye_sub(i,j) = eye(i,j)
       end do
    end do
    Q_sub = 0
    R_sub = 0

    mu = A_sub(size(A_sub,1),size(A_sub,2))
    Amu_sub = A_sub - mu*eye_sub
    call QR_houseH(Amu_sub,Q_sub,eye_sub)
    R_sub = Amu_sub
    
    A_sub = matmul(R_sub,Q_sub)
    A_sub = A_sub + mu*eye_sub
    ew_n = A_sub(size(A_sub,1),size(A_sub,2))

    do while (abs(ew_diff) > tol)
       mu = A_sub(size(A_sub,1),size(A_sub,2))
       Amu_sub = A_sub - mu*eye_sub
       call QR_houseH(Amu_sub,Q_sub,eye_sub)
       R_sub = Amu_sub

       A_sub = matmul(R_sub,Q_sub)
       A_sub = A_sub + mu*eye_sub

       ew_diff = abs(ew_n - A_sub(size(A_sub,1),size(A_sub,2)))
       ew_n = A_sub(size(A_sub,1),size(A_sub,2))
       iters = iters + 1
    end do

    do i = 1, size(A_sub,1)
       do j = 1, size(A_sub,1)
          A(i,j) = A_sub(i,j)
       end do
    end do


  end subroutine QR_ews_Shift


  subroutine g_jacobi(A,b,x,tol,count,upper_tol)
    !subroutine: g_jacobi                                                                                                  
    !purpose: solves the equation Ax = b using the Gauss-Jacobi algorithm                                                 
    !inputs/outputs: A -- m x m input matrix written as A = R+D                                                           
    !      b -- m long vector                                                                                             
    !      x -- m long vector holding the solution upon output                                                            
    !      tol -- specified tolerance to determine successful program termination                                         
    !      count -- counts number of iterations until convergence
    !      upper_tol -- upper limit on error, ||b - Ax||_2 specified by user

    implicit none

    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: b(:), x(:)
    real(dp), intent(in) :: tol
    integer, intent(inout) :: count
    real(dp), intent(inout) :: upper_tol

    real(dp) :: R(size(A,1),size(A,2)), D(size(x))
    real(dp) :: y(size(x)), error(size(x))
  
    character(len=40) :: fileName
    character(len=5) :: num
    integer :: i
    
    y = 0

    R = A
    do i = 1, size(A,1)
       R(i,i) = 0.d0     !form the R matrix, which is A with zeros on diagonal
       D(i) = A(i,i)     !form D vector, which is 1 / all diagonal elements of A
    end do

    write(num,101) int(D(1)) !prepare file to write error to
    101 format(I5.5)
    fileName = 'errorJ_1' // num // '.dat'
    error = 1
    count = 1

    do while (norm2(error) > tol) !run loop until desired accuracy is achieved                                           
       y = (b - matmul(R,x))         !y = b - Rx
       error = (b - matmul(A,x))     !loop termination criteria: ||b - Ax||_2 < tol
       open(1,file = fileName)       !write error to file
       write(1,*) norm2(error)
      
       y = y/D  !algorithm is that we replace x with D^-1(b - Rx)
       x = y    !repace old x with new
       count = count + 1 !update iteration count

       if (norm2(error) > upper_tol) then
          print 600, upper_tol
          600 format('Error upper limit reached: ', e10.4)
          print*, "Algorithm divergence..."
          exit                                                                                                              
       end if

    end do
    count = count - 1 !account for final iteration which terminated loop

    close(1)

  end subroutine g_jacobi




  subroutine g_seidel(A,b,x,tol,count,upper_tol)
    !subroutine: g_seidel
    !purpose: Solves the equation Ax = b using the Gauss-Seidel algorithm
    !inputs/outputs: A -- m x m input matrix
    !      b -- m long vector
    !      x -- m long vector holding the solution on output
    !      tol -- specified tolerance to determine successful program termination  
    !      count -- counts number of iterations until convergence
    !      upper_tol -- upper limit on error, ||b - Ax||_2 specified by user   

    implicit none

    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: b(:), x(:)
    real(dp), intent(in) :: tol
    integer, intent(inout) :: count
    real(dp), intent(inout) :: upper_tol

    real(dp) :: L(size(A,1),size(A,2))
    real(dp) :: U(size(A,1),size(A,2))
    real(dp) :: D(size(A,1),size(A,2))
    real(dp) :: DL(size(A,1),size(A,2))
    real(dp) :: bUx(size(x)), error(size(x))
    real(dp) :: y(size(x))
    real(dp) :: summ
 
    character(len=40) :: fileName
    character(len=5) :: num

    integer :: i,j
    
    D = 0
    L = 0
    U = 0
    DL = 0
    bUx = 0

    do i = 1, size(A,1) !forming U, L, and D matrices
       D(i,i) = A(i,i)
       do j = 1, size(A,2)
          if (j > i) then
             U(i,j) = A(i,j)
          else if (j < i) then
             L(i,j) = A(i,j)
          
          end if
       end do
    end do
    
    DL = D + L             !form new variable (D + L)
    bUx = b - matmul(U,x)  !form new variable (b - Ux)

!    call matPrint(DL,shape(DL))
!    call matPrint(U,shape(U))


    write(num,101) int(D(1,1))  !preparing file to write error to
    101 format(I5.5)
    fileName = 'errorS_1' // num // '.dat'
    
    error = 1
    summ = 0
    y = bUx
    count = 1

    do while (norm2(error) > tol)
       
       !perform a forward substitution to solve the lower triangular system (D+L)*y = (b-Ux) for y
       do i = 1, size(x,1)
          summ = bUx(i)
          do j = 1, i-1
             summ = summ - y(j)*DL(i,j)
          end do
          y(i) = summ/DL(i,i)
       end do

       !calculate error, write error to file, replace new x with old
       error = (b - matmul(A,x))
       open(1,file = fileName)
       write(1,*) norm2(error)

       x = y
       bUx = b - matmul(U,x)   !reform (b - Ux) since x updated
       count = count + 1 !update iteration count

       if (norm2(error) > upper_tol) then
          print 600, upper_tol
          600 format('Error upper limit reached: ', e10.4)
          print*, "Algorithm divergence..."
          exit
       end if

    end do

    count = count - 1 !account for final iteration which terminated loop   

    close(1)


  end subroutine g_seidel

  subroutine CG_smart(A,b,x,tol,count)
    !subroutine: CG_smart
    !purpose: Smart conjugate gradient algorithm without pre-conditioning to solve Ax = b
    !inputs/outputs: A -- m x m input matrix
    !      b -- m long input vector
    !      x -- m long solution output vector
    !      tol -- desired accuracy of algorithm
    !      count -- stores total number of iterations upon output

    implicit none
 
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: b(:), x(:)
    real(dp), intent(in) :: tol
    integer, intent(inout) :: count
    real(dp) :: r(size(x))
    real(dp) :: p(size(x))
    real(dp) :: y(size(x))
    real(dp) :: alpha, beta, ee, ee_new
    real(dp) :: sym_tol

    !Check if A is symmetric
    sym_tol = 1.0d-12

    if (norm2(A - transpose(A)) > sym_tol) then
       print*, "Matrix A is not symmetric!"
       stop
    end if

    r = b - matmul(A,x)
    p = r
    ee = (norm2(r))**2
    count = 1
    
    do while (sqrt(ee) > tol)
       y = matmul(A,p)    !calculate temporary vector to save time
       
       alpha = ee/dot_product(p,y)
       x = x + alpha*p
       r = r - alpha*y  !this formula for r uses y, avoids calculating Ax
       ee_new = (norm2(r))**2
       beta = ee_new/ee
       p = r + beta*p
       ee = ee_new
       count = count + 1
    end do

    count = count - 1 !account for final iteration which terminated loop   

  end subroutine CG_smart

  subroutine CG_precond(A,b,x,M,tol,count)
    !subroutine: CG_precond
    !purpose: Conjugate gradient algorithm with preconditioning to solve Ax = b
    !inputs/outputs: A -- m x m input matrix
    !      b -- m long input vector
    !      x -- m long solution output vector
    !      M -- preconditioner matrix
    !      tol -- desired accuracy of algorithm
    !      count -- stores total number of iterations upon output
    
    implicit none
    
    real(dp), intent(inout) :: A(:,:), M(:,:)
    real(dp), intent(inout) :: b(:), x(:)
    real(dp), intent(in) :: tol
    integer, intent(inout) :: count
    real(dp) :: r(size(x))
    real(dp) :: p(size(x))
    real(dp) :: y(size(x)), z(size(x))
    integer(dp) :: s(size(A,1))
    real(dp) :: A_t(size(A,1),size(A,2))
    real(dp) :: alpha, beta, ee, ee_new
    real(dp) :: sym_tol
    Logical :: sTest
    real(dp) :: r2(size(x,1),1)
    real(dp) :: z2(size(x,1),1)
    real(dp) :: m_vect(size(A,1))
    integer :: i
    
    r2 = 0
    z2 = 1
    s = 0
    sTest = .false.
    
    sym_tol = 1.0d-12
    if (norm2(A - transpose(A)) > sym_tol) then
       print*, "Matrix A is not symmetric!"
       stop
    end if

    do i = 1,size(M,1)
       m_vect(i) = M(i,i)
    end do

    r = b - matmul(A,x)

!    call gausElim_pp(M,r2,shape(M),shape(r2),sTest)
!    call backSub(M,r2,z2)  !z = M^-1(r), M diagonal matrix                                                    
!    call matPrint(z2,shape(z2))
!    r = r2(:,1)
!    z = z2(:,1)

    z = r/m_vect  !z = M^-1(r), since we know M is diagonal matrix
    p = z
    ee = dot_product(r,z)
    
    count = 1
    do while (sqrt(ee) > tol)
       y = matmul(A,p)
       alpha = ee/dot_product(p,y)
       x = x + alpha*p
       r = r - alpha*y
       
      ! do i = 1,size(r)
      !    r2(i,1) = r(i)
      !    z2(i,1) = z(i)
      ! end do
       
      ! call LU_backSub(M,size(M,1),r2,s,z2)  !z = M^-1(r), 
       z = r/m_vect    !M diagonal matrix, already decomposed from routine call above
      ! r = r2(:,1)
      ! z = z2(:,1)

       ee_new = dot_product(r,z)
       beta = ee_new/ee
       p = z + beta*p
       ee = ee_new
       count = count + 1
    end do
      
    count = count - 1
 
  end subroutine CG_precond

  subroutine readMat(filename,mat)

    implicit none
    character(len=*) :: filename
    real(dp) :: mat(:,:)
    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! Note that entries must be separated by a tab.

    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) i,j

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)
    
  end subroutine readMat




end module LinAl
