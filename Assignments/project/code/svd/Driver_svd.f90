!Main program. To compile and run, use command make LinAl, and then ./LinAl which will print
!all results to the screen including the specified singular values after the SVD has been complete

!This program reads in a .dat image as input, and uses LAPACK's dgesvd() to compute the SVD of the image
!Compressed images are also computed using the specified amount of singular values, k.


Program Driver_svd

  use utility, only: dp

  implicit none

  real(dp), allocatable :: A(:,:),A_s(:,:),A_im(:,:),A_small(:,:)
  real(dp), allocatable :: U(:,:),VT(:,:)
  real(dp), allocatable :: S(:),WORK(:)
  integer, allocatable :: ks(:)

  character(len=10) :: JOBU
  character(len=10) :: JOBVT
  integer :: m,n,lda,ldu,ldvt,lwork,lworkmax,info

  character(len=40) :: Data, fileName
  character(len=5) :: num
  
  integer, dimension (1:2) :: order1 = (/1, 2/)
  integer, dimension (1:2) :: order2 = (/2, 1/)  

  integer :: i,j


  JOBU = 'A' !specifies options for computing all or part of matrix U
             !A = all M columns returned in U
             !S = first min(m,n) columns of U returned in U
             !O = first min(m,n) columns of U overwritten in A
             !N = no columns of U (no left singular vectors) computed

  JOBVT = 'A' !specifies options for computing all or part of matrix transpose(V)
              !A = all
              !S = first min(m,n) rows of V^T returned in array VT
              !O = first min(m,n) rows of V^T overwritten on array A
              !N = no rows of V^T (no right singular vectors) are computed

  !JOBU and JOBVT cannot both be 'O'

  m = 3355    ! number of rows in matrix A
  n = 5295    ! number of columns in matrix A

!  m = 6
!  n = 5


  lda = m  !leading dimension of array A
  ldu = m  !leading dimension of U
  ldvt = n !leading dimension of VT
  lworkmax = 20000 !dimension of WORK array
  
  
  allocate(A(m,n))
  allocate(A_s(size(A,1),size(A,2)))
  A = 0
  A_s = 0
  ! allocate(A_small(lda,n))
 


  !-------- Read in dog_bw_data.dat file as input --------                                                                                              
  Data = 'dog_bw_data.dat'
                                                                                                                                      
  open(10,file=Data)

  do i = 1, size(A,1)
     read(10,*) A(i,:)
  end do

  close(10)

  print*, 'shape of data matrix A:', shape(A)
  A_s = A !store saved copy of data matrix


 ! A_small(:,:) = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/),order=order2)
  
!  A_small(:,:) = reshape((/8.79,9.93,9.83,5.45,3.16, 6.11,6.91,5.04,-0.27,7.98,&
!                         -9.15,-7.93,4.86,4.85,3.01, 9.57,1.64,8.83,0.74,5.80, &
!                         -3.49,4.02,9.80,10.0,4.27, 9.84,0.15,-8.99,-6.02,-5.31/),(/6,5/), order = order2)
!  call matPrint(A_small,shape(A_small))

  allocate(S(min(m,n)))        !singular values of A, sorted so that S(i) >= S(i+1)
  allocate(U(ldu,m))    !shape(ldu,m) if JOBU = 'A', or shape(ldu,min(m,n)) if JOBU = 'S'
  allocate(VT(ldvt,n))  !shape(n,n) if JOBVT = 'A', shape(min(m,n),n) if JOBVT = 'S'

  lwork = 30000
  allocate(WORK(lwork))

!  S = 0
!  U = 0
!  VT = 0
!  WORK = 0

!  -------- we can first query the optimal workspace if needed ------
!  -------- to query, allocate WORK(lworkmax), and use lwork=-1 ----
!  lwork = -1


 ! call dgesvd(JOBU,JOBVT,m,n,A,lda,S,U,ldu,VT,ldvt,WORK,lwork,info)
  
 ! lwork = nint(WORK(1))
 ! print*, 'new lwork =',lwork

 ! deallocate(WORK)
 ! allocate(WORK(lwork))

! -------- above code to query workspace --------
! -------- below we now form SVD of matrix A --------

  call dgesvd(JOBU,JOBVT,m,n,A,lda,S,U,ldu,VT,ldvt,WORK,lwork,info)

  print*, "Info is:", info

!  print*, "Singular values:", S
!  print*, "U matrix:", U

!  S = 0
!  U = 0
!  VT = 0
!  WORK = 0



  !-----computing compressed images -----

  allocate(ks(9))
  
  ks = (/20,40,80,160,320,640,1280,2560,3355/)
!  ks = (/20/)
  allocate(A_im(size(A,1),size(A,2)))
  A_im = 0

  !first loop to add rank 1 updates to existing matrix A_im
  do i = 1, ks(1)
     A_im = A_im + S(i)*matmul(U(:,i:i),VT(i:i,:))
  end do

!------- write compressed matrix data to a data file --------
  write(num,100) ks(1)
  100 format(I5.5)
  fileName = 'Image_appn_1' // num // '.dat'
  open(1, file = fileName)
  do i = 1, size(A_im,1)
     write(1,*) A_im(i,:)
  end do
  close(1)

!----- write error value to data file ------
  open(5,file = 'svd_errors.dat')
  write(5,*) (norm2(A_s - A_im)/(m*n))


  !now loop over remaining values of k
  do j = 2, size(ks)
     do i = ks(j-1)+1,ks(j)
        A_im = A_im + S(i)*matmul(U(:,i:i),VT(i:i,:))
     end do
     
     write(num,101) ks(j) 
     101  format(I5.5)
     fileName = 'Image_appn_1' // num // '.dat'
     open(1, file = fileName)
     do i = 1, size(A_im,1)
        write(1,*) A_im(i,:)
     end do
     close(1)

     !write error value to data file
     write(5,*) (norm2(A_s - A_im)/(m*n))
  end do

  close(5)
  
  !---- display to screen the first 10 singular values and singular values for each k, S(k) ----
  print*,""
  print*,"First 10 singular values:" 
  do i = 1,10
     print*,'sigma',i,'=',S(i)
  end do

  !---- display to screen rest singular values for ks array ----
  print*,""
  print*,"Rest singular values:"
  do i = 1, size(ks)
     print*,'sigma',ks(i),'=',S(ks(i))
  end do
  


  
  deallocate(A_im)
  deallocate(A)
!  deallocate(A_small)
  deallocate(S)
  deallocate(U)
  deallocate(VT)
  deallocate(WORK)



end program Driver_svd
