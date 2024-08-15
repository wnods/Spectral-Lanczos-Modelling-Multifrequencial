Module sparsemul

! Este programa carrega duas matrizes esparsas, e as converte em matrizes
! no formato  'compressed row format' também conhecido como Harwell-Boeing format,

! As[]   -> armazena os valores não-nulos da matriz A
! asub[] -> armazena os indices da linha a qual os valores pertencem
! xa[]   -> armazena os indices que indicam o inicio de cada coluna
!           considerando o posicionamento desses valores no vetor As[]


! Em seguida realiza a multiplicação entre elas.

!na,ma,nnza --> rows A, columns A, noz zeroes A

! Autor Diego C. Miranda     28 jan 2020

Implicit None
Save
Contains

subroutine matmul_sparse(matA,na,ma,matb,nb,mb,matC)
Implicit None

! matA - matriz A cheia
! na   - numero de linhas de A
! ma   - numero de colunas de A
! matB - matriz ou vetor B
! nb   - numero de linhas de B
! mb   - numero de colunas de B
!-----------------------------------------------------

Integer,Intent(in)  :: na,ma,nb,mb
Real,Intent(in)     :: matA(na,ma),matB(nb,mb)

Integer             :: nnza,nnzb
Real,Allocatable    :: valA(:),valB(:),c(:)
Integer,Allocatable :: rowA(:),rowB(:),colPtrA(:),colPtrB(:),jc(:),ic(:)
Integer,Allocatable :: iwa(:),iwb(:),ndegr(:)
Integer		    :: job,nzmax
Integer	            :: flag, ierr

Integer             :: i,j,k

Real,Intent(out)    :: matC(na,mb)


nnza = count(matA /= 0)
nnzb = count(matB /= 0)


!compressed row storage -----------
Allocate(valA(nnza),valB(nnzb),rowA(nnza),rowB(nnzb),colPtrA(na+1),colPtrB(nb+1))

!convertendo Matriz A para crs
colPtrA = 0
k = 1
Do i = 1,na
	flag = 0
	Do j = 1,ma
		If (matA(i,j) /= 0) Then
			valA(k) = matA(i,j)
			rowA(k) = j     !   linha variando mais rapido
			If (flag == 0) Then
				colPtrA(i) = k				
				flag = 1
			End If
			k = k + 1
		End If
	End Do
End Do
colPtrA(na+1) = nnza + 1

Do i = 1,na
	If (colPtrA(i) == 0) colPtrA(i) = colPtrA(i+1)
End Do

!convertendo Matriz B para crs
colPtrB = 0
k = 1
Do i = 1,nb
	flag = 0
	Do j = 1,mb
		If (matB(i,j) /= 0) Then
			valB(k) = matB(i,j)
			rowB(k) = j     !   linha variando mais rapido
			If (flag == 0) Then
				colPtrB(i) = k				
				flag = 1
			End If
			k = k + 1
		End If
	End Do
End Do
colPtrB(nb+1) = nnzb + 1

Do i = 1,nb
	If (colPtrb(i) == 0) colPtrB(i) = colPtrB(i+1)
End Do


Allocate(ndegr(na),iwb(mb),iwa(ma))

 Call amubdg(na,ma,mb,rowA,colPtrA,rowB,colPtrB,ndegr,nzmax,iwb,nnza,nnzb)  

Allocate(c(nzmax),jc(nzmax),ic(na+1))   !ic(ma+1))  !
job = 1

! multiplicacao matricial esparsa
 
 Call amub(na, ma, job, valA, rowA, colPtrA, valB, rowB, colPtrB, c, jc, ic, nzmax,  iwa, ierr,nnza,nnzb)
!Call cpu_time(finish)


 Call csrdns ( na, mb, c, jc, ic, nzmax, matC, ierr )



end subroutine
!=========================================================



subroutine amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr ,nnza,nnzb)

!*****************************************************************************80
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, job indicator.  When JOB = 0, only the
!    structure is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer ( kind = 4 ) NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer ( kind = 4 ) work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ::  ncol
  integer :: nrow
  integer :: nzmax

  integer :: nnza,nnzb

  real::  a(nnza)
  real :: b(nnzb)
  real :: c(nzmax)
  integer :: ia(nrow+1)
  integer :: ib(ncol+1)
  integer :: ic(nrow+1)!ic(ncol+1)!
  integer::  ierr
  integer :: ii
  integer :: iw(ncol)
  integer :: ja(nnza)
  integer :: jb(nnzb)
  integer :: jc(nzmax)
  integer :: jcol
  integer :: jj
  integer :: job
  integer::  jpos
  integer::  k
  integer::  ka
  integer::  kb
  integer::  len
  real::  scal
  logical:: values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              if ( values ) then
                c(jpos) = c(jpos) + scal * b(kb)
              end if
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
end subroutine

!============================================================================================

subroutine amubdg ( nrow, ncol, ncolb, ja, ia, jb, ib, ndegr, nnz, iw, nnza, nnzb )

!*****************************************************************************80
!
!! AMUBDG gets the number of nonzero elements in each row of A * B.
!
!  Discussion:
!
!    The routine also computes the total number of nonzero elements in A * B.
!
!    Method: A' * A = sum [over i = 1, nrow]  a(i)^T a(i)
!    where a(i) = i-th row of  A.  We must be careful not to add  the
!    elements already accounted for.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix A.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix A,
!    (and the row dimension of B).
!
!    Input, integer ( kind = 4 ) NCOLB, the column dimension of the matrix B.
!
!    Input, ja, ia= row structure of input matrix A: ja = column indices of
!    the nonzero elements of A stored by rows.
!    ia = pointer to beginning of each row in ja.
!
!    Input, jb, ib, the row structure of input matrix B: jb = column indices of
!    the nonzero elements of A stored by rows.
!    ib is a pointer to beginning of each row in jb.
!
!    Output, integer ( kind = 4 ) NDEGR(NROW), contains the degrees (the number
!    of nonzeros in each row of the matrix A * B.
!
!    Output, integer ( kind = 4 ) NNZ, the number of nonzero elements 
!    found in A * B.
!
!    Workspace, integer ( kind = 4 ) IW(NCOLB).
!
  implicit none

  integer :: ncol
  integer :: ncolb
  integer :: nrow

  integer :: nnza,nnzb


  integer :: ia(nrow+1)
  integer :: ib(ncol+1)
  integer :: ii
  integer :: iw(ncolb)
  integer :: j
  integer :: ja(nnza)
  integer :: jb(nnzb)
  integer :: jc
  integer :: jr
  integer :: k
  integer :: last
  integer :: ldg
  integer :: ndegr(nrow)
  integer :: nnz

  iw(1:ncolb) = 0
  ndegr(1:nrow) = 0

  do ii = 1, nrow
!
!  For each row of A.
!
    ldg = 0
!
!  End-of-linked list.
!
    last = -1

    do j = ia(ii), ia(ii+1)-1
!
!  Row number to be added.
!
        jr = ja(j)

        do k = ib(jr), ib(jr+1)-1
           jc = jb(k)
!
!  Add one element to the linked list.
!
           if ( iw(jc) == 0 ) then
              ldg = ldg + 1
              iw(jc) = last
              last = jc
           end if

         end do

    end do

    ndegr(ii) = ldg
!
!  Reset IW to zero.
!
    do k = 1, ldg
      j = iw(last)
      iw(last) = 0
      last = j
     end do

  end do

  nnz = sum ( ndegr(1:nrow) )

  return
end subroutine




subroutine csrdns ( nrow, ncol, a, ja, ia, nzmax, dns, ierr )

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the dimension of the DNS array.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  implicit none

  integer :: ncol
  integer :: nrow
  integer :: nzmax

  real :: a(nzmax)
  real  ::dns(nrow,ncol)
  integer :: i
  integer :: ia(nrow+1)
  integer :: ierr
  integer :: j
  integer :: ja(nzmax)
  integer :: k
  
  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end
















end module







































