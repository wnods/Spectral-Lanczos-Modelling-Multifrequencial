!include 'mkl_pardiso.f90'
include '/opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/include/mkl_pardiso.f90'
module mod_mtTM
!--
!    Modulo com as subrotinas de Elementos Finitos 2D para o metodo MT.
!Autor: Anderson Almeida. 19/05/2020
!--
use variaveis_globais
use arrayspds
!--
contains        

subroutine montagem_matriz_vetorfonte_globais(n,nnz,G,F)
implicit none
integer,intent(in):: n,nnz
complex(db),intent(out) :: G(nnz), F(n) 

integer :: i,j,k,ie,m
real(db),allocatable :: a(:),b(:),c(:)
real(db),allocatable :: x(:),z(:)
integer,allocatable :: inel(:)
real(db),allocatable :: mass(:,:),stiff(:,:)
complex(db),allocatable :: Ke(:,:), Fe(:), Ex(:)
real(db),allocatable :: Kema(:,:), Kest(:,:), massa(:,:), rigidez(:,:)
real(db):: Ae
real(db):: etae, sige, w, delet, etap
complex(db):: zetae, Hy

allocate(a(3),b(3),c(3))
allocate(x(3),z(3))
allocate(inel(3))
allocate(mass(3,3),stiff(3,3))
allocate(Ke(3,3), Fe(3), Ex(3))
allocate(Kema(3,3), Kest(3,3))
allocate(massa(3,3),rigidez(3,3))


G = 0.d0
F = 0.d0

w = 2.d0*pi*freq

do i = 1,nelm

   inel = mel(i,:)
   x = mnod(inel(:),1)
   z = mnod(inel(:),2)

   a(1) = x(2)*z(3) - x(3)*z(2) 
   a(2) = x(3)*z(1) - x(1)*z(3) 
   a(3) = x(1)*z(2) - x(2)*z(1)
     
   b(1) = z(2) - z(3)
   b(2) = z(3) - z(1)
   b(3) = z(1) - z(2)

   c(1) = x(3) - x(2)
   c(2) = x(1) - x(3)
   c(3) = x(2) - x(1)
   
   Ae = 0.5d0*sum(a)
   sige  = (1.d0/prop(i))

   etae  = sige   
   zetae = im*w*muo

   mass(1,1) = b(1)**2 + c(1)**2; mass(1,2) = b(1)*b(2)+c(1)*c(2) ; mass(1,3) = b(1)*b(3)+c(1)*c(3)
   mass(2,1) = mass(1,2); mass(2,2) = b(2)**2 + c(2)**2 ; mass(2,3) = b(2)*b(3)+c(2)*c(3)
   mass(3,1) = mass(1,3) ; mass(3,2) = mass(2,3) ; mass(3,3) = b(3)**2 + c(3)**2

   mass = (1.d0/(4.d0*Ae)) * mass

   rigidez = mass;

   stiff(1,1) = 2.d0; stiff(1,2) = 1.d0 ; stiff(1,3) = 1.d0
   stiff(2,1) = stiff(1,2); stiff(2,2) = 2.d0 ; stiff(2,3) = 1.d0
   stiff(3,1) = stiff(1,3) ; stiff(3,2) = stiff(2,3) ; stiff(3,3) = 2.d0

   stiff = (Ae/12.d0) * stiff
   
   massa = stiff;

   Ke = (rigidez/etae) + zetae*massa  ! stiff + mass...

   Kest = rigidez/etae
   Kema = massa

   etap  = 1.d0/rhob 
   delet = etae - etap

   if( dabs(delet) > 1.d-7 )then
       do j = 1,3
          call primario(z(j),etap,w,Hy,Ex(j))
       end do
       Fe = -(delet/(6.d0*etae))*sum(Ex)*c; 
   else    
       Fe = 0.d0
   end if
   
   call set_arrays_pds(cdim,ngl,nnod,nnz,ja,inel,Ke,ia,G)

   call set_arrays_pds_r(cdim,ngl,nnod,nnz,ja,inel,Kest,ia,mstiff)
   call set_arrays_pds_r(cdim,ngl,nnod,nnz,ja,inel,Kema,ia,mmass)

   do j = 1,3
      F(inel(j)) = F(inel(j)) + Fe(j)
   end do

!-------------
   do m = 1,3
      do k = 1,3
         bmass(inel(m),inel(k)) = bmass(inel(m),inel(k)) + massa(m,k)
         bstiff(inel(m),inel(k)) = bstiff(inel(m),inel(k)) + (rigidez(m,k)/etae)
      end do 
   end do
!-------------

end do

end subroutine

!--

subroutine primario(z,sigmp,w,Hy,Ex)
!-
implicit none
   
REAL(8),INTENT(IN)     :: Z,SIGMP,w
COMPLEX(8),INTENT(OUT) :: Ex,Hy
   
real(8),parameter      :: pi  = 3.141592653589793d0 
real(8),parameter      :: muo = 4.d0*pi*1.d-7, eo = 8.85418*1.d-12
complex(8),parameter   :: ic  = (0.d0,1.d0) 
COMPLEX(8)             :: k_1,u_1,ZET_0,ZET_1,H_0,H_1,Rtm_0
   
   k_1   = dsqrt(w*muo*SIGMP/2.D0)*( 1.d0 - ic )
   u_1   = ic*k_1
   
   ZET_0 = DSQRT(muo/eo)
   ZET_1 = u_1/SIGMP
   
   Rtm_0 = (ZET_0 - ZET_1)/(ZET_0 + ZET_1)
   
   H_0   = 1.d0/(1.D0 + Rtm_0)
   
   H_1   = H_0*( 1.d0 + Rtm_0 )
   
   Hy    = H_1*CDEXP(-u_1*z)
   Ex    = ZET_1*Hy
   
end subroutine
!--

subroutine solver_pardiso ( ngl , nglob , nnz , Vnn , ia , ja , vfonte )

!   use mkl_pardiso
   
   implicit none
   integer,intent(in)     :: ngl,nglob,nnz
   integer,intent(in)     :: ia(ngl*nglob+1) , ja(nnz)
   complex(db),intent(in)    :: Vnn(nnz)
   complex(db),intent(inout) :: vfonte(ngl*nglob)

   integer(8), ALLOCATABLE  :: pt(:)   
!   TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
   INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
   INTEGER error1
   INTEGER, ALLOCATABLE :: iparm( : )
   INTEGER i, idum(1)
   COMPLEX(KIND=db) ddum(1)
   
   complex(db),allocatable :: x(:)
   Real(db) :: t1,t2
   
   !%%%%%%% Pardiso %%%%%%
   
   allocate(x(ngl*nglob))
   
   ! Initialize the solver.
   
   call cpu_time(t1)
   
   !%%%%
   
   nrhs = 1 
   maxfct = 1 
   mnum = 1
   
   ALLOCATE( iparm ( 64 ) )
   
   do i = 1, 64
      iparm(i) = 0
   end do 
   
   iparm(1) = 1 ! no solver default
   iparm(2) = 2 ! fill-in reordering from METIS
   iparm(4) = 0 ! no iterative-direct algorithm
   iparm(5) = 0 ! no user fill-in reducing permutation
   iparm(6) = 0 ! =0 solution on the first n compoments of x
   iparm(8) = 9 ! numbers of iterative refinement steps
   iparm(10) = 13 ! perturbe the pivot elements with 1E-13
   iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
   iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
   iparm(14) = 0 ! Output: number of perturbed pivots
   iparm(18) = -1 ! Output: number of nonzeros in the factor LU
   iparm(19) = -1 ! Output: Mflops for LU factorization
   iparm(20) = 0 ! Output: Numbers of CG Iterations
   
   error  = 0 ! initialize error flag
   msglvl = 0 ! print statistical information
   mtype  = 6 ! symmetric, indefinite
   
   ALLOCATE ( pt ( 64 ) )
   do i = 1, 64
      pt( i ) =  0 
!      pt( i )%DUMMY =  0 
   end do
   
   phase = 11 ! only reordering and symbolic factorization
   
   CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, vnn, ia, ja, &
   idum, nrhs, iparm, msglvl, ddum, ddum, error)    
   WRITE(*,*) 'Reordering completed ... '
   
   !.. Factorization.
   phase = 22 ! only factorization
   CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, vnn, ia, ja, &
   idum, nrhs, iparm, msglvl, ddum, ddum, error)
   WRITE(*,*) 'Factorization completed ... '
   
   !.. Back substitution and iterative refinement
   iparm(8) = 2 ! max numbers of iterative refinement steps
   phase = 33 ! only factorization
   CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, vnn, ia, ja, &
   idum, nrhs, iparm, msglvl, vfonte, x, error)
   WRITE(*,*) 'Solve completed ... '
   
   !.. Termination and release of memory
   phase = -1 ! release internal memory
   CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, ddum, idum, idum, &
   idum, nrhs, iparm, msglvl, ddum, ddum, error1)
   
   IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )
   
   !%%%%
   
   call cpu_time(t2)
   
   print*,'Tempo para resolver o sistema linear :',(t2-t1)/60.,'minutos'
   print*,'  '
   
   vfonte = x
   
   !%%%%%%%%%%%%%%%%%%%
   
   deallocate(x)
   
end subroutine

!--

subroutine derivada_MFB( ng, nel, npd, vnd, Uap, matel , v_props , matc, Grad ) 

implicit none

integer,intent(in)  :: ng , nel , npd 
integer,intent(in)  :: matel(nel,3)
real(8),intent(in)  :: matc(ng,2), v_props(nel)
complex(8),intent(in)  :: Uap(ng)
integer,intent(in)  :: vnd(npd) 
complex(8),intent(out) :: Grad(npd,2)

integer :: i,j,k,l,cont,nelv
integer,allocatable :: vel(:)
integer :: aux(20),no(3)
real(8) :: sumsig,a(3),b(3),c(3),x(3),z(3),Ael
complex(8) :: sumx,sumz,sumex,sumez
real(8) :: Sig(npd)

Grad=0.d0;

do i = 1,npd

   cont = 0
   do j = 1,nel
        do l = 1,3
            if( vnd(i).eq.int(matel(j,l)) )then
              cont = cont + 1
              aux(cont) = j
            end if
        end do  
   end do
   nelv = cont 
   allocate(vel(nelv))

   do j = 1,nelv
       vel(j) = aux(j) 
   end do

!== set sigma ====
   sumsig = 0.d0
   do j = 1,nelv
      sumsig = sumsig + ( 1.d0/( v_props(vel(j)) ) ) 
   end do
   Sig(i) = sumsig/nelv 
!=================

   sumx = 0.
   sumz = 0.
   do j = 1,nelv
        
     do l = 1,3
        x(l)  = matc( int( matel(vel(j),l) ) , 1 )
        z(l)  = matc( int( matel(vel(j),l) ) , 2 )
        no(l) = int( matel( vel(j),l ) ) 
     end do 

     a(1) = x(2)*z(3)-x(3)*z(2)
     a(2) = x(3)*z(1)-x(1)*z(3)
     a(3) = x(1)*z(2)-x(2)*z(1)

     b(1) = z(2)-z(3)
     b(2) = z(3)-z(1)
     b(3) = z(1)-z(2)
 
     c(1) = x(3)-x(2)
     c(2) = x(1)-x(3)
     c(3) = x(2)-x(1)

     Ael = dabs(dabs(sum(a))/2.d0)
     
     sumex = 0. ; sumez = 0.
     do l = 1,3
            sumex = sumex + b(l)*Uap(no(l))
            sumez = sumez + c(l)*Uap(no(l))
     end do 
     sumex = sumex/(2.*Ael)
     sumez = sumez/(2.*Ael)

     sumx = sumx + sumex
     sumz = sumz + sumez
                         
   end do

   Grad(i,1) = (sumx/nelv)/sig(i)  ! Campo Ez^s =  1/sig * dHy^s/dx
   Grad(i,2) = -(sumz/nelv)/sig(i) ! Campo Ex^s = -1/sig * dHy^s/dz
 
   deallocate(vel)

end do


end subroutine

!--

subroutine set_arrays_pds_r(cdim,ngl,nnos,nnz,ic,el,ml,npr,vl)

! Rotina que monta o vetor com os valores nao nulos da matriz global

! cdim  = Escolha da dimensao do problema 2 -> 2D ou 3 -> 3D;
! nel   = Numero de elementos da malha de EF
! nnos  = Numero de nos da malha
! ngl   = Nuero de graus de liberdade
! nnz   = Numero de valores nao nulos da matriz global, inicialmente.
! ic    = indice de coluna dos elementos nao nulos de cada linha da matriz global.
! el    = Vetor com os indices do elemento
! ml    = Matriz local
! npr   = Vetor com a informacao da quantidade de elementos nao nulos em cada linha da matriz global
! vl    = Valores nao nulos da matriz global

implicit none

integer,intent(in) :: cdim,ngl,nnos,nnz,el(:)
integer,intent(in) :: ic(nnz),npr(ngl*nnos+1)

!integer,intent(inout) :: vl(nnz),ml(:,:)
!complex(dpcs),intent(inout) :: vl(nnz),ml(:,:)
real(dpcs),intent(inout) :: vl(nnz),ml(:,:)

integer :: i,j,k,i_f,j_f,ii,jj,pos

if(cdim==2)then
    i_f=ngl*3 ; j_f=i_f
elseif(cdim==3)then
    i_f=ngl*4 ; j_f=i_f
else
    print*,'Escolha a dimensao correta !';stop
end if

do i = 1,i_f
    ii = el(i)
       do j = 1,j_f
          jj = el(j)
              if(ii.le.jj)then
                 do k = npr(ii),npr(ii+1)-1               
                    if( jj.eq.ic(k) )then
                      vl(k) = vl(k) + ml(i,j)   
                      exit
                    end if
                 end do 
              end if 
       end do
end do

end subroutine

!--

SUBROUTINE LU_Gauss_compacto2( A, b, n)                               
!====================================================================
!     Fatora A em LU e armazena em A. resolver o sistema (U*x) = b e 
!     guardar a solução em b                                         
!     A    - Matriz do sistema          =============================
!     b    - Vetor fonte                - Versão 1.0   01 / 01 / 2012
!     n    - Número de variáveis        - Autor: C.M. Barriga Nunes  
!====================================================================
IMPLICIT NONE
INTEGER,INTENT(in):: n
!REAL(8), DIMENSION(:,:),INTENT(inout):: A(n,n)
!REAL(8), DIMENSION(:),INTENT(inout):: b(n)
COMPLEX(8), DIMENSION(:,:),INTENT(inout):: A(n,n)
COMPLEX(8), DIMENSION(:),INTENT(inout):: b(n)
INTEGER:: i, j, k, cont

A(2:n,1) = A(2:n,1) / A(1,1)
cont=1
DO i=2, n
 cont = cont + 1
   do k=1, i-1
   b(i) = b(i) - b(k) * A(i,k)
   enddo
   DO j=cont, n
     DO k=1, i-1
       A(i,j) = A(i,j) - A(k,j) * A(i,k)
     END DO
     IF( j == n )EXIT
     IF( i == n )EXIT
     DO k=1, i-1
       A(j+1,i) = A(j+1,i) - A(j+1,k) * A(k,i)
     END DO
     A(j+1,i) =A(j+1,i) / A(i,i)
  END DO
END DO
b(n) = b(n) / A(n,n)
DO i=n-1, 1, -1
  b(i) = (b(i) - sum(A(i,i+1:n) * b(i+1:n))) / A(i,i)
END DO
END SUBROUTINE LU_Gauss_compacto2

!--

end module
