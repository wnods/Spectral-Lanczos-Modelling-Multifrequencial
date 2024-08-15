program mt2D

!----------------------------
!   Programa para o calculo do campo eletrico e magnetico
!para o modo TM do metodo Magnetotelurico, em meios 2D.
!Autor: Anderson Almeida.   19/05/2020        
!----------------------------

use variaveis_globais
use mod_mtTM 

implicit none
integer  :: i,j,k,l,cnt
real(db) :: ti,tf,auxr,w,sigb,cx,cy,cz
character(len=100) :: files, auxc
complex(db),allocatable :: G(:), F(:), vcc(:)
integer :: nnz, erro, iaux
integer,allocatable :: bounds(:), nbd(:)
integer :: Nc, flagb = 33, flagrx = 22, inb, qnbd, nRx
complex(db) :: Hyp, Exp
integer,allocatable :: indRx(:)
complex(db),allocatable :: ExEz(:,:), Glob(:,:)
real(db) :: tii, tff, sm
complex(db) :: Imp

character(len=20)      :: frq, pasta

call getarg(1,frq)
call getarg(2,pasta)
print*,trim(frq),' ',trim(pasta)

print*,''
print*,'Gerando a malha para o modelo...'
print*,''
call execute_command_line('gfortran gera_poly.f90 -o gpoly.x')
call execute_command_line('./gpoly.x')
call execute_command_line('triangle -pjq23aA modelo2D.poly')
print*,''
print*,'feito.'
print*,''

files = 'modelo2D'

!-- Lendo o arquivo de elementos
open(10,file=trim(files)//'.1.ele',status='old',action='read')
read(10,*)nelm ; allocate(mel(nelm,3),prop(nelm),flag(nelm));
do i = 1,nelm
   read(10,*)iaux,mel(i,:),flag(i)
end do; close(10)
!-- Lendo os arquivos dos nos
open(10,file=trim(files)//'.1.node',status='old',action='read')
read(10,*)nnod ; allocate(mnod(nnod,2),bounds(nnod));
do i = 1,nnod
   read(10,*)iaux,mnod(i,:),bounds(i)
end do; close(10)
!--

!-- lendo as resistividades do modelo

open(10,file='input_malha.in',status='old',action='read')
read(10,*)auxc
read(10,*)auxr
read(10,*)auxr
read(10,*)auxc
read(10,*)rhob
read(10,*)auxc
read(10,*)Nc
allocate(rhoc(Nc))
do i = 1,Nc
   read(10,*)rhoc(i)
   read(10,*)auxr
   read(10,*)auxr
end do

read(frq,*)freq
   
!-- Impondo as proriedades ao modelo

do i = 1,nelm

   if(flag(i)==10)then
      prop(i) = rhob
   else   
     do j = 1,Nc
      if(flag(i)==((j+1)*10))then
         prop(i) = rhoc(j)
      end if
     end do    
   end if

end do

!-- Nos da borda --

inb = 0
do i = 1,nnod
   if(bounds(i)==flagb) inb = inb + 1 ! Contador para os nos da borda
end do
qnbd = inb
allocate(nbd(qnbd), vcc(qnbd))
vcc = 0.d0
inb = 0
do i = 1,nnod
   if(bounds(i)==flagb)then
      inb = inb + 1
      nbd(inb) = i
   end if   
end do   

close(10)

!--

Nv = nnod ! Numero de variaveis

!-- Iniciando as matrizes CSR

call initial_arrays(cdim,nelm,nnod,ngl,nnz,mel,ia,ja)

allocate(G(nnz),F(Nv))
allocate(mmass(nnz),mstiff(nnz))

allocate(bmass(nnod,nnod),bstiff(nnod,nnod)); bmass = 0.d0; bstiff = bmass;

print*,'nnz,nnod',nnz,nnod

mmass = 0.d0 ; mstiff = 0.d0

!-- Montagem da matriz e vetor fontes globais

call cpu_time(ti)
call montagem_matriz_vetorfonte_globais(Nv,nnz,G,F)
call cpu_time(tf)
print*,''
print*,'Tempo para montar as matrizes:',tf-ti,'seg.'
print*,''

!call recount_arrays(ngl,nnod,ja,ia,G,nnz)

call condicao_contorno( cdim, ngl, qnbd , nbd , vcc , ia , nnz , G , nnod , F )

open(12,file='MassStif.dat',status='replace',action='write')
open(13,file='ja.dat',status='replace',action='write')
open(14,file='ia.dat',status='replace',action='write')
open(15,file='vecB.dat',status='replace',action='write')

do i = 1,nnz
   write(12,'(2(E20.10E3,2x))')mmass(i), mstiff(i)
   write(13,*)ja(i)
end do
do i = 1,nnod+1
   write(14,*)ia(i)
   if(i<=nnod) write(15,'(2(E20.10E3,2x))')dreal(F(i)),dimag(F(i))
end do

do i = 1,qnbd
   bmass(nbd(i),nbd(i)) = 1.d20
   bstiff(nbd(i),nbd(i)) = 1.d20
end do

call execute_command_line("mkdir "//trim(pasta))

!--
print*,'escrevendo as matrizes'
!open(100,file='mass.dat',status='replace',action='write' )
!open(110,file='stiff.dat',status='replace',action='write')
open(120,file='realB.dat',status='replace',action='write')
open(121,file='imagB.dat',status='replace',action='write')

print*,nnod

do i = 1,nnod
   sm = sum(bmass(i,:)) ; if(sm == 0.d0) print*,i
!   write(100,'(<nnod>(E20.10E3,2x))')(bmass(i,j),j=1,nnod)   
!   write(110,'(<nnod>(E20.10E3,2x))')(bstiff(i,j),j=1,nnod)
!   write(100,'((E20.10E3,2x))')(bmass(i,j),j=1,nnod)   
!   write(110,'((E20.10E3,2x))')(bstiff(i,j),j=1,nnod)

   write(120,'(1(E20.10E3,2x))')dreal(F(i))
   write(121,'(1(E20.10E3,2x))')dimag(F(i))
end do
! close(100);
! close(110); 
 close(120); 
 close(121)
print*,'feito'
!--

allocate(Glob(nnod,nnod))

print*,'Resolvendo o sisitema linear'

call cpu_time(ti)

!Glob = bstiff + dcmplx(0.d0,1.d0)*(2.d0*pi*freq)*muo*bmass
!call LU_Gauss_compacto2( Glob, F, nnod)
call solver_pardiso ( ngl , nnod , nnz , G , ia , ja , F )

call cpu_time(tf)
print*,''
print*,'Tempo para resolver o sistema:',tf-ti,'seg.'
print*,''

open(40,file='Hyfield.dat',status='replace',action='write')
open(44,file='coordenadas_malha.dat',status='replace',action='write')
open(50,file='coordenadas_receptores.in',status='old',action='read')
open(41,file='REAL.dat',status='replace',action='write')
open(42,file='IMAG.dat',status='replace',action='write')

read(50,*)nRx
read(50,*)cx,cy,cz

w = 2.d0*pi*freq
sigb = 1.d0/rhob

 cz = 0.d0
call primario(cz,sigb,w,Hyp,Exp)

allocate(indRx(nRx))

 cnt = 0
do i = 1,nnod
   if(bounds(i)==flagrx)then
      cnt=cnt+1
      indRx(cnt) = i
!      write(40,'(E20.10E3,2x,E20.10E3)')dreal(Hyp+F(i)), dimag(Hyp+F(i))
      write(40,'(E20.10E3,2x,E20.10E3)')dreal(F(i)), dimag(F(i))
   end if  
      write(44,'(f20.10,2x,f20.10)')(mnod(i,1)), (mnod(i,2)) 
      write(41,'(E20.10E3,2x)')dreal(F(i)) 
      write(42,'(E20.10E3,2x)')dimag(F(i)) 
end do   

!--

allocate(ExEz(nRx,2)) ; ExEz = 0.d0

call derivada_MFB( nnod, nelm, nRx, indRx, F, mel , prop , mnod, ExEz ) 

open(60,file='ExEzfield.dat',status='replace',action='write')
open(61,file='indices_receptores.dat',status='replace',action='write')

do i = 1,nRx
!     write(60,'(4(E20.10E3,2x))')dreal(Exp+ExEz(i,2)),dimag(Exp+ExEz(i,2)),dreal(ExEz(i,1)),dimag(ExEz(i,1))
     write(60,'(4(E20.10E3,2x))')dreal(ExEz(i,2)),dimag(ExEz(i,2)),dreal(ExEz(i,1)),dimag(ExEz(i,1))
     write(61,*) indRx(i)
end do

open(70,file='Res_apar.dat',status='replace',action='write')

do i = 1,nRx
   Imp = (Exp+ExEz(i,2))/(Hyp+F(indRx(i)))
   write(70,'(2(E20.10E3,2x))')(1.d0/(w*muo))*abs( Imp )**2 , (180.d0/pi) * datan2(dimag(Imp),real(Imp))
!   write(70,'(1(E20.10E3,2x))')(1.d0/(w*muo))*abs((ExEz(i,2))/(F(indRx(i))))**2
end do

call execute_command_line(" cp realB.dat imagB.dat REAL.dat IMAG.dat ./"//trim(pasta))
call execute_command_line(" cp indices_receptores.dat coordenadas_receptores.in ./matrizes")

end program
