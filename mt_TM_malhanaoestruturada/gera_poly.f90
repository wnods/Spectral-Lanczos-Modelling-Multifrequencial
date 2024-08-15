program gpoly
implicit none

integer,parameter :: db = 8
integer :: i,j,k
real(db),allocatable :: cbd(:,:), rhoc(:), Rx(:,:)
real(db) :: xim, xfm, zim, zfm
real(db) :: rhob, aux, Area, quantel
integer  :: Nc, Nr, cno, lbbd = 33
character(len=20) :: auxc
real(db) :: inc = 2.d0

open(10,file='input_malha.in',status='old',action='read')

read(10,*)auxc
read(10,*)xim, xfm
read(10,*)zim, zfm

read(10,*)auxc
read(10,*)rhob

read(10,*)auxc
read(10,*)Nc, quantel

allocate(cbd(2*Nc,2),rhoc(Nc))
j = 1
do i = 1,Nc
   read(10,*)rhoc(i)
   read(10,*)cbd(j,:)
   read(10,*)cbd(j+1,:)
   j = j + 2
end do    

open(20,file='coordenadas_receptores.in',status='old',action='read')
read(20,*)Nr
allocate(Rx(Nr,2))
do i = 1,Nr
   read(20,*)Rx(i,1),aux,Rx(i,2)
end do

open(30,file='modelo2D.poly',status='replace',action='write')

write(30,'(4I5)') 5*Nr + 4*Nc + 4, 2,  0, 1
write(30,*)''
write(30,*)'# Regiao dos nos'
write(30,*)''
write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')1,xim,zim,lbbd
write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')2,xfm,zim,lbbd
write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')3,xfm,zfm,lbbd
write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')4,xim,zfm,lbbd

write(30,*)''
cno = 4
j   = 1
do i = 1,Nc
   write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+j, cbd(j,1), cbd(j,2), 0
   write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+j+1, cbd(j+1,1), cbd(j,2), 0
   write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+j+2, cbd(j+1,1), cbd(j+1,2), 0
   write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+j+3, cbd(j,1), cbd(j+1,2), 0
   j = j + 4      
end do
write(30,*)''
cno = Nc*4 + 4
do i = 1,Nr
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+i, Rx(i,1), Rx(i,2), 22
end do
write(30,*)''
cno = Nc*4 + 4 + Nr
do i = 1,Nr
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+i, Rx(i,1)+inc, Rx(i,2), 0
end do
write(30,*)''
cno = Nc*4 + 4 + 2*Nr
do i = 1,Nr
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+i, Rx(i,1)-inc, Rx(i,2), 0
end do
write(30,*)''
cno = Nc*4 + 4 + 3*Nr
do i = 1,Nr
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+i, Rx(i,1), Rx(i,2)+inc, 0
end do
write(30,*)''
cno = Nc*4 + 4 + 4*Nr
do i = 1,Nr
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5)')cno+i, Rx(i,1), Rx(i,2)-inc, 0
end do


write(30,*)''
write(30,*)'# Regiao das faces'
write(30,*)''
write(30,'(2I5)') 4 + 4*Nc, 1
write(30,*)''

write(30,'(I6,2x,I6,2x,I6,2x,I5)')1,1,2,lbbd
write(30,'(I6,2x,I6,2x,I6,2x,I5)')2,2,3,lbbd
write(30,'(I6,2x,I6,2x,I6,2x,I5)')3,3,4,lbbd
write(30,'(I6,2x,I6,2x,I6,2x,I5)')4,4,1,lbbd

write(30,*)''
cno = 4
j   = 1 ; k = 4;
do i = 1,Nc
   write(30,'(I6,2x,I6,2x,I6,2x,I5)')cno+j, k+1, k+2, 0
   write(30,'(I6,2x,I6,2x,I6,2x,I5)')cno+j+1, k+2, k+3, 0
   write(30,'(I6,2x,I6,2x,I6,2x,I5)')cno+j+2, k+3, k+4, 0
   write(30,'(I6,2x,I6,2x,I6,2x,I5)')cno+j+3, k+4, k+1, 0
   j = j + 4 ; k = k + 4
   write(30,*)''
end do

write(30,*)'Buracos'
write(30,*) 0
write(30,*)''
write(30,*)'Propriedades das regioes'
write(30,*)''
write(30,'(I6)')1+Nc
write(30,*)''

Area = (xfm-xim)*(zfm-zim)/quantel

!write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5,2x,F15.8)')1, 0.d0, zfm - 100.d0, 10, -1.0
write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5,2x,F15.8)')1, 0.d0, zfm - 100.d0, 10, Area/2.d0

cno = 1
j   = 1
do i = 1,Nc
    Area = (cbd(j+1,1)-cbd(j,1))*(cbd(j+1,2)-cbd(j,2))/quantel
    write(30,'(I6,2x,F15.8,2x,F15.8,2x,I5,2x,F15.8)')cno+i, (cbd(j,1)+cbd(j+1,1))/2.d0, &
         (cbd(j,2)+cbd(j+1,2))/2.d0, (1+i)*10, Area
    j = j + 2
end do

end program
