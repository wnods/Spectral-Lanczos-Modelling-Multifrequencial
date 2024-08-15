Program test_sparsemul
use sparsemul

Implicit None

Integer             :: na,ma,nnza,nnzb,nb,mb
Real,Allocatable    :: matA(:,:),matB(:,:),matC(:,:),matC2(:,:)
Integer             :: i,j
Real		    :: start,finish

Open(unit = 1, file = 'dimA.txt', status = 'old', ACTION = 'read')
Open(unit = 2, file = 'dimB.txt', status = 'old', ACTION = 'read')
Open(unit = 3, file = 'matA.txt', status = 'old', ACTION = 'read')
Open(unit = 4, file = 'matB.txt', status = 'old', ACTION = 'read')

Read(1,*) na,ma,nnza
Read(2,*)nb,mb,nnzb


Allocate(matA(na,ma),matB(nb,mb),matC(na,mb),matC2(na,mb))

Do i = 1,na
	Read(3,*) (matA(i,j),j = 1,ma)
End Do

Do i = 1,nb
	Read(4,*) (matB(i,j),j = 1,mb)
End Do


 
 Call matmul_sparse(matA,na,ma,matB,nb,mb,matC)



Print*,'Programa Terminado!'
End Program
