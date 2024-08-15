module variaveis_globais

integer,parameter :: db = 8
integer :: nelm,nnod,Nv
integer,allocatable  :: mel(:,:), flag(:)
real(db),allocatable :: mnod(:,:)
real(db),parameter   :: pi = 3.14159265358979323846264338327950288419716d0
real(db),parameter   :: muo = 4.d0*pi*1.d-7, epo = 8.8541878176d-12
real(db) :: freq
complex(db),parameter :: im = dcmplx(0.d0,1.d0)
real(db),allocatable  :: prop(:), rhoc(:)
real(db) :: rhob
integer :: cdim = 2, ngl = 1

integer,allocatable :: ia(:), ja(:)
complex(db),allocatable :: vnz(:), vft(:)
real(db),allocatable :: mmass(:), mstiff(:)

real(db),allocatable :: bmass(:,:), bstiff(:,:) 

end module
