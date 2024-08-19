PROGRAM SDLM
  

  IMPLICIT NONE
  REAL(8),PARAMETER:: PI      = 3.141592653589793238462643383279502884197d0
  REAL(8),ALLOCATABLE,DIMENSION(:,:):: Stiff_Mat, Diag, Mass_Mat,Imat, MatA,Diag_inv,Hess_Mat
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:):: V, V_1, Mat_temp, Mat_1, Mat_2, Mat_3, Mat_4,eig_vec,eig_vec_transp, Mat_eig,Inv_Mat_eig,Conj_Mat_eig
  REAL(8),ALLOCATABLE,DIMENSION(:)::Source_Vec_real, Source_Vec_imag,mod_Mat_eig,vec_e
  REAL(8),ALLOCATABLE,DIMENSION(:)::alpha,alpha_tmp,beta,mod_diag, beta_tmp,rnd
  REAL(4),ALLOCATABLE,DIMENSION(:)::eig
  COMPLEX*16,ALLOCATABLE,DIMENSION(:)::Source_Vec, VecB , Vec_temp, Vec_res,V_0,vec_1,vec_2,field
  COMPLEX(8)::imag, t1, t2,t3
  REAL(8):: freq,beta_0,mod_VecB,sum_VecB,eig_0,omega,ni,mod_field
  INTEGER(4)::nlin,ntemp,ncol,ierr,j,j1,j2,k1,k2,k3,k4,k10,right_index,k,left_index,opt

  CHARACTER*64::fname,fname1, fname2, fname3, fname4,fname5,fname6,fname7,fname8
  REAL time_begin, time_end


  CALL CPU_TIME(time_begin)


  imag=(0.d0,1.d0, KIND=8);  !i imaginario puro

!==================================================
!Entrada da matriz
!==================================================
 write(*,*)'Entre com a opcao para a matriz: calculada(1), lida de de arquivo(2), matrizes de EM (3)'
 read(*,*)opt

elseif(opt==1)then
    ! Set dimensions for the matrix
    nlin = 10
    ncol = 10

    ! Allocate memory for MatA
    ALLOCATE(MatA(nlin, ncol), STAT=ierr)
    IF (ierr /= 0) THEN
        STOP 'Nao foi possivel allocar memoria para MatA'
    END IF

    ! Initialize MatA as a diagonal matrix with incrementing values
    MatA = 0.0d0
    DO j = 1, nlin
        MatA(j, j) = 0.01 * j
    END DO
    MatA(10, 10) = 1000.0d0  ! Set a large value at the bottom right

    ! You can add more logic here as required by your program
    
elseif(opt==2)then
	write(*,*)'Entre com o nome do arquivo da matriz'
	read(*,*)fname7
	open(111, FILE=fname7,status ='old')
 	read(111,*)nlin,ncol

  ALLOCATE(MatA(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
    STOP 'Nao foi possivel allocar memoria para MatA'
  END IF

  do j1=1,nlin
    read(111,*)(MatA(j1,j2),j2=1,ncol)
  end do

 

elseif(opt==3)then
!MAT Stiffness
!write(*,*)'Entre com a dimensão das matrizes de Stiffness e Mass'
!read(*,*)
nlin=1089
ncol=1089


  ALLOCATE(Stiff_Mat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Stiff_Mat'
  END IF


	!write(*,*)'Entre com arquivo da parte real da matriz de Stiffness'
	!read(*,*) 
    fname1 =  'Stiff.dat'!'real_Stiff10.dat' !'matC_real.dat' 

	!write(*,*)'Entre com arquivo da parte imaginaria da matriz de Stiffness'
	!read(*,*) 

  	open(10, FILE=fname1,status ='old')
  	
   do j1=1,nlin
        read(10,*)(Stiff_Mat(j1,j2),j2=1,ncol)  
   end do

 !Entrada da matriz Mass 
 ! write(*,*)'Entre com arquivo da parte real da matriz de Massa'
 ! read(*,*) 
fname3 =  'Mass.dat'!'real_Mass10.dat' !'matT_real.dat'
 ! write(*,*)'Entre com arquivo da parte imaginaria da matriz de Massa'
 ! read(*,*) 

  open(11, FILE=fname3,status ='old')


  ALLOCATE(Mass_Mat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mass_Mat'
  END IF
    


  do j1=1,nlin
        read(11,*)(Mass_Mat(j1,j2),j2=1,ncol)    
   end do

!Entrada do vetor b fonte 
 
 !write(*,*)'Entre com arquivo d a parte real do vetor fonte'
 !read(*,*) 
fname5 =   'realB.dat'!'real_B10.dat' !'vecB_real.dat'


! write(*,*)'Entre com arquivo d a parte imag do vetor fonte'
! read(*,*) 
fname6 = 'imagB.dat'!'imag_B10.dat' !'vecB_imag.dat'

  open(15, FILE=fname5,status ='old')
  open(25, FILE=fname6,status ='old')


 
  ALLOCATE(Source_Vec_real(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec_real'
  END IF

  ALLOCATE(Source_Vec_imag(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec_imag'
  END IF


   do j1=1,nlin
        read(15,*)(Source_Vec_real(j1))  
        read(25,*)(Source_Vec_imag(j1))
   end do



  ALLOCATE(Source_Vec(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec'
  END IF

   Source_Vec = Source_Vec_real + imag*Source_Vec_imag

  ALLOCATE(MatA(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para MatA'
  END IF

 end if


!============================================================
!Outras alocações de memoria
ntemp=nlin-1
  ALLOCATE(alpha_tmp(nlin))

  Allocate(beta_tmp(ntemp))

  ALLOCATE(Diag(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Diag'
  END IF

  ALLOCATE(Imat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Imat'
  END IF

  ALLOCATE(Mat_eig(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_eig'
  END IF

  ALLOCATE(Inv_Mat_eig(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_eig'
  END IF


  ALLOCATE(Conj_Mat_eig(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_eig'
  END IF



  ALLOCATE(VecB(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

  ALLOCATE(Vec_1(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

  ALLOCATE(Vec_2(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF



  ALLOCATE(mod_Mat_eig(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF



  ALLOCATE(field(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF




  ALLOCATE(vec_e(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

 ALLOCATE(Vec_res(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

 ALLOCATE(eig(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


 ALLOCATE(eig_vec(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

ALLOCATE(eig_vec_transp(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(Diag_inv(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(Vec_temp(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(V(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF


ALLOCATE(V_0(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

ALLOCATE(V_1(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

ALLOCATE(Mat_temp(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

ALLOCATE(alpha(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para alpha'
  END IF

ALLOCATE(mod_diag(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para alpha'
  END IF

ALLOCATE(beta(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para beta'
  END IF

 ALLOCATE(Mat_1(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF

 ALLOCATE(Mat_2(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF


ALLOCATE(Mat_3(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF

ALLOCATE(Mat_4(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF

ALLOCATE(Hess_Mat(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Hess_Mat'
  END IF



if(opt==3)then
	Diag = 0.d0;
        Diag_inv = 0.d0;
	VecB = 0.d0;
	Imat = 0.d0;
	mod_VecB=0.d0;
	Mat_temp=0.d0;
        MatA = 0.d0;
        V=0.d0;

	do j1=1,nlin
		Diag(j1,j1) = SUM(Mass_Mat(j1,:));
                Imat(j1,j1) = 1;      
 		Diag_inv(j1,j1) = 1/sqrt(Diag(j1,j1));  !D^-1/2        
                VecB(j1)= Diag_inv(j1,j1)*Source_Vec(j1);
	end do



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPEN(41,FILE = 'VecB1.dat')
do j1=1,nlin
 WRITE(41,'(<nlin>(x,f17.4))') (VecB(j1)) 
end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           do j1=1,nlin
	             do j2=1,ncol
	                     MatA (j1,j2) =  Diag_inv(j1,j1)*Stiff_Mat(j1,j2)*Diag_inv(j2,j2);
		     end do
             end do



	mod_VecB = sqrt(dot_product(VecB,VecB));
         write(*,*)'mod',mod_VecB

	V(:,1) =  VecB/(mod_VecB);


elseif(opt==1)then

!falta
elseif(opt==2)then
	!MatA = MatA_real + imag*MatA_imag
	CALL random_number(rnd)
	beta_0 = 0;
	V_0 = 0;
	V(:,1) =  rnd/NORM2(rnd) 
endif


OPEN(44,FILE = 'matA.dat')
do j1=1,nlin
 WRITE(44,'(<nlin>(x,f17.4))') (matA(j1,j2),j2=1,nlin) 
end do


OPEN(5, file="eigvn.txt")
WRITE(5,*) "#Approximated eigenvalues"


!=============================================
!Recursao de lanczos
!m=numero de iterações
!=============================================
k1=1;
k2=0;
k3=0;
k4=0;
beta= 0.d0; 
alpha=0.d0;
beta_0 = 0.d0;
V_0 = 0.d0;

j=1;
 

k10=0



 DO WHILE (j.ne.nlin)
   if (k1==1)then
       
 	 Vec_temp= matmul(MatA,V(:,j));  !Avj
      
         Mat_temp(1,:) = V(:,j);              !transp(vj)            
         
         alpha(j) = DOT_PRODUCT(Mat_temp(1,:),Vec_temp(:));  
         
         Vec_res(:) = Vec_temp(:) - beta_0*V_0(:) - alpha(j)*V(:,j);

	 beta(j) = SQRT(DOT_PRODUCT(Vec_res(:),Vec_res(:))); 
        
         k1=j+1

         V(:,k1) = Vec_res(:)/beta(j)
   		
  else
    j=j+1
        Vec_temp(:) = matmul(MatA(:,:),V(:,j));  !Avj
      
        Mat_temp(1,:) = V(:,j);             !transp (vj)   
	alpha(j) = DOT_PRODUCT(Mat_temp(1,:),Vec_temp(:)); 

           k2=j-1
 	
        Vec_res(:) = Vec_temp(:) - beta(k2)*V(:,k2) - alpha(j)*V(:,j); !recursao lanczos     
       
    	beta(j) = SQRT(DOT_PRODUCT(Vec_res(:),Vec_res(:))); 
        k3=j+1


        if(beta(j)==0)then
          write(*,*)'beta nulo'
         stop
	else
		V(:,k3) = Vec_res(:)/beta(j)
        end if 

  end if
      	k10=j  !numero de iteraçoes para a determinacao da matriz tridiagonal
	write(*,*)'iter',k10
!write(*,*)j, alpha(j)
!pause
       
END DO

!%%%%%%%%%%%% Fim da iteração de Lanczos, produtos: base ortornolmal e mat tridiag%%%%%%%%%%%%%%%%%%%%%%%%%%


Mat_eig = 0.d0; !matriz diagonal com os autovalores. Matriz da formula 29
eig_vec = (0.d0,0.d0)

omega = 2*pi*10 !2*pi*f
ni = 4e-7*pi

write(5,*)'Iteracoes',j
write(*,*)'ultimo valor de beta',beta(j),j

  alpha_tmp = alpha
  beta_tmp((/1:ntemp/))  = beta((/1:ntemp/))

write(*,*)'ok1'

 CALL eigenvalues_complex(alpha_tmp,beta_tmp,eig_vec,k10) !eig_vect autovec compl e alpha_mim autoval reias da mat tridiag
 	
	eig=alpha_tmp

write(*,*)'ok2'

do j1=1,nlin
   write(111,*) eig(j1)
 	     Mat_eig(j1,j1) = eig(j1) + imag*omega*ni     !ja calculado (eig+j*omega*ni*identidade)**-1
 	     mod_Mat_eig(j1) = (real(Mat_eig(j1,j1))**2 + aimag(Mat_eig(j1,j1))**2);
             Conj_Mat_eig(j1,j1) = conjg(Mat_eig(j1,j1))
             Inv_Mat_eig(j1,j1)  = Conj_Mat_eig(j1,j1)/mod_mat_eig(j1)
            
enddo
write(*,*)'ok3'

j1=1

!do j1=2,nlin

 !	if(abs(abs(eig(j1))-abs(eig_0))<0.00001)then
 !	else

 !	     eig_0=eig(j1)
   !          write(5,*) eig_0
 !	     Mat_eig(j1,j1) = eig(j1)+imag*omega*ni     !ja calculado (eig+j*omega*ni*identidade)**-1
 !	     mod_Mat_eig(j1) = (real(Mat_eig(j1,j1))**2 + aimag(Mat_eig(j1,j1))**2);
  !           Conj_Mat_eig(j1,j1) = conjg(Mat_eig(j1,j1))
   !          Inv_Mat_eig(j1,j1)  = Conj_Mat_eig(j1,j1)/mod_mat_eig(j1)
     
 !   end if
               
!end do
  

write(*,*)'Calculo dos auto-valores --OK'

	open(31,file = 'Alphas.dat')

do j1=1,nlin
write(31,*)alpha(j1)
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Hess_Mat = 0


	do j1=1,nlin
    		Hess_Mat(j1,j1) = alpha(j1);
    		if(j1>1.and.j1<nlin-1)then
    
    			Hess_Mat(j1,j1+1) = beta(j1);
    			Hess_Mat(j1+1,j1) = beta(j1);
    		elseif(j1==1)then
         		Hess_Mat(j1,j1+1) = beta(j1);
         		Hess_Mat(j1+1,j1) = beta(j1);
    		elseif(j1==nlin)then  
        		Hess_Mat(nlin,nlin-1)=beta(nlin-1);
    		end if
	end do

	open(36,file = 'mat_tridiag.dat')

	do j1=1,nlin
 		write(36,'(<nlin>(x,f17.4))') (Hess_Mat(j1,j2),j2=1,nlin) 
	end do


         open(40,FILE = 'Base_ort.dat')

	do j1=1,nlin
 		write(40,'(<nlin>(x,f20.10))') real(V(j1,:)) 
	end do

write(*,*)'Calculo da mat tridiag--OK'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!SDLM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do j1=1,nlin
do j2=1,nlin
  eig_vec_transp(j1,j2) = eig_vec(j2,j1)
enddo
enddo

vec_e = 0.d0;
vec_e(1) = 1.d0;

!mult matrizes
vec_1 = 0.d0;
vec_2 = 0.d0;

CALL MATMAT(nlin,ncol,V,&
            nlin,ncol,eig_vec,&
            mat_1)

write(*,*)'ok4'
       
CALL MATMAT(nlin,ncol,mat_1,&
            nlin,ncol,Inv_Mat_eig,&
            mat_2)

write(*,*)'ok5'

vec_1 = 0.d0;

CALL MATMAT(nlin,ncol,Mat_2,&
            nlin,ncol,eig_vec_transp,&
            mat_3)

write(*,*)'ok6'
!multiplacacao de matriz or vetor

CALL MATVEC(nlin,ncol,Mat_3,vec_e,vec_1)
 write(*,*)'ok7'

do j1=1,nlin
field(j1) = mod_VecB*vec_1(j1)*Diag_inv(j1,j1);
end do 

mod_field = sqrt(dot_product(field,field))

!field = mod_VecB*V.*eig_vec.*Mat_eig.*eig_vec_transp.*vec_e 
OPEN(47,FILE = 'Campo_new.dat')
	do j1=1,nlin
 		write(47,'(<nlin>(x,E20.10E3))') (field(j1)) 
	end do



 write(*,*)'Acabou1',mod_field


    CALL CPU_TIME (time_end)

  WRITE(*,*) 'Time of operation was ', time_end - time_begin, ' seconds'


CONTAINS
SUBROUTINE eigenvalues_complex(d,u,eig_vec,n)
   
 !   Use nag_library, Only: nag_wp, x04dbf, zhbtrd, zstedc
 implicit none
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    CHARACTER :: COMPZ = 'I'                   !opcao I calcula os autovetores de H tambem
    INTEGER :: n,lgn			       !dimensão da matriz tridiagonal
    REAL(dp), dimension(n)::d                  !array with diagonal
    REAL(dp), dimension(n - 1):: u             !E upper diagonal
    COMPLEX*16,dimension(:,:), allocatable :: eig_vec !Z autovetores da mat tridiag
    INTEGER::LDZ                                !=n dimension of the array Z
    COMPLEX*16,dimension(:), allocatable :: WORK
    INTEGER :: LWORK, LIWORK                  !assumir -1
    REAL(dp), dimension(:), allocatable::RWORK
  !  Real (Kind=nag_wp)               :: rdum(1)
   ! Integer                          :: idum(1)
    INTEGER, dimension(:), allocatable :: IWORK
    INTEGER::LRWORK,INFO

!lgn = ceiling(log(real(n,kind=nag_wp))/log(2.0_nag_wp))
lgn = ceiling(log(real(n,kind=dp))/log(2.0_dp))
LDZ = N

ALLOCATE(WORK(1), IWORK(1),rwork(1))

 lwork  = -1
 lrwork = -1
 liwork = -1

CALL zstedc (COMPZ,N,d,u,eig_vec,n,WORK,lwork,RWORK,lrwork,IWORK,liwork,INFO ) 

!Call zstedc('V',n,d,e,z,ldz,cdum,lwork,rdum,lrwork,idum,liwork,info)

      lwork = max(n*n, nint(real(work(1))))
      lrwork = max(1+3*n+2*n*lgn+4*n*n, nint(rwork(1)))
      liwork = max(6+6*n+5*n*lgn, iwork(1))

DEALLOCATE(WORK, IWORK,RWORK)		

Allocate (work(lwork), rwork(lrwork), iwork(liwork))

CALL zstedc (COMPZ,N,d,u,eig_vec,n,WORK,lwork,RWORK,lrwork,IWORK,liwork,INFO ) 

DEALLOCATE(WORK, IWORK,RWORK)		


!DSTEVR(JOBZ,RANGE,N,D,E,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO ) 		
!DSTEVR(JOBZ,RANGE,n,d,u,VL,VU,IL,IU,ABSTOL,n,eig,Z,n,ISUPPZ,WORK, -1, IWORK, -1,INFO)



END SUBROUTINE eigenvalues_complex


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SUBROUTINE MATMAT(Anrow,Ancol,Amatrix,Bnrow,Bncol,Bmatrix,Cmatrix)
    IMPLICIT NONE
    INTEGER(4)::Anrow,Ancol,Bnrow,Bncol,i,j
    COMPLEX(8)   ::Amatrix(Anrow,Ancol),Bmatrix(Bnrow,Bncol),Cmatrix(Anrow,Bncol)

    DO j=1,Bncol
       DO i=1,Anrow
          Cmatrix(i,j) = 0
       END DO
    END DO
    DO j=1,Bncol
       CALL GAXPY(Anrow,Ancol,Amatrix,Bmatrix(:,j),Cmatrix(:,j))
    END DO
END SUBROUTINE MATMAT



  SUBROUTINE GAXPY(Nrow,Ncol,Amatrix,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,Ncol,i,j
    COMPLEX(8)   ::Amatrix(Nrow,Ncol),xvector(Ncol),yvector(Nrow),x_i

    DO i=1,Ncol
       x_i = xvector(i)
       CALL SAXPY(Nrow,x_i,Amatrix(:,i),yvector)
    END DO
  END SUBROUTINE GAXPY
  
  
  SUBROUTINE SAXPY(Nrow,alpha,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,i
    COMPLEX(8)   ::alpha,xvector(Nrow),yvector(Nrow)

    DO i=1,Nrow
       yvector(i) = yvector(i) + alpha*xvector(i)
    END DO

  END SUBROUTINE SAXPY

 
  SUBROUTINE MATVEC(Nrow,Ncol,Amatrix,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,Ncol,i,j
    COMPLEX(8)   ::Amatrix(Nrow,Ncol),yvector(Nrow),x_i
    REAL(8):: xvector(Ncol)

    DO i=1,Nrow
       yvector(i) = 0.D0
    END DO

    CALL GAXPY2(Nrow,Ncol,Amatrix,xvector,yvector)

  END SUBROUTINE MATVEC
 SUBROUTINE GAXPY2(Nrow,Ncol,Amatrix,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,Ncol,i,j
    COMPLEX(8)   ::Amatrix(Nrow,Ncol),yvector(Nrow)
    REAL(8):: xvector(Ncol),x_i

    DO i=1,Ncol
       x_i = xvector(i)
       CALL SAXPY2(Nrow,x_i,Amatrix(:,i),yvector)
    END DO
  END SUBROUTINE GAXPY2


  SUBROUTINE SAXPY2(Nrow,alpha,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,i
    COMPLEX(8)   ::xvector(Nrow),yvector(Nrow)
    REAL(8)::alpha

    DO i=1,Nrow
       yvector(i) = yvector(i) + alpha*xvector(i)
    END DO

  END SUBROUTINE SAXPY2

  
END PROGRAM SDLM

