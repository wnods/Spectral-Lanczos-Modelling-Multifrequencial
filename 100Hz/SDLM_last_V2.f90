PROGRAM SDLM_last
  !======================================================================================================
  ! Programa para calculo do Spectral Decomposition Lanczos Method
  ! programado para o discretizador de diferenças finitas
  !Ellen Gomes
  ! 28/11/2019
  ! 16/12/2019
  ! 24/12/2019
  ! Os testes no Matlab indicam que a matriz A é simetrica!!!!! 
  !27/12/2019
  !28/12/2019
  ! 07/01/2020
  !15/01/2020
  !20/01/2020 Tentativa de separação da contribuicao da fonte no vetor B
!27/01/2020 Otimização
!11/02/2020 numero de iteracoes do alg de lanczos
!01/04/2020 Versao do prgrama que roda os dados so anderson. Esta versão e a otimizada
!Utilização de subrotinas para matrizes esparcas - 06/04/2020
  !=====================================================================================================
  IMPLICIT NONE
  REAL(8),PARAMETER:: PI      = 3.141592653589793238462643383279502884197d0
  REAL(8),ALLOCATABLE,DIMENSION(:,:):: Stiff_Mat, Diag, Mass_Mat,Imat, Diag_inv,Hess_Mat
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:):: V, Mat_temp,eig_vec,eig_vec_transp,Mat_eig,Inv_Mat_eig,Conj_Mat_eig,Mat_source,&
                                          Mat_1,Mat_2,Mat_3, MatB,MatC,MatA,MAtA1
  REAL(8),ALLOCATABLE,DIMENSION(:)::Source_Vec_real, Source_Vec_imag, mod_Mat_eig
  REAL(8),ALLOCATABLE,DIMENSION(:)::alpha,alpha_tmp,beta, beta_tmp
  REAL(4),ALLOCATABLE,DIMENSION(:)::eig
  COMPLEX(8),ALLOCATABLE,DIMENSION(:)::Source_Vec, VecB , Vec_temp, Vec_res,V_0,vec_1,field,vec_e
  COMPLEX(8)::imag, t1, t2,t3
  REAL(8):: freq,beta_0,mod_VecB,sum_VecB,eig_0,omega,ni,mod_field,soma
  INTEGER(4)::nlin,niter,ntemp,ncol,ierr,j,j1,j2,j4,j5,k1,k2,k3,k4,k10,right_index,k,left_index,opt
!Matriz Esparsa
  COMPLEX(8),ALLOCATABLE,DIMENSION(:)::Val_MatA
  INTEGER(4),ALLOCATABLE,DIMENSION(:)::row_MatA,colPtr_MatA
  INTEGER(4)::nnza

  CHARACTER*64::fname,fname1, fname2, fname3, fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11
  REAL time_begin, time_end


 CALL CPU_TIME(time_begin)
 


  imag=(0.d0,1.d0);  !i imaginario puro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Entrada de Dados
!==================================================
!1 - Matrizes
!==================================================

!MAT Stiffness
!write(*,*)'Entre com a dimensão das matrizes de Stiffness e Mass'
!read(*,*)
nlin=5041
ncol=nlin


  ALLOCATE(Stiff_Mat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Stiff_Mat'
  END IF


  ALLOCATE(Mass_Mat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mass_Mat'
  END IF
 

  ALLOCATE(Source_Vec_real(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec_real'
  END IF

  ALLOCATE(Source_Vec_imag(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec_imag'
  END IF

ALLOCATE(Source_Vec(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Source_Vec'
  END IF




	!write(*,*)'Entre com arquivo da parte real da matriz de Stiffness'
	!read(*,*) 
    fname1 =  'Stiff.dat'!'real_Stiff10.dat' !'matC_real.dat' 

	!write(*,*)'Entre com arquivo da parte imaginaria da matriz de Stiffness'
	!read(*,*) 

    fname3 =  'Mass.dat'!'real_Mass10.dat' !'matT_real.dat'
 ! write(*,*)'Entre com arquivo da parte imaginaria da matriz de Massa'
 ! read(*,*) 


 !write(*,*)'Entre com arquivo d a parte real do vetor fonte'
 !read(*,*) 
fname5 = 'realB100.dat'!'real_B10.dat' !'vecB_real.dat'


fname6 = 'imagB100.dat'
! write(*,*)'Entre com arquivo d a parte imag do vetor fonte'
! read(*,*) 


!fname9 = 'Source100.dat'!'imag_B10.dat' !'vecB_imag.dat'



open(15, FILE=fname5,status ='old')

  do j1=1,nlin
        read(15,*)Source_Vec_real(j1)
 end do
open(25, FILE=fname6,status ='old')
  do j1=1,nlin
        read(25,*)Source_Vec_imag(j1) 
  end do   

 open(11, FILE=fname3,status ='old')
  do j1=1,nlin
        read(11,*)(Mass_Mat(j1,(/1:ncol/))) 
  end do
 open(10, FILE=fname1,status ='old')
  	
   do j1=1,nlin
        read(10,*)(Stiff_Mat(j1,(/1:ncol/)))  
   end do

 !Entrada da matriz Mass 
 ! write(*,*)'Entre com arquivo da parte real da matriz de Massa'
 ! read(*,*) 


     CALL CPU_TIME (time_end)

 write(*,*) 'Time of operation was ', time_end - time_begin, ' seconds'


 ! do j1=1,nlin
 !       read(11,*)(Mass_Mat(j1,j2),j2=1,ncol)  
 !  end do
!====================================================
!2 - Vetor b fonte 
!====================================================




 !  do j1=1,nlin
  !      read(15,*)Source_Vec_real(j1)
   !     read(25,*)Source_Vec_imag(j1)
       
  ! end do

   Source_Vec = Source_Vec_real + imag*Source_Vec_imag


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Alocações de memória
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Preparacao do sistema para uso na recurssao de Lanczos
  ALLOCATE(Diag(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Diag'
  END IF


  ALLOCATE(Diag_inv(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Diag_inv'
  END IF

  ALLOCATE(VecB(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(Imat(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Imat'
  END IF


  ALLOCATE(MatA(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para MatA'
  END IF

  ALLOCATE(MatA1(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para MatA'
  END IF
 
!-----------------------------------------------------------------------
!Alocação de memória para a recurssão de Laczos
!V base ortornormal, dim=niterxniter;
!V_0 base ortornormal de indice 0;
!Vec_temp, Vec_res variaveis usadas durante a recurssao, dim=nlin;
!alpha e beta, produtos do algoritmo de recurssao, forma a matriz tridiagonal, dim=niter;
!Hess_Mat, produto da recurssao, dim=niterxniter;
!-----------------------------------------------------------------------

   niter=1335 !numero de iterações calculado do nuremo de autovetores nao iguais (difereca entre eles <0.0001)



  ALLOCATE(V(nlin,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para V'
  END IF

  
  ALLOCATE(V_0(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF

  
  ALLOCATE(Vec_temp(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(Mat_temp(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_temp'
  END IF


ALLOCATE(MatB(nlin,ncol),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para MatB'
  END IF



ALLOCATE(MatC(nlin,1),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para MatC'
  END IF


  ALLOCATE(alpha(niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para alpha'
  END IF

  ALLOCATE(beta(niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para beta'
  END IF


  ALLOCATE(Vec_res(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para VecB'
  END IF


  ALLOCATE(Hess_Mat(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Hess_Mat'
  END IF
!-----------------------------------------------------------------------
!Calculo dos autovalores da Hessiana
!alpha_temp, beta_temp varlores de alpha e beta a serem usados no calculo dos autovalores; 
!Os auto-valoes seram alocados em alpha_temp, dim = niter;
!eig_vec os auto-vetores, dim=niterxniter;
! eig, os auto-valores de alpha_temp realocados;
! Mat_eig, matriz dada por: = eig(j1) + imag * omega * ni. Dim = niterxniter;
! mod_Mat_eig, vetor com os modolos da Mat_eig, dim = niter;
! Conj_Mat_eig, conjugado da matriz Mat_eig, dim = niterxniter;
! Inv_Mat_eig, inversa da matriz Mat_eig, dim = niterxniter;
!-----------------------------------------------------------------------


  ntemp=niter-1

  ALLOCATE(alpha_tmp(niter))
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para alpha_tmp'
  END IF

  Allocate(beta_tmp(ntemp))
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para beta_tmp'
  END IF


  ALLOCATE(eig_vec(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para eig_vec'
  END IF


  ALLOCATE(eig(niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para eig'
  END IF


  ALLOCATE(Mat_eig(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_eig'
  END IF

  ALLOCATE(mod_Mat_eig(niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para mod_Mat_eig'
  END IF

  ALLOCATE(Conj_Mat_eig(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Conj_Mat_eig'
  END IF

  ALLOCATE(Inv_Mat_eig(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Inv_Mat_eig'
  END IF
  
!-----------------------------------------------------------------------
!Variáveis envolvidas no calculo da decomposição espectral de Lanczos
! eig_vec_transp transposta da matriz eig_vec, dim = niterxniter;
! Mat_1 matriz usada no calculo de matrizes, dim=nlinxniter;
! Mat_2 matriz usada no calculo de matrizes, dim=nlinxniter;
! Mat_3 matriz usada no calculo de matrizes, dim=nlinxniter;
! vec_e vetor usado no calculo, dim = niter;
! Vec_1 vetor usado no calculo;
!field campo calculado, dim = nlin;
!-----------------------------------------------------------------------


  ALLOCATE(eig_vec_transp(niter,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para eig_vec_transp'
  END IF


  ALLOCATE(Mat_1(nlin,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_1'
  END IF


  ALLOCATE(Mat_2(nlin,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_1'
  END IF


  ALLOCATE(Mat_3(nlin,niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_1'
  END IF


  ALLOCATE(vec_e(niter),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para vec_e'
  END IF


  ALLOCATE(field(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para field'
  END IF


  ALLOCATE(Vec_1(nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Vec_1'
  END IF


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculos
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	Diag = 0.d0;
        Diag_inv = 0.d0;
	VecB = 0.d0;
	Imat = 0.d0;
	mod_VecB=0.d0;
	Mat_temp=0.d0;
        MatA = 0.d0;
        V=0.d0;
      

	do j1=1,nlin
		Diag(j1,j1) =  SUM(Mass_Mat(j1,:));                    
                Imat(j1,j1) = 1.d0;      
 		Diag_inv(j1,j1) = 1.d0/sqrt(Diag(j1,j1));  !D^-1/2
                VecB(j1)= Diag_inv(j1,j1)*Source_Vec(j1);
                 if(Diag(j1,j1)==0.d0)then
                   write(*,*)'Valores Nulos',j1,Diag(j1,j1)
			pause
		endif
	end do




           do j1=1,nlin
	             do j2=1,ncol
	                     MatA (j1,j2) =  Diag_inv(j1,j1)*Stiff_Mat(j1,j2)*Diag_inv(j2,j2);
		     end do
             end do



!Reescrevendo a matriz esparsa
  nnza = count(matA /= 0)

  ALLOCATE(Val_MatA(nnza),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Val_MatA'
  END IF

 ALLOCATE(row_MatA(nnza),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para row_MatA'
  END IF

 ALLOCATE(colPtr_MatA(nlin+1),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para colPtr_MatA'
  END IF


  CALL mat_sparse(MatA,nlin,ncol,row_MatA,ColPtr_MatA,Val_MatA,nnza)


!==========================================================
!==========================================================
!==========================================================


	mod_VecB = sqrt(dot_product(VecB,VecB));

	V(:,1) =  VecB/(mod_VecB);


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
 
MatB = 0.d0;
k10=0
        
 DO WHILE (j.ne.niter)
   if (k1==1)then

        matB(((/1:nlin/)),1)=V(((/1:nlin/)),j);
         

        Call matmul_sparse2(nlin,ncol,row_MatA,ColPtr_MatA,Val_MatA,nnza,V(((/1:nlin/)),j),nlin,1,MatC)
        

        vec_temp((/1:nlin/)) = MatC(((/1:nlin/)),1)
    
        
        Mat_temp(1,((/1:nlin/))) = V(((/1:nlin/)),j);              !transp(vj)            
         
         alpha(j) = DOT_PRODUCT(Mat_temp(1,((/1:nlin/))),Vec_temp((/1:nlin/)));  
         
         Vec_res((/1:nlin/)) = Vec_temp((/1:nlin/)) - beta_0*V_0((/1:nlin/)) - alpha(j)*V(((/1:nlin/)),j);

	 beta(j) = SQRT(DOT_PRODUCT(Vec_res((/1:nlin/)),Vec_res((/1:nlin/)))); 
        
         k1=j+1

         V(((/1:nlin/)),k1) = Vec_res(((/1:nlin/)))/beta(j)
 
   !write(*,*)'alpha,beta1',k1,j,alpha(j),beta(j)
   !pause
  else
    j=j+1

      matB(((/1:nlin/)),1)=V(((/1:nlin/)),j);

    Call matmul_sparse2(nlin,ncol,row_MatA,ColPtr_MatA,Val_MatA,nnza,V(((/1:nlin/)),j),nlin,1,MatC)
        
        
       vec_temp((/1:nlin/)) = MatC(((/1:nlin/)),1)
        
     
        Mat_temp(1,((/1:nlin/))) = V(((/1:nlin/)),j);             !transp (vj)   
	
        
       alpha(j) = DOT_PRODUCT(Mat_temp(1,((/1:nlin/))),Vec_temp(((/1:nlin/)))); 

           k2=j-1

        Vec_res((/1:nlin/)) = Vec_temp((/1:nlin/)) - beta(k2)*V(((/1:nlin/)),k2) - alpha(j)*V(((/1:nlin/)),j); !recursao lanczos     
       

    	beta(j) = SQRT(DOT_PRODUCT(Vec_res((/1:nlin/)),Vec_res((/1:nlin/)))); 

        k3=j+1
	!write(*,*)'alpha,beta',j,alpha(j),beta(j)
        !pause

        if(beta(j)==0)then
          write(*,*)'VAlores de beta nulo'
         stop
	else
		V(((/1:nlin/)),k3) = Vec_res((/1:nlin/))/beta(j)
        end if 
      
  end if
      	k10=j  !numero de iteraçoes para a determinacao da matriz tridiagonal

       
END DO

!%% Fim da iteração de Lanczos, produtos: base ortornormal e mat tridiag%%%%


!========================
!Calculo dos autovalores
!========================

Mat_eig = 0.d0; !matriz diagonal com os autovalores. Matriz da formula 29
eig_vec = (0.d0,0.d0)


write(*,*)'ultimo valor de beta',beta(j),j,k10

  alpha_tmp = alpha
  beta_tmp((/1:ntemp/))  = beta((/1:ntemp/))

 CALL eigenvalues_complex(alpha_tmp,beta_tmp,eig_vec,k10) !eig_vect autovec compl e alpha_mim autoval reias da mat tridiag

 	
	eig=alpha_tmp;



write(*,*)'Calculo dos auto-valores --OK'

write(*,*)'Razao', minval(eig)/maxval(eig)
!OPEN(71,FILE = 'autovalores.dat')

!do j1=1,niter
! WRITE(71,*) eig(j1) 
!end do



!OPEN(72,FILE = 'autovalores_new.dat')
  
!eig_0 = eig(1)
!WRITE(72,*) eig_0

!do j1=2,niter
 
 !  if( eig(j1) - eig_0>0.0001)then
 !     eig_0 = eig(j1)
 !     WRITE(72,*) eig_0    
 !   end if    
!end do


!%%%%%%%%Fim do calculo dos autovalores%%%%%%%

!==========================================
!Decomposição espctral de Lanczos
!==========================================

omega = 2*pi*100 !2*pi*f
ni = 4e-7*pi

do j1=1,niter
 !           write(111,*) eig(j1)
 
 	     Mat_eig(j1,j1) = eig(j1) + imag*omega*ni     !aqui e calculado (eig+j*omega*ni*identidade)**-1
 	 !    mod_Mat_eig(j1) = (real(Mat_eig(j1,j1))**2 + aimag(Mat_eig(j1,j1))**2);
         !    Conj_Mat_eig(j1,j1) = conjg(Mat_eig(j1,j1))
          !   Inv_Mat_eig(j1,j1)  = Conj_Mat_eig(j1,j1)/mod_mat_eig(j1)
           Inv_Mat_eig(j1,j1)  = conjg(Mat_eig(j1,j1))/((real(Mat_eig(j1,j1))**2 + aimag(Mat_eig(j1,j1))**2))
            
enddo

write(*,*)'depois da mult'
do j1=1,niter
do j2=1,niter
  eig_vec_transp(j1,j2) = eig_vec(j2,j1)
enddo
enddo

!eig_vec_transp = transpose(eig_vec)


call matmul_sparse3(V,nlin,niter,eig_vec,Inv_Mat_eig,eig_vec_transp,vec_e,vec_1)

do j1=1,nlin
field(j1) = vec_1(j1)*Diag_inv(j1,j1);
end do 


field = field*mod_VecB;

OPEN(47,DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Campos_5041/100Hz',FILE = 'Campo_100_teste.dat')
	do j1=1,nlin
 		write(47,'(<nlin>(x,E20.10E3))') (field(j1)) 
	end do

 write(*,*)'Fim dos calculos'






CONTAINS

!==========================================================================================
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
    COMPLEX*8,dimension(:), allocatable :: WORK
    INTEGER :: LWORK, LIWORK                  !assumir -1
    REAL(dp), dimension(:), allocatable::RWORK
  !  Real (Kind=nag_wp)               :: rdum(1)
   ! Integer                          :: idum(1)
    INTEGER, dimension(:), allocatable :: IWORK
    INTEGER::LRWORK,INFO


lgn = ceiling(log(real(n,kind=dp))/log(2.0_dp))
LDZ = N


ALLOCATE(WORK(1), IWORK(1),rwork(1))


 lwork  = -1
 lrwork = -1
 liwork = -1

CALL zstedc (COMPZ,N,d,u,eig_vec,n,WORK,lwork,RWORK,lrwork,IWORK,liwork,INFO ) 
      lwork = max(n*n, nint(real(work(1))))
      lrwork = max(1+3*n+2*n*lgn+4*n*n, nint(rwork(1)))
      liwork = max(6+6*n+5*n*lgn, iwork(1))

DEALLOCATE(WORK, IWORK,RWORK)		

Allocate (work(lwork), rwork(lrwork), iwork(liwork))

CALL zstedc (COMPZ,N,d,u,eig_vec,n,WORK,lwork,RWORK,lrwork,IWORK,liwork,INFO ) 

DEALLOCATE(WORK, IWORK,RWORK)		

END SUBROUTINE eigenvalues_complex
!==========================================================================================

!==========================================================================================
SUBROUTINE MATMAT(Anrow,Ancol,Amatrix,Bnrow,Bncol,Bmatrix,Cmatrix)
    IMPLICIT NONE
    INTEGER(4)::Anrow,Ancol,Bnrow,Bncol,i,j
    COMPLEX(8)   ::Amatrix(Anrow,Ancol),Bmatrix(Bnrow,Bncol),Cmatrix(Anrow,Bncol)

  !  DO j=1,Bncol
  !     DO i=1,Anrow
          Cmatrix = 0.d0
  !     END DO
  !  END DO
    DO j=1,Bncol
       CALL GAXPY(Anrow,Ancol,Amatrix,Bmatrix(:,j),Cmatrix(:,j))
    END DO
END SUBROUTINE MATMAT
!==========================================================================================

!==========================================================================================
SUBROUTINE MATMAT2(Anrow,Ancol,Amatrix,Bnrow,Bncol,Bmatrix,Cmatrix)
    IMPLICIT NONE
    INTEGER(4)::Anrow,Ancol,Bnrow,Bncol,i,j
    COMPLEX(8)   ::Amatrix(Anrow,Ancol),Bmatrix(Bnrow,Bncol),Cmatrix(Anrow,Bncol)

          Cmatrix = 0.d0  
    DO j=1,Bncol
         DO i=1,Ancol
          DO k=1,Anrow
             Cmatrix(k,j) = Cmatrix(k,j) + Bmatrix(i,j)*Amatrix(k,i)
          END DO
         END DO
    END DO
END SUBROUTINE MATMAT2
!==========================================================================================

!==========================================================================================
  SUBROUTINE GAXPY(Nrow,Ncol,Amatrix,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,Ncol,i,j
    COMPLEX(8)   ::Amatrix(Nrow,Ncol),xvector(Ncol),yvector(Nrow),x_i

    DO i=1,Ncol
       x_i = xvector(i)
       CALL SAXPY(Nrow,x_i,Amatrix(:,i),yvector)
    END DO
  END SUBROUTINE GAXPY
!==========================================================================================
  
!==========================================================================================  
  SUBROUTINE SAXPY(Nrow,alpha,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,i
    COMPLEX(8)   ::alpha,xvector(Nrow),yvector(Nrow)

    DO i=1,Nrow
       yvector(i) = yvector(i) + alpha*xvector(i)
    END DO

  END SUBROUTINE SAXPY
!==========================================================================================

!========================================================================================== 
  SUBROUTINE MATVEC(Nrow,Ncol,Amatrix,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,Ncol,i,j
    COMPLEX(8)::Amatrix(Nrow,Ncol),yvector(Nrow),x_i
    REAL(8):: xvector(Ncol)

  !  DO i=1,Nrow
       yvector = 0.D0
  !  END DO

    CALL GAXPY2(Nrow,Ncol,Amatrix,xvector,yvector)

  END SUBROUTINE MATVEC
!==========================================================================================

!==========================================================================================
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
!==========================================================================================

!==========================================================================================
  SUBROUTINE SAXPY2(Nrow,alpha,xvector,yvector)
    IMPLICIT NONE
    INTEGER(4)::Nrow,i
    COMPLEX(8)   ::xvector(Nrow),yvector(Nrow)
    REAL(8)::alpha

    DO i=1,Nrow
       yvector(i) = yvector(i) + alpha*xvector(i)
    END DO

  END SUBROUTINE SAXPY2
!==========================================================================================


!================================================================================
! Subrotinas para utilização de matrizes esparsas
! As subrotinas abaixo multiplica duas matrizes esparsas, e as converte em matrizes
! no formato  'compressed row format' também conhecido como Harwell-Boeing format,

! As[]   -> armazena os valores não-nulos da matriz A
! asub[] -> armazena os indices da linha a qual os valores pertencem
! xa[]   -> armazena os indices que indicam o inicio de cada coluna
!           considerando o posicionamento desses valores no vetor As[]
! Em seguida realiza a multiplicação entre elas.

!na,ma,nnza --> rows A, columns A, noz zeroes A
! Autor Diego C. Miranda     28 jan 2020



!==========================================================================================
subroutine mat_sparse(matA,na,ma,rowA,ColPtrA,ValA,nnz)
Implicit None
!Reescreve a matriz esparsa
! matA(in)    - matriz A cheia
! na(in)      - numero de linhas de A
! ma(in)      - numero de colunas de A
!rowA(out)    -
!ColPtrA(out) - 
!ValA(out)    - 
!nnzA(out)    - 
! !-----------------------------------------------------

Integer(4),Intent(in)  :: na,ma
COMPLEX(8),Intent(in):: matA(na,ma)

Integer(4),Intent(in):: nnz
COMPLEX(8),Intent(out),Allocatable:: valA(:)
Integer(4),Intent(out),Allocatable :: rowA(:),colPtrA(:)
Integer(4)	            :: flag, ierr

Integer(4)             :: i,j,k,nnza



nnza=nnz
!compressed row storage -----------
Allocate(valA(nnza),rowA(nnza),colPtrA(na+1))

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

end subroutine
!==========================================================================================


!=======================================================================================
subroutine matmul_sparse2(na,ma,rowA,ColPtrA,ValA,nnz,matb,nb,mb,matC)
Implicit None

! matA - matriz A cheia
! na   - numero de linhas de A
! ma   - numero de colunas de A
! matB - matriz ou vetor B
! nb   - numero de linhas de B
! mb   - numero de colunas de B
!-----------------------------------------------------
Integer(4),Intent(in)  :: na,ma,nb,mb

Integer(4),Intent(in):: nnz
COMPLEX(8),Intent(in):: valA(nnz)
Integer(4),Intent(in):: rowA(nnz),colPtrA(na+1)
 
COMPLEX(8),Intent(in)::matB(nb,mb)

Integer(4)                 :: nnza,nnzb
COMPLEX(8),Allocatable    ::valB(:),c(:)
Integer(4),Allocatable     :: rowB(:),colPtrB(:),jc(:),ic(:)
Integer(4),Allocatable     :: iwa(:),iwb(:),ndegr(:)
Integer(4)		   :: job,nzmax
Integer(4)	           :: flag, ierr

Integer(4)                 :: i,j,k,j4

COMPLEX(8),Intent(out)    :: matC(na,mb)

nnzb = count(matB /= 0)
nnza = nnz


MatC = 0.d0;

!compressed row storage -----------
Allocate(valB(nnzb),rowB(nnzb),colPtrB(nb+1))


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

 Allocate(c(nzmax),jc(nzmax),ic(na+1))   
job = 1


 Call amub(na, ma, job, valA, rowA, colPtrA, valB, rowB, colPtrB, c, jc, ic, nzmax,  iwa, ierr,nnza,nnzb)




 Call csrdns ( na, mb, c, jc, ic, nzmax, matC, ierr )


end subroutine
!==========================================================================================
subroutine amub (nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
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

  Complex(8)::  a(nnza)
  Complex(8) :: b(nnzb)
  Complex(8) :: c(nzmax)
  integer :: ia(nrow+1)
  integer :: ib(ncol+1)
  integer :: ic(nrow+1)
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
  integer::  k,k2,k3,k4
  integer::  ka
  integer::  kb
  integer::  len
  Complex(8)::  scal
  logical:: values

  values = ( job /= 0 )
  len   = 0
  ic(1) = 1
  ierr  = 0
  k4    = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow

  		if(ia(ii+1)-1==-1)then
      			k2=ia(ii)
        
		else
      			k2=ia(ii+1)-1
 		endif
!  Row I.

		if(ia(ii).ne.0)then
 				do ka = ia(ii), k2
       					
        				if ( values ) then
              					scal = a(ka)
                                                 
      	        			end if
        				jj = ja(ka)
                                           
	   				if(ib(jj+1)-1==-1)then
      	       					k3=ib(jj)
   	   				else
      	      					k3=ib(jj+1)-1
 	   				endif
                                               
           				if(ib(jj).ne.0)then
      	      					do kb = ib(jj), k3
                                                            
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

               							k4 =k4+1
               							
           				        end do
    					endif
   			     end do
	
    			do k = ic(ii), len
        			iw(jc(k)) = 0
   		        end do

    			ic(ii+1) = len + 1
	endif
end do
  return
end subroutine

!==========================================================================================

!==========================================================================================

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
  integer :: k,k2,k3
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
if(ia(ii+1)-1==-1)then
      k2=ia(ii)
else
      k2=ia(ii+1)-1
 endif

  if(ia(ii).ne.0)then

  	do j = ia(ii), k2
			!
			!  Row number to be added.
			!
        		jr = ja(j)
                        if(ib(jr+1)-1==-1)then
      				k3=ib(jr)
   			else
      				k3=ib(jr+1)-1
 			endif
			if(ib(jr).ne.0)then
        			do k = ib(jr), k3
           					jc = jb(k)
          					!  Add one element to the linked list.
          					if ( iw(jc) == 0 ) then
              						ldg = ldg + 1
              						iw(jc) = last
              						last = jc
            					end if
         			end do
			end if
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


  endif

end do

  nnz = sum ( ndegr(1:nrow) )

  return
end subroutine


!==========================================================================================


!==========================================================================================
subroutine csrdns ( nrow, ncol, a, ja, ia, nzmax, dns, ierr2 )

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

  complex(8) :: a(nzmax)
  complex(8) ::dns(nrow,ncol)
  integer :: i
  integer :: ia(nrow+1)
  integer :: ierr2
  integer :: j
  integer :: ja(nzmax)
  integer :: k
  


  ierr2 = 0
  dns(1:nrow,1:ncol) = 0.d0

  do i = 1, nrow


    do k = ia(i), ia(i+1)-1
     j = ja(k)
      if ( ncol < j ) then
        ierr2 = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end

!==========================================================================================

!==========================================================================================
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
Complex(8),Intent(in)     :: matA(na,ma)
Complex(8), Intent(in)    :: matB(nb,mb)

Integer             :: nnza,nnzb
Complex(8),Allocatable  :: valA(:)
Complex(8),Allocatable:: valB(:),c(:)
Integer,Allocatable :: rowA(:),rowB(:),colPtrA(:),colPtrB(:),jc(:),ic(:)
Integer,Allocatable :: iwa(:),iwb(:),ndegr(:)
Integer		    :: job,nzmax
Integer	            :: flag, ierr

Integer             :: i,j,k

Complex(8),Intent(out)    :: matC(na,mb)


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


Do i = 1,na+1
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


Do i = 1,nb+1
    If (colPtrb(i) == 0) colPtrB(i) = colPtrB(i+1)
End Do



Allocate(ndegr(na),iwb(mb),iwa(ma))

 Call amubdg(na,ma,mb,rowA,colPtrA,rowB,colPtrB,ndegr,nzmax,iwb,nnza,nnzb)  


Allocate(c(nzmax),jc(nzmax),ic(na+1))   !ic(ma+1))  !
job = 1


! multiplicacao matricial esparsa

 Call amub(na, ma, job, valA, rowA, colPtrA, valB, rowB, colPtrB, c, jc, ic, nzmax,  iwa, ierr,nnza,nnzb)

Do i = 1,na+1
    If (ic(i) == 0) then
	      if (ic(i+1)==0)then
	      
			ic(i) = ic(i-1)
	      else
                          ic(i) = ic(i+1)
	      endif
     endif
End Do


 Call csrdns ( na, mb, c, jc, ic, nzmax, matC, ierr )

end subroutine
!==========================================================================================


!============================================================================================
!=============================================================================================

subroutine matmul_sparse3(matA,na,ma,matB,matC,matD,matE,matR)
Implicit None

! matA - V
!matmul_sparse3(V,nlin,niter,eig_vec,Inv_Mat_eig,eig_vec_transp,vec_e,vec_1)
!MatB - eig_vec
!MatC - Inv_mat_eig
!MatD - eig_vec_trans
!MatE - vece

! na   - nlin
! ma   - ninter
!Mat_temp1 - Mat1

! matB - matriz ou vetor B
! nb   - numero de linhas de B
! mb   - numero de colunas de B
!-----------------------------------------------------

Integer,Intent(in)  :: na,ma
Complex(8),Intent(in)     :: matA(na,ma),matB(ma,ma),MatC(ma,ma),MatD(ma,ma),MatE(na,1)
!Real(8),Intent(in)        :: MatE(na,1)

Integer             :: nnza,nnzb,nnzc,nnzd,nnze
Complex(8),Allocatable  :: valA(:),valB(:),valC(:),valD(:),valE(:),valR(:),valtemp1(:),valtemp2(:),valtemp3(:),valtemp4(:) 


Integer,Allocatable :: rowA(:),rowB(:),colPtrA(:),colPtrB(:),rowC(:),rowD(:),colPtrC(:),colPtrD(:),rowE(:),colPtrE(:)
Integer,Allocatable :: rowR(:),colPtrR(:),jtemp1(:),itemp1(:),jtemp2(:),itemp2(:),jtemp3(:),itemp3(:),jtemp4(:),itemp4(:)
Integer,Allocatable :: iwa(:),iwb(:),ndegr(:),iwb1(:)
Integer		    :: job,nzmax,nb,mb,nc,mc,nd,md,ne,me,nnz3,nnz1,nnz2,nnztemp
Integer	            :: flag, ierr

Integer             :: i,j,k

Complex(8),Intent(out)    :: matR(na,1)

nnza = count(matA /= 0)
nnzb = count(matB /= 0)
nnzc = count(matC /= 0)
nnzd = count(matD /= 0)
nnze = count(matE /= 0)


!compressed row storage -----------

Allocate(valA(nnza),valB(nnzb),rowA(nnza),rowB(nnzb),colPtrA(na+1),colPtrB(ma+1))
Allocate(valC(nnzc),valD(nnzd),rowC(nnzc),rowD(nnzd),colPtrC(ma+1),colPtrD(ma+1))
Allocate(valE(nnze),rowE(nnz3),colPtrE(ma+1))

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


Do i = 1,na+1
    If (colPtrA(i) == 0) colPtrA(i) = colPtrA(i+1)
End Do


!convertendo Matriz B para crs
nb = ma
mb = ma
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


Do i = 1,nb+1
    If (colPtrb(i) == 0) colPtrB(i) = colPtrB(i+1)
End Do


!Convertend a matriz C para crs

colPtrC = 0
nc=ma
mc=ma
k = 1
Do i = 1,nc
	flag = 0
	Do j = 1,mc
		If (matC(i,j) /= 0) Then
			valC(k) = matC(i,j)
			rowC(k) = j     !   linha variando mais rapido
			If (flag == 0) Then
				colPtrC(i) = k				
				flag = 1
			End If
			k = k + 1
		End If
	End Do
End Do
colPtrC(nc+1) = nnzc + 1


Do i = 1,nc+1
    If (colPtrC(i) == 0) colPtrC(i) = colPtrC(i+1)
End Do



!convertendo Matriz D para crs
colPtrD = 0
k = 1
nd=ma
md=ma
Do i = 1,nd
	flag = 0
	Do j = 1,md
		If (matD(i,j) /= 0) Then
			valD(k) = matD(i,j)
			rowD(k) = j     !   linha variando mais rapido
			If (flag == 0) Then
				colPtrD(i) = k				
				flag = 1
			End If
			k = k + 1
		End If
	End Do
End Do
colPtrD(nd+1) = nnzd + 1


Do i = 1,nd+1
    If (colPtrD(i) == 0) colPtrD(i) = colPtrD(i+1)
End Do



!convertendo Matriz E para crs (e vetor)
colPtrE = 0
k = 1
ne=ma
Do i = 1,ne
	flag = 0
	!Do j = 1,ma
		If (matE(i,1) /= 0) Then
			valE(k) = matE(i,1)
			rowE(k) = 1     !   linha variando mais rapido
			If (flag == 0) Then
				colPtrE(i) = k				
				flag = 1
			End If
			k = k + 1
		End If
	!End Do
End Do
colPtrE(ne+1) = nnze + 1


Do i = 1,ne+1
    If (colPtrE(i) == 0) colPtrE(i) = colPtrE(i+1)
End Do


!Alocacoes de memorias para cada multiplicacao
!______________________________
!______________________________
!Primeira multiplicação
!______________________________
!______________________________

Allocate(ndegr(na),iwb(mb),iwa(ma))

 Call amubdg(na,ma,mb,rowA,colPtrA,rowB,colPtrB,ndegr,nzmax,iwb,nnza,nnzb)  


Allocate(valtemp1(nzmax),jtemp1(nzmax),itemp1(na+1))   !aqui c é mat_temp1(na.ma)
job = 1

! multiplicacao matricial esparsa

!amub(na, ma, job, valA, rowA, colPtrA, valB, rowB, colPtrB, c, jc, ic, nzmax,  iwa, ierr,nnza,nnzb)



 Call amub(na, ma, job, valA, rowA, colPtrA, valB, rowB, colPtrB, valtemp1, jtemp1, itemp1, nzmax,  iwa, ierr,nnza,nnzb)

Do i = 1,na+1
    If (itemp1(i) == 0) then
	      if (itemp1(i+1)==0)then
	      
			itemp1(i) = itemp1(i-1)
	      else
                          itemp1(i) = itemp1(i+1)
	      endif
     endif
End Do

!______________________________
!______________________________
!Segunda multiplicação                 !c=Mat_temp1
!______________________________
!______________________________
nnztemp=0
nnztemp =  count(valtemp1 /= 0)

!Allocate(ndegr(na),iwb(mc),iwa(ma))

 Call amubdg(na,ma,mc,jtemp1,itemp1,rowC,colPtrC,ndegr,nzmax,iwb,nnztemp,nnzc)  


Allocate(valtemp2(nzmax),jtemp2(nzmax),itemp2(na+1))   !aqui c2 é mat_temp2(na,ma)
job = 1

! multiplicacao matricial esparsa

 Call amub(na, ma, job, valtemp1, jtemp1, itemp1, valC, rowC, colPtrC, valtemp2, jtemp2, itemp2, nzmax,  iwa, ierr,nnztemp,nnzc)

Do i = 1,na+1
    If (itemp2(i) == 0) then
	      if (itemp2(i+1)==0)then
	      
			itemp2(i) = itemp2(i-1)
	      else
                          itemp2(i) = itemp2(i+1)
	      endif
     endif
End Do

!______________________________
!______________________________
!Terceira multiplicação                 !c1=Mat_temp2
!______________________________
!______________________________

nnztemp=0

nnztemp =  count(valtemp2 /= 0)


!Allocate(ndegr(na),iwb(md),iwa(ma))

 Call amubdg(na,ma,md,jtemp2,itemp2,rowD,colPtrD,ndegr,nzmax,iwb,nnztemp,nnzd)  


Allocate(valtemp3(nzmax),jtemp3(nzmax),itemp3(na+1))   !aqui c2 é mat_temp2(na,ma)
job = 1

! multiplicacao matricial esparsa

 Call amub(na, ma, job, valtemp2, jtemp2, itemp2, valD, rowD, colPtrD, valtemp3, jtemp3, itemp3, nzmax,  iwa, ierr,nnztemp,nnzd)

Do i = 1,na+1
    If (itemp3(i) == 0) then
	      if (itemp3(i+1)==0)then
	      
			itemp3(i) = itemp3(i-1)
	      else
                          itemp3(i) = itemp3(i+1)
	      endif
     endif
End Do

!______________________________
!______________________________
!Quarta multiplicação                 !c1=Mat_temp2
!______________________________
!______________________________

nnztemp=0
nnztemp =  count(valtemp3 /= 0)

me=1

Allocate(iwb1(me))

 Call amubdg(na,ma,me,jtemp3,itemp3,rowE,colPtrE,ndegr,nzmax,iwb1,nnztemp,nnze)  
 


Allocate(valtemp4(nzmax),jtemp4(nzmax),itemp4(na+1))   !aqui c2 é mat_temp2(na,ma)
job = 1


! multiplicacao matricial esparsa

 Call amub(na, ma, job, valtemp3, jtemp3, itemp3, valE, rowE, colPtrE, valtemp4, jtemp4, itemp4, nzmax,  iwa, ierr,nnztemp,nnze)

Do i = 1,na+1
    If (itemp4(i) == 0) then
	      if (itemp4(i+1)==0)then
	      
			itemp4(i) = itemp4(i-1)
	      else
                          itemp4(i) = itemp4(i+1)
	      endif
     endif
End Do

nnztemp =  count(valtemp4 /= 0)


 Call csrdns ( na, 1, valtemp4, jtemp4, itemp4, nzmax, matR, ierr )

end subroutine
!==========================================================================================

!==========================================================================================
subroutine amub2(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
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

  Complex(8)::  a(nnza)
  Complex(8) :: b(nnzb)
  Complex(8) :: c(nzmax)
  integer :: ia(nrow+1)
  integer :: ib(ncol+1)
  integer :: ic(nrow+1)
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
  integer::  k,k2,k3,k4
  integer::  ka
  integer::  kb
  integer::  len
  Complex(8)::  scal
  logical:: values

  values = ( job /= 0 )
  len   = 0
  ic(1) = 1
  ierr  = 0
  k4    = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow

  		if(ia(ii+1)-1==-1)then
      			k2=ia(ii)
        
		else
      			k2=ia(ii+1)-1
 		endif
!  Row I.

		if(ia(ii).ne.0)then
 				do ka = ia(ii), k2
       					
        				if ( values ) then
              					scal = a(ka)
                                                 
      	        			end if
        				jj = ja(ka)
                                           
	   				if(ib(jj+1)-1==-1)then
      	       					k3=ib(jj)
   	   				else
      	      					k3=ib(jj+1)-1
 	   				endif
                                               
           				if(ib(jj).ne.0)then
      	      					do kb = ib(jj), k3
                                                            
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

               							k4 =k4+1
               							
           				        end do
    					endif
   			     end do
	
    			do k = ic(ii), len
        			iw(jc(k)) = 0
   		        end do

    			ic(ii+1) = len + 1
	endif

end do
  return
end subroutine


  
END PROGRAM SDLM_last

