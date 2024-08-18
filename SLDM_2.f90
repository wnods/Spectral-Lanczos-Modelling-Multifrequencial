PROGRAM SDLM_2
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
!01/04/2020 Versao do prgrama que roda os dados so anderson. Esta versão não esta otimizada
!10/12/20 Ultima atualização
  !=====================================================================================================
  IMPLICIT NONE
  REAL(8),PARAMETER:: PI      = 3.141592653589793238462643383279502884197d0
  REAL(8),ALLOCATABLE,DIMENSION(:,:):: Stiff_Mat, Diag, Mass_Mat,Imat, MatA,Diag_inv,Hess_Mat
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:):: V, Mat_temp,eig_vec,eig_vec_transp,Mat_eig,Inv_Mat_eig,Conj_Mat_eig,Mat_source,&
                                          Mat_1,Mat_2,Mat_3
  REAL(8),ALLOCATABLE,DIMENSION(:)::Source_Vec_real, Source_Vec_imag, mod_Mat_eig,vec_e
  REAL(8),ALLOCATABLE,DIMENSION(:)::alpha,alpha_tmp,beta, beta_tmp
  REAL(4),ALLOCATABLE,DIMENSION(:)::eig,freq
  COMPLEX*16,ALLOCATABLE,DIMENSION(:)::Source_Vec, VecB , Vec_temp, Vec_res,V_0,vec_1,field
  COMPLEX(8)::imag, t1, t2,t3
  REAL(8):: freq_base,beta_0,mod_VecB,sum_VecB,eig_0,omega,ni,mod_field,soma
  INTEGER(4)::nlin,niter,ntemp,ncol,ierr,j,j1,j2,k1,k2,k3,k4,k10,right_index,k,left_index,opt,nfreq

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

!Mat Stiffness (rigidez)
write(*,*)'Entre com a dimensão das matrizes de Stiffness e Mass'
read(*,*)nlin
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

!Faixa de frequências

write(*,*)'Entre com número de frequências a ser utilizado na decomposição de lanczos'
read(*,*)nfreq
write(*,*)'Entre com as frequências a serem utilizadas na decomposição de lanczos'

do j1=1,nfreq
read(*,*)freq(j2)
end do

write(*,*)'Entre com a frequência base'
read(*,*)freq_base


   write(*,*)'Entre com nome do arquivo da matriz de Stiffness'
   read(*,*)  fname1       !  'Stiff.dat'    

   write(*,*)'Entre com nome do arquivo da matriz de massa'
   read(*,*) fname3   !'Mass.dat'


 write(*,*)'Entre com arquivo da parte real do vetor fonte'
 read(*,*) fname5 

 write(*,*)'Entre com arquivo da parte imag do vetor fonte'
 read(*,*) fname6 



!Endereço das pastas onde estão or arquivos. Essa parte deve ser atualizada!!

open(15, DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Artigo1/Modelo2_f001_3625/1000hz', FILE=fname5,status ='old')

open(25, DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Artigo1/Modelo2_f001_3625/1000hz', FILE=fname6,status ='old')
   

 open(11, DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Artigo1/Modelo2_f001_3625', FILE=fname3,status ='old')

 open(10, DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Artigo1/Modelo2_f001_3625', FILE=fname1,status ='old')
  	
   do j1=1,nlin
       read(10,*)(Stiff_Mat(j1,j2),j2=1,ncol)  
       read(11,*)(Mass_Mat(j1,j2),j2=1,ncol) 
       read(15,*)Source_Vec_real(j1)
       read(25,*)Source_Vec_imag(j1) 
 end do



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



!-----------------------------------------------------------------------
!Alocação de memório para a recurssão de Laczos
!V base ortornormal, dim=niterxniter;
!V_0 base ortornormal de indice 0;
!Vec_temp, Vec_res variaveis usadas durante a recurssao, dim=nlin;
!alpha e beta, produtos do algoritmo de recurssao, forma a matriz tridiagonal, dim=niter;
!Hess_Mat, produto da recurssao, dim=niterxniter;
!-----------------------------------------------------------------------

   Write(*,*)'Entre com o número de iterações'
   read(*,*) niter    != 2000!793!numero de iterações calculado do nuremo de autovetores nao iguais (difereca entre eles <0.0001)



  ALLOCATE(V(nlin,niter+1),STAT=ierr)
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


  ALLOCATE(Mat_temp(nlin,nlin),STAT=ierr)
  IF(ierr)THEN
  STOP 'Nao foi possivel allocar memoria para Mat_temp'
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
              write(500,*)sqrt(Diag(j1,j1))
 		Diag_inv(j1,j1) = 1.d0/sqrt(Diag(j1,j1));  !D^-1/2
                VecB(j1)= Diag_inv(j1,j1)*Source_Vec(j1);
                 if(Diag(j1,j1)==0.d0)then
                   write(*,*)'linha nula',j1,Diag(j1,j1)
		endif
                
	end do



           do j1=1,nlin
	             do j2=1,ncol
	                     MatA (j1,j2) =  Diag_inv(j1,j1)*Stiff_Mat(j1,j2)*Diag_inv(j2,j2);
		     end do
             end do



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
 

k10=0

 DO WHILE (j.ne.niter)
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
        write(100,*)j,beta(j),alpha(j)
  end if
      	k10=j  !numero de iteraçoes para a determinacao da matriz tridiagonal
!	write(*,*)'iter',k10

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

 !WRITE(71,*) eig(j1) 
!end do



!OPEN(72,FILE = 'autovalores_new.dat')
  
!eig_0 = eig(1)
!WRITE(72,*) eig_0

!do j1=2,niter
 
 !  if( eig(j1) - eig_0>0.0001)then
  !    eig_0 = eig(j1)
   !   WRITE(72,*) eig_0    
   ! end if    
!end do


!%%%%%%%%Fim do calculo dos autovalores%%%%%%%

!==========================================
!Decomposição espctral de Lanczos
!==========================================

do j2=1,nfreq
omega = 2*pi*freg(j2)!2*pi*f
ni = 4e-7*pi

do j1=1,niter
          
 	     Mat_eig(j1,j1) = eig(j1) + imag*omega*ni     !aqui e calculado (eig+j*omega*ni*identidade)**-1
 	     mod_Mat_eig(j1) = (real(Mat_eig(j1,j1))**2 + aimag(Mat_eig(j1,j1))**2);
             Conj_Mat_eig(j1,j1) = conjg(Mat_eig(j1,j1))
             Inv_Mat_eig(j1,j1)  = Conj_Mat_eig(j1,j1)/mod_mat_eig(j1)
            
enddo

do j1=1,niter
do j2=1,niter
  eig_vec_transp(j1,j2) = eig_vec(j2,j1)
enddo
enddo

vec_e = 0.d0;
vec_e(1) = 1.d0;

!mult matrizes
vec_1 = 0.d0;

CALL MATMAT2(nlin,niter,V,&
            niter,niter,eig_vec,&
            mat_1)


write(*,*)'ok1'
       
CALL MATMAT2(nlin,niter,mat_1,&
            niter,niter,Inv_Mat_eig,&
            mat_2)

write(*,*)'ok2'

vec_1 = 0.d0;

CALL MATMAT2(nlin,niter,Mat_2,&
            niter,niter,eig_vec_transp,&
            mat_3)

write(*,*)'ok3'
!multiplacacao de matriz or vetor

CALL MATVEC(nlin,niter,Mat_3,vec_e,vec_1)
 write(*,*)'ok4'


do j1=1,nlin
field(j1) = vec_1(j1)*Diag_inv(j1,j1);
end do 
field = field*mod_VecB;

write(*,*)'Entre com o nome do arquivo do campo calculado', freq(j2)

read(*,*)fname2   !Tem que ajeitar aqui Campo_1000_001new.dat
!Caminho da pasta onde o campo caluculado está sendo guardado. Deve ser atualizado!!!
OPEN(47,DEFAULTFILE='~/Área de Trabalho/PEROBA/estudo_19_04-19/Programs/Dados_Anderson/Artigo1/Modelo2_f001_3625/plot',FILE = fname2)
	do j1=1,nlin
 		write(47,'(<nlin>(x,E20.10E3))') (field(j1)) 
	end do
end do

 write(*,*)'Fim dos calculos'

 CALL CPU_TIME (time_end)

 write(*,*) 'Time of operation was ', time_end - time_begin, ' seconds'




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

  !  DO i=1,Nrow
       yvector = 0.D0
  !  END DO

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

  
END PROGRAM SDLM_2

