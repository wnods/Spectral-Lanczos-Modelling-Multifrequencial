PROGRAM SDLM_last
  IMPLICIT NONE
  
  ! Definição de constantes e variáveis
  REAL(8), PARAMETER :: PI = 3.141592653589793238462643383279502884197d0
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Stiff_Mat, Diag, Mass_Mat, Imat, MatA, Diag_inv, Hess_Mat
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: V, Mat_temp, eig_vec, eig_vec_transp, Mat_eig, Inv_Mat_eig, Conj_Mat_eig, &
                                             Mat_1, Mat_2, Mat_3
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Source_Vec_real, Source_Vec_imag, mod_Mat_eig, vec_e
  REAL(8), ALLOCATABLE, DIMENSION(:) :: alpha, alpha_tmp, beta, beta_tmp
  REAL(4), ALLOCATABLE, DIMENSION(:) :: eig
  COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: Source_Vec, VecB, Vec_temp, Vec_res, V_0, vec_1, field
  COMPLEX(8) :: imag, t1, t2, t3
  REAL(8) :: freq, beta_0, mod_VecB, sum_VecB, eig_0, omega, ni, mod_field
  INTEGER(4) :: nlin, niter, ntemp, ncol, ierr, j, j1, j2, k1, k2, k3, k4, k10, right_index, k, left_index, opt
  CHARACTER(len=64) :: fname, fname1, fname3, fname5, fname6
  REAL :: time_begin, time_end, rnd
  
  ! Novas variáveis ALLOCATABLE
  COMPLEX(8), ALLOCATABLE :: Vec_2(:), V_1(:,:), Mat_4(:,:)
  REAL(8), ALLOCATABLE :: mod_diag(:)
  
  ! Início do cálculo do tempo de execução
  CALL CPU_TIME(time_begin)
  
  ! Definindo a constante imaginária
  imag = (0.d0, 1.d0)

  ! Configuração do tamanho das matrizes
  nlin = 1521
  ncol = nlin
  
  ! Alocação e leitura das matrizes Stiffness e Mass
  ALLOCATE(Stiff_Mat(nlin, ncol), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Stiff_Mat'
  
  fname1 = 'Stiff.dat'
  OPEN(10, FILE=fname1, STATUS='old')
  DO j1 = 1, nlin
      READ(10, *) Stiff_Mat(j1, :)
  END DO
  CLOSE(10)
  
  ALLOCATE(Mass_Mat(nlin, ncol), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mass_Mat'
  
  fname3 = 'Mass.dat'
  OPEN(11, FILE=fname3, STATUS='old')
  DO j1 = 1, nlin
      READ(11, *) Mass_Mat(j1, :)
  END DO
  CLOSE(11)
  
  ! Alocação e leitura do vetor de fonte
  ALLOCATE(Source_Vec_real(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Source_Vec_real'
  
  ALLOCATE(Source_Vec_imag(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Source_Vec_imag'
  
  fname5 = 'realB100.dat'
  fname6 = 'imagB100.dat'
  OPEN(15, FILE=fname5, STATUS='old')
  OPEN(25, FILE=fname6, STATUS='old')
  DO j1 = 1, nlin
      READ(15, *) Source_Vec_real(j1)
      READ(25, *) Source_Vec_imag(j1)
  END DO
  CLOSE(15)
  CLOSE(25)
  
  ALLOCATE(Source_Vec(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Source_Vec'
  
  Source_Vec = Source_Vec_real + imag * Source_Vec_imag
  
  ! Alocação de memória para variáveis usadas na recursão de Lanczos
  ALLOCATE(Diag(nlin, nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Diag'
  
  ALLOCATE(Diag_inv(nlin, nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Diag_inv'
  
  ALLOCATE(VecB(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para VecB'
  
  ALLOCATE(Imat(nlin, ncol), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Imat'
  
  ALLOCATE(MatA(nlin, ncol), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para MatA'
  
  niter = nlin
  
  ALLOCATE(V(nlin, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para V'
  
  ALLOCATE(V_0(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para V_0'
  
  ALLOCATE(Vec_temp(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Vec_temp'
  
  ALLOCATE(Mat_temp(nlin, nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mat_temp'
  
  ALLOCATE(alpha(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para alpha'
  
  ALLOCATE(beta(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para beta'
  
  ALLOCATE(Vec_res(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Vec_res'
  
  ALLOCATE(Hess_Mat(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Hess_Mat'
  
  ntemp = niter - 1
  ALLOCATE(alpha_tmp(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para alpha_tmp'
  
  ALLOCATE(beta_tmp(ntemp), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para beta_tmp'
  
  ALLOCATE(eig_vec(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para eig_vec'
  
  ALLOCATE(eig(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para eig'
  
  ALLOCATE(Mat_eig(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mat_eig'
  
  ALLOCATE(mod_Mat_eig(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para mod_Mat_eig'
  
  ALLOCATE(Conj_Mat_eig(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Conj_Mat_eig'
  
  ALLOCATE(Inv_Mat_eig(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Inv_Mat_eig'
  
  ALLOCATE(eig_vec_transp(niter, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para eig_vec_transp'
  
  ALLOCATE(Mat_1(nlin, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mat_1'
  
  ALLOCATE(Mat_2(nlin, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mat_2'
  
  ALLOCATE(Mat_3(nlin, niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para Mat_3'
  
  ALLOCATE(vec_e(niter), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para vec_e'
  
  ALLOCATE(field(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para field'
  
  ALLOCATE(vec_1(nlin), STAT=ierr)
  IF (ierr /= 0) STOP 'Nao foi possivel allocar memoria para vec_1'
  
  ! Cálculos principais
  Diag = 0.d0
  Diag_inv = 0.d0
  VecB = 0.d0
  Imat = 0.d0
  mod_VecB = 0.d0
  Mat_temp = 0.d0
  MatA = 0.d0
  V = 0.d0

  DO j1 = 1, nlin
      Diag(j1, j1) = SUM(Mass_Mat(j1, :))
      Imat(j1, j1) = 1.d0      
      Diag_inv(j1, j1) = 1.d0 / SQRT(Diag(j1, j1))
      IF (Diag(j1, j1) == 0.d0) THEN
          WRITE(*, *) j1, Diag_inv(j1, j1)
          
      END IF
      VecB(j1) = Diag_inv(j1, j1) * Source_Vec(j1)
  END DO

  DO j1 = 1, nlin
      DO j2 = 1, ncol
          MatA(j1, j2) = Diag_inv(j1, j1) * Stiff_Mat(j1, j2) * Diag_inv(j2, j2)
      END DO
  END DO

  mod_VecB = SQRT(DOT_PRODUCT(VecB, VecB))
  WRITE(*, *) 'mod', mod_VecB
  
  V(:, 1) = VecB / mod_VecB

  ! Recursão de Lanczos
  k1 = 1
  beta = 0.d0
  alpha = 0.d0
  beta_0 = 0.d0
  V_0 = 0.d0
  j = 1

  DO WHILE (j /= niter)
      IF (k1 == 1) THEN
          Vec_temp = MATMUL(MatA, V(:, j))
          Mat_temp(1, :) = V(:, j)
          alpha(j) = DOT_PRODUCT(Mat_temp(1, :), Vec_temp(:))
          Vec_res(:) = Vec_temp(:) - beta_0 * V_0(:) - alpha(j) * V(:, j)
          beta(j) = SQRT(DOT_PRODUCT(Vec_res(:), Vec_res(:)))
          k1 = j + 1
          V(:, k1) = Vec_res(:) / beta(j)
      ELSE
          j = j + 1
          Vec_temp(:) = MATMUL(MatA(:,:), V(:, j))
          Mat_temp(1, :) = V(:, j)
          alpha(j) = DOT_PRODUCT(Mat_temp(1, :), Vec_temp(:))
          k2 = j - 1
          Vec_res(:) = Vec_temp(:) - beta(k2) * V(:, k2) - alpha(j) * V(:, j)
          beta(j) = SQRT(DOT_PRODUCT(Vec_res(:), Vec_res(:)))
          k3 = j + 1

          IF (beta(j) == 0) THEN
              WRITE(*, *) 'beta nulo'
              STOP
          ELSE
              V(:, k3) = Vec_res(:) / beta(j)
          END IF
          WRITE(100, *) j, beta(j), alpha(j)
      END IF
      k10 = j  ! número de iterações para a determinação da matriz tridiagonal
      WRITE(*, *) 'iter', k10
  END DO

  ! Cálculo dos autovalores
  Mat_eig = 0.d0
  eig_vec = (0.d0, 0.d0)

  WRITE(*, *) 'último valor de beta', beta(j), j, k10

  alpha_tmp = alpha
  beta_tmp(1:ntemp) = beta(1:ntemp)

  CALL eigenvalues_complex(alpha_tmp, beta_tmp, eig_vec, k10)

  eig = alpha_tmp

  WRITE(*, *) 'Cálculo dos autovalores --OK'
  WRITE(*, *) 'Razão', MINVAL(eig) / MAXVAL(eig)
  OPEN(71, FILE = 'autovalores.dat')

  DO j1 = 1, niter
      WRITE(71, *) eig(j1)
  END DO
  CLOSE(71)

  OPEN(72, FILE = 'autovalores_new.dat')

  eig_0 = eig(1)
  WRITE(72, *) eig_0

  DO j1 = 2, niter
      IF (eig(j1) - eig_0 > 0.0001) THEN
          eig_0 = eig(j1)
          WRITE(72, *) eig_0
      END IF    
  END DO
  CLOSE(72)

  ! Decomposição espectral de Lanczos
  omega = 2 * PI * 100
  ni = 4e-7 * PI

  DO j1 = 1, niter
      Mat_eig(j1, j1) = eig(j1) + imag * omega * ni
      mod_Mat_eig(j1) = (REAL(Mat_eig(j1, j1)) ** 2 + AIMAG(Mat_eig(j1, j1)) ** 2)
      Conj_Mat_eig(j1, j1) = CONJG(Mat_eig(j1, j1))
      Inv_Mat_eig(j1, j1) = Conj_Mat_eig(j1, j1) / mod_Mat_eig(j1)
  END DO

  DO j1 = 1, niter
      DO j2 = 1, niter
          eig_vec_transp(j1, j2) = eig_vec(j2, j1)
      END DO
  END DO

  vec_e = 0.d0
  vec_e(1) = 1.d0

  CALL MATMAT(nlin, niter, V, &
              niter, niter, eig_vec, &
              Mat_1)

  WRITE(*, *) 'ok1'

  CALL MATMAT(nlin, niter, Mat_1, &
              niter, niter, Inv_Mat_eig, &
              Mat_2)

  WRITE(*, *) 'ok2'

  CALL MATMAT(nlin, niter, Mat_2, &
              niter, niter, eig_vec_transp, &
              Mat_3)

  WRITE(*, *) 'ok3'

  CALL MATVEC(nlin, niter, Mat_3, vec_e, vec_1)
  WRITE(*, *) 'ok4'

  DO j1 = 1, nlin
      field(j1) = mod_VecB * vec_1(j1) * Diag_inv(j1, j1)
  END DO

  mod_field = SQRT(DOT_PRODUCT(field, field))

  OPEN(47, FILE = 'Campo_100.dat')
  DO j1 = 1, nlin
      WRITE(47, '(E20.10E3)') field(j1)
  END DO
  CLOSE(47)

  WRITE(*, *) 'Fim dos cálculos'

  CALL CPU_TIME(time_end)
  WRITE(*, *) 'Tempo de execução foi', time_end - time_begin, 'segundos'

CONTAINS

  SUBROUTINE eigenvalues_complex(d, u, eig_vec, n)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(0.d0)
    CHARACTER :: COMPZ = 'I'
    INTEGER :: n, lgn
    REAL(dp), DIMENSION(n) :: d
    REAL(dp), DIMENSION(n - 1) :: u
    COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: eig_vec
    INTEGER :: LDZ
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: LWORK, LIWORK
    REAL(dp), DIMENSION(:), ALLOCATABLE :: RWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: LRWORK, INFO

    lgn = CEILING(LOG(REAL(n, KIND=dp)) / LOG(2.0_dp))
    LDZ = N

    ALLOCATE(WORK(1), IWORK(1), RWORK(1))
    LWORK = -1
    LRWORK = -1
    LIWORK = -1

    ! Chamada fictícia para obter tamanhos de trabalho necessários
    ! CALL zstedc(COMPZ, N, d, u, eig_vec, n, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

    LWORK = MAX(n * n, NINT(REAL(WORK(1))))
    LRWORK = MAX(1 + 3 * n + 2 * n * lgn + 4 * n * n, NINT(RWORK(1)))
    LIWORK = MAX(6 + 6 * n + 5 * n * lgn, IWORK(1))

    DEALLOCATE(WORK, IWORK, RWORK)		
    ALLOCATE(WORK(LWORK), RWORK(LRWORK), IWORK(LIWORK))

    ! Chamada para calcular autovalores e autovetores
    ! CALL zstedc(COMPZ, N, d, u, eig_vec, n, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

    DEALLOCATE(WORK, IWORK, RWORK)		
  END SUBROUTINE eigenvalues_complex

  SUBROUTINE MATMAT(Anrow, Ancol, Amatrix, Bnrow, Bncol, Bmatrix, Cmatrix)
    IMPLICIT NONE
    INTEGER(4) :: Anrow, Ancol, Bnrow, Bncol, i, j
    COMPLEX(8) :: Amatrix(Anrow, Ancol), Bmatrix(Bnrow, Bncol), Cmatrix(Anrow, Bncol)

    DO j = 1, Bncol
       DO i = 1, Anrow
          Cmatrix(i, j) = 0
       END DO
    END DO
    DO j = 1, Bncol
       CALL GAXPY(Anrow, Ancol, Amatrix, Bmatrix(:, j), Cmatrix(:, j))
    END DO
  END SUBROUTINE MATMAT

  SUBROUTINE GAXPY(Nrow, Ncol, Amatrix, xvector, yvector)
    IMPLICIT NONE
    INTEGER(4) :: Nrow, Ncol, i, j
    COMPLEX(8) :: Amatrix(Nrow, Ncol), xvector(Ncol), yvector(Nrow), x_i

    DO i = 1, Ncol
       x_i = xvector(i)
       CALL SAXPY(Nrow, x_i, Amatrix(:, i), yvector)
    END DO
  END SUBROUTINE GAXPY
  
  SUBROUTINE SAXPY(Nrow, alpha, xvector, yvector)
    IMPLICIT NONE
    INTEGER(4) :: Nrow, i
    COMPLEX(8) :: alpha, xvector(Nrow), yvector(Nrow)

    DO i = 1, Nrow
       yvector(i) = yvector(i) + alpha * xvector(i)
    END DO
  END SUBROUTINE SAXPY

  SUBROUTINE MATVEC(Nrow, Ncol, Amatrix, xvector, yvector)
    IMPLICIT NONE
    INTEGER(4) :: Nrow, Ncol, i, j
    COMPLEX(8) :: Amatrix(Nrow, Ncol), yvector(Nrow)
    REAL(8) :: xvector(Ncol)

    DO i = 1, Nrow
       yvector(i) = 0.d0
    END DO

    CALL GAXPY2(Nrow, Ncol, Amatrix, xvector, yvector)
  END SUBROUTINE MATVEC

  SUBROUTINE GAXPY2(Nrow, Ncol, Amatrix, xvector, yvector)
    IMPLICIT NONE
    INTEGER(4) :: Nrow, Ncol, i, j
    COMPLEX(8) :: Amatrix(Nrow, Ncol), yvector(Nrow)
    REAL(8) :: xvector(Ncol), x_i

    DO i = 1, Ncol
       x_i = xvector(i)
       CALL SAXPY2(Nrow, x_i, Amatrix(:, i), yvector)
    END DO
  END SUBROUTINE GAXPY2

  SUBROUTINE SAXPY2(Nrow, alpha, xvector, yvector)
    IMPLICIT NONE
    INTEGER(4) :: Nrow, i
    COMPLEX(8) :: xvector(Nrow), yvector(Nrow)
    REAL(8) :: alpha

    DO i = 1, Nrow
       yvector(i) = yvector(i) + alpha * xvector(i)
    END DO
  END SUBROUTINE SAXPY2
  
END PROGRAM SDLM_last
