MODULE module_2body_THO
  !
  !
  USE nrtype
  USE vardef_ri
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: f_x_mat
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: f_x2_mat
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: xthomatel
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: x2thomatel
  !
CONTAINS
  !
  SUBROUTINE comp_2body_E1(apar, aval_2b, avec_2b, dim_2b, dim_eig2body, f_x_mat)
    !     
    !     Compute matrix elements < (2B) psi_i | x_1 + x_2 | (2B) psi_A > in the HO case
    !     with i = 1, dim_eig2body FOR ALL psi_A
    !
    !
    !     INPUT 
    !           apar        ---> THO length
    !           aval_2b     ---> 2-BODY EIGENVALUES
    !           avec_2b     ---> 2-BODY EIGENVECTORS 
    !           dim_2b      ---> DIMENSION OF THE 2-BODY SYSTEM
    !           f_x_mat     ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
    !
    !     by Currix TM
    !
    !
    ! Arguments
    REAL(KIND = DP), INTENT(IN) ::  apar
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  f_x_mat, avec_2b
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b, dim_eig2body
    !
    !
    ! Local variables
    INTEGER(KIND = I4B) :: index_I, index_J
    REAL(KIND = DP) :: e1_result, e1_total
    REAL(KIND = DP), DIMENSION(1:dim_2b) :: avec_i, avec_j
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "tho_1D_2b"
    !
    file = 'E1'
    DO index_I = 1, dim_eig2body
       !
       e1_total = 0.0_DP
       !
       avec_i(:) = avec_2b(:,index_I)
       !
       WRITE(filename, '(A, "_",A,"_",I1, ".dat")') TRIM(prog), TRIM(file), index_I
       OPEN(UNIT = 200+ index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(200+index_I,*) "E1"
       !
       DO index_J = 1, dim_2b
          !
          avec_j(:) = avec_2b(:,index_J)
          !
          CALL comp_2body_E1_val(avec_i, avec_j, dim_2b, f_x_mat, e1_result)
          !
          e1_result = e1_result/apar
          !
          e1_total = e1_total + e1_result**2
          !
          WRITE(UNIT = 200 + index_I, FMT = *) index_J, aval_2b(index_J), e1_result, e1_result**2
          !
       ENDDO
       CLOSE(UNIT = 200 + index_I)
       !
       WRITE(*,*) 'Total ', index_I,' E1 transition strength = ', e1_total
       !
    ENDDO
    !
  END SUBROUTINE comp_2body_E1
  !
  SUBROUTINE comp_2body_E1_val(eigvec_a, eigvec_b, dim_2b, f_x_mat, e1_result)
    !     
    !     Compute matrix elements < (2B) psi_B | x_1 + x_2 | (2B) psi_A > in the HO case
    !
    !     INPUT 
    !           eigvec_a   ---> 2-BODY EIGENVECTOR A 
    !           eigvec_b   ---> 2-BODY EIGENVECTOR B
    !           dim_2b     ---> DIMENSION OF THE 2-BODY SYSTEM
    !           f_x_mat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
    !
    !     OUTPUT
    !           e1_result  ---> < (2B) psi_B | x_1 + x_2 | (2B) psi_A >
    !
    !     by Currix TM
    !
    !
    ! Arguments
    REAL(KIND = DP), DIMENSION(1:dim_2b), INTENT(IN) ::  eigvec_a, eigvec_b
    REAL(KIND = DP), DIMENSION(1:dim_2b,1:dim_2b), INTENT(IN) ::  f_x_mat
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b
    !
    REAL(KIND = DP), INTENT(OUT) ::  e1_result
    !
    e1_result = 0.0_DP
    !
    e1_result = DOT_PRODUCT(eigvec_a, MATMUL(f_x_mat,eigvec_b))
    !
  END SUBROUTINE comp_2body_E1_val
  !
  SUBROUTINE comp_x_tho_matel(dim_1b_THO, THar_1b_Bas, der_S_x)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_1b_THO
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: der_S_x 
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: THar_1b_Bas
    !
    REAL(KIND = DP) :: integ, error
    REAL(KIND = DP), DIMENSION(1:dim_X) :: matel
    INTEGER(KIND = I4B) :: Ifail, j, k
    !
    !
    xthomatel = 0.0_DP
    matel = 0.0_DP
    DO j = 1, dim_1b_THO
       DO k = j, dim_1b_THO
          !
          integ = 0.0_DP
          Ifail = 0
          !
          IF (j /= k) THEN
             matel = X_Grid(:)*THar_1b_Bas(:,k)*THar_1b_Bas(:,j)*der_S_x(:)   !E' qui il problema per k = 3 !!!!!!!!!
             !
             CALL D01GAF(X_Grid, matel, dim_X, integ, error, Ifail)
             IF(Ifail /= 0) print*, "Integration error"
          ENDIF
          xthomatel(k,j) = integ
          xthomatel(j,k) = xthomatel(k,j)
       ENDDO
    ENDDO
    !
  END SUBROUTINE comp_x_tho_matel
  !
  SUBROUTINE comp_tho_mat_e1(dim_core, dim_1b, dim_2b, avec_1b, two_body_bas, f_x_mat, xthomatel)
    !     
    !     Compute matrix elements < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 >
    !
    !     INPUT 
    !           dim_core     ---> NUMBER OF 1-BODY BOUND STATES IN THE CORE
    !           dim_1b       ---> 1-BODY HARMONIC BASIS DIMENSION. 
    !           avec_1b      ---> 1-BODY BASIS : HARMONIC COMPONENTS
    !           dim_2b       ---> DIMENSION OF THE 2-BODY SYSTEM
    !           two_body_bas ---> LABELS OF THE 2-BODY |(s) n_1 n_2> BASIS
    !
    !     OUTPUT
    !           f_x_mat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
    !
    !     by Currix TM
    !
    !    
    ! ARGUMENTS 
    INTEGER(KIND = I4B), INTENT(IN) :: dim_core, dim_1b, dim_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: xthomatel
    !
    !     
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: f_x_mat
    !
    INTEGER(KIND = I4B) :: index_I, index_J
    !    
    !     
    DO index_I = 1, dim_2b
       DO index_J = index_I, dim_2b
          !  
          f_x_mat(index_J, index_I) = fxho(index_I, index_J, dim_core, dim_1b, avec_1b, two_body_bas, xthomatel) 
          f_x_mat(index_I, index_J) = f_x_mat(index_J, index_I)
          !
       ENDDO
    ENDDO
    !
  CONTAINS
    !
    FUNCTION fxho(index_I, index_J, dim_core, dim_1b, avec_1b, two_body_bas, xthomatel)
      !
      INTEGER(KIND = I4B), INTENT(IN) ::  index_I, index_J, dim_core, dim_1b
      INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
      REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b
      REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: xthomatel
      !
      REAL(KIND = DP) :: fxho
      !
      ! Local variables
      INTEGER(KIND = I4B) ::  n1_i, n2_i, n1_j, n2_j
      INTEGER(KIND = I4B) ::  ki, kj
      REAL(KIND = DP) :: fxho1, fxho2, fxho3, fxho4
      !
      REAL(KIND = DP) :: ddeltak    
      EXTERNAL ddeltak
      !
      n1_i = two_body_bas(index_I,1)
      n1_j = two_body_bas(index_J,1)
      n2_i = two_body_bas(index_I,2)
      n2_j = two_body_bas(index_J,2)
      !     
      fxho = 0.5_DP*SQRT((2.0_DP-ddeltak(n1_i,n2_i))*(2.0_DP-ddeltak(N1_J,N2_J)))
      !
      fxho1 = 0.0_DP
      !
      IF (INT(ddeltak(n2_i,n2_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fxho1 = fxho1 + avec_1b(ki,dim_core+n1_i)*avec_1b(kj,dim_core+n1_j)*xthomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fxho2 = 0.0_DP
      IF (INT(ddeltak(n1_i,n1_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fxho2 = fxho2 + avec_1b(ki,dim_core+n2_i)*avec_1b(kj,dim_core+n2_j)*xthomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fxho3 = 0.0_DP
      IF (INT(ddeltak(n1_i,n2_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fxho3 = fxho3 + avec_1b(ki,dim_core+n2_i)*avec_1b(kj,dim_core+n1_j)*xthomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fxho4 = 0.0_DP
      IF (INT(ddeltak(n2_i,n1_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fxho4 = fxho4 + avec_1b(ki,dim_core+n1_i)*avec_1b(kj,dim_core+n2_j)*xthomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fxho = fxho*(fxho1+fxho2+fxho3+fxho4)
      !
    END FUNCTION fxho
    !
  END SUBROUTINE comp_tho_mat_e1
  !
  !
  SUBROUTINE comp_2body_E2(apar, aval_2b, avec_2b, dim_2b, dim_eig2body, f_x2_mat)
    !     
    !     Compute matrix elements < (2B) psi_i | x^2_1 + x^2_2 | (2B) psi_A > in the HO case
    !     with i = 1, dim_eig2body FOR ALL psi_A
    !
    !
    !     INPUT 
    !           apar        ---> HO length
    !           aval_2b     ---> 2-BODY EIGENVALUES
    !           avec_2b     ---> 2-BODY EIGENVECTORS 
    !           dim_2b      ---> DIMENSION OF THE 2-BODY SYSTEM
    !           f_x2_mat     ---> < (s) n'_1 n'_2 | x^2_1 + x^2_2 | (s) n_1 n_2 > matrix
    !
    !     by Currix TM
    !
    !
    ! Arguments
    REAL(KIND = DP), INTENT(IN) ::  apar
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  f_x2_mat, avec_2b
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b, dim_eig2body
    !
    !
    ! Local variables
    INTEGER(KIND = I4B) :: index_I, index_J
    REAL(KIND = DP) :: e2_result, e2_total
    REAL(KIND = DP), DIMENSION(1:dim_2b) :: avec_i, avec_j
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "tho_1D_2b"
    !
    file = 'E2'
    DO index_I = 1, dim_eig2body
       !
       e2_total = 0.0_DP
       !
       avec_i(:) = avec_2b(:,index_I)
       !
       WRITE(filename, '(A, "_",A,"_",I1, ".dat")') TRIM(prog), TRIM(file), index_I
       OPEN(UNIT = 300+ index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(300+index_I,*) "E2"
       !
       DO index_J = 1, dim_2b
          !
          avec_j(:) = avec_2b(:,index_J)
          !
          CALL comp_2body_E2_val(avec_i, avec_j, dim_2b, f_x2_mat, e2_result)
          !
          e2_result = e2_result/(apar**2)
          !
          e2_total = e2_total + e2_result**2
          !
          WRITE(UNIT = 300 + index_I, FMT = *) index_J, aval_2b(index_J), e2_result, e2_result**2
          !
       ENDDO
       CLOSE(UNIT = 300 + index_I)
       !
       WRITE(*,*) 'Total ', index_I,' E2 transition strength = ', e2_total
       !
    ENDDO
    !
  END SUBROUTINE comp_2body_E2
  !
  SUBROUTINE comp_2body_E2_val(eigvec_a, eigvec_b, dim_2b, f_x2_mat, e2_result)
    !     
    !     Compute matrix elements < (2B) psi_B | x^2_1 + x^2_2 | (2B) psi_A > in the HO case
    !
    !     INPUT 
    !           eigvec_a   ---> 2-BODY EIGENVECTOR A 
    !           eigvec_b   ---> 2-BODY EIGENVECTOR B
    !           dim_2b     ---> DIMENSION OF THE 2-BODY SYSTEM
    !           f_x2_mat      ---> < (s) n'_1 n'_2 | x^2_1 + x^2_2 | (s) n_1 n_2 > matrix
    !
    !     OUTPUT
    !           e2_result  ---> < (2B) psi_B | x^2_1 + x^2_2 | (2B) psi_A >
    !
    !     by Currix TM
    !
    !
    ! Arguments
    REAL(KIND = DP), DIMENSION(1:dim_2b), INTENT(IN) ::  eigvec_a, eigvec_b
    REAL(KIND = DP), DIMENSION(1:dim_2b,1:dim_2b), INTENT(IN) ::  f_x2_mat
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b
    !
    REAL(KIND = DP), INTENT(OUT) ::  e2_result
    !
    e2_result = 0.0_DP
    !
    e2_result = DOT_PRODUCT(eigvec_a, MATMUL(f_x2_mat,eigvec_b))
    !
  END SUBROUTINE comp_2body_E2_val
  !
  SUBROUTINE comp_x2_tho_matel(x2thomatel, dim_1b_THO, dim_X, X_Grid, THar_1b_Bas, der_S_x)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_1b_THO, dim_X
    REAL(KIND = DP), DIMENSION(1:dim_X), INTENT(IN) :: der_S_x, X_Grid
    REAL(KIND = DP), DIMENSION(1:dim_X,1:dim_1b_THO), INTENT(IN) :: THar_1b_Bas
    REAL(KIND = DP), DIMENSION(1:dim_1b_THO,1:dim_1b_THO), INTENT(OUT) :: x2thomatel
    !
    REAL(KIND = DP) :: integ, error
    REAL(KIND = DP), DIMENSION(1:dim_X) :: matel
    INTEGER(KIND = I4B) :: Ifail, j, k
    !
    x2thomatel = 0.0_DP
    DO j = 1, dim_1b_THO
       DO k = j, dim_1b_THO
          !
          integ = 0.0_DP
          Ifail = 0
          !
          IF (j /= k) THEN
             matel = (X_Grid(:)**2)*THar_1b_Bas(:,k)*THar_1b_Bas(:,j)*der_S_x(:)
             !
             CALL D01GAF(X_Grid, matel, dim_X, integ, error, Ifail)
          ENDIF
          x2thomatel(k,j) = integ
          x2thomatel(j,k) = x2thomatel(k,j)
       ENDDO
    ENDDO
    !
  END SUBROUTINE comp_x2_tho_matel
  !
  SUBROUTINE comp_tho_mat_e2(dim_core, dim_1b, dim_2b, avec_1b, two_body_bas, f_x2_mat, x2thomatel)
    !     
    !     Compute matrix elements < (s) n'_1 n'_2 | x^2_1 + x^2_2 | (s) n_1 n_2 >
    !
    !     INPUT 
    !           dim_core     ---> NUMBER OF 1-BODY BOUND STATES IN THE CORE
    !           dim_1b       ---> 1-BODY HARMONIC BASIS DIMENSION. 
    !           avec_1b      ---> 1-BODY BASIS : HARMONIC COMPONENTS
    !           dim_2b       ---> DIMENSION OF THE 2-BODY SYSTEM
    !           two_body_bas ---> LABELS OF THE 2-BODY |(s) n_1 n_2> BASIS
    !
    !     OUTPUT
    !           f_x2_mat      ---> < (s) n'_1 n'_2 | x^2_1 + x^2_2 | (s) n_1 n_2 > matrix
    !
    !     by Currix TM
    !
    !    
    ! ARGUMENTS 
    INTEGER(KIND = I4B), INTENT(IN) :: dim_core, dim_1b, dim_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  x2thomatel
    !     
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: f_x2_mat
    !
    INTEGER(KIND = I4B) :: index_I, index_J
    !    
    !     
    DO index_I = 1, dim_2b
       DO index_J = index_I, dim_2b
          !  
          f_x2_mat(index_J, index_I) = fx2ho(index_I, index_J, dim_core, dim_1b, avec_1b, two_body_bas) 
          f_x2_mat(index_I, index_J) = f_x2_mat(index_J, index_I)
          !
       ENDDO
    ENDDO
    !
  CONTAINS
    !
    FUNCTION fx2ho(index_I, index_J, dim_core, dim_1b, avec_1b, two_body_bas)
      !
      INTEGER(KIND = I4B), INTENT(IN) ::  index_I, index_J, dim_core, dim_1b
      INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
      REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b
      !
      REAL(KIND = DP) :: fx2ho
      !
      ! Local variables
      INTEGER(KIND = I4B) ::  n1_i, n2_i, n1_j, n2_j
      INTEGER(KIND = I4B) ::  ki, kj
      REAL(KIND = DP) :: fx2ho1, fx2ho2, fx2ho3, fx2ho4
      !
      REAL(KIND = DP) :: ddeltak  
      EXTERNAL ddeltak
      !
      n1_i = two_body_bas(index_I,1)
      n1_j = two_body_bas(index_J,1)
      n2_i = two_body_bas(index_I,2)
      n2_j = two_body_bas(index_J,2)
      !     
      fx2ho = 0.5_DP*SQRT((2.0_DP-ddeltak(n1_i,n2_i))*(2.0_DP-ddeltak(N1_J,N2_J)))
      !
      fx2ho1 = 0.0_DP
      !
      IF (INT(ddeltak(n2_i,n2_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fx2ho1 = fx2ho1 + avec_1b(ki,dim_core+n1_i)*avec_1b(kj,dim_core+n1_j)*x2thomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fx2ho2 = 0.0_DP
      IF (INT(ddeltak(n1_i,n1_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fx2ho2 = fx2ho2 + avec_1b(ki,dim_core+n2_i)*avec_1b(kj,dim_core+n2_j)*x2thomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fx2ho3 = 0.0_DP
      IF (INT(ddeltak(n1_i,n2_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fx2ho3 = fx2ho3 + avec_1b(ki,dim_core+n2_i)*avec_1b(kj,dim_core+n1_j)*x2thomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fx2ho4 = 0.0_DP
      IF (INT(ddeltak(n2_i,n1_j)) == 1) THEN
         DO ki = 1, dim_1b
            DO kj = 1, dim_1b
               fx2ho4 = fx2ho4 + avec_1b(ki,dim_core+n1_i)*avec_1b(kj,dim_core+n2_j)*x2thomatel(ki-1, kj-1)
            ENDDO
         ENDDO
      ENDIF
      !
      fx2ho = fx2ho*(fx2ho1+fx2ho2+fx2ho3+fx2ho4)
      !
    END FUNCTION fx2ho
    !
  END SUBROUTINE comp_tho_mat_e2
  !
END MODULE module_2body_THO
