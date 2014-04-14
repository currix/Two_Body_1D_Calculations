MODULE module_2body_ISQW
  !
  !
  USE nrtype
  USE vardef_ri
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: f_x_mat
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: f_x2_mat
  !
CONTAINS
  !
  ! E1 CALCULATION
  !
  SUBROUTINE comp_2body_E1(aval_2b, avec_2b, dim_2b, dim_eig2body, f_x_mat)
    !     
    !     Compute matrix elements < (2B) psi_i | x_1 + x_2 | (2B) psi_A > in the ISQW case
    !     with i = 1, dim_eig2body FOR ALL psi_A
    !
    !
    !     INPUT 
    !           apar        ---> HO length
    !           aval_2b     ---> 2-BODY EIGENVALUES
    !           avec_2b     ---> 2-BODY EIGENVECTORS 
    !           dim_2b      ---> DIMENSION OF THE 2-BODY SYSTEM
    !           f_x_mat     ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
    !
    !     by Currix TM
    !
    !
    ! Arguments
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  f_x_mat, avec_2b
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b, dim_eig2body
    !
    !calcola e1_total e scrive risultati
    ! Local variables
    INTEGER(KIND = I4B) :: index_I, index_J
    REAL(KIND = DP) :: e1_result, e1_total
    REAL(KIND = DP), DIMENSION(1:dim_2b) :: avec_i, avec_j
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "isqw_1D_2b"
    !
    WRITE(*,*) 'E1 RESULTS '
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
       WRITE(200+index_I,*) "#E1"
       !
       DO index_J = 1, dim_2b
          !
          avec_j(:) = avec_2b(:,index_J)
          !
          CALL comp_2body_E1_val(avec_i, avec_j, dim_2b, f_x_mat, e1_result)
          !
          e1_result = e1_result*((4.0_dp*X_max)/(PI_D**2)) ! constant factor
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
    IMPLICIT NONE
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
    !DO I = 1, dim_2b
    !   DO J = 1, dim_2b
          !     
    !      e1_result = e1_result  + EIGVEC_A(J)*EIGVEC_B(I)*F_X_MAT(I,J) 
          !
    !   ENDDO
    !ENDDO
    !
  END SUBROUTINE comp_2body_E1_val
  !
  SUBROUTINE comp_isqw_mat_e1(dim_bound, dim_1b, dim_2b, avec_1b_Box, two_body_bas, f_x_mat)
    !     
    !     Compute matrix elements < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 >
    !
    !     INPUT 
    !           dim_bound     ---> NUMBER OF 1-BODY BOUND STATES IN THE CORE
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
    IMPLICIT NONE
    !    
    ! ARGUMENTS 
    INTEGER(KIND = I4B), INTENT(IN) :: dim_bound, dim_1b, dim_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b_Box
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    !     
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: f_x_mat
    !
    INTEGER(KIND = I4B) :: I, J
    DO I = 1, DIM_2B
       DO J = I, DIM_2B
          !
          F_X_MAT(J,I) =  FXSQW(I, J, DIM_bound, DIM_1B, AVEC_1b_Box, two_Body_BAS)
          F_X_MAT(I,J) =  F_X_MAT(J,I)
          !     
       ENDDO
    ENDDO
    !
  CONTAINS
    !
    FUNCTION fxsqw(I, J, dim_bound, dim_1b, avec_1b_Box, two_body_bas)
      INTEGER(KIND = I4B), INTENT(IN) :: I, J, DIM_bound, DIM_1B
      INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
      !
      REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: AVEc_1b_Box
      !
      REAL(KIND = DP) :: fxsqw
      !
      INTEGER(KIND = I4B) :: N1I, N2I, N1J, N2J
      INTEGER(KIND = I4B) :: KI, KJ
      REAL(KIND = DP) :: FXSQW1, FXSQW2, FXSQW3, FXSQW4
      !
      REAL(KIND = DP) :: DDELTAK, X_BOX_MATEL      
      EXTERNAL DDELTAK, X_BOX_MATEL
      !
      N1I = two_body_BAS(I,1)
      N1J = two_body_BAS(J,1)
      N2I = two_body_BAS(I,2)
      N2J = two_body_BAS(J,2)
      !     
      FXSQW = 0.5_dp*SQRT((2.0_dp-DDELTAK(N1I,N2I))*&
           (2.0_dp-DDELTAK(N1J,N2J)))
      !    
      FXSQW1 = 0.0_dp
      IF (INT(DDELTAK(N2I,N2J)) == 1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FXSQW1 = FXSQW1 + &
                    AVEC_1b_Box(KI,DIM_bound+N1I)*AVEC_1b_Box(KJ,DIM_bound+N1J)*X_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !
      FXSQW2 = 0.0_dp
      IF (INT(DDELTAK(N1I,N1J))==1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FXSQW2 = FXSQW2 + &
                    AVEC_1b_Box(KI,DIM_bound+N2I)*AVEC_1b_Box(KJ,DIM_bound+N2J)*X_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !
      FXSQW3 = 0.0_dp
      IF (INT(DDELTAK(N1I,N2J))==1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FXSQW3 = FXSQW3 + & 
                    AVEC_1b_Box(KI,DIM_bound+N2I)*AVEC_1b_Box(KJ,DIM_bound+N1J)*X_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !     
      FXSQW4 = 0.0_dp
      IF (INT(DDELTAK(N2I,N1J))==1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FXSQW4 = FXSQW4 + &
                    AVEC_1b_Box(KI,DIM_bound+N1I)*AVEC_1b_Box(KJ,DIM_bound+N2J)*X_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !
      FXSQW = FXSQW*(FXSQW1+FXSQW2+FXSQW3+FXSQW4)
    END FUNCTION fxsqw
    !
  END SUBROUTINE comp_isqw_mat_e1
  !
  ! E2 CALCULATION
  !
  SUBROUTINE comp_2body_E2(aval_2b, avec_2b, dim_2b, dim_eig2body, f_x2_mat)
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
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  f_x2_mat, avec_2b
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b, dim_eig2body
    !
    !calcola e1_total e scrive risultati
    ! Local variables
    INTEGER(KIND = I4B) :: index_I, index_J
    REAL(KIND = DP) :: e2_result, e2_total
    REAL(KIND = DP), DIMENSION(1:dim_2b) :: avec_i, avec_j
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "isqw_1D_2b"
    !
    file = 'E2'
    !
    WRITE(*,*) 'E2 RESULTS '
    !
    DO index_I = 1, dim_eig2body
       !
       e2_total = 0.0_DP
       !
       avec_i(:) = avec_2b(:,index_I)
       !
       WRITE(filename, '(A, "_",A,"_",I1, ".dat")') TRIM(prog), TRIM(file), index_I
       OPEN(UNIT = 300+ index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(300+index_I,*) "#E2"
       !
       DO index_J = 1, dim_2b
          !
          avec_j(:) = avec_2b(:,index_J)
          !
          !print*, avec_i, avec_j
          CALL comp_2body_E2_val(avec_i, avec_j, dim_2b, f_x2_mat, e2_result)
          !
          e2_result = e2_result*(X_max**2) ! constant factor
          !
          !print*, e2_result 
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

  END SUBROUTINE comp_2body_E2
  !
  SUBROUTINE comp_2body_E2_val(eigvec_a, eigvec_b, dim_2b, f_x2_mat, e2_result)
    !     
    !     Compute matrix elements < (2B) psi_B | x^2_1 + x^2_2 | (2B) psi_A > 
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
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b
    REAL(KIND = DP), DIMENSION(1:dim_2b), INTENT(IN) ::  eigvec_a, eigvec_b
    REAL(KIND = DP), DIMENSION(1:dim_2b,1:dim_2b), INTENT(IN) ::  f_x2_mat
    !
    REAL(KIND = DP), INTENT(OUT) ::  e2_result
    !
    E2_RESULT = 0.0_DP
    !
    !print*, eigvec_a, eigvec_b, f_x2_mat
    e2_result = DOT_PRODUCT(eigvec_a, MATMUL(f_x2_mat,eigvec_b))
    !
  END SUBROUTINE comp_2body_E2_val
  !
  SUBROUTINE comp_isqw_mat_e2(dim_BOUND, dim_1b, dim_2b, avec_1b_box, two_body_bas, f_x2_mat)
    !     
    !     Compute matrix elements < (s) n'_1 n'_2 | x^2_1 + x^2_2 | (s) n_1 n_2 >
    !
    !     INPUT 
    !           dim_bound     ---> NUMBER OF 1-BODY BOUND STATES IN THE CORE
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
    INTEGER(KIND = I4B), INTENT(IN) :: dim_bound, dim_1b, dim_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_1b_box
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    !    
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: f_x2_mat 
    !
    INTEGER(KIND = I4B) :: I, J
    !     
    !     
    DO I = 1, DIM_2B
       DO J = I, DIM_2B
          !     
          F_X2_MAT(J,I) = FX2SQW(I, J) 
          F_X2_MAT(I,J) =  F_X2_MAT(J,I)
          !     simmetria
       ENDDO
       !print*, f_x2_mat
    ENDDO
    !
    RETURN
    !
  CONTAINS
    !
    FUNCTION FX2SQW(I, J)
      !
      IMPLICIT NONE
      !
      INTEGER(KIND = I4B), INTENT(IN) :: I, J
      REAL(KIND = DP) :: FX2SQW
      !
      !REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: AVECSQW
      !
      INTEGER(KIND = I4B) :: N1I, N2I, N1J, N2J
      INTEGER(KIND = I4B) :: KI, KJ
      REAL(KIND = DP) :: FX2SQW1, FX2SQW2, FX2SQW3, FX2SQW4
      !
      REAL(KIND = DP) :: DDELTAK, X2_BOX_MATEL      
      EXTERNAL DDELTAK, X2_BOX_MATEL
      !
      N1I = two_body_bas(I,1)
      N1J = two_body_bas(J,1)
      N2I = two_body_bas(I,2)
      N2J = two_body_bas(J,2)
      !     
      !     
      FX2SQW = 0.5_DP*SQRT((2.0_DP-DDELTAK(N1I,N2I))*&
           (2.0D0-DDELTAK(N1J,N2J)))
      !     
      FX2SQW1 = 0.0_DP
      IF ( INT(DDELTAK(N2I,N2J)) == 1 ) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FX2SQW1 = FX2SQW1 + &
                    AVEC_1b_box(KI,DIM_bound+N1I)*AVEC_1b_box(KJ,DIM_bound+N1J)*&
                    X2_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !     
      FX2SQW2 = 0.0_DP
      IF (INT(DDELTAK(N1I,N1J)) == 1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FX2SQW2 = FX2SQW2 + &
                    AVEC_1b_box(KI,DIM_bound+N2I)*AVEC_1b_box(KJ,DIM_bound+N2J)*&
                    X2_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !     
      FX2SQW3 = 0.0_DP
      IF (INT(DDELTAK(N1I,N2J)) == 1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FX2SQW3 = FX2SQW3 + &
                    AVEC_1b_box(KI,DIM_bound+N2I)*AVEC_1b_box(KJ,DIM_bound+N1J)* &
                    X2_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !     
      FX2SQW4 = 0.0_DP
      IF (INT(DDELTAK(N2I,N1J))== 1) THEN
         DO KI = 1, DIM_1B
            DO KJ = 1, DIM_1B
               FX2SQW4 = FX2SQW4 + &
                    AVEC_1b_box(KI,DIM_bound+N1I)*AVEC_1b_box(KJ,DIM_bound+N2J)* &
                    X2_BOX_MATEL(KI, KJ)
            ENDDO
         ENDDO
      ENDIF
      !
      !
      FX2SQW = FX2SQW*(FX2SQW1+FX2SQW2+FX2SQW3+FX2SQW4)
    END FUNCTION FX2SQW
    !
  END SUBROUTINE comp_isqw_mat_e2
  !
END MODULE module_2body_ISQW
