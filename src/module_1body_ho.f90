MODULE module_1body_ho
  !
  USE nrtype
  USE constants
  USE vardef_ri
  USE pot_param
  !
  IMPLICIT NONE
  !
  !
  ! One body variables
  !
  ! DIMENSION OF THE 1body HO BASIS
  INTEGER(KIND = I4B) :: dim_1b_HO
  !
  ! HARMONIC BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_1b_Har
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Har_1b_Bas, Avec_1b_Har, Avec_1b_Har_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_1b_Har_Der_X
  !
  !
CONTAINS
  !
  !
  SUBROUTINE HO_1D_BASIS(apar, NdimH, Iprint)
    !     
    !     COMPUTES A 1D HARMONIC BASIS
    !
    !     INPUT  :: apar    --> LENGTH SCALE OF THE PROBLEM
    !               NdimH   --> DIMENSION + 1 OF THE HARMONIC BASIS
    !
    !
    !     OUTPUT (Module) :: HAR_1b_BAS  --> MATRIX WITH HARMONIC BASIS
    !
    !     FORMAT :: Iprint  --> VERBOSITY CONTROL
    !
    !
    !     by Currix TM.
    !
    !
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    INTEGER(KIND = I4B), INTENT(IN) :: NdimH, Iprint
    REAL(KIND = DP), INTENT(IN) :: apar 
    !
    !
    ! OTHER VARIABLES
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: HO_norm_test
    REAL(KIND = DP) :: PI14, apar2, ainteg, Error
    INTEGER(KIND = I4B) :: kx, Ifail, Ierr
    !
    IF (Iprint > 1) PRINT*, "BUILDING HARMONIC BASIS"
    !
    PI14 = SQRT(SQRT(PI_D))
    apar2 = apar*apar
    !
    !
    !     HO n = 0
    Har_1b_Bas(:,1) = SQRT(apar)/(PI14)*EXP(-apar2*X_Grid*X_Grid/2.0_DP)
    !
    !     HO n = 1
    IF (NdimH > 1) Har_1b_Bas(:,2) = (SQRT(apar)/(PI14))*SQRT(2.0_DP)*apar*X_Grid*EXP(-apar2*X_Grid*X_Grid/2.0_DP)
    !
    !    RECURRENCE RELATION (WATCH OUT THE INDEXES :: MATRIX START AT 1, NOT 0)
    DO kx = 2, NdimH-1
       Har_1b_Bas(:,kx+1) = &
            SQRT(2.0_DP/(1.0_DP*kx))*apar*X_Grid*Har_1b_Bas(:,kx) - &
            SQRT((1.0_DP*(kx-1))/(1.0_DP*kx))*Har_1b_Bas(:,kx-1)
    ENDDO
    !
    !     TESTING NORMALIZATION 
    IF (Iprint > 1) THEN
       !
       ! DEFINE HO_norm_test
       ALLOCATE(HO_norm_test(1:dim_X), STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "HO_norm_test allocation request denied."
          STOP
       ENDIF
       !
       DO kx = 1, NdimH
          HO_norm_test = Har_1b_Bas(:,kx)*Har_1b_Bas(:,kx)
          Ifail = 0
          !
          CALL D01GAF(X_Grid, HO_norm_test, dim_X, ainteg, Error, Ifail)
          PRINT*, "HARMONIC FUNCTION ", kx, " NORMALIZATION", ainteg
       ENDDO
       DEALLOCATE(HO_norm_test, STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "HO_norm_test deallocation request denied."
          STOP
       ENDIF
       !
       WRITE(*,*) "DONE"
       !
    ENDIF
    !     
    !     
  END SUBROUTINE HO_1D_BASIS
  !
  !
  SUBROUTINE HARDIAG(apt, Iprint)
    !     
    !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A HARMONIC BASIS
    !
    !     INPUT  :: 
    !               apt     --> PROBLEM LENGTH SCALE (fm^-1)
    !
    !     FORMAT :: IPRINT  --> VERBOSITY CONTROL
    !
    !     by Currix TM.
    !
    !
    !USE nrtype
    !USE constants
    !USE pot_param
    !USE vardef_ri
    !
    ! Lapack 95
    USE LA_PRECISION, ONLY: WP => DP
    USE F95_LAPACK, ONLY: LA_SYEVR
    !
    IMPLICIT NONE
    !
    !
    ! ARGUMENTS
    INTEGER(KIND = I4B), INTENT(IN) :: Iprint
    REAL(KIND = DP), INTENT(IN) ::  apt
    !
    ! POTENTIAL FUNCTION
    INTERFACE Potf
       !
       ELEMENTAL FUNCTION Potf(X)
         !
         !     WOODS-SAXON 1D POTENTIAL
         !
         USE nrtype
         USE constants
         USE pot_param
         !
         IMPLICIT NONE
         !
         ! ARGUMENTS
         REAL(KIND = DP), INTENT(IN) :: X
         REAL(KIND = DP) :: Potf
       END FUNCTION Potf
       !
    END INTERFACE Potf
    !
    !     HAMILTONIAN MATRIX
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  Ham_Mat
    !     POTENTIAL VECTOR
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE ::  Pot_vec
    !
    !
    INTEGER(KIND = I4B) ::  I, J, KX, IFAIL, Ierr
    REAL(KIND = DP) ::  apt2, AJ, AHAM, ERROR
    !
    apt2 = apt*apt*h_sq_over_m/2.0_DP
    !
    IF (Iprint > 1) PRINT*, "BUILDING HAMILTONIAN MATRIX"
    !     
    !   
    ALLOCATE(Ham_Mat(1:dim_1b_ho,1:dim_1b_ho), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Ham_Mat allocation request denied."
       STOP
    ENDIF
    ! 
    ALLOCATE(Pot_vec(1:dim_X), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Pot_vec allocation request denied."
       STOP
    ENDIF
    ! 
    Ham_Mat = 0.0_DP
    !
    !
    !     HAMILTONIAN MATRIX
    !     
    DO I = 1, dim_1b_ho
       !
       DO J = I, dim_1b_ho
          AJ = REAL(J - 1, DP)
          !     KINETIC ENERGY (ADIMENSIONAL UNITS)
          IF (I.EQ.J-2) Ham_Mat(J,I) = Ham_Mat(J,I) -  apt2*SQRT(AJ*(AJ-1.0_DP))/2.0_DP
          IF (I.EQ.J) Ham_Mat(J,I)   = Ham_Mat(J,I) +  apt2*(2.0_DP*AJ+1.0_DP)/2.0_DP
          !
          !     POTENTIAL ENERGY
          Pot_vec = Potf(X_Grid)*har_1b_bas(:,I)*har_1b_bas(:,J)
          !
          !     INTEGRATION
          IFAIL = 0
          CALL D01GAF(X_Grid,Pot_vec,dim_X,AHAM,ERROR,IFAIL)
          Ham_Mat(J,I) = Ham_Mat(J,I) + AHAM
          !
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE(Pot_vec, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Pot_vec deallocation request denied."
       STOP
    ENDIF
    !     DIAGONALIZATION USING LAPACK
    CALL LA_SYEVR(A=Ham_Mat, W=aval_1b_har, JOBZ='V', UPLO='L')
    !
    !     EIGENVECTOR MATRIX
    !
    AVEC_1B_HAR = Ham_Mat
    !
    AVEC_1B_HAR_X = 0.0_DP
    !
    DO KX = 1, dim_X
       DO I = 1, dim_1b_ho
          DO J = 1, dim_1b_ho
             AVEC_1b_HAR_X(KX,I) = AVEC_1B_HAR_X(KX,I) + Ham_Mat(J,I)*HAR_1b_BAS(KX,J)
          ENDDO
       ENDDO
    ENDDO
    !
    DEALLOCATE(Ham_Mat, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Ham_Mat deallocation request denied."
       STOP
    ENDIF
    !    
  END SUBROUTINE HARDIAG
  !
  SUBROUTINE output_1b_results(Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS)
    !
    !USE nrtype
    !USE vardef_ri
    !
    IMPLICIT NONE
    !
    !
    INTERFACE Potf
       !
       ELEMENTAL FUNCTION Potf(x)
         !
         USE nrtype
         USE constants
         USE pot_param
         !
         IMPLICIT NONE
         !
         ! ARGUMENTS
         REAL(KIND = DP), INTENT(IN) :: x
         REAL(KIND = DP) :: Potf
       END FUNCTION Potf
       !
    END INTERFACE Potf
    !
    ! Argument definition
    INTEGER, INTENT(IN) :: Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS
    !
    !
    ! Local variables
    INTEGER :: I_index, X_index
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "ho_1D_1b"
    !
    IF (Iprint > 0) THEN
       PRINT*, "EIGENVALUES IN A ", dim_1b_HO, " DIM HARMONIC BASIS"
       DO I_index = 1, dim_1b_HO
          PRINT*, I_index, Aval_1b_Har(I_index)
       ENDDO
    ENDIF
    !
    ! SAVING HARMONIC BASIS
    IF (i_save_1b_bas == 1) THEN
       file = 'basis'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       IF ( dim_1b_HO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       IF ( dim_1b_HO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(70,*) "# HO  dim_1b_HO = ", dim_1b_HO, " Box radius = ", X_max, " fm"
       WRITE(70,*) "#Grid     Harmonic basis basis"
       !
       DO X_index = 1, dim_X
          WRITE(70,11) X_grid(X_index), har_1b_bas(X_index,1:dim_1b_HO)
       ENDDO
       CLOSE(UNIT = 70)
    ENDIF
    !
    ! SAVING EIGENVECTORS
    file = 'eigenvectors'
    WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
    IF ( dim_1b_HO < 10) THEN !to avoid spaces
       WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
    ENDIF
    IF ( dim_1b_HO > 99) THEN 
       WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
    ENDIF
    OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    WRITE(71,*) "# HO  dim_1b_HO = ", dim_1b_HO, " Box radius = ", X_max, " fm"
    WRITE(71,*) "#Grid    Eigenvectors"
    !
    IF (i_save_1b_wf == 1) THEN
       DO X_index = 1, dim_X
          WRITE(71,11) X_grid(X_index), avec_1b_har_x(X_index,1:dim_1b_HO)
       ENDDO
       CLOSE(UNIT = 71)
    ENDIF
    !
    ! SAVING ENERGIES
    IF (i_save_1b_en == 1) THEN
       !
       file = 'eigenvalues'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       IF ( dim_1b_HO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       IF ( dim_1b_HO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(73,*) "# HO  dim_1b_HO = ", dim_1b_HO, " Box radius = ", X_max, " fm" 
       WRITE(73,*) "# Eigenvalues"
       DO I_index = 1, dim_1b_HO 
          WRITE(73,10) I_index, Aval_1b_har(I_index)
       ENDDO
       !
       file = 'eigenvec_der'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       IF ( dim_1b_HO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       IF ( dim_1b_HO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       OPEN(UNIT = 72, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(72,*) "# HO  dim_1b_HO = ", dim_1b_HO, " Box radius = ", X_max, " fm"
       WRITE(72,*) "# Eigenvectors derivatives"
       do X_index = 1, dim_X
          WRITE(72,11) X_grid(X_index), Potf(X_grid(X_index)), &
               (10.0_DP*avec_1b_har_x(X_index,I_index)+aval_1b_har(I_index), I_index = 1, dim_1b_ho)
       ENDDO
       !
       file = 'pot_eigvec'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       IF ( dim_1b_HO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       IF ( dim_1b_HO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_HO
       ENDIF
       OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(74,*) "# HO  dim_1b_HO = ", dim_1b_HO, " Box radius = ", X_max, " fm"
       WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
       !
       DO X_index = 1, dim_X
          WRITE(74,11) X_grid(X_index), Potf(X_grid(X_index)), &
               (10.0_DP*avec_1b_har_x(X_index,I_index)**2+aval_1b_har(I_index), I_index = 1, dim_1b_ho)
       ENDDO
       !
       CLOSE(UNIT = 72)
       CLOSE(UNIT = 73)
       CLOSE(UNIT = 74)
    ENDIF
    !
10  FORMAT (1X,I6,1X,E16.8)
11  FORMAT (1X,E14.6,1X,150E16.8)
    !
  END SUBROUTINE output_1b_results
  !
  !
  !
END MODULE module_1body_ho
  
