MODULE module_1body_tho
  !
  USE nrtype
  USE constants
  USE vardef_ri
  USE pot_param
  !
  IMPLICIT NONE
  !
  !THO CONSTANTS
  REAL(KIND = DP) ::  m = 4.0_dp, gamma, b, K_eff, ratio
  !
  !SCALING FUNCTION
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: S_x, der_S_x, der2_S_x, S_over_x
  !
  !
  ! One body variables
  !
  ! DIMENSION OF THE 1body THO BASIS
  INTEGER(KIND = I4B) :: dim_1b_THO
  !
  ! HARMONIC BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_1b_THar
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: THar_1b_Bas, Avec_1b_THar, Avec_1b_THar_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_1b_THar_Der_X
  !
  !
  !
CONTAINS
  !
  !
  ELEMENTAL FUNCTION lst(x)
    !
    !     ANALYTICAL 1D THO LST
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  lst
    !
    REAL(KIND = DP) :: aux
    !
    aux = ABS(x)**(-m)
    lst = ( aux + SQRT(aux) * (gamma**(-m)) )**(-1.0_DP/REAL(m,DP))
    IF (x < 0.0_DP) lst = - lst
    !
  END FUNCTION lst
  !
  !
  ELEMENTAL FUNCTION sox(x)
    !
    !     ANALYTICAL 1D THO LST s(x)/x
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  sox
    !
    sox = ( 1.0_DP + (gamma**(-m))*SQRT(ABS(x))**m )**(-1.0_DP/REAL(m,DP))
    !
  END FUNCTION sox
  !
  ELEMENTAL FUNCTION d_lst(x)
    !
    !     ANALYTICAL 1D THO LST (1st derivative)
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  d_lst
    !
    d_lst = (sox(x)/2.0_DP)*( 1.0_DP + sox(x)**m )
    !
  END FUNCTION d_lst
  !
  ELEMENTAL FUNCTION d2_lst(x)
    !
    !     ANALYTICAL 1D THO LST (2nd derivative)
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  d2_lst
    !
    REAL(KIND = DP) :: aux1, aux2, vx
    !
    vx = ABS(x)
    aux1 = SQRT(vx)**(m-2)
    aux2 = SQRT(vx)**m
    !
    d2_lst = -(sox(x)/4.0_DP)*(aux1)/(gamma**m + aux2)*( 1.0_DP + (m+1.0_DP)*sox(x)**m )
    IF (x < 0.0_DP) d2_lst = - d2_lst
    !
  END FUNCTION d2_lst
  !
  !
  SUBROUTINE THO_1D_BASIS(apar, NdimH, Iprint)
    !     
    !     COMPUTES A 1D ANALYTIC THO BASIS
    !
    !     INPUT  :: 
    !               apar    --> LENGTH SCALE OF THE PROBLEM
    !               NdimH   --> DIMENSION + 1 OF THE TRANSFORMED HARMONIC BASIS
    !               THO constants (INCLUDE THEM AS ARGUMENTS)
    !
    !     OUTPUT :: THar_1b_Bas  --> MATRIX WITH TRANSFORMED HARMONIC BASIS
    !
    !     FORMAT :: Iprint  --> VERBOSITY CONTROL
    !
    !     $Id: module_1body_tho.f90,v 1.1 2013/07/01 15:49:24 curro Exp $
    !
    !     by Currix TM.
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
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: THO_norm_test
    REAL(KIND = DP) :: PI14, apar2, ainteg, Error
    INTEGER(KIND = I4B) :: kx, Ifail, Ierr, X_dim
    !
    !
    ! GAMMA PARAMETER
    gamma = ratio/apar
    !
    IF (Iprint > 1) PRINT*, "BUILDING THO BASIS"
    !
    ! DEFINING SCALING FUNCTION 
    S_x = 0.0_dp
    der_S_x = 0.0_dp
    der2_S_x = 0.0_dp
    !
    S_x = lst(X_grid)
    !
    S_over_x = sox(X_grid)
    !
    der_S_x = d_lst(X_grid)
    !
    der2_S_x = d2_lst(X_grid)
    !
    PI14 = SQRT(SQRT(PI_D))
    apar2 = apar*apar
    !
    !     HO n = 0
    THar_1b_Bas(:,1) = SQRT(apar)/(PI14)*EXP(-apar2*S_x*S_x/2.0_DP)
    !
    !     HO n = 1
    IF (NdimH > 1) THar_1b_BAS(:,2) = (SQRT(apar)/(PI14))*SQRT(2.0_DP)*apar*S_x*EXP(-apar2*S_x*S_x/2.0_DP)
    !
    !    RECURRENCE RELATION (WATCH OUT THE INDEXES :: MATRIX START AT 1, NOT 0)
    DO kx = 2, NdimH-1
       THar_1b_Bas(:,kx+1) = &
            SQRT(2.0_DP/(1.0_DP*kx))*apar*S_x*THar_1b_Bas(:,kx) - &
            SQRT((1.0_DP*(kx-1))/(1.0_DP*kx))*THar_1b_Bas(:,kx-1)
    ENDDO
    !
    !
    !     TESTING NORMALIZATION (DESTROYS HARBAS(X_GRID,DIMENSION + 2 )!!!)
    IF (Iprint > 1) THEN
       !
       ! DEFINE THO_norm_test
       X_dim = SIZE(X_GRID)
       ALLOCATE(THO_norm_test(1:X_dim), STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "THO_norm_test allocation request denied."
          STOP
       ENDIF
       !
       DO kx = 1, NdimH
          THO_norm_test = THar_1b_Bas(:,kx)*THar_1b_Bas(:,kx)*der_S_x(:)
          !
          Ifail = 0
          !
          CALL D01GAF(X_GRID, THO_norm_test, X_dim, ainteg, Error, Ifail)
          PRINT*, "THO FUNCTION ", kx, " NORMALIZATION", ainteg
       ENDDO
       DEALLOCATE(THO_norm_test, STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "THO_norm_test deallocation request denied."
          STOP
       ENDIF
       !
       WRITE(*,*) "DONE"
    ENDIF
    !     
    !     
  END SUBROUTINE THO_1D_BASIS
  !
  !
  SUBROUTINE THARDIAG(apt, Iprint)
    !     
    !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A HARMONIC BASIS
    !
    !     INPUT  :: 
    !               apt     --> PROBLEM LENGTH SCALE (fm^-1)
    !
    !     FORMAT :: IPRINT  --> VERBOSITY CONTROL
    !
    !  $Id: module_1body_tho.f90,v 1.1 2013/07/01 15:49:24 curro Exp $
    !
    !     by Currix TM.
    !
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
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE ::  Pot_vec, Kin_vec_1, Kin_vec_2, Kin_vec_3, Kin_vec_4
    !
    !
    INTEGER(KIND = I4B) ::  I, J, KX, IFAIL, Ierr
    REAL(KIND = DP) ::  AI, AJ, AHAM, ERROR, K1, K2, K3, K4, Kinetic_En
    !
    !
    IF (Iprint > 1) PRINT*, "BUILDING HAMILTONIAN MATRIX"
    !     
    !   
    ALLOCATE(Ham_Mat(1:dim_1b_THO,1:dim_1b_THO), STAT = Ierr)
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
    ALLOCATE(Kin_vec_1(1:dim_X), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_1 allocation request denied."
       STOP
    ENDIF
    !
    ALLOCATE(Kin_vec_2(1:dim_X), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_2 allocation request denied."
       STOP
    ENDIF
    !
    ALLOCATE(Kin_vec_3(1:dim_X), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_3 allocation request denied."
       STOP
    ENDIF
    !
    ALLOCATE(Kin_vec_4(1:dim_X), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_4 allocation request denied."
       STOP
    ENDIF
    ! 
    Ham_Mat = 0.0_DP
    !
    !
    !     HAMILTONIAN MATRIX
    ! 
    DO J = 1, dim_1b_THO
       AJ = (J-1)*1.0_dp
       !    
       DO I = 1, dim_1b_THO
          AI = (I-1)*1.0_dp
          !     KINETIC ENERGY (ADIMENSIONAL UNITS)
          Kin_vec_1 = 0.0_dp
          Kin_vec_2 = 0.0_dp
          Kin_vec_3 = 0.0_dp
          Kin_vec_4 = 0.0_dp
          !
          !same parity case
          parity_if : IF (mod(I,2) == mod(J,2)) THEN
             !
             ! Four Terms
             !
             ! Term 1
             Kin_vec_1 = ((der_S_x(:))**(-1))*((der2_S_x(:))**2)*THar_1b_Bas(:,I)*THar_1b_Bas(:,J)
             !
             ! Term 2
             IF ( J == 1 ) THEN
                Kin_vec_2 = der_S_x(:)*der2_S_x(:)*THar_1b_Bas(:,I)*(-SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
             ELSE
                Kin_vec_2 = der_S_x(:)*der2_S_x(:)*THar_1b_Bas(:,I)*(SQRT(AJ/2.0_dp)*THar_1b_Bas(:,J-1) -&
                     SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
             ENDIF
             !
             ! Term 3
             IF ( I == 1) THEN
                Kin_vec_3 = der_S_x(:)*der2_S_x(:) * THar_1b_Bas(:,J) * (-SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1))
             ELSE 
                Kin_vec_3 = der_S_x(:)*der2_S_x(:) * THar_1b_Bas(:,J) * (SQRT(AI/2.0_dp)*THar_1b_Bas(:,I-1) -&
                     SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1))
             ENDIF
             !
             ! Term 4
             IF (I==1) THEN
                IF(J==1) THEN
                   Kin_vec_4 = ((der_S_x(:))**3) * &
                        ( - SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1)) * &
                        ( - SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
                ELSE
                   Kin_vec_4 = ((der_S_x(:))**3) * &
                        ( - SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1)) * &
                        (SQRT(AJ/2.0_dp)*THar_1b_Bas(:,J-1) - SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
                ENDIF
             ELSE
                IF(J==1) THEN
                   Kin_vec_4 = ((der_S_x(:))**3) * &
                        (SQRT(AI/2.0_dp)*THar_1b_Bas(:,I-1) - SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1)) * &
                        ( - SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
                ELSE
                   Kin_vec_4 = ((der_S_x(:))**3) * &
                        (SQRT(AI/2.0_dp)*THar_1b_Bas(:,I-1) - SQRT((AI+1)/2.0_dp)*THar_1b_Bas(:,I+1)) * &
                        (SQRT(AJ/2.0_dp)*THar_1b_Bas(:,J-1) - SQRT((AJ+1)/2.0_dp)*THar_1b_Bas(:,J+1))
                ENDIF
             ENDIF
             !
          ENDIF parity_if
          !
          !     INTEGRATION 
          IFAIL = 0
          CALL D01GAF(X_Grid,Kin_vec_1,dim_X,K1,ERROR,IFAIL)
          IFAIL = 0
          CALL D01GAF(X_Grid,Kin_vec_2,dim_X,K2,ERROR,IFAIL)
          IFAIL = 0
          CALL D01GAF(X_Grid,Kin_vec_3,dim_X,K3,ERROR,IFAIL)
          IFAIL = 0
          CALL D01GAF(X_Grid,Kin_vec_4,dim_X,K4,ERROR,IFAIL)
          !
          Kinetic_en =((1.0_DP/8.0_DP)*h_sq_over_m)*K1 +((apt/4.0_DP)*h_sq_over_m)*K2 &
               +((apt/4.0_DP)*h_sq_over_m)*K3 + ( ((apt*apt)/2.0_DP )*h_sq_over_m)*K4
          !
          !     POTENTIAL ENERGY
          Pot_vec = Potf(X_Grid)*THar_1b_BAS(:,I)*THar_1b_BAS(:,J)*der_S_x(:)
          !
          !     INTEGRATION
          IFAIL = 0
          CALL D01GAF(X_Grid,Pot_vec,dim_X,AHAM,ERROR,IFAIL)
          !
          Ham_Mat(I,J) = Kinetic_en + AHAM
          !
       ENDDO
       !
    ENDDO
    !
    !
    DEALLOCATE(Pot_vec, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Pot_vec deallocation request denied."
       STOP
    ENDIF
    !
    DEALLOCATE(Kin_vec_1, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_1 deallocation request denied."
       STOP
    ENDIF
    !
    DEALLOCATE(Kin_vec_2, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_2 deallocation request denied."
       STOP
    ENDIF
    !
    DEALLOCATE(Kin_vec_3, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_3 deallocation request denied."
       STOP
    ENDIF
    !
    DEALLOCATE(Kin_vec_4, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Kin_vec_4 deallocation request denied."
       STOP
    ENDIF
    !     DIAGONALIZATION USING LAPACK
    CALL LA_SYEVR(A=Ham_Mat, W=AVAL_1b_THar, JOBZ='V', UPLO='L')
    !
    !
    !  EIGENVECTOR MATRIX
    !
    AVEC_1b_THar = Ham_Mat
    !
    AVEC_1b_THar_X = 0.0_DP
    !
    DO KX = 1, dim_X
       DO I = 1, dim_1b_THO
          DO J = 1, dim_1b_THO
             AVEC_1b_THar_X(KX,I) = AVEC_1b_THar_X(KX,I) + Ham_Mat(J,I)*THar_1b_BAS(KX,J)*SQRT(der_S_x(KX))
          ENDDO
       ENDDO
    ENDDO
    !
    !
    DEALLOCATE(Ham_Mat, STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "Ham_Mat deallocation request denied."
       STOP
    ENDIF
    !    
  END SUBROUTINE THARDIAG
  !
  SUBROUTINE output_1b_results(Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS)
    !
    IMPLICIT NONE
    !
    INTERFACE  Potf
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
    INTEGER(KIND = I4B), INTENT(IN) :: Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS
    !
    !
    ! Local variables
    INTEGER(KIND = I4B) :: i, kx
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "tho_1D_1b"
    !
    !
    IF (Iprint > 0) THEN
       PRINT*, "EIGENVALUES IN A ", dim_1b_THO, " DIM T-HARMONIC BASIS"
       DO I = 1, dim_1b_THO
          PRINT*, I, AVAL_1b_THar(I)
       ENDDO
    ENDIF
    !
    ! SAVING GROUND STATE, XGRID, AND POTENTIAL INFORMATION
    !IF (i_GS == 1) THEN
    !   !
    !   OPEN(UNIT=9,FILE='gs_wavefunction.dat')
    !   !         
    !   WRITE(9,*) "# ", reduced_mass, "         # Reduced Mass"
    !   WRITE(9,*) "# ", Param_pot, "        # Potential parameters"
    !   WRITE(9,*) "# ", X_min, X_max, Delta_X, "      # xmin xmax dx"
    !   WRITE(9,*) "# ", Aval_1b_THar(1), "     # G.S. energy"
    !   !
    !   DO kx = 1, dim_X
    !      WRITE(9,*) X_grid(kx), Avec_1b_THar(kx,1)
    !   ENDDO
    !   !  
    !   CLOSE(9)
    !   !
    !ENDIF
    !
    ! SAVING T-HARMONIC BASIS
    IF (i_save_1b_BAS == 1) THEN
       file = 'basis'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       IF ( dim_1b_THO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       IF ( dim_1b_THO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(70,*) "# THO  dim_1b_THO = ", dim_1b_THO, " Box radius = ", X_max, " fm"
       WRITE(70,*) "#Grid    T-Harmonic basis"
       DO kx = 1, dim_X
          WRITE(70,11) X_grid(kx), THar_1b_Bas(kx,1:dim_1b_THO)*SQRT(der_S_x(kx))
       ENDDO
       CLOSE(UNIT = 70)
    ENDIF
    !
    ! SAVING EIGENVECTORS
    IF (i_save_1b_WF == 1) THEN
       file = 'eigenvectors'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       IF ( dim_1b_THO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       IF ( dim_1b_THO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(71,*) "# THO  dim_1b_THO = ", dim_1b_THO, " Box radius = ", X_max, " fm"
       WRITE(71,*) "#Grid  Eigenvectors"
       DO kx = 1, dim_X
          WRITE(71,11) X_grid(kx), Avec_1b_THar_X(kx,1:dim_1b_THO)
       ENDDO
       CLOSE(UNIT = 71)
    ENDIF
    !
    ! SAVING EIGENVECTOR DERIVATIVES
    !IF (i_save_WF == 1) THEN
    !   file = 'eigvec_der'
    !   WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
    !   IF ( dim_THO < 10) THEN !to avoid spaces
    !      WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
    !   ENDIF
    !   IF ( dim_THO > 99) THEN 
    !      WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
    !   ENDIF
    !   OPEN(UNIT = 72, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !   WRITE(72,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
    !   WRITE(72,*) "#Grid  Eigenvector derivatives"
    !   DO kx = 1, dim_X
    !      WRITE(72,11) X_grid(kx), Avec_THO_Der_X(kx,1:dim_THO)
    !   ENDDO
    !   CLOSE(UNIT = 72)
    !ENDIF
    !
    !
    ! SAVING ENERGIES
    IF (i_save_1b_EN == 1) THEN
       file = 'eigenvalues'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       IF ( dim_1b_THO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       IF ( dim_1b_THO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(73,*) "# THO  dim_1b_THO = ", dim_1b_THO, " Box radius = ", X_max, " fm"
       WRITE(73,*) "# Eigenvalues"
       DO I = 1, dim_1b_THO
          WRITE(73,10) I, Aval_1b_THar(I)
       ENDDO
       CLOSE(UNIT = 73)
       !
       file = 'pot_eigvec'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       IF ( dim_1b_THO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       IF ( dim_1b_THO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(74,*) "# THO  dim_1b_THO = ", dim_1b_THO, " Box radius = ", X_max, " fm"
       WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
       !
       file = 'pot_eigvec2'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       IF ( dim_1b_THO < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       IF ( dim_1b_THO > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_THO
       ENDIF
       OPEN(UNIT = 75, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(75,*) "# THO  dim_1b_THO = ", dim_1b_THO, " Box radius = ", X_max, " fm"
       WRITE(75,*) "#Grid    Potential    10*Eigenfunctions^2+eigenvalue"
       DO kx = 1, dim_X
          WRITE(74,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_1b_THar_X(kx,1:dim_1b_THO) + Aval_1b_THar(1:dim_1b_THO)
          WRITE(75,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_1b_THar_X(kx,1:dim_1b_THO)**2 + Aval_1b_THar(1:dim_1b_THO)
       ENDDO
       CLOSE(UNIT = 74)
       CLOSE(UNIT = 75)
    ENDIF
    !
    !
10  FORMAT (1X,I6,1X,E16.8)
11  FORMAT (1X,E14.6,1X,300E16.8)
    !
  END SUBROUTINE output_1b_results
  !
  !
  !
END MODULE module_1body_tho
  
