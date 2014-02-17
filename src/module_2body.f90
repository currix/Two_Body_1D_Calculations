MODULE module_2body
  !
  USE nrtype
  USE vardef_ri
  !
  IMPLICIT NONE
  !
  ! Residual Interaction 
  !
  ! M_ri_param -> MAXIMUM NUMBER OF RESIDUAL INTERACTION PARAMETERS
  INTEGER(KIND = I4B), PARAMETER :: M_ri_param = 5
  ! POTENTIAL PARAMETER VALUES AND UNITS
  REAL(KIND = DP), DIMENSION(1:M_ri_param) :: param_ri
  !
  !
  ! Two body variables
  !
  INTEGER(KIND = I4B) :: n_bound_bas ! Number of bound states in the basis
  REAL(KIND = DP) :: E_threshold     ! Threshold energy
  INTEGER(KIND = I4B) :: dim_bound   ! number of bound states in the core
  INTEGER(KIND = I4B) :: dim_2body   ! 2 body basis dimension
  INTEGER(KIND = I4B) :: dim_eig2body! # of 2-body eigenvectors considered
  INTEGER(KIND = I4B) :: i_pri, i_flag_graph, np3d ! Format variables
  INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE :: two_body_bas
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: rho_c
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: aval_2b
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: avec_2b
  REAL(KIND = DP), DIMENSION(:,:,:), ALLOCATABLE :: avec_2b_X
  !
CONTAINS
  !
  SUBROUTINE bas_2_body(dim_1b, E_threshold, Aval_1b, dim_bound, dim_2b, two_body_bas, Iprint)
    !
    ! Two-body basis building subroutine
    !
    ! Input  : dim_1b, E_threshold, Aval_1b, Iprint
    ! Output : dim_bound, dim_2b, two_body_bas
    !
    !  by Currix TM
    !
    ! ARGUMENTS
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_1b, Iprint
    REAL(KIND = DP), INTENT(IN) :: E_threshold
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Aval_1b
    !
    INTEGER(KIND = I4B), INTENT(OUT) :: dim_bound, dim_2b
    INTEGER(KIND = I4B), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: two_body_bas
    !
    ! LOCAL VARIABLES
    !
    INTEGER(KIND = I4B) :: N_bound, I_index, J_index, I_count, dim_c, Ierr
    !
    IF (Iprint > 1) WRITE(*,*) "Building 2-BODY BASIS indexes"
    !  
    ! Number of bound states
    N_bound = 0
    DO I_index = 1, dim_1b
       IF (Aval_1b(I_index) > 0.0_DP) THEN
          N_bound = I_index - 1
          EXIT
       END IF
    END DO
    !     Core dimension
    dim_bound = N_bound - n_bound_bas
    !
    !     Continuum dimension and basis
    !
    ! Include states under E_threshold
    dim_c = 0
    IF (E_threshold == 0.0_DP) THEN 
       dim_c = dim_1b - dim_bound
    ELSE
       DO I_index = dim_bound + 1, dim_1b
          IF (Aval_1b(I_index) > E_threshold) EXIT
          dim_c = dim_c + 1
       ENDDO
    ENDIF
    !
    !
    ! Continuum dimension and basis (two step process)
    !
    ! Continuum dimension
    dim_2b = 0
    IF (E_threshold == 0.0D0) THEN ! No threshold energy
       dim_2b = ( dim_c*(dim_c + 1) ) / 2
    ELSE
       DO I_index = 1, dim_c
          DO J_index = I_index, dim_c
             IF (Aval_1b(dim_bound + I_index)+Aval_1b(dim_bound + J_index) <= E_threshold) dim_2b = dim_2b + 1
          ENDDO
       ENDDO
    ENDIF
    !
    ! Continuum basis 
    ALLOCATE(two_body_bas(1:dim_2b,2), STAT = Ierr)
    IF (Ierr /= 0) THEN
       PRINT*, "X_grid allocation request denied."
       STOP
    ENDIF
    ! Include 2body states that accomplish E_1 + E_2 <= E_c
    I_count = 0
    DO I_index = 1, dim_c
       DO J_index = I_index, dim_c
          IF (E_threshold == 0.0D0) THEN ! No threshold energy
             I_count = I_count + 1
             two_body_bas(I_count,1) = I_index
             two_body_bas(I_count,2) = J_index
          ELSE
             IF (Aval_1b(dim_bound + I_index)+Aval_1b(dim_bound + J_index) <= E_threshold) THEN
                I_count = I_count + 1
                two_body_bas(I_count,1) = I_index
                two_body_bas(I_count,2) = J_index
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    IF (Iprint > 0) THEN 
       WRITE(*,*) "Core dimension = ", dim_bound
       WRITE(*,*) "1-body basis dimension E < E_th = ", dim_c
       WRITE(*,*) "2-body basis dimension = ", dim_2b
    ENDIF
    !
    !
  END SUBROUTINE bas_2_body
  !
  !
  !
  SUBROUTINE core_matt_den(X_grid, Avec_1b_x, dim_bound, rho_c, Iprint)
    !     
    ! COMPUTE CORE DENSITY MATTER
    !
    !     INPUT  :: 
    !               X_grid      --> VECTOR WITH X GRID
    !               Avec_1b_v   --> ARRAY WITH THE 1BODY BOUND EIGENVECTORS (CORE)
    !               dim_bound    --> 1 BODY BOUND STATES DIMENSION (CORE)
    !
    !     OUTPUT :: 
    !               RHOC    --> VECTOR WITH CORE MATTER DENSITY
    !
    !     FORMAT :: IPRINT  --> VERBOSITY CONTROL
    !
    !     by Currix TM.
    !
    !
    !
    !     ARGUMENTS
    INTEGER, INTENT(IN) :: dim_bound, Iprint
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) ::  X_grid
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  Avec_1b_X
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(OUT) :: rho_c
    !
    !
    INTEGER :: I_index
    INTEGER :: KX
    !
    REAL(KIND = DP) :: FNOP
    !
    !
    IF (IPRINT > 1) WRITE(*,*) "BUILDING 1-BODY CORE MATTER DENSITY ... "
    !
    rho_c = 0.0_DP
    !     
    DO I_index = 1, dim_bound
       !
       DO KX = 1, dim_X
          !
          FNOP = 2.0_DP       ! TAKES INTO ACCOUNT TWO PARTICLES PER BOUND LEVEL
          rho_c(KX) = rho_c(KX) + FNOP*Avec_1b_X(KX, I_index)**2
          !
       ENDDO
       !     
    ENDDO
    !     
    IF (IPRINT > 2) THEN
       !
       WRITE(*,*) "CORE MATTER DENSITY"
       !
       DO KX = 1, dim_X
          !     
          WRITE(*,*) X_grid(KX), rho_c(KX)
          !     
       ENDDO
    ENDIF
    !
    IF (IPRINT > 1) WRITE(*,*) "BUILDING 1-BODY CORE MATTER DENSITY ... DONE"
    !
  END SUBROUTINE core_matt_den
  !
  SUBROUTINE diag_2body_ham(X_grid, rho_c, avec_1b_X, dim_core, two_body_bas, &
       Aval_1b, dim_2b, Aval_2b, avec_2b, Ipri, Iprint, Iflag)
    !     
    !     Compute and diagonalize 2-body 1d hamiltonian using 
    !     a contact residual interaction
    !
    !     INPUT  :: 
    !
    !       X_grid       --> VECTOR WITH X GRID
    !       rho_c        --> VECTOR WITH CORE MATTER DENSITY
    !       avec_1b_X    --> ARRAY WITH THE 1BODY EIGENVECTORS (X dependence)
    !       dim_core     --> 1-BODY STATES IN THE CORE
    !       two_body_bas --> ARRAY WITH n1 n2 VALUE FOR EACH 2BODY BASIS ELEMENT
    !       Aval_1b      --> VECTOR WITH 1-BODY ENERGIES
    !       dim_2b       --> 2-BODY BASIS DIMENSION 
    !
    !     OUTPUT :: 
    !
    !       Avec_2b    --> ARRAY WITH EIGENVECTORS (1-BODY DEPENDENCE)
    !       Aval_2b    --> VECTOR WITH 2-BODY ENERGIES
    !
    !     FORMAT :: 
    !
    !       Ipri  --> DISPLAY RESIDUAL INTERACTION DIAGONAL MATRIX ELEMENTS
    !       Iprint  --> VERBOSITY CONTROL
    !       Iflag        --> INTEGER :: 0 EIGENVALUES 
    !                              1 EIGENVALUES AND GS EIGENVECTOR (R1 = R2)
    !                              2 EIGENVALUES AND GS EIGENVECTOR
    !
    !     by Currix TM.
    !
    !
    !     LAPACK 95
    USE LA_PRECISION, ONLY: WP => DP
    USE F95_LAPACK, ONLY: LA_SYEVR
    !
    !     ARGUMENTS
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) ::  X_grid, rho_c
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  avec_1b_X
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) ::  Aval_1b
    INTEGER, INTENT(IN) :: dim_core, dim_2b
    INTEGER, INTENT(IN) :: Ipri, Iprint, Iflag
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(OUT)   :: aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: avec_2b
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index_I, index_J, N_1, N_2, N_1_P, N_2_P
    INTEGER(KIND = I4B) :: index_X, IFAIL
    REAL(KIND = DP) :: WFI, WFJ, SQFACI, SQFACJ, AIN, ERROR
    REAL(KIND = DP), DIMENSION(dim_X) :: ri_1, ri_2
    !
    IF (IPRINT > 1) WRITE(*,*) "BUILDING 2-BODY HAMILTONIAN MATRIX ..."
    !     
    avec_2b = 0.0_DP
    !
    !     HAMILTONIAN MATRIX BUILDING
    !
    !
    ! MAIN LOOP
    indexI: DO index_I = 1, dim_2b
       !
       N_1_P = two_body_bas(index_I,1)
       N_2_P = two_body_bas(index_I,2)
       !
       IF (N_1_P == N_2_P) THEN
          SQFACI = 1.0_DP
       ELSE
          SQFACI = SQRT(2.0_DP)
       ENDIF
       !     
       indexJ: DO index_J = index_I, dim_2b
          !     
          N_1 = two_body_bas(index_J,1)
          N_2 = two_body_bas(index_J,2)
          !     
          IF (N_1 == N_2) THEN
             SQFACJ = 1.0_DP
          ELSE
             SQFACJ = SQRT(2.0_DP)
          ENDIF
          !
          !     DIAGONAL PART (SUM OF 1-BODY ENERGIES)
          IF (N_1_P == N_1 .AND. N_2_P == N_2) &
               avec_2b(index_J,index_I) = avec_2b(index_J,index_I) + &
               SQFACI*SQFACJ*(aval_1b(dim_core + N_1) + aval_1b(dim_core + N_2))/2.0_DP
          IF (N_1_P == N_2 .AND. N_2_P == N_1) & 
               avec_2b(index_J, index_I) = avec_2b(index_J, index_I) + &
               SQFACI*SQFACJ*(aval_1b(dim_core + N_1) + aval_1b(dim_core + N_2))/2.0_DP
          !
          !     SELECTION RULE N1 + N2 + N1P + N2P = EVEN NUMBER
          IF (MOD(N_1 + N_2 + N_1_P + N_2_P,2) == 0) THEN
             !
             !ri_1(:) = avec_1b_X(:,N_1_P)*avec_1b_X(:,N_2_P)* &
             !     avec_1b_X(:,N_1)*avec_1b_X(:,N_2)
             DO index_X = 1, dim_X
                !
                WFI = avec_1b_X(index_X, dim_core + N_1_P)*avec_1b_X(index_X, dim_core + N_2_P)
                WFJ = avec_1b_X(index_X, dim_core + N_1)*avec_1b_X(index_X, dim_core + N_2)
                !     
                !     RESIDUAL INTERACTION
                !     1ST TERM
                ri_1(index_X) = WFI*WFJ
                !     
                !     RESIDUAL INTERACTION
                !     2ND TERM
                ri_2(index_X) =  WFI*WFJ*((rho_c(index_X)/param_ri(3))**param_ri(4))
                !     
             ENDDO
             !
             !     INTEGRATION (1ST TERM)
             IFAIL = 0
             CALL D01GAF(X_grid,ri_1,dim_X,AIN,ERROR,IFAIL)
             !     
             IF (Iprint > 3) THEN 
                WRITE(*,*) "1ST TERM RESIDUAL INTERACTION"
                WRITE(*,*) "N1P, N2P = ", N_1_P, N_2_P, " N1, N2 = ,", &
                     N_1, ",", N_2, " I, J = ", &
                     index_I, ",", index_J
                WRITE(*,*) param_ri(1), AIN, param_ri(1)*AIN
             ENDIF
             !
             avec_2b(index_J,index_I) = avec_2b(index_J,index_I) + &
                  SQFACI*SQFACJ*param_ri(1)*AIN
             !
             IF (ipri == 1 .AND. index_J == index_I) &
                  WRITE(*,*) index_I,"-TH ELEMENT VOLUME TERM RES. INT. = ", &
                  SQFACI*SQFACJ*param_ri(1)*AIN
             !
             !    INTEGRATION (2ND TERM)
             IFAIL = 0
             CALL D01GAF(X_grid,ri_2,dim_X,AIN,ERROR,IFAIL)
             !
             !     
             IF (IPRINT > 3) THEN 
                WRITE(*,*) "2ND TERM RESIDUAL INTERACTION"
                WRITE(*,*) param_ri(2), AIN, param_ri(2)*AIN
             ENDIF
             !     
             avec_2b(index_J,index_I) = avec_2b(index_J,index_I) + & 
                  SQFACI*SQFACJ*param_ri(2)*AIN
             !
             IF (ipri == 1 .AND. index_J == index_I) &
                  WRITE(*,*) index_I,"-TH ELEMENT DENSITY TERM RES.INT. = ", &
                  SQFACI*SQFACJ*param_ri(2)*AIN
             !
          ELSE
             !
             avec_2b(index_J,index_I) = 0.0_DP
             !
          ENDIF
          !
       ENDDO indexJ
       !  
    ENDDO indexI
    !     
    !
    IF (IPRINT >= 1) WRITE(*,*) "2body Hamiltonian diagonalization"
    !     
    !     DIAGONALIZATION USING LAPACK
    IF (IFLAG >= 1) THEN
       CALL LA_SYEVR(A=avec_2b, W=aval_2b, JOBZ='V', UPLO='L')
    ELSE
       CALL LA_SYEVR(A=avec_2b, W=aval_2b, JOBZ='N', UPLO='L')
    ENDIF
    !
    IF (IPRINT > 1) THEN
       WRITE(*,*) "2P EIGENVALUES, ", dim_2b, " DIM BASIS"
       DO index_I = 1, dim_2b
          WRITE(*,*) index_I, aval_2b(index_I)
       ENDDO
    ENDIF
    !
  END SUBROUTINE diag_2body_ham
  !
  !
  SUBROUTINE two_body_results(X_grid, avec_1b_X, aval_1b, aval_2b, avec_2b, &
       dim_core, dim_2b, two_body_bas, num_2b_states, avec_2b_X, &
       Iflag, Np3d, Iprint)
    !
    !     2-BODY OUTPUT AND BUILDING SPATIAL DEPENDENCE OF 2-BODY EIGENVECTORS
    !
    !     INPUT  :: 
    !           X_grid      --> VECTOR WITH X GRID
    !           avec_1b_X   --> ARRAY WITH 1BODY EIGENVECTORS
    !           aval_1b     --> VECTOR WITH 1-BODY ENERGIES
    !           aval_2b     --> VECTOR WITH 2-BODY ENERGIES
    !           avec_2b     --> ARRAY WITH EIGENVECTORS (2-BODY DEPENDENCE)
    !           dim_core     --> 1-BODY STATES IN THE CORE
    !           dim_2b       --> 2-BODY BASIS DIMENSION 
    !           num_2b_states--> Number of 2-body eigenstates whose X1, X2 dependence is computed.
    !           two_body_bas --> ARRAY WITH n1 n2 VALUE FOR EACH 2BODY BASIS ELEMENT
    !
    !
    !     OUTPUT :: 
    !               avec_2b_X  --> ARRAY WITH EIGENVECTORS (X DEPENDENCE)
    !
    !     FORMAT :: 
    !
    !       Ipri    --> DISPLAY RESIDUAL INTERACTION DIAGONAL MATRIX ELEMENTS
    !       Iprint  --> VERBOSITY CONTROL
    !       Iflag   --> INTEGER :: 0 SAVE NOTHING 
    !                              1 ALL EIGENVALUES (FILES 4X) 
    !                              2 EIGENVALUES AND GS EIGENVECTOR (R1 = R2)
    !                              3 EIGENVALUES AND GS EIGENVECTOR
    !       NP3D    --> ONLY EACH NP3D POINTS ARE SAVED IN THE 3D 
    !                       EIGENVECTOR REPRESENTATION
    !
    !
    !     OUTPUT UNITS
    !
    !     by Currix TM & Moschini GmbH
    !
    !
    !
    !
    !     ARGUMENTS
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) ::  X_grid
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) ::  avec_1b_X, avec_2b
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) ::  Aval_1b, Aval_2b
    INTEGER, INTENT(IN) :: dim_core, dim_2b, num_2b_states
    INTEGER, INTENT(IN) :: Iprint, Iflag, Np3d
    !
    REAL(KIND = DP), DIMENSION(:,:,:), INTENT(OUT) :: avec_2b_X
    !
    ! Local variables
    INTEGER index_I, index_J, N_1, N_2
    INTEGER index_X1, index_X2
    DOUBLE PRECISION SQFACJ
    CHARACTER(LEN=100) :: format_1, format_2
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "1D_2body_"
    !
    !
    IF (Iprint > 1) WRITE(*,*) "BUILDING 2-BODY EIGENVECTORS"
    !
    ! Prepare FORMAT statements
110 FORMAT(1X, I7, 1X, F20.8)
    WRITE(format_1, '("(1X, ", I5, "F20.8)")'), num_2b_states
    WRITE(format_2, '("(1X, I7, 1X, ", I7, "E20.10)")'), num_2b_states + 1
    !
    ! Standard output :: Eigenvalues
    WRITE(UNIT=*,FMT='(1X, I6, 2X, F20.8)') dim_2b, 2.0_DP*Aval_1b(dim_core+1)
    WRITE(UNIT=*, FMT=format_1) (Aval_2b(index_I), index_I = 1, num_2b_states)
    !
    !
    IF (Iflag >= 1) THEN
       ! Write 2-BODY ENERGIES
       file = 'eigenvalues'
       WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
       OPEN(UNIT = 45, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(45,*) "#1D 2body eigenenergies"
       DO index_I = 1, dim_2b
          WRITE(UNIT=45,FMT=110) index_I, aval_2b(index_I)
       ENDDO
       CLOSE(UNIT = 45)
    ENDIF
    !
    !
    IF (Iflag >= 2) THEN
       !
       ! Write num_2b_states states data for correlation plots, with components as a function of uncorrelated energy
       file = 'corr_plots'
       WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
       OPEN(UNIT = 100, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(100,*) "# Write num_2b_states states data for correlation plots, with components as a function of uncorrelated energy"
       DO index_I = 1, dim_2b
          WRITE(UNIT=100, FMT=format_2) index_I,  &
               aval_1b(dim_core + two_body_bas(index_I,1))+aval_1b(dim_core + two_body_bas(index_I,2)), &
               (avec_2b(index_I,index_J), index_J = 1, num_2b_states)
       ENDDO
       CLOSE(UNIT = 100)
       !
       !     X1 = X2 CASE (PLANE BISECTOR)
       file = 'uncorr_eigenvec'
       WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
       OPEN(UNIT = 55, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(55,*) "# X1 = X2 CASE (PLANE BISECTOR)"
       DO index_X1 = 1, dim_X
          !
          N_1 = two_body_bas(1,1)
          N_2 = two_body_bas(1,2)
          !
          !     55  ::  X1 Z Z^2 uncorrelated minimum energy state
          WRITE(UNIT=55,FMT=*) X_grid(index_X1), avec_1b_X(index_X1, dim_core + N_1)*avec_1b_X(index_X1, dim_core + N_2)
          !     
       ENDDO
       CLOSE(UNIT = 55)
       !
       !     X1 = X2 CASE (PLANE BISECTOR)
       !     55 + I ::  X1 Z Z^2
       !     num_2b_states :: Maximum number of 2-body eigenvectors to deal with.
       !
       ! Initialize array     
       avec_2b_X = 0.0_DP
       !     
       file = 'ith_eigenvector'
       !
       DO index_I = 1, num_2b_states
          !
          WRITE(filename, '(A, "_",A,"_",i1, ".dat")') TRIM(prog), TRIM(file), index_I
          OPEN(UNIT = 55+index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
          WRITE(55+index_I,*) "#  55 + i (i = 1, ... , )   x_k  i-th eigenvector(x_k)", &
               " #i-th eigenvector(x_k)**2   XX(KX1), AVEC2B(KX1,1,I), AVEC2B(KX1,1,I)**2 "
          !
          DO index_X1 = 1, dim_X
             !     
             DO index_J = 1, dim_2b
                !
                N_1 = two_body_bas(index_J,1)
                N_2 = two_body_bas(index_J,2)
                !
                IF (N_1.EQ.N_2) THEN
                   SQFACJ = 0.5_DP
                ELSE
                   SQFACJ = 1.0_DP/SQRT(2.0_DP)
                ENDIF
                !
                avec_2b_x(index_X1,1,index_I) = avec_2b_x(index_X1,1,index_I) + &
                     avec_2b(index_J,index_I)*SQFACJ*2.0_DP* &
                     avec_1b_X(index_X1,dim_core + N_1)*avec_1B_X(index_X1,dim_core + N_2)
                !
             ENDDO
             !
             WRITE(UNIT=55+index_I,FMT=*) X_grid(index_X1),avec_2b_X(index_X1,1,index_I),avec_2b_X(index_X1,1,index_I)**2
             !
          ENDDO
          CLOSE(UNIT = 55+index_I)
       ENDDO
    ENDIF
    !     
    IF (Iflag == 3) THEN
       !
       !     75 uncorrelated minimum energy state
       file = 'vector_u77'
       WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
       OPEN(UNIT = 77, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(77,*) "#  x_k   x_l uncorrelated eigenvector(x_k,x_l)", &
            "# XX(KX1),    XX(KX2), BASIS(KX1,1)*BASIS(KX2,1)"
       !
       file = 'vector2_u78'
       WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
       OPEN(UNIT = 78, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(78,*) "# x_k   x_l   uncorrelated eigenvector(x_k,x_l)**2", &
            "#XX(KX1), XX(KX2), (BASIS(KX1,I)*BASIS(KX2,I))**2"
       !
       DO index_X1 = 1, dim_X, NP3D
          !     
          N_1 = two_body_bas(1,1)
          N_2 = two_body_bas(1,2)
          !
          WRITE(77,*) (avec_1b_X(index_X1,dim_core + N_1)*avec_1b_X(index_X2,dim_core + N_2), index_X2=1, dim_X, NP3D)  
          WRITE(78,*) ((avec_1b_X(index_X1,dim_core + N_1)*avec_1b_X(index_X2,dim_core + N_2))**2, index_X2=1, dim_X, NP3D)  
          !
       ENDDO
       !
       CLOSE(UNIT = 77)
       CLOSE(UNIT = 78)
       !
       !
       avec_2b_X = 0.0_DP
       !
       !
       DO index_I = 1, num_2b_states
          !
          file = 'vector_u77+i'
          WRITE(filename, '(A, "_",A,"_",i1, ".dat")') TRIM(prog), TRIM(file), index_I
          OPEN(UNIT = 77+index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
          WRITE(77+index_I,*) "#x_k  x_l            i-th eigenvector(x_k,x_l)", &
               " # XX(KX1),    XX(KX2),               AVEC2B(KX1,KX2,I)"    
          !
          file = 'vector2_u78+i'
          WRITE(filename, '(A, "_",A,"_",i1, ".dat")') TRIM(prog), TRIM(file), index_I
          OPEN(UNIT = 78+index_I, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
          WRITE(78+index_I,*) "#x_k    x_l   i-th eigenvector(x_k,x_l)**2", &
               " #XX(KX1),    XX(KX2),               AVEC2B(KX1,KX2,I)**2 "        
          !
          DO index_X1 = 1, dim_X, NP3D
             !
             DO index_X2 = 1, dim_X, NP3D
                !
                DO index_J = 1, dim_2b
                   !
                   N_1 = two_body_bas(index_J,1)
                   N_2 = two_body_bas(index_J,2)
                   !
                   IF (N_1 == N_2) THEN
                      SQFACJ = 0.5_DP
                   ELSE
                      SQFACJ = 1.0_DP/SQRT(2.0_DP)
                   ENDIF
                   !
                   avec_2b_X(index_X1,index_X2,index_I) = avec_2b_X(index_X1,index_X2,index_I) + &
                        avec_2b(index_J,index_I)*SQFACJ* &
                        ( avec_1b_X(index_X1,dim_core + N_1)*avec_1b_X(index_X2,dim_core + N_2) + & 
                        avec_1b_X(index_X1,dim_core + N_2)*avec_1b_X(index_X2,dim_core + N_1) )
                ENDDO
                !
                !
             ENDDO
             !
             !print*, "avec_2b_X(:,:,1)", avec_2b_X(:,:,1)
             !
             WRITE(77+index_I,*) (avec_2b_X(index_X1,index_X2,index_I), index_X2=1, dim_X, NP3D)  
             WRITE(78+index_I,*) (avec_2b_X(index_X1,index_X2,index_I)**2, index_X2=1, dim_X, NP3D)  
             !
             !
          ENDDO
          !
          CLOSE(UNIT = 77+index_I)
          CLOSE(UNIT = 78+index_I)
          !
       ENDDO
       !     
    ENDIF
    !
    IF (IPRINT.GT.1) WRITE(*,*) "BUILDING 2-BODY EIGENVECTORS ... DONE"
    !
  END SUBROUTINE two_body_results
  !
  !
  SUBROUTINE anom_density(Eigvector, dim_2b, two_body_bas, a_d_value)
    !     
    ! Calculation of the anomalous density
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: Eigvector
    INTEGER(KIND = I4B), INTENT(IN) :: dim_2b
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    !     
    REAL(KIND = DP), INTENT(OUT) :: a_d_value
    !
    !
    !    Local Variables
    INTEGER(KIND = I4B) :: index_I
    !
    a_d_value = 0.0_DP
    !
    DO index_I = 1, dim_2b
       !     
       IF ( two_body_bas(index_I,1) == two_body_bas(index_I,2) ) a_d_value = a_d_value + Eigvector(index_I)
       !
    ENDDO
    !
    a_d_value = ABS(a_d_value)
    !
  END SUBROUTINE anom_density
  !
  SUBROUTINE compute_anom_density(X_grid, aval_2b, avec_2b, avec_2b_X, dim_X, dim_2b, num_2b_states, two_body_bas)
    !     
    !
    ! Two-method calculation of the anomalous density
    !
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, aval_2b
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_2b
    REAL(KIND = DP), DIMENSION(:,:,:), INTENT(IN) :: avec_2b_X
    INTEGER(KIND = I4B), INTENT(IN) :: dim_X, dim_2b, num_2b_states
    INTEGER(KIND = I4B), DIMENSION(:,:), INTENT(IN) :: two_body_bas
    !
    !     
    ! Local variables
    REAL(KIND = DP), DIMENSION(1:dim_X) :: eigvec_2b_X 
    REAL(KIND = DP), DIMENSION(1:dim_2b) :: eigvec_2b
    REAL(KIND = DP) :: tot_anom_density, a_d_value_1, a_d_value_2, ERROR
    INTEGER(KIND = I4B) :: index_I, IFAIL
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "1D_2body_"
    !
    !
    WRITE(UNIT=*,FMT=11) 
11  FORMAT(1X, "ANOMALOUS DENSITY RESULTS")
    !
    DO index_I = 1, num_2b_states
       !
       eigvec_2b(:) =  avec_2b(:,index_I)
       eigvec_2b_X(:)  = avec_2b_X(:,1,index_I)
       !
       CALL D01GAF(X_grid, eigvec_2b_X, dim_X, a_d_value_1, ERROR, IFAIL)
       !
       CALL anom_density(eigvec_2b, dim_2b, two_body_bas, a_d_value_2)
       !
       !                             a_d_value_2 : algebraic result Eq. (35) notes \label{T0bis}
       !                             a_d_value_1 : integration result Eq. (31) notes \label{T0}
       WRITE(UNIT=*,FMT=12) index_I, a_d_value_2, a_d_value_1
       !
    ENDDO
    !
12  FORMAT(1X,I6,2ES22.10)
    !
    !
    !
    tot_anom_density = 0.0_DP
    !
    file = 'anom_density'
    WRITE(filename, '(A, "_",A, ".dat")') TRIM(prog), TRIM(file)
    OPEN(UNIT = 101, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    WRITE(101,*) "# anomalous density "
    !
    DO index_I = 1, dim_2b
       eigvec_2b(:) =  avec_2b(:,index_I)
       CALL anom_density(eigvec_2b, dim_2b, two_body_bas, a_d_value_2)
       WRITE(101,*) index_I, aval_2b(index_I), a_d_value_2, a_d_value_2**2
       tot_anom_density = tot_anom_density + a_d_value_2**2
    ENDDO
    !
    CLOSE (UNIT = 101)
    !
    WRITE(*,*) 'TOTAL ANOM. DENSITY = ', tot_anom_density
    !
  END SUBROUTINE compute_anom_density

END MODULE module_2body
