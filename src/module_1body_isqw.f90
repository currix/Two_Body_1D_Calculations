MODULE module_1body_isqw
  !
  USE nrtype
  USE constants
  USE vardef_ri
  USE pot_param
  !
  IMPLICIT NONE
  !
  ! One body variables
  !
  ! DIMENSION OF THE 1body BOX BASIS
  INTEGER(KIND = I4B) :: dim_1b_BOX
  !
  ! ISQW BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_1b_Box
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Box_1b_Bas, Avec_1b_Box, Avec_1b_Box_X
  ! EIGENVECTOR DERIVATIVES
  !REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_1b_Box_Der_X
  !
  !
  !
CONTAINS
  !
  !
  SUBROUTINE ISQW_1D_BASIS(XMpar, Iprint)
    !     
    !     COMPUTES THE 1D BOX BASIS
    !
    !     INPUT     :: XMPAR  --> PI/(2 X_M)
    !    
    !     VARIABLES :: DIM_X  --> XGRID DIMENSION
    !                  X_GRID --> VECTOR WITH X GRID
    !                  X_MAX  --> BOX RADIUS
    !                  NDIM   --> DIMENSION OF THE BOX BASIS
    !
    !     OUTPUT    :: BOX_1b_BAS --> MATRIX WITH BOX BASIS
    !
    !     FORMAT    :: IPRINT --> VERBOSITY CONTROL
    !
    !     by Currix TM.
    !
    !
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    INTEGER(KIND = I4B), INTENT(IN) :: Iprint
    REAL(KIND = DP), INTENT(IN) ::  XMpar
    !
    !
    ! OTHER VARIABLES
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: ISQW_norm_test
    REAL(KIND = DP) :: xmi, error, ainteg
    INTEGER(KIND = I4B) :: kx, ind, Ifail, Ierr
    !
    IF (Iprint > 1) PRINT*, "BUILDING BOX BASIS"
    !
    xmi = 1.0_dp/SQRT(X_max)
    !
    DO kx = 1, dim_X
       DO ind = 1, dim_1b_BOX, 2
          Box_1b_Bas(kx,ind) = xmi*COS(ind*XMpar*X_grid(kx))
       ENDDO
       DO ind = 2, dim_1b_BOX, 2
          Box_1b_Bas(kx,ind) = xmi*SIN(ind*XMpar*X_grid(kx))
       ENDDO
    ENDDO
    !
    !    TESTING NORMALIZATION 
    IF (Iprint > 1) THEN
       !
       ! DEFINE ISQW_norm_test
       ALLOCATE(ISQW_norm_test(1:dim_X), STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "ISQW_norm_test allocation request denied."
          STOP
       ENDIF
       !
       DO kx = 1, dim_1b_BOX
          ISQW_norm_test = Box_1b_Bas(:,kx)*Box_1b_Bas(:,kx)
          Ifail = 0
          !
          CALL D01GAF(X_grid, ISQW_norm_test, dim_X, ainteg, error, Ifail)
          PRINT*, "BOX FUNCTION ", kx, " NORMALIZATION", ainteg
       ENDDO
       !
       DEALLOCATE(ISQW_norm_test, STAT = Ierr)    
       IF (Ierr /= 0) THEN
          PRINT*, "ISQW_norm_test deallocation request denied."
          STOP
       ENDIF
       !
       WRITE(*,*) "DONE"
    ENDIF
    !     
    !     
  END SUBROUTINE ISQW_1D_BASIS
  !
  !
  SUBROUTINE BOXDIAG(XMpar, Iprint)
    !
    !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A BOXBASIS
    !
    !     INPUT      :: XMPAR    --> PI/(2 X_MAX)
    !
    !     FORMAT     :: IPRINT  --> VERBOSITY CONTROL
    !
    !     by Currix TM.
    !
    ! Lapack 95
    USE LA_PRECISION, ONLY: WP => DP
    USE F95_LAPACK, ONLY: LA_SYEVR
    !
    IMPLICIT NONE
    !
    !
    !
    !ARGUMENTS
    INTEGER(KIND = I4B), INTENT(IN) :: Iprint
    REAL(KIND = DP), INTENT(IN) :: XMpar
    !
    INTEGER(KIND = I4B) :: Ierr
    REAL(KIND = DP) :: XMpar2
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
    REAL(KIND = DP) :: ai, aham, error
    INTEGER(KIND = I4B) :: i, j, kx, Ifail
    !
    XMpar2 = XMpar*XMpar*h_sq_over_m/2.0_dp
    !
    IF (IPRINT>1) PRINT*, "BUILDING BOX HAMILTONIAN MATRIX"
    !
    !   
    ALLOCATE(Ham_Mat(1:dim_1b_BOX,1:dim_1b_BOX), STAT = Ierr)
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
    !    HAMILTONIAN MATRIX
    !
    DO i = 1, dim_1b_BOX
       ai = 1.0_dp*i
       !     KINETIC ENERGY 
       Ham_Mat(i,i) = Ham_Mat(i,i) + XMpar2*ai*ai
       !     
       DO j = i, dim_1b_BOX
          !     POTENTIAL ENERGY
          Pot_vec = Potf(X_Grid)*Box_1b_Bas(:,i)*Box_1b_Bas(:,j)
          !
          !     INTEGRATION
          Ifail = 0
          CALL D01GAF(X_grid, Pot_vec,dim_X,aham,error,Ifail)
          Ham_Mat(j,i) = Ham_Mat(j,i) + aham
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
    !
    !
    !     DIAGONALIZATION USING LAPACK 95
    CALL LA_SYEVR(A=Ham_Mat, W=Aval_1b_Box, JOBZ='V', UPLO='L')
    !     
    !
    !     EIGENVECTOR MATRIX
    Avec_1b_Box = Ham_Mat
    Avec_1b_Box_x = 0.0_dp
    !
    DO kx = 1, dim_X
       DO i = 1, dim_1b_BOX
          DO j = 1, dim_1b_BOX
             Avec_1b_Box_x(kx,i) = Avec_1b_Box_x(kx,i) + Ham_Mat(j,i)*Box_1b_Bas(kx,j)
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
  END SUBROUTINE BOXDIAG
  !
  ! subroutine output results in temp_output_isqw.f90
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
    INTEGER(KIND = I4B) :: i, kx, kkx
    CHARACTER(LEN=56) :: prog, file, filename
    prog = "isqw_1D_1b"
    !
    IF (Iprint > 0) THEN
       WRITE(*,*) "EIGENVALUES IN A ", dim_1b_BOX, " DIM BOX BASIS"
       DO i = 1, dim_1b_BOX
          WRITE(*,*) i, Aval_1b_Box(i)
       ENDDO
    ENDIF
    !
    !    SAVING BOX BASIS
    IF (I_save_1b_BAS==1) THEN
       file = 'basis'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       IF ( dim_1b_BOX < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       IF ( dim_1b_BOX > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(70,*) "# INFSQ  dim_1b_BOX = ", dim_1b_BOX, " Box radius = ", X_max, " fm"
       WRITE(70,*) "#Grid     Infinite squared well basis"
       DO kx = 1, dim_X
          WRITE(70,11) X_grid(kx), (Box_1b_bas(kx,kkx), kkx=1, dim_1b_BOX)
       ENDDO
       CLOSE(UNIT = 70)
    ENDIF
    !
    !     SAVING EIGENVECTORS
    IF (I_save_1b_WF==1) THEN
       file = 'eigenvectors'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       IF ( dim_1b_BOX < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       IF ( dim_1b_BOX > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(71,*) "# BOX  dim_1b_BOX = ", dim_1b_BOX, " Box radius = ", X_max, " fm"
       WRITE(71,*) "#Grid    Eigenvectors"
       DO kx = 1, dim_X
          WRITE(71,11) X_grid(kx), (Avec_1b_Box_X(kx,kkx), kkx=1, dim_1b_BOX)
       ENDDO
       CLOSE(UNIT = 71)
    ENDIF
    !
    !
    !     SAVING ENERGIES
    IF (I_save_1b_EN==1) THEN
       file = 'eigenvalues'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       IF ( dim_1b_BOX < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       IF ( dim_1b_BOX > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(73,*) "# BOX  dim_1b_BOX = ", dim_1b_BOX, " Box radius = ", X_max, " fm"
       WRITE(73,*) "# Eigenvalues"
       DO i = 1, dim_1b_BOX
          WRITE(73,10) i, Aval_1b_Box(i)
       ENDDO
       CLOSE(UNIT = 73)
       !
       file = 'pot_eigvec'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       IF ( dim_1b_BOX < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       IF ( dim_1b_BOX > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(74,*) "# BOX  dim_1b_BOX = ", dim_1b_BOX, " Box radius = ", X_max, " fm"
       WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
       !
       file = 'pot_eigvec2'
       WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       IF ( dim_1b_BOX < 10) THEN !to avoid spaces
          WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       IF ( dim_1b_BOX > 99) THEN 
          WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_1b_BOX
       ENDIF
       OPEN(UNIT = 75, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
       WRITE(75,*) "# BOX  dim_1b_BOX = ", dim_1b_BOX, " Box radius = ", X_max, " fm"
       WRITE(75,*) "#Grid    Potential    10*Eigenfunctions^2+eigenvalue"
       DO kx = 1, dim_X
          WRITE(74,11) X_grid(kx), Potf(X_grid(kx)), &
               (10.0_dp*Avec_1b_Box_X(kx,kkx)+Aval_1b_Box(kkx), kkx=1, dim_1b_BOX)
       ENDDO
       DO kx = 1, dim_X
          WRITE(75,11) X_grid(kx), Potf(X_grid(kx)), &
               (10.0_dp*Avec_1b_Box_X(kx,kkx)*Avec_1b_Box_X(kx,kkx)+Aval_1b_Box(kkx), &
               kkx=1, dim_1b_BOX)
       ENDDO
       CLOSE(UNIT = 74)
       CLOSE(UNIT = 75)
    ENDIF
    !
10  FORMAT (1X,I6,1X,E16.8)
11  FORMAT (1X,E14.6,1X,150E16.8)
    !
  END SUBROUTINE output_1b_results
  !
  !
  !
END MODULE module_1body_isqw
  
