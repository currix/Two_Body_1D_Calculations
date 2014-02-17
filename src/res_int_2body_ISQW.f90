PROGRAM res_int_ISQW_basis
  !
  ! 
  !     PROGRAM TO DIAGONALIZE A MATTER DENSITY DEPENDENT RESIDUAL 
  !     INTERACTION FOR A 2BODY PROBLEM USING A CONTINUUM DISCRETIZATION 
  !     BASED ON A TRUNCATED 2 BODY INFINITE SQUARE WELL BASIS
  !
  !     SINGLET CASE -> S = 0, Symmetrised spatial components
  !
  !     $Id: res_int_2body_ISQW.f90,v 1.1 2013/07/01 15:52:33 curro Exp $
  !
  !     by Currix TM
  !
  !
  USE nrtype
  USE constants
  USE pot_param
  USE vardef_ri
  USE module_1body_isqw
  USE module_2body
  USE module_2body_isqw
  !
  IMPLICIT NONE
  !
  ! Output control (1-body)
  INTEGER(KIND = I4B) :: Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS 
  ! Output control (2-body)
  INTEGER(KIND = I4B) :: i_comp_anom_den, i_comp_e1, i_comp_e2
  !
  ! AUXILIARY VARIABLES
  INTEGER(KIND = I4B) :: Index = 0, Ierr
  REAL(KIND = DP) :: XMpar
  !
  !
  !    INTERFACE BLOCKS
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
  !
  ! DATA INPUT
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_X/     X_min, X_max
  NAMELIST/INP_DIM/   dim_X, dim_1b_BOX
  NAMELIST/INP_MASS/  iad, reduced_mass
  NAMELIST/INP_POT/   Param_pot
  NAMELIST/INP_RI/    Param_ri
  NAMELIST/INP_2B/    n_bound_bas, E_threshold, dim_eig2body
  NAMELIST/INP_GRAPH/ i_pri, i_flag_graph, np3d
  NAMELIST/INP_AUX1/  Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS 
  NAMELIST/INP_AUX2/  i_comp_anom_den, i_comp_e1, i_comp_e2
  !
  !
  !
  ! NAMELIST FILE
  ! READING INPUT
  READ(UNIT=*,NML=INP_X)
  !
  READ(UNIT=*,NML=INP_DIM)
  dim_X = dim_X + 2 ! TO ACCOMODATE X_min AND X_max
  !
  READ(UNIT=*,NML=INP_MASS)
  !
  READ(UNIT=*,NML=INP_POT)
  !
  READ(UNIT=*,NML=INP_RI)
  !
  READ(UNIT=*,NML=INP_2B)
  !
  READ(UNIT=*,NML=INP_GRAPH)
  !
  READ(UNIT=*,NML=INP_AUX1)
  !  
  READ(UNIT=*,NML=INP_AUX2)
  !  
  ! PROGRAM VERSION
  IF (Iprint > 1) PRINT*, "Version! $Id: res_int_2body_ISQW.f90,v 1.1 2013/07/01 15:52:33 curro Exp $"
  !
  ! DEFINE PROBLEM UNITS
  IF (iad == 1) THEN 
     h_sq_over_m = 1.0_DP/reduced_mass
  ELSE
     h_sq_over_m = H2OM/reduced_mass
  ENDIF
  !
  MASS_MEV = UM0*reduced_mass
  !
  !
  ! (1) DEFINE X GRID
  !
  ALLOCATE(X_grid(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid allocation request denied."
     STOP
  ENDIF
  !
  Delta_X = (X_max-X_min)/REAL(dim_X - 1, DP)  
  X_grid = X_min + Delta_X*REAL((/ (Index, Index = 0, dim_X-1) /),DP)
  !
  IF (Iprint > 1) PRINT*, ' X grid step = ', Delta_X, 'fm'
  IF (Iprint > 5) PRINT*, ' X grid = ', X_grid, 'fm'
  !
  !
  ! ONE - BODY    SECTION
  !
  !  (2)   LENGTH SCALE PARAMETER
  !
  XMpar = PI_D/(2.0_dp*X_max)
  !
  IF (Iprint>1) & 
       WRITE(*,*) "BOX BASIS CALCULATION XMpar = ", XMpar, &
       " DIMENSION ", dim_1b_BOX     
  !
  !
  ! (3) Build Basis
  !
  ALLOCATE(Box_1b_Bas(1:dim_X, 1:dim_1b_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Box_1b_Bas allocation request denied."
     STOP
  ENDIF
  !
  CALL ISQW_1D_BASIS(XMpar, Iprint)
  !
  !
  ! (4) Build and diagonalize 1-body Hamiltonian
  !
  ALLOCATE(Aval_1b_Box(1:dim_1b_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_Box allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_Box(1:dim_1b_BOX, 1:dim_1b_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_BOX allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_Box_X(1:dim_X, 1:dim_1b_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_Box_X allocation request denied."
     STOP
  ENDIF
  !
  !ALLOCATE(Avec_1b_Box_Der_X(1:dim_X, 1:dim_1b_BOX), STAT = Ierr)    
  !IF (Ierr /= 0) THEN
  !   PRINT*, "Avec_1b_Box_Der_X allocation request denied."
  !   STOP
  !ENDIF
  !
  CALL BOXDIAG(XMpar, Iprint)
  !
  ! (5) Output 1-body results 
  CALL output_1b_results(Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS)
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Two-body section
  !
  !
  ! (6) Build 2-body basis
  CALL bas_2_body(dim_1b_BOX, E_threshold, Aval_1b_Box, dim_bound, dim_2body, two_body_bas, Iprint)
  !
  !
  ! (7) Compute core matter density
  ALLOCATE(rho_c(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "rho_c allocation request denied."
     STOP
  ENDIF
  !
  CALL core_matt_den(X_grid, Avec_1b_Box_x, dim_bound, rho_c, Iprint)
  !
  !
  ! (8) Build and Diagonalize 2-body Hamiltonian
  ALLOCATE(aval_2b(1:dim_2body), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "aval_2b allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(avec_2b(1:dim_2body,1:dim_2body), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "avec_2b allocation request denied."
     STOP
  ENDIF
  !
  !
  CALL diag_2body_ham(X_grid, rho_c, avec_1b_Box_X, dim_bound, &
       two_body_bas, Aval_1b_Box, dim_2body, aval_2b, avec_2b, &
       i_pri, Iprint, i_flag_graph)
  !
  !
  ALLOCATE(avec_2b_X(1:dim_X,1:dim_X,1:dim_eig2body), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "avec_2b_X allocation request denied."
     STOP
  ENDIF
  !
  !
  ! (9) Data output
  CALL two_body_results(X_grid, avec_1b_Box_X, aval_1b_Box, aval_2b, avec_2b, &
       dim_bound, dim_2body, two_body_bas, dim_eig2body, avec_2b_X, &
       I_flag_graph, Np3d, Iprint)
  !
  !
  ! (10) Anomalous density
  IF (i_comp_anom_den == 1) &
       CALL compute_anom_density(X_grid, aval_2b, avec_2b, avec_2b_X, dim_X, dim_2body, dim_eig2body, two_body_bas)
  !
  ! (11) E1 Transition intensity
  IF (i_comp_e1 /= 0) THEN
     !     COMPUTE AUXILIARY MATRIX
     !           fxmat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
     ALLOCATE(f_x_mat(1:dim_2body,1:dim_2body), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x_mat allocation request denied."
        STOP
     ENDIF
     !
     CALL comp_isqw_mat_e1(dim_bound, dim_1b_BOX, dim_2body, avec_1b_BOX, two_body_bas, f_x_mat)
     !
     CALL comp_2body_E1(aval_2b, avec_2b, dim_2body, dim_eig2body, f_x_mat)
     !
     DEALLOCATE(f_x_mat, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x_mat deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  ! (12) E2 Transition intensity
  IF (i_comp_e2 /= 0) THEN
     !     COMPUTE AUXILIARY MATRIX
     !           fxmat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
     ALLOCATE(f_x2_mat(1:dim_2body,1:dim_2body), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x2_mat allocation request denied."
        STOP
     ENDIF
     !
     CALL comp_isqw_mat_e2(dim_bound, dim_1b_BOX, dim_2body, avec_1b_Box, two_body_bas, f_x2_mat)
     !
     CALL comp_2body_E2(aval_2b, avec_2b, dim_2body, dim_eig2body, f_x2_mat)
     !
     DEALLOCATE(f_x2_mat, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x2_mat deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  !
  DEALLOCATE(aval_2b, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "aval_2b deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(avec_2b, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "avec_2b deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(avec_2b_X, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "avec_2b_X deallocation request denied."
     STOP
  ENDIF
  !
  STOP 'Sayonara baby...'
  !
END PROGRAM res_int_ISQW_basis
