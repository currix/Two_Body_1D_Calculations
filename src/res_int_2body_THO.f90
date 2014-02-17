PROGRAM res_int_THO_basis
  !
  ! 
  !     PROGRAM TO DIAGONALIZE A MATTER DENSITY DEPENDENT RESIDUAL 
  !     INTERACTION FOR A 2BODY PROBLEM USING A CONTINUUM DISCRETIZATION 
  !     BASED ON A TRUNCATED 2 BODY TRANSFORMED HARMONIC BASIS
  !
  !     SINGLET CASE -> S = 0, Symmetrised spatial components
  !
  !     $Id: res_int_2body_THO.f90,v 1.1 2013/07/01 15:52:18 curro Exp $
  !
  !     by Currix TM
  !
  !
  USE nrtype
  USE constants
  USE pot_param
  USE vardef_ri
  USE module_1body_tho
  USE module_2body
  USE module_2body_tho
  !
  IMPLICIT NONE
  !
  ! force constant kmin
  REAL(KIND = DP) :: amin
  ! Inverse oscillator length a = (\nu K/\hbar^2)^(1/4)  (fm^{-1})
  REAL(KIND = DP) :: apar
  !
  ! Output control (1-body)
  INTEGER(KIND = I4B) :: Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS 
  ! Output control (2-body)
  INTEGER(KIND = I4B) :: i_comp_anom_den, i_comp_e1, i_comp_e2
  !
  ! AUXILIARY VARIABLES
  INTEGER(KIND = I4B) :: Index = 0, Ierr
  !
  REAL(KIND = DP) :: compute_amin_tho
  EXTERNAL compute_amin_tho
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
  NAMELIST/INP_X/     X_min, X_max, ratio
  NAMELIST/INP_DIM/   dim_X, dim_1b_THO
  NAMELIST/INP_MASS/  iad, reduced_mass
  NAMELIST/INP_POT/   Param_pot
  NAMELIST/INP_RI/    Param_ri
  NAMELIST/INP_2B/    n_bound_bas, E_threshold, dim_eig2body
  NAMELIST/INP_GRAPH/ i_pri, i_flag_graph, np3d
  NAMELIST/INP_AUX1/  Iprint, i_save_1b_EN, i_save_1b_WF, i_save_1b_BAS 
  NAMELIST/INP_AUX2/  i_comp_anom_den, i_comp_e1, i_comp_e2
  !
  ! PROGRAM VERSION
  IF (Iprint > 1) PRINT*, "$Id: res_int_2body_THO.f90,v 1.1 2013/07/01 15:52:18 curro Exp $"
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
  ! ONE - BODY    SECTION
  !
  ! (2) COMPUTE OSCILLATOR LENGTH
  !
  IF (IAD == 1) THEN
     !
     b = (SQRT(SQRT(h_sq_over_m)))
     apar = 1.0_DP/b
     !
  ELSE
     !
     ! Compute optimum apar value
     apar = compute_amin_tho(Iprint)
     b =  1.0_DP/apar
     !
  ENDIF
  !
  !
  ! (3) Build Basis ALLOCTE SCALING FUNCTIONS!!!!!%%%%%%%%% E AGGIUNGI LST_MOD!!
  !
  IF (Iprint > 1) PRINT*,  "TRANSFORMED HARMONIC BASIS CALCULATION apar = ", &
       apar, "fm^-1 bmin = ", 1/amin, " fm;  DIMENSION = ", dim_1b_THO     
  !
  !  Add one for the calculation of derivatives in the wfp subroutine
  ALLOCATE(THar_1b_Bas(1:dim_X, 1:dim_1b_THO + 1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "THar_1b_Bas allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(der_S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der_S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(der2_S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der2_S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(S_over_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_over_x allocation request denied."
     STOP
  ENDIF
  !
  CALL THO_1D_BASIS(apar, dim_1b_THO + 1, Iprint)
  !
  !
  ! (4) Build and diagonalize 1-body Hamiltonian
  !
  ALLOCATE(Aval_1b_THar(1:dim_1b_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_THar allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_THar(1:dim_1b_THO, 1:dim_1b_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_THar_X(1:dim_X, 1:dim_1b_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar_X allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_THar_Der_X(1:dim_X, 1:dim_1b_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar_Der_X allocation request denied."
     STOP
  ENDIF
  !
  CALL Thardiag(apar, Iprint)
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
  CALL bas_2_body(dim_1b_THO, E_threshold, Aval_1b_THar, dim_bound, dim_2body, two_body_bas, Iprint)
  !
  !
  ! (7) Compute core matter density
  ALLOCATE(rho_c(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "rho_c allocation request denied."
     STOP
  ENDIF
  !
  CALL core_matt_den(X_grid, Avec_1b_Thar_x, dim_bound, rho_c, Iprint)
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
  CALL diag_2body_ham(X_grid, rho_c, avec_1b_Thar_X, dim_bound, &
       two_body_bas, Aval_1b_Thar, dim_2body, aval_2b, avec_2b, &
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
  CALL two_body_results(X_grid, avec_1b_Thar_X, aval_1b_thar, aval_2b, avec_2b, &
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
     !
     IF (Iprint > 1) WRITE(UNIT=*,FMT=*) "Computing E2 transition intensity"
     !
     !     COMPUTE AUXILIARY MATRIX
     !           fxmat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
     ALLOCATE(f_x_mat(1:dim_2body,1:dim_2body), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x_mat allocation request denied."
        STOP
     ENDIF
     !
     print*, dim_1b_THO
     ALLOCATE(xthomatel(1:dim_1b_THO,1:dim_1b_THO), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "xthomatel allocation request denied."
        STOP
     ENDIF
     !
     CALL comp_x_tho_matel(dim_1b_THO, THar_1b_Bas, der_S_x)
     !
     CALL comp_tho_mat_e1(dim_bound, dim_1b_THO, dim_2body, avec_1b_Thar, two_body_bas, f_x_mat, xthomatel)
     !
     CALL comp_2body_E1(apar, aval_2b, avec_2b, dim_2body, dim_eig2body, f_x_mat)
     !
     DEALLOCATE(xthomatel, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "xthomatel deallocation request denied."
        STOP
     ENDIF
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
     !
     IF (Iprint > 1) WRITE(UNIT=*,FMT=*) "Computing E2 transition intensity"
     !     COMPUTE AUXILIARY MATRIX
     !           fxmat      ---> < (s) n'_1 n'_2 | x_1 + x_2 | (s) n_1 n_2 > matrix
     ALLOCATE(f_x2_mat(1:dim_2body,1:dim_2body), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "f_x2_mat allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(x2thomatel(1:dim_1b_THO,1:dim_1b_THO), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "x2thomatel allocation request denied."
        STOP
     ENDIF
     !
     CALL comp_x2_tho_matel(x2thomatel, dim_1b_THO, dim_X, X_Grid, THar_1b_Bas, der_S_x)
     !
     CALL comp_tho_mat_e2(dim_bound, dim_1b_THO, dim_2body, avec_1b_Thar, two_body_bas, f_x2_mat, x2thomatel)
     !
     CALL comp_2body_E2(apar, aval_2b, avec_2b, dim_2body, dim_eig2body, f_x2_mat)
     !
     DEALLOCATE(x2thomatel, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "x2thomatel deallocation request denied."
        STOP
     ENDIF
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
  DEALLOCATE(S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(der_S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der_S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(der2_S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der2_S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(S_over_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_over_x deallocation request denied."
     STOP
  ENDIF
  !
  !
  STOP 'Sayonara baby...'
  !
END PROGRAM res_int_THO_basis
