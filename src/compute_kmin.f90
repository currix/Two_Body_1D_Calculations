FUNCTION Compute_kmin(Iprint)
  !
  USE nrtype
  USE constants
  USE vardef_ri
  USE pot_param
  USE module_1body_ho 
  USE module_egs_kmin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: Iprint
  REAL(KIND = DP) :: Compute_kmin
  !
  !
  ! LOCAL VARIABLES
  INTEGER :: dim_HO_temp, Ierr, Ifail, nev
  REAL(KIND = DP) :: fmult, k0, k1, k10, kmin, egs0, tol
  !
  !
  dim_HO_temp = dim_1b_HO
  dim_1b_HO = 1
  !
  fmult = 5.0_DP/1.5_DP
  !
  !
  ALLOCATE(Har_1b_Bas(1:dim_X, 1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Har_1b_Bas allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Aval_1b_Har(1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_Har allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_Har(1:1, 1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_Har allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_Har_X(1:dim_X, 1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_Har_X allocation request denied."
     STOP
  ENDIF
  !
  !
  DO
     fmult = fmult*1.5_DP
     ! LENGTH SCALE PARAMETER INITIAL EVALUATION
     k0 = 0.0_DP
     !     PARAM(1) -> Well depth  &  PARAM(2)  -> Well length unit 
     IF (Param_pot(2) /= 0.0_DP) THEN
        k1 = ABS(fmult*Param_pot(1)/(Param_pot(2)*Param_pot(2)))
     ELSE
        k1 = ABS(fmult*Param_pot(1))
     ENDIF
     k10 = k1
     !     
     !
     Ifail = 0
     tol = 0.0_DP
     nev = 50
     !
     CALL E04ABF(EGS, tol, tol, k0, k1, nev, kmin, egs0, Ifail)
     !
     IF (Iprint > 1) PRINT*, "kmin = ", kmin, " kmin - k10 = ", kmin-k10, & 
          " a = ", SQRT(sqrt(kmin/h_sq_over_m)), &
          " EGS = ", egs0
     !     
     IF (ABS(kmin-k10) > 1.0D-6) EXIT
     !
  ENDDO
  !
  dim_1b_HO = dim_HO_temp
  !
  DEALLOCATE(Har_1b_Bas, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Har_1b_Bas deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Aval_1b_Har, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_Har deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_1b_Har, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_Har deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_1b_Har_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_Har_X deallocation request denied."
     STOP
  ENDIF
  !
  IF (Iprint > 1) PRINT*, "MINIMUM:: kmin = ", kmin, & 
       " a = ", SQRT(sqrt(kmin/h_sq_over_m)), &
       " EGS = ", egs0
  !
  !
  Compute_kmin = kmin
  !
END FUNCTION Compute_kmin
