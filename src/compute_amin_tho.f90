FUNCTION Compute_amin_tho(Iprint)
  !
  USE nrtype
  USE constants
  USE vardef_ri
  USE pot_param
  USE module_1body_tho 
  USE module_egs_amin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: Iprint
  REAL(KIND = DP) :: Compute_amin_tho
  !
  !
  ! LOCAL VARIABLES
  INTEGER :: dim_THO_temp, Ierr, Ifail, nev
  REAL(KIND = DP) :: fmult, a0, a1, a10, amin, egs0, tol
  !
  !
  !     OSCILLATOR LENGTH b = (\nu K/\hbar^2)^(-1/4)  (fm)
  dim_THO_temp = dim_1b_THO
  dim_1b_THO = 1
  !
  !
  ! Basis dimension = dim + 1 = 2 for Hamiltonian diagonalization
  ALLOCATE(THar_1b_Bas(1:dim_X, 1:2), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "THO_Bas allocation request denied."
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
  ALLOCATE(Aval_1b_THar(1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_THar allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_THar(1:1, 1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_1b_THar_X(1:dim_X, 1:1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar_X allocation request denied."
     STOP
  ENDIF
  !
  !
  DO
     fmult = 1.0_DP
     ! LENGTH SCALE PARAMETER INITIAL EVALUATION
     a0 = 0.0_DP
     !
     ! DEFINE THO CONSTANTS
     IF (Param_pot(2) /= 0.0_DP) THEN
        a1 = 1.0_DP/SQRT(SQRT(ABS(fmult*Param_pot(1)/(Param_pot(2)*Param_pot(2)))/h_sq_over_m))
     ELSE
        a1 = 1.0_DP/SQRT(SQRT(ABS(fmult*Param_pot(1))/h_sq_over_m))
     ENDIF
     a1 = fmult*0.5_DP
     a10 = a1
     !     
     !
     !
     Ifail = 0
     tol = 0.0_DP ! se la fissi a zero come defoult è la sqrt(epsilon)
     nev = 50
     !
     IF (Iprint > 2) print*, "a1 = ", a1
     CALL EGS_THO(a1,egs0)
     IF (Iprint > 2) print*, "energy = ", egs0
     a1 = 2.0_DP*a1
     IF (Iprint > 2) print*, "a1 = ", a1
     CALL EGS_THO(a1,egs0)
     IF (Iprint > 2) print*, "energy = ", egs0
     !
     CALL E04ABF(EGS_THO, tol, tol, a0, a1, nev, amin, egs0, Ifail)
     !
     IF (Iprint > 2) PRINT*, "amin = ", amin, " amin - a10 = ", amin-a10, & 
          " bmin = ", 1/amin," fm; EGS = ", egs0
     !     
     IF (ABS(amin-a10) > 1.0D-6) EXIT
     EXIT
     !
  ENDDO
  !
  dim_1b_THO = dim_THO_temp
  !
  DEALLOCATE(THar_1b_Bas, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "THar_1b_Bas deallocation request denied."
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
  DEALLOCATE(Aval_1b_THar, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_1b_THar deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_1b_THar, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_1b_THar_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_1b_THar_X deallocation request denied."
     STOP
  ENDIF
  !
  Compute_amin_tho = amin
  !
END FUNCTION Compute_amin_tho

