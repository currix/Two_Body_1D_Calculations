MODULE module_egs_amin
  !
  ! $Id: e_gstate_tho.f90,v 1.1 2013/07/01 15:47:33 curro Exp $
  !
  USE nrtype
  USE constants
  USE pot_param
  USE vardef_ri
  USE module_1body_tho 
  !
  IMPLICIT NONE
  !
CONTAINS
  SUBROUTINE EGS_THO(KVAL,FVAL)
    !
    REAL(KIND = DP), INTENT(IN)  :: KVAL
    REAL(KIND = DP), INTENT(OUT) :: FVAL
    !
    REAL(KIND = DP) :: APAR
    !
    !
    APAR = SQRT(SQRT(KVAL/H_SQ_OVER_M))
    !
    CALL THO_1D_BASIS(APAR, 1, 1)
    !
    CALL THARDIAG(APAR, 1)
    !
    FVAL = aval_1b_thar(1)
    !
    RETURN
  END SUBROUTINE EGS_THO
  !
END MODULE module_egs_amin
