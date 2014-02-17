MODULE module_egs_kmin
  !
  ! $Id: e_gstate_ho.f90,v 1.1 2013/06/07 07:23:37 curro Exp $
  !
  USE nrtype
  USE constants
  USE pot_param
  USE vardef_ri
  USE module_1body_ho 
  !
  IMPLICIT NONE
  !
CONTAINS
  SUBROUTINE EGS(KVAL,FVAL)
    !
    REAL(KIND = DP), INTENT(IN)  :: KVAL
    REAL(KIND = DP), INTENT(OUT) :: FVAL
    !
    REAL(KIND = DP) :: APAR
    !
    !
    APAR = SQRT(SQRT(KVAL/H_SQ_OVER_M))
    !
    CALL HO_1D_BASIS(APAR, 1, 1)
    !
    CALL HARDIAG(APAR, 1)
    !
    FVAL = aval_1b_har(1)
    !
    RETURN
  END SUBROUTINE EGS
  !
END MODULE module_egs_kmin
