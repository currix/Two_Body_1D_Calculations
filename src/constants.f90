MODULE constants
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! PHYSICAL CONSTANTS
  REAL(KIND = DP), PARAMETER :: HC = 197.32858E0_DP  ! MeV fm
  REAL(KIND = DP), PARAMETER :: UM0 = 938.92635E0_DP ! MeV/c^2 / amu
  REAL(KIND = DP), PARAMETER :: H2OM = 41.4713768_DP ! MeV/fm^2 (hbar^2/1 amu)
  !
  ! REDUCED MASS (amu), H_SQ_OVER_M = h^2/M MeV/fm^2, MASS_MEV = MASS (MeV/c^2)
  REAL(KIND = DP) :: reduced_mass, h_sq_over_m, mass_mev
  !
  ! FLAG FOR NATURAL UNITS (H2OM = 1)
  INTEGER(KIND = I4B) ::  iad  
END MODULE constants
