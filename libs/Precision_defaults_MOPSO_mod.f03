! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  MOPSOL = MOPSO_Light algorithm
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! *** MOPSO_precision_defaults_mod.f03 (Version 1.0)  
! ***************************************************************************************
! *** Module Application:
! *     Set the default precision to real kind in MOPSO lib-real-variables 
! *     
! *** Auxiliary libraries:
!
! *** Initialization:
! *     use MOPSO_Precision_defaults_mod
! *
! ***  Internal Functions / Subroutines (Privet)
! *    
! ***************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!****************************************************************************
module Precision_defaults_MOPSO_mod
  implicit none
! ** Module real type constant 
  integer(kind=2), parameter :: rp=8            
end module Precision_defaults_MOPSO_mod
