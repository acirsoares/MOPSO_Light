! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  Pareto Front - Class (general)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! *** ParetoFront_mod - Class, abstract (Version 1.0)
! ***************************************************************************************
! *** Object Description:
! * An objetc to store a set of Single Pareto Elements (SPE) 
! * Contains specific managing procedures to get, compare, add, 
! * and eliminate its elements.
! *
! *** Variables and Parameters: 
! *- Nmaxsse = The maximum number of storable solution elements
! *- nX = the Number of parameters 
! *- nY = the Number of objective functions
! *- X = an array of parameters 
! *- Y = an array of parameters 
! *- P = an SPE_ParetoFront_type  ! see module ParetoFront_Structures_mod
! *- NCSP = The current number of stored solutions in Pareto Front
! *
! *** External Functions / Subroutines (Public)
! *- Function - NonDominatedByPFS (Y) -> Logical  : tests if an element Y is nondominated by the Pareto Front Elements 
! *- Function - PFScontains (Y) -> Logical   : tests if the Pareto Front Set contains the element Y
! *- Subroutine - Try_Solution_PFS (P)    :  If a point P(X,Y) is nondominated by the Pareto Front then the Pareto 
! *    Front receives it and may eliminate one or more points whether the new entrance dominates it.    
! *- OrdenatePFS (i)  :  orders its elements according to the objective function "i" values.
! *- Function get_X (k) -> Objct%P(k)%X 
! *- Function get_Y (k) -> Objct%P(k)%Y
! *- Function get_P (k) -> Objct%P(k)
! *- Function get_P_Y (k,j) -> Objct%P(k)%Y(j)  
! *- Function get_NCSP -> Objct%NCSP  
! *
! * OBS:  Pareto Front is expected to hold the lower values of the function F(X) = Y
! ***************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ** Module Structures: sets the basic type for Paretro front elements
! ** TYPE Single Pareto Element (SPE) 
module ParetoFront_Structures_mod
  use Precision_defaults_MOPSO_mod

  implicit none

  private 

  public :: SPE_ParetoFront_type

  type SPE_ParetoFront_type      
    real(kind=rp), allocatable ::  X(:)    ! Parameters (X) 
    real(kind=rp), allocatable ::  Y(:)    ! Function Values (Y)
  end type SPE_ParetoFront_type

end module ParetoFront_Structures_mod

! ** Module ParetoFront_mod: sets the abstract class ParetoFront
module ParetoFront_mod
  use Precision_defaults_MOPSO_mod ! real(kind=rp)
  use ParetoFront_Structures_mod  ! SPE_ParetoFront_type

  implicit none

  private 

  public :: ParetoFront_class, SPE_ParetoFront_type

  ! ** Class Definition
  type, abstract :: ParetoFront_class

  contains
    ! ** Procedures 
    procedure (Sort_PFS_abstract) , deferred :: Sort_PFS
    procedure (NonDominatedByPFS_abstract) , deferred :: NonDominatedByPFS
    procedure (PFScontains_abstract) , deferred :: PFScontains
    procedure (TrySolution_PFS_abstract) , deferred :: TrySolution_PFS
    procedure (get_X_abstract) , deferred :: get_X 
    procedure (get_Y_abstract) , deferred :: get_Y 
    procedure (get_P_abstract) , deferred :: get_P 
    procedure (get_P_y_abstract) , deferred :: get_P_y 
    procedure (get_NCSP_abstract) , deferred :: get_NCSP
  end type ParetoFront_class  

  ! ** Interfaces 
  interface

    function get_X_abstract(PF,k) result(g_X)
      use Precision_defaults_MOPSO_mod
      import :: ParetoFront_class
      implicit none 
      class(ParetoFront_class) :: PF
      integer, intent(in) :: k
      real(kind=rp), allocatable ::  g_X(:)    
    end function get_X_abstract

    function get_Y_abstract(PF,k) result(g_Y)
      use Precision_defaults_MOPSO_mod
      import :: ParetoFront_class
      implicit none 
      class(ParetoFront_class) :: PF
      integer, intent(in) :: k
      real(kind=rp), allocatable ::  g_Y(:)     
      end function get_Y_abstract

    function get_P_abstract(PF,k) result(g_P)
      use ParetoFront_Structures_mod
      import :: ParetoFront_class
      implicit none 
      class(ParetoFront_class) :: PF
      integer, intent(in) :: k
      type (SPE_ParetoFront_type) ::  g_P    
    end function get_P_abstract

    function get_NCSP_abstract(PF) 
      import :: ParetoFront_class
      implicit none 
      class(ParetoFront_class) :: PF
      integer :: get_NCSP_abstract
    end function get_NCSP_abstract

    function get_P_Y_abstract(PF,k,j) 
      use Precision_defaults_MOPSO_mod
      import :: ParetoFront_class
      implicit none 
      class(ParetoFront_class) :: PF
      integer, intent(in) :: k,j
      real(kind=rp) :: get_P_Y_abstract
    end function get_P_Y_abstract

    subroutine Sort_PFS_abstract(PF,j)
      import :: ParetoFront_class
      implicit none 
      class (ParetoFront_class),intent (inout) :: PF
      integer, intent(in) :: j
    end subroutine Sort_PFS_abstract

    logical function NonDominatedByPFS_abstract (PF,Y)       
      use Precision_defaults_MOPSO_mod
      import :: ParetoFront_class
      implicit none 
      class (ParetoFront_class) :: PF
      real(kind=rp), dimension(:), intent(in) :: Y  
!      real(kind=rp),optional, dimension(:), intent(in) :: S 
    end function NonDominatedByPFS_abstract

    logical function PFScontains_abstract (PF,Y)       
      use Precision_defaults_MOPSO_mod
      import :: ParetoFront_class
      implicit none 
      class (ParetoFront_class) :: PF
      real(kind=rp), dimension(:), intent(in) :: Y 
    end function PFScontains_abstract

    subroutine TrySolution_PFS_abstract (PF,P)       
      use Precision_defaults_MOPSO_mod
      use ParetoFront_Structures_mod
      import :: ParetoFront_class
      implicit none 
      class (ParetoFront_class) :: PF
      type(SPE_ParetoFront_type), intent(in) :: P  
    end subroutine TrySolution_PFS_abstract

  end interface  

contains

end module ParetoFront_mod
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
