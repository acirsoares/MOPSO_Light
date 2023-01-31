! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  ObjectiveFunction_MOPSOL - Class (abstract)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************************
!*** ObjectiveFunction_MOPSOL_mod - Class(abstract) (Version 1.0)
!*****************************************************************************************
! *** Object Description:
! *     A parent abstract object to set an objective function R(n) in R(n).
! *     It recieves a vector with Real variable values (X) with a set of limit 
! *     values (min,max) and returs a vector with the Real functions values (Y)
! *
! *** Initialization:
! *     use ObjectiveFunction_MOPSOL_mod              ! Parent Type mode
! *     implicit none
! *     type, public, extends(ObjectiveFunction_MOPSO_Class)::::"OF_NAME"    ! Child Type
! *     contains
! *         procedure ::Function_X
! *     end type OF_test1_Class
! *** Variables and Parameters: 
! *    -     limsupX(:) , liminfX(:) ! Domain limits to Objective Function
! *    integer nX ! Number of variables
! *    integer nY ! Number of functions
! *
! *** External Functions / Subroutines (Public)
! *    - subroutine Function_X(Obj,X,Y)           
! *
! ***************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module ObjectiveFunction_MOPSO_mod
  use Precision_defaults_MOPSO_mod ! real(kind=rp)
  implicit none
! ** Module Objects / Variables / Procedures 
  private 
  public ObjectiveFunction_MOPSO_Class
! ** Object_Class definition
  type, abstract :: ObjectiveFunction_MOPSO_Class
    real(kind=rp), allocatable :: limsupX(:) , liminfX(:) ! Domain limits to Objective Function
  contains
    procedure (Function_X_implemented),deferred :: Function_X
    procedure (get_nX_implemented),deferred :: get_nX
    procedure (get_nY_implemented),deferred :: get_nY
  end type ObjectiveFunction_MOPSO_Class

  interface
    subroutine Function_X_implemented(ObjFunc,X,Y) 
      use Precision_defaults_MOPSO_mod
      import :: ObjectiveFunction_MOPSO_Class
      implicit none   
      class(ObjectiveFunction_MOPSO_Class) :: ObjFunc
      real(kind=rp), intent(in), dimension (:) :: X   
      real(kind=rp), intent(out), dimension (:) :: Y  
    end subroutine Function_X_implemented

    integer function get_nX_implemented(ObjFunc)
      import :: ObjectiveFunction_MOPSO_Class
      implicit none   
      class(ObjectiveFunction_MOPSO_Class) :: ObjFunc
    end function get_nX_implemented     

    integer function get_nY_implemented(ObjFunc)
      import :: ObjectiveFunction_MOPSO_Class
      implicit none   
      class(ObjectiveFunction_MOPSO_Class) :: ObjFunc
    end function get_nY_implemented     
  end interface
contains

end module ObjectiveFunction_MOPSO_mod
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
