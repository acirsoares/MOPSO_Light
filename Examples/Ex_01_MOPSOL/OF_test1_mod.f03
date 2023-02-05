! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  OF_test1_Class - An implementation example of the object :
!         ObjectiveFunction_MOPSOL_Class (abstract) for MOPSO algorithms
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************************
!*** OF_test1_mod - (Version 1.0)
!*****************************************************************************************
! *** Object Description:
! *     A chield object that implements an objective function test R(2) in R(2).
! *     It receives a vector with Real variable values [x1,x2] with a set of limits 
! *     values [(min1, max1),(min2,max2)] and returns a vector with the Real 
! *     functions values [y1,y2].
! *   Objective Function:
! *       Y(1)=X(1)
! *       g=2.0_rp-DEXP(-((X(2)-0.2_rp)/0.004_rp)**2)-0.8_rp*DEXP(-((X(2)-0.6_rp)/0.4_rp)**2)
! *       Y(2)=g/X(1)
! *   Domain limits:
! *       liminfX(1)=0.1
! *       limsupX(1)=1.0
! *       liminfX(2)=0.1
! *       limsupX(2)=1.0
! *   Pareto Front :
! *       local minimum x(2)=0.6  x(1)= [0.1,1.0] 
! *       global minimum x(2)=0.2 x(1)= [0.1,1.0]
! *
! *** Initialization:
! *     use OF_test1_mod              ! Chield Type mode
! *     implicit none
! *     type(OF_test1_Class) :: "OF_NAME" 
! *     "OF_NAME" = OF_test1_Class()
! *
! ***************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module OF_test1_mod
  use Precision_defaults_MOPSO_mod   ! real(kind=rp)
  use ObjectiveFunction_MOPSO_mod    ! Abstract Parent Class mode
  implicit none

  private

  public OF_test1_Class

! ** Module Objects / Variables / Procedures 
  type, public, extends(ObjectiveFunction_MOPSO_Class)::OF_test1_Class  ! Child Type
    integer,private :: nFunc      ! Function selecting number
    integer,private :: nX
    integer,private :: nY
  contains
    procedure ::Function_X=>ObjectiveFunction_test1
    procedure ::get_nX=>get_nX_test1
    procedure ::get_nY=>get_nY_test1
    final::DESTRUCTOR_OF_test1
  end type OF_test1_Class

  interface OF_test1_Class
    module procedure CONSTRUCTOR_OF_test1
  end interface OF_test1_Class  

contains
  ! ** Constructor
  !  Initializes Objective Function (ObjFunc)
  function CONSTRUCTOR_OF_test1() result (ObjFunc)
    implicit none

    type(OF_test1_Class) :: ObjFunc
    integer :: i

    ! Number of parameters (nx) and objective functions (ny)
    ObjFunc%nX=2
    ObjFunc%nY=2
    ! Upper and lower limits to domain variables 
    allocate  (ObjFunc%limsupX(ObjFunc%nX),ObjFunc%liminfX(ObjFunc%nX))
    ObjFunc%liminfX(1)=0.1
    ObjFunc%limsupX(1)=1.0
    ObjFunc%liminfX(2)=0.1
    ObjFunc%limsupX(2)=1.0

    return
  end function CONSTRUCTOR_OF_test1  

  ! ** Destructor
  !  Destruct the Pareto Front (VPS)
  subroutine DESTRUCTOR_OF_test1(ObjFunc)
    type(OF_test1_Class) :: ObjFunc

    return
  end subroutine DESTRUCTOR_OF_test1

  ! ** Objective Function 
  !   (Function test 4 (eq. 5 in Proceeding Series of the Brazilian Society of Computational and Applied Mathematics, v. 7, n. 1, 2020.))
  !   https://proceedings.sbmac.org.br/sbmac/article/view/2947
  subroutine ObjectiveFunction_test1(ObjFunc,X,Y) 
    implicit none

    class(OF_test1_Class)::ObjFunc
    real(kind=rp), intent(in), dimension (:) :: X
    real(kind=rp), intent(out), dimension (:) :: Y
    real(kind=rp) ::  g
 
    Y(1)=X(1)
    g=2.0_rp-DEXP(-((X(2)-0.2_rp)/0.004_rp)**2)-0.8_rp*DEXP(-((X(2)-0.6_rp)/0.4_rp)**2)
    Y(2)=g/X(1)

    return
  end subroutine ObjectiveFunction_test1

  ! ** get_nX_test1
  !  Returns the number of parameter (independent variables): nX
  integer function get_nX_Test1(ObjFunc) 
    class(OF_test1_Class)::ObjFunc

    get_nX_Test1 = ObjFunc%nX

    return
  end function get_nX_Test1

  ! ** get_nY_test1
  !  Returns the number of functions (dependent variables): nY
  integer function get_nY_Test1(ObjFunc) 
    class(OF_test1_Class)::ObjFunc

    get_nY_Test1 = ObjFunc%nY

    return
  end function get_nY_Test1

end module OF_test1_mod
