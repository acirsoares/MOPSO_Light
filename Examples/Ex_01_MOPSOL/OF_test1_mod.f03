!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C*****************************************************************************************
!C             =============               ====================
!C             "  Module   "               " OFtestFortran_mod "
!C             =============               ====================
!C
!C    It calculates the 
!C*****************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module OF_test1_mod
  use Precision_defaults_MOPSO_mod ! real(kind=rp)
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
!C************************************************************************
!C             ==============    
!C             Constructor
!C             ==============              
!C************************************************************************
!C  ** Initializes Objective Function (ObjFunc)
  function CONSTRUCTOR_OF_test1() result (ObjFunc)
    type(OF_test1_Class) :: ObjFunc
    integer :: i
    ObjFunc%nX=2
    ObjFunc%nY=2
  !***********  Upper and lower limits  to domain variables (MOPSO use these limits to set X(i)
    allocate  (ObjFunc%limsupX(ObjFunc%nX),ObjFunc%liminfX(ObjFunc%nX))
    ObjFunc%liminfX(1)=0.1
    ObjFunc%limsupX(1)=1.0
    ObjFunc%liminfX(2)=0.1
    ObjFunc%limsupX(2)=1.0
 return
  end function CONSTRUCTOR_OF_test1  
 
!C************************************************************************
!C             ===========    
!C             Destructor
!C             ===========              
!C************************************************************************
!C  ** Destruct the Pareto Front (VPS)
  subroutine DESTRUCTOR_OF_test1(ObjFunc)
    type(OF_test1_Class) :: ObjFunc
!   Depois de pronto, excluir a linha a seguir.
    write(*,*)"DESTRUCTOR of OF_test1_Class WORKS OK"
  end subroutine DESTRUCTOR_OF_test1

!  Objective Function (Function test 4 (eq. 5 in Proceeding Series of the Brazilian Society of Computational and Applied Mathematics, v. 7, n. 1, 2020.))
!  **  https://proceedings.sbmac.org.br/sbmac/article/view/2947
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

! get_nX_test1
! ** Returns the number of parameter (independent variables): nX
  integer function get_nX_Test1(ObjFunc) 
    class(OF_test1_Class)::ObjFunc
    get_nX_Test1 = ObjFunc%nX
    return
  end function get_nX_Test1

! get_nY_test1
! ** Returns the number of functions (dependent variables): nY
  integer function get_nY_Test1(ObjFunc) 
    class(OF_test1_Class)::ObjFunc
    get_nY_Test1 = ObjFunc%nY
    return
  end function get_nY_Test1

end module OF_test1_mod
