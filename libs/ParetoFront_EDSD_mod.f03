! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  ParetoFront_EDSD_mod - Class (general)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! ***ParetoFront_EDSD_mod - Class (general) (Version 1.0)
! ***************************************************************************************
! *** Object Description:
! *  A set of Pareto solution elements stored in an allocatable vector "PFS(:)"
! * contains specific managing procedures to get, compare, add, and eliminate
! * its elements. It uses the Euclidian Distance (ED) criteria to eliminate
! * the excess points and the S-Dominance criteria (SD) to assess dominance 
! * between two points.  
! *
! *** Initialization:
! * USE ParetoFront_EDSD_mod
! * TYPE(ParetoFront_EDSD_class) :: "NAME"
! * "NAME" = ParetoFront_EDSD_class(Nmaxsse,Nx,Ny)
! *
! *** Variables and Parameters: 
! *- PFS(:) (X,Y) = A set of non-dominated solutions  
! *- Nmaxsse = The maximum number of storable solution elements
! *- Nx = The Number of parameters 
! *- Ny = The Number of functions
! *- NCSP = The current number of stored solutions in Pareto Front
! *
! *** External Functions / Subroutines (Public)
! *- Constructor
! *- Destructor
! *- Function - NonDominatedByPFS (Y,s) -> Logical  : tests if an element Y is nondominated by the Pareto Front Elements 
! *- Function - PFScontains (Y) -> Logical   : tests if the Pareto Front Set contains the element Y
! *- Subroutine - Try_Solution_PFS (P)    :  If a point P(X,Y) is nondominated by the Pareto Front then the Pareto 
! *    Front receives it and may eliminate one or more points whether the new entrance dominates it.    
! *- OrdenatePFS (i)  :  orders its elements according to the objective function "i" values.
! *- Function get_X (k) -> Objct%P(k)%X
! *- get_Y (k) -> Objct%P(k)%Y
! *- get_P (k) -> Objct%P(k)
! *- get_P_Y (k,j) -> Objct%P(k)%Y(j)
! *- get_NCSP -> Objct%NCSP
! *
! *** Internal Functions / Subroutines (Privet)
! *- Subroutine - InsertPointIntoPFS_EDSD (P) 
! *- Subroutine - ExchangePositionPFS_EDSD(i,j)
! *    
! * OBS:  Pareto Front is set to store the lower values of the function F(X) = Y
! ***************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ** Module ParetoFront_EDSD_mod: sets the class ParetoFront_EDSD
module ParetoFront_EDSD_mod
  use Precision_defaults_MOPSO_mod ! real(kind=rp)
  use ParetoFront_mod ! Contains the abstract class ParetoFront

  implicit none

  private 

  public ParetoFront_EDSD_class, SPE_ParetoFront_type 

  ! ** STRUCTURES
  type FPE_ParetoFront_EDSD_type     !  
    type(SPE_ParetoFront_type):: P         !  Pareto Element (domain X e image Y)
    real(kind=rp) , allocatable ::  d (:) !  Euclidian distance
  end type FPE_ParetoFront_EDSD_type

  type Dist_ParetoFront_EDSD_type      
    integer  L,C   !Distance Element (Line L and column C from the Distance Matrix)               
  end type Dist_ParetoFront_EDSD_type 

  ! ** Class Definition
  type, public, extends(ParetoFront_class) :: ParetoFront_EDSD_class

    ! ** Variables 
    type(Dist_ParetoFront_EDSD_type),private  :: minD  !  Minimum distance element (minD) 
    type(FPE_ParetoFront_EDSD_type),private, allocatable :: PFS_EDSD(:)  !  Vetor PFS contains PFSsize elements (= Nmaxsse + 1)
    integer, private :: NCSP_EDSD   ! Current Number of Solutions in Pareto's Front (NCSP)  
    integer, private :: Nmaxsse_EDSD ! The maximum storable solution elements (Nmaxsse) 
    integer, private :: PFSsize_EDSD ! Size of the structure (PFSsize = Nmaxsse +1)
    integer, private :: Nx_EDSD  !  Number of parameters
    integer, private :: NY_EDSD  !  Number of functions
    real(kind=rp), private, allocatable :: S(:) 


  contains
    ! ** Class procedures
    procedure :: Sort_PFS => Sort_PFS_EDSD
    procedure :: NonDominatedByPFS => NonDominatedByPFS_EDSD 
    procedure :: PFScontains => PFScontains_EDSD
    procedure :: TrySolution_PFS => TrySolution_PFS_EDSD
    procedure :: get_X => get_X_EDSD
    procedure :: get_Y => get_Y_EDSD
    procedure :: get_P => get_P_EDSD
    procedure :: get_P_y => get_P_y_EDSD
    procedure :: get_NCSP => get_NCSP_EDSD
    procedure :: Set_S_dominanceCriteria
    procedure, private :: ExchangePositionPFS_EDSD 
    procedure, private :: InsertPointIntoPFS_EDSD
    final :: DESTRUCTOR_ParetoFront_EDSD_class 
  end type ParetoFront_EDSD_class  

  ! ** Module Interfaces 
  interface ParetoFront_EDSD_class
    module procedure CONSTRUCTOR_ParetoFront_EDSD_class
  end interface ParetoFront_EDSD_class  

  ! ** Module procedures
contains  

  !C  Constructor
  !C  ** Initializes Pareto Front (PF)
  function CONSTRUCTOR_ParetoFront_EDSD_class(Nmaxsse,Nx,Ny,S) result(PF)
    implicit none
    type(ParetoFront_EDSD_class) :: PF
    integer, intent(IN) :: Nmaxsse,Nx,Ny
    real(kind=rp), optional, dimension(:), intent(in) :: S 
    integer k
    PF%Nx_EDSD=Nx
    PF%Ny_EDSD=Ny
    PF%Nmaxsse_EDSD=Nmaxsse
    PF%PFSsize_EDSD=Nmaxsse+1
    allocate (PF%S(Ny))
    PF%S=0.0_rp
    if(present(S)) PF%S=S
    allocate  (PF%PFS_EDSD(PF%PFSsize_EDSD))
    do k=1,Nmaxsse
      allocate  (PF%PFS_EDSD(k)%P%X(Nx))
      allocate  (PF%PFS_EDSD(k)%P%Y(Ny))
      allocate  (PF%PFS_EDSD(k+1)%d(k))
    end do
    allocate  (PF%PFS_EDSD(PF%PFSsize_EDSD)%P%X(Nx))
    allocate  (PF%PFS_EDSD(PF%PFSsize_EDSD)%P%Y(Ny))
    PF%NCSP_EDSD=0

    return
  end function CONSTRUCTOR_ParetoFront_EDSD_class

  !C  Destructor
  !C  ** Destructs the Pareto Front (PF)
  subroutine DESTRUCTOR_ParetoFront_EDSD_class(PF)
    type(ParetoFront_EDSD_class) :: PF
    deallocate  (PF%PFS_EDSD)
    write(*,*) "DESTRUCTOR of ParetoFront_EDSD_class WORKS OK"
  end subroutine DESTRUCTOR_ParetoFront_EDSD_class

  !C  get_X
  !C  ** Returns vector X from element "k": PF%PFS_EDSD(k)%P%X
  function get_X_EDSD(PF,k) RESULT(g_X)
    class(ParetoFront_EDSD_class) :: PF
    integer, INTENT(IN) :: k
    real(kind=rp), allocatable ::  g_X(:)    
    allocate  (g_X(PF%Nx_EDSD))
    g_X = PF%PFS_EDSD(k)%P%X
    return
  end function get_X_EDSD

  !C  get_Y
  !C  ** Returns vector Y from element "k": PF%PFS_EDSD(i)%P%Y
  function get_Y_EDSD(PF,k) RESULT(g_Y)
    class(ParetoFront_EDSD_class) :: PF
    integer, INTENT(IN) :: k
    real(kind=rp), allocatable ::  g_Y(:)     
    allocate  (g_Y(PF%Ny_EDSD))
    g_Y = PF%PFS_EDSD(k)%P%Y
    return
  end function get_Y_EDSD

  !C  get_P
  !C  ** Returns the single Pareto element "k": PF%PFS_EDSD(k)%P
  function get_P_EDSD(PF,k) RESULT(g_P)
    class(ParetoFront_EDSD_class) :: PF
    integer, INTENT(IN) :: k
    type (SPE_ParetoFront_type) ::  g_P    
    allocate  (g_P%X(PF%Nx_EDSD))
    allocate  (g_P%Y(PF%Ny_EDSD))
    g_P = PF%PFS_EDSD(k)%P
    return
  end function get_P_EDSD

  !C  get_NCSP
  !C  ** Returns the current number of elements in Pareto front: PF%NCSP_EDSD
  function get_NCSP_EDSD(PF) 
    class(ParetoFront_EDSD_class) :: PF
    integer get_NCSP_EDSD
    get_NCSP_EDSD = PF%NCSP_EDSD
    return
  end function get_NCSP_EDSD

  !C  get_P_Y
  !C  ** Returns vector Y from element "k": PF%PFS_EDSD(i)%P%Y 
    function get_P_Y_EDSD(PF,k,j) 
    class(ParetoFront_EDSD_class) :: PF
    integer, INTENT(IN) :: k,j
    real(kind=rp) get_P_Y_EDSD
    get_P_Y_EDSD = PF%PFS_EDSD(k)%P%Y(j)
    return
  end function get_P_Y_EDSD

  !C  Set_S_dominanceCriteria value
  !C  ** Sets the S_dominance criteria
  subroutine Set_S_dominanceCriteria(PF,S)
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    real(kind=rp), dimension(:), intent(in) :: S 
    PF%S=S
    return
  end subroutine Set_S_dominanceCriteria

  !C  ExchangePositionPFS_EDSD
  !C  ** Exchange the position i with j of Pareto Front Elements
  subroutine ExchangePositionPFS_EDSD(PF,i,j)
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    integer, INTENT (IN) :: i,j
    integer k

    PF%PFS_EDSD(PF%NCSP_EDSD+1)%P=PF%PFS_EDSD(i)%P

    do k = 1,i-1     ! updates the "i" line from "d" matrix
      PF%PFS_EDSD(PF%NCSP_EDSD+1)%d(k)=PF%PFS_EDSD(i)%d(k)
    end do
    do k = i+1,PF%NCSP_EDSD   ! updates the "i" column from "d" matrix
      PF%PFS_EDSD(PF%NCSP_EDSD+1)%d(K)=PF%PFS_EDSD(k)%d(i)
    end do

    PF%PFS_EDSD(i)%P=PF%PFS_EDSD(j)%P

    do k = 1,i-1     ! updates the "i" line from "d" matrix
      PF%PFS_EDSD(i)%d(k)=PF%PFS_EDSD(j)%d(k)
    end do
    do k = i+1,j-1   ! updates the "i" column from "d" matrix
      PF%PFS_EDSD(k)%d(i)=PF%PFS_EDSD(j)%d(K)
    end do

    PF%PFS_EDSD(j)%P=PF%PFS_EDSD(PF%NCSP_EDSD+1)%P     ! replace the point
    do k = 1,j-1     ! updates the "i" line from "d" matrix
      PF%PFS_EDSD(j)%d(k)=PF%PFS_EDSD(PF%NCSP_EDSD+1)%d(k)
    end do
    do k = j+1,PF%NCSP_EDSD-1   ! updates the "i" column from "d" matrix
      PF%PFS_EDSD(k)%d(j)=PF%PFS_EDSD(PF%NCSP_EDSD+1)%d(K)
    end do
    return
  end subroutine ExchangePositionPFS_EDSD

  !C  OrdenatePFS
  !C  ** Sort the Pareto front elements in ascending order (SelectionSort) 
  !C  according to the function "j" value.
  subroutine Sort_PFS_EDSD(PF,j)
    implicit none
    class (ParetoFront_EDSD_class), intent (inout) :: PF
    integer, INTENT(IN) :: j
    integer i,k,lower
    do i = 1,PF%NCSP_EDSD-1
      lower=i
      do k=i+1,PF%NCSP_EDSD
        if (PF%PFS_EDSD(k)%P%Y(j).LT.PF%PFS_EDSD(lower)%P%Y(j)) lower=k
      end do
      if (lower.NE.i) CALL PF%ExchangePositionPFS_EDSD(i,lower)      
    end do
    return
  end subroutine Sort_PFS_EDSD

  !C  NonDominatedByPFS
  !C  ** Test if a particle Y is non_S-dominated by the Pareto Front (PF)
  logical function NonDominatedByPFS_EDSD (PF,Y)       
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    real(kind=rp), dimension(:), intent(in) :: Y  
    integer i,j
    NonDominatedByPFS_EDSD=.TRUE.                       
    do i= 1,PF%NCSP_EDSD
      do j = 1,PF%NY_EDSD
        IF (Y(j).LT.(PF%PFS_EDSD(i)%P%Y(j)-PF%S(j))) GOTO 40
      end do
      NonDominatedByPFS_EDSD=.FALSE.  
      return     
40    continue
    end do
    return
  end function NonDominatedByPFS_EDSD

  !C  PFScontains
  !C  ** Tests if PFS contains the particle Y
  logical function PFScontains_EDSD (PF,Y)       
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    real(kind=rp), DIMENSION(:), intent(in) :: Y 
    integer i,j                                    
    PFScontains_EDSD=.FALSE.                  
    do i= 1,PF%NCSP_EDSD
      do j = 1,PF%NY_EDSD
        IF (Y(j).NE.(PF%PFS_EDSD(i)%P%Y(j))) GOTO 40
      end do
      PFScontains_EDSD=.TRUE.   
      return              
40    continue
    end do
    return
  end function PFScontains_EDSD

  !C  InsertParticleIntoPFS
  !C  ** Insert the point "P" into Pareto Front with a "S" dominance criteria
  !C  ** If the number of Points exceeds Nmaxsse, it eliminates the   
  !C  ** particle with the lowest Euclidian Distance.            
  subroutine InsertPointIntoPFS_EDSD (PF,P)
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    type(SPE_ParetoFront_type), intent(in) :: P  ! the new point
    integer i,ia,ie,j,k,je
    real(kind=rp) smd !  

    !  Tests if the entering point dominates any stored point at 
    !  Pareto Front (PFS). If true, then it eliminates the old position  
    !  by overlapping it with the values of the PFS end position.
    ie=0
    do i = 1,PF%NCSP_EDSD
      ia=i-ie           ! ia corresponds to the current i value
      do j = 1,PF%NY_EDSD
        IF (P%Y(j)-PF%S(j) .GT.PF%PFS_EDSD(ia)%P%Y(j)) GOTO 50
      end do
      !  Remove the dominated point
      PF%PFS_EDSD(ia)%P=PF%PFS_EDSD(PF%NCSP_EDSD)%P     ! replace the point
      do k = 1,ia-1     ! updates the "i" line from "d" matrix
        PF%PFS_EDSD(ia)%d(k)=PF%PFS_EDSD(PF%NCSP_EDSD)%d(k)
      end do
      do k = ia+1,PF%NCSP_EDSD-1   ! updates the "i" column from "d" matrix
        PF%PFS_EDSD(k)%d(ia)=PF%PFS_EDSD(PF%NCSP_EDSD)%d(K)
      end do
      PF%NCSP_EDSD=PF%NCSP_EDSD-1
      ie=ie+1            ! ie = number of eliminated points
50    continue                                    
    end do
    !  Set a new position and insert the new point P
    PF%NCSP_EDSD=PF%NCSP_EDSD+1                
    PF%PFS_EDSD(PF%NCSP_EDSD)%P=P
    !  Calculates the Euclidian distances from the new point to
    !  each point in the Pareto Front
    do i = 1,PF%NCSP_EDSD-1
!          write(*,*)"Ok i=",i,PF%NCSP_EDSD        
      PF%PFS_EDSD(PF%NCSP_EDSD)%d(i)=0.0d0
!          write(*,*)"Ok i=",i  
      do j = 1,PF%NY_EDSD
        PF%PFS_EDSD(PF%NCSP_EDSD)%d(i)=PF%PFS_EDSD(PF%NCSP_EDSD)%d(i)+((PF%PFS_EDSD(PF%NCSP_EDSD)%P%Y(j)-PF%PFS_EDSD(i)%P%Y(j)))**2
      end do
      PF%PFS_EDSD(PF%NCSP_EDSD)%d(i)=DSQRT(PF%PFS_EDSD(PF%NCSP_EDSD)%d(i))
    end do
    ! Tests if the Pareto Front is full. If true it eliminates  
    ! the point with the small Euclidian Distance.
    if (PF%NCSP_EDSD.EQ.PF%PFSsize_EDSD) then
    ! Finds the minimum distance
    PF%minD%L=2
    PF%minD%C=1
    do i=3,PF%NCSP_EDSD
      do j=1,i-1
        if(PF%PFS_EDSD(i)%d(j).LT.PF%PFS_EDSD(PF%minD%L)%d(PF%minD%C))THEN
          PF%minD%L=i
          PF%minD%C=j
        end if
      end do
    end do
    !  Selects the particle with the second lowest distance
    je=PF%minD%L    ! je corresponds to the point to be eliminated
    if(PF%minD%C.EQ.1)THEN
      if(PF%minD%L.EQ.2)THEN
        smd=PF%PFS_EDSD(3)%d(1)      
        do i=4,PF%NCSP_EDSD          
          if (PF%PFS_EDSD(i)%d(1).LT.smd) smd=PF%PFS_EDSD(i)%d(1)
        end do
        do i=3,PF%NCSP_EDSD          
          if (PF%PFS_EDSD(i)%d(2).LT.smd) GOTO 60
        end do
      else
        smd=PF%PFS_EDSD(2)%d(1)      
        do i=3,PF%minD%L-1      
          if (PF%PFS_EDSD(i)%d(1).LT.smd) smd=PF%PFS_EDSD(i)%d(1)
        end do
        do i=PF%minD%L+1,PF%NCSP_EDSD   
          if (PF%PFS_EDSD(i)%d(1).LT.smd) smd=PF%PFS_EDSD(i)%d(1)
        end do
        do i=2,PF%minD%L-1      
          if (PF%PFS_EDSD(PF%minD%L)%d(i).LT.smd) GOTO 60
        end do
        do i=PF%minD%L+1,PF%NCSP_EDSD
          if (PF%PFS_EDSD(i)%d(PF%minD%L).LT.smd) GOTO 60
        end do
      end if
    else
      smd=PF%PFS_EDSD(PF%minD%C)%d(1)    
      do i=2,PF%minD%C-1         
        if (PF%PFS_EDSD(PF%minD%C)%d(i).LT.smd) smd=PF%PFS_EDSD(PF%minD%C)%d(i)
      end do
      do i=PF%minD%C+1,PF%NCSP_EDSD     
        if (PF%PFS_EDSD(i)%d(PF%minD%C).LT.smd) smd=PF%PFS_EDSD(i)%d(PF%minD%C)
      end do
      do i=1,PF%minD%L-1         
        if (PF%PFS_EDSD(PF%minD%L)%d(i).LT.smd) GOTO 60
      end do
      do i=PF%minD%L+1,PF%NCSP_EDSD     
        if (PF%PFS_EDSD(i)%d(PF%minD%L).LT.smd) GOTO 60
      end do
    end if
    je=PF%minD%C     
60  continue
    !Eliminates the selected point and updates the matrix "d"
    PF%PFS_EDSD(je)%P=PF%PFS_EDSD(PF%NCSP_EDSD)%P                     
    do k = 1,je-1                               
      PF%PFS_EDSD(je)%d(k)=PF%PFS_EDSD(PF%NCSP_EDSD)%d(k)
    end do
    do k = je+1,PF%NCSP_EDSD-1                        
      PF%PFS_EDSD(k)%d(je)=PF%PFS_EDSD(PF%NCSP_EDSD)%d(k)
    end do
      PF%NCSP_EDSD=PF%NCSP_EDSD-1
    end if
    return
  end subroutine InsertPointIntoPFS_EDSD

  !C  TrySolution_PFS
  !C  ** Tests if a point with solution Y is Non-Dominated By PFS,
  !C  and if .true., it tries to insert the point into Parteto Front.
  subroutine TrySolution_PFS_EDSD (PF,P)       
    implicit none
    class (ParetoFront_EDSD_class) :: PF
    type(SPE_ParetoFront_type), intent(in) :: P  
    if ( PF%NonDominatedByPFS(P%Y) ) THEN 
      call PF%InsertPointIntoPFS_EDSD(P)               
    end if
    return
  end subroutine TrySolution_PFS_EDSD

end module ParetoFront_EDSD_mod
