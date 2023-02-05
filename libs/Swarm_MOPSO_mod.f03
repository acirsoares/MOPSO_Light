! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  Swarm = Swarm structure and function for MOPSO algorithms
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! *** Swarm_MOPSOL_mod - Class, (Version 1.0)
! ***************************************************************************************
! *** Module Application:
! * - A public structure to store a set of particles (Swarm) with X,Y,bestLocalPosition, 
! *    bestGlobalPosition, and particle velocity. 
! * - Contains a specific initialization procedure. 
! *
! *** Auxiliary libraries:
! * - Precision_defaults_MOPSO_mod.f03 ;  
! * - ParetoFront_EDSD_mod.f03 ;  
! * - ObjectiveFunction_MOPSOL_mod.f03 ;
! * - SobolSequence_mod.f03 .
! *
! *** Initialization:
! *   use(Swarm_MOPSOL_mod)
! *   type(Swarm_MOPSO_type),allocatable :: "NAME" (:) 
! *   SW = New_Swarm(OF,NPS,nx,ny,NIM,Random) 
! *       OF = Objective_Function
! *       NPS = Number of Particles in Swarm
! *       nX = Number of variables/Parameters
! *       nY = Number of functions
! *       NIM = Initialization method
! *       RANDOM = .True.(clock dependent)  .False.(clock independent)
! *    
! *** Variables and Parameters: 
! *- type(SPE_ParetoFront_type):: P
! *- type(SPE_ParetoFront_type):: B 
! *- allocatable :: VX (:)
! *- allocatable :: LX (:) 
! *
! *** External Subroutines (Public)
! *- function - New_Swarm
! *
! ***************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ** Module ParetoFront_mod: sets the abstract class ParetoFront
module Swarm_MOPSO_mod
  use Precision_defaults_MOPSO_mod ! real(kind=rp)
  use ParetoFront_mod  
  use ObjectiveFunction_MOPSO_mod
  use SobolSequence_mod
  implicit none

  private :: Local_RAND

  public :: Swarm_MOPSO_type, New_Swarm

  type :: Swarm_MOPSO_type
    type(SPE_ParetoFront_type):: P,Lb        ! Single Particle Element (P= Particle Position; Lb=Particle Best Position)
    real(kind=rp), allocatable :: Vx (:)     ! Particle Velocity (nX dimensions)
    real(kind=rp), allocatable :: Gb (:)     ! Global best (leader) position (nX dimensions)
  end type Swarm_MOPSO_type  

contains

  function New_Swarm(OF,PF,nSP,nX,nY,ncase,Random) result (SW)
    implicit none

    class(ObjectiveFunction_MOPSO_Class):: OF
    class(ParetoFront_class):: PF
    integer, intent (in) :: nSP,nX,nY,ncase
    logical,intent (in) :: Random
    type(Swarm_MOPSO_type), allocatable :: SW(:)
    type(SobolSequenceGenerator_CLASS) :: SSG
    integer i,j
    real , allocatable, dimension ( : ) :: SobolX
    real :: R1

    allocate  (SW(nSP))
    do i=1,nSP
      allocate  (SW(i)%VX(nX))
      allocate  (SW(i)%Gb(nX))
      allocate  (SW(i)%P%X(nX))
      allocate  (SW(i)%Lb%X(nX))
      allocate  (SW(i)%P%Y(nY))
      allocate  (SW(i)%Lb%Y(nY))
    end do
    select case (ncase)

    case (1)
      do i=1,nSP
        do j=1,nX
          R1=Local_RAND(Random)
          SW(i)%P%X(j)=OF%liminfX(j)+R1*(OF%limsupX(j)-OF%liminfX(j))
          SW(i)%Vx(j)=0.0d0   
        end do
        call OF%Function_X(SW(i)%P%X,SW(i)%P%Y)
        SW(i)%Lb=SW(i)%P
      end do

    case (2)
      allocate (SobolX(nX))
      SSG=SobolSequenceGenerator_CLASS(nX)
      do i=1,nSP
        SobolX=SSG%i4_sobol()
        do j=1,nX
          SW(i)%P%X(j)=OF%liminfX(j)+SobolX(j)*(OF%limsupX(j)-OF%liminfX(j))
          SW(i)%Vx(j)=0.0d0   
        end do
        call OF%Function_X(SW(i)%P%X,SW(i)%P%Y)
        SW(i)%Lb=SW(i)%P
      end do
      deallocate (SobolX)

    case (3)
      do i=1,nSP
        do j=1,nX
          R1=Local_RAND(Random)
          SW(i)%P%X(j)=OF%liminfX(j)+R1*(OF%limsupX(j)-OF%liminfX(j))
          R1=Local_RAND(Random)
          SW(i)%Vx(j)=OF%liminfX(j)+R1*(OF%limsupX(j)-OF%liminfX(j))-SW(i)%P%X(j)    ! Inicializa as velocidades das partículas do SWame
        end do
        call OF%Function_X(SW(i)%P%X,SW(i)%P%Y)
        SW(i)%Lb=SW(i)%P
      end do

    case (4)
      allocate (SobolX(nX))
      SSG=SobolSequenceGenerator_CLASS(nX)
      do i=1,nSP
        SobolX=SSG%i4_sobol()
        do j=1,nX
          SW(i)%P%X(j)=OF%liminfX(j)+SobolX(j)*(OF%limsupX(j)-OF%liminfX(j))
          R1=Local_RAND(Random)
          SW(i)%Vx(j)=OF%liminfX(j)+R1*(OF%limsupX(j)-OF%liminfX(j))-SW(i)%P%X(j)    ! Inicializa as velocidades das partículas do SWame
        end do
        call OF%Function_X(SW(i)%P%X,SW(i)%P%Y)
        SW(i)%Lb=SW(i)%P
      end do
      deallocate (SobolX)

    case DEFAULT
      write(*,*) "ERROR - Particle Swarm INITIALIZATION not SET "

    end select
    ! Try to insert the initial swarm particles into Pareto Front.
    do i=1,nSP
      call PF%TrySolution_PFS(SW(i)%P)
    end do

    return
  end function New_Swarm

  real function Local_RAND(RANDOM)
    implicit none

    logical,INTENT(IN) :: RANDOM

    if (RANDOM) then
      call RANDOM_NUMBER (Local_RAND)
      return
    end if  
    Local_RAND=RAND()

    return
  end function Local_RAND  

end module Swarm_MOPSO_mod
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
