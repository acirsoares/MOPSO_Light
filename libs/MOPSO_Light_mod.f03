! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  MOPSO_Light: multi-objective optimization algorithm based on PSO (Particle Swarm Optimization)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! *** MOPSOLight_mod.f03 (Version 1.1)  
! ***************************************************************************************
! *** Algorithm Application:
! *     MOPSO_Light returns a Pareto Multidimensional Front applying the PSO 
! *      algorithm adapted for multi-objective optimization problems.
! *** Main MOPSO_Light features                                                       
! *      - Uses ParetoFront_EDSD_mod based on Euclidean distance to store and manage
! *         the Pareto Front                          
! *      - Global Leadership fixed for N=3 to 8 iterations
! *      - Randomly chosen global leadership for each particle
! *      - Dominance criteria to select local leadership                                                    
! *      - Swarm initialization with Random or Sobol distribution position 
! *      - Swarm initialization with Random or null initial velocity 
! *      OBS:  Pareto Front is set to deal with minimization problems.
! *
! *** Auxiliary libraries:
! *     ParetoFront_EDSD_mod.f03 ;  
! *     ObjectiveFunction_MOPSOL_mod.f03;
! *     Swarm_MOPSOL_mod.f03;
! *
! *** Initialization:
! *     
! *     "Parameters"%NPFS = ...              ! Number of storage positions for particles from Pareto's Front
! *     "Parameters"%nSP  = ...              ! Number of particles at the swarm
! *     "Parameters"%NITMOPSO = ...          ! Method Iterations
! *     "Parameters"%NCompleteRuns = ...     ! Number of MOPSO_L iterations
! *     "Parameters"%Random = ...            ! Logical .True. or .False. to set random function
! *     "Parameters"NIM = ...                ! Code number to initialize swarm particles POSITION (P) and VELOCITY (V)  (1) P=RANDOM V=0; (2) P=SOBOL(max of NX=6 dimensions) V=0;(3) P=RANDOM V=RANDOM;(4) P=SOBOL V=RANDOM.
! *     "Parameters"%NFI1_4 = ...            ! The number of fixed leadership iterations in the first quarter
! *     "Parameters"%NFI2_4 = ...            ! The number of fixed leadership iterations in the second quarter 
! *     "Parameters"%NFI2_2 = ...            ! The number of fixed leadership iterations in the second half
! *
! *     "OBJ_FUNC" = OFtestFortran_type()
! *
! ***  External Structures / Functions / Subroutines (Public)
! *    - MOPSO_Light_Parameters_type 
! *    - MOPSO_Light
! *    
! ***  Internal Structures / Functions / Subroutines (Private)
! *    - Local_RAND
! *    
! ***********************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!____________________________________________________________________________________
!|REFERENCES                                                                        |
!|REF#[1] Proceeding Series of the Brazilian Society of Computational and Applied   |
!|          Mathematics, v. 7, n. 1, 2020.                                          |
!|REF#[2] R. Eberhart and J. Kennedy. A new optimizer using particle swarm theory.  |
!|          In Proceedings of the Sixth International Symposium on Micro Machine    |
!|          and Human Science, 39-43, 1995. DOI: 10.1109/MHS.1995.494215.           |
!|REF#[3] A. Chatterjee and P. Siarry. Nonlinear inertia weight variation for       |
!|          dynamic adaptation in particle swarm optimization. Computers &          |
!|          operations research, 33:859-871, 2006. DOI:10.1016/j.cor.2004.08.012    |
!|REF#[4] Original FORTRAN77 version by Bennett Fox.                                |
!|        FORTRAN90 version by John Burkardt.                                       |
!|        I4_SOBOL generates a new quasirandom Sobol vectors                        |
!|        at https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html            |          |
!|__________________________________________________________________________________|

! ********** MOPSOLight
module MOPSO_Light_mod
  use Precision_defaults_MOPSO_mod
  use ParetoFront_EDSD_mod
  use ObjectiveFunction_MOPSO_mod
  use Swarm_MOPSO_mod
  implicit none

  private Local_RAND

  public MOPSO_Light_Parameters_type, MOPSO_Light, Set_MOPSO_Light_Parameters

  ! ** STRUCTURES 
  type,public :: MOPSO_Light_Parameters_type
    integer :: NPFS                       ! The maximum number of storable solution elements
    integer :: nSP                        ! The number of particles in SWARM
    integer :: NITMOPSO                   ! The number of MOPSO cycles 
    logical :: Random                     ! Random procedure: .true. = Random_Number() .false. = RAND()
    integer :: NIM                        ! Code number to initialize the swarm
    integer :: NFI1_4                     ! number of fixed iterations in the first quarter.
    integer :: NFI2_4                     ! number of fixed iterations in the second quarter.
    integer :: NFI2_2                     ! number of fixed iterations in the second half.
    real(kind=rp), allocatable :: S_dominance(:)           ! Dominance Criteria
  end type MOPSO_Light_Parameters_type

contains

  subroutine Set_MOPSO_Light_Parameters(nY,MP,ParametersFileName)
      implicit none

  integer,intent(in):: nY
  type(MOPSO_Light_Parameters_type) :: MP      !     The structure that contains all optimization parameters 
  character (len=*), intent(in) :: ParametersFileName     !     The file name that contains the MOPSO parameters.
  integer i
    open (UNIT=30,FILE=ParametersFileName)
    read(30,*)MP%NPFS              ! Number of storage positions for particles from Pareto's Front
    read(30,*)MP%nSP               ! Number of particles at the swarm
    read(30,*)MP%NITMOPSO          ! Method Iterations
    read(30,*)MP%Random            ! Logical .True. or .False. to set random function
    read(30,*)MP%NIM               ! Code number to initialize swarm particles POSITION (P) and VELOCITY (V)  (1) P=RANDOM V=0; (2) P=SOBOL(max of NX=6 dimensions) V=0;(3) P=RANDOM V=RANDOM;(4) P=SOBOL V=RANDOM.
    read(30,*)MP%NFI1_4            ! Number of fixed interations on the first quarter
    read(30,*)MP%NFI2_4            ! Number of fixed interations on the second quarter
    read(30,*)MP%NFI2_2            ! Number of fixed interations on the second half
    allocate (MP%S_dominance(ny))
    do i = 1 , nY
      read(30,*)MP%S_dominance(i)
    end do
    close (UNIT=30, STATUS='KEEP')

    return
  end subroutine Set_MOPSO_Light_Parameters



  function MOPSO_Light (OF,MLP,PFin) result(PF)
    implicit none

    class(ObjectiveFunction_MOPSO_Class), intent(in):: OF     
    type(MOPSO_Light_Parameters_type), intent(in):: MLP 
    type(ParetoFront_EDSD_class),optional, intent(in)  :: PFin
    type(ParetoFront_EDSD_class) :: PF
    ! ** MOPSO_Light algorithm variables
    type(Swarm_MOPSO_type), allocatable :: SW (:)  ! Particle Swarm  with nSP particles
    integer :: i,j,k,L,m
    integer :: Nat,iNAT                            ! NAT , Number of iterations to keep the global leader fixed, iNAT auxiliar index
    integer :: IqNIT,Iq1,Iq2                       ! First and second quarters of iterative process
    integer :: nX,nY                               
    real(kind=rp) :: dpL                           ! An auxiliary variable that measures the distance to the domain bound
    real(kind=rp) :: omegamin,omegamax,q,C1,C2     ! PSO Parameters
    real(kind=rp) :: C_omega                       ! The auxiliar PSO Parameter
    parameter (omegamin=0.4_rp)                    ! REF #Chatterjee and Siarry, 2006 [3]
    parameter (omegamax=0.9_rp)                    ! REF #Chatterjee and Siarry, 2006 [3]
    parameter (q=1.2_rp)                           ! REF #Chatterjee and Siarry, 2006 [3]
    parameter (C1=2.05_rp)                         ! REF #Eberhart and Kennedy, 1995 [2]
    parameter (C2=2.05_rp)                         ! REF #Eberhart and Kennedy, 1995 [2]
    real(kind=rp) :: omega1                        ! PSO Parameter
    real R1,R2                                     ! The Random auxiliar variables
    real(kind=rp) :: maVA,VA                       ! Auxiliary variables
    logical EXCHANGE/.false./                      ! Auxiliary variable
    logical OUTDOM/.true./                         ! Auxiliary variable

    ! Initialize dimensions variables
    nX=OF%get_nX()      
    nY=OF%get_nY()      

    ! Initializes the ParetoFront_EDSD Object
    if (present(PFin))then
      PF=PFin
    else
      if ((allocated(MLP%S_dominance))) then
        PF=ParetoFront_EDSD_class(MLP%NPFS,nX,nY,MLP%S_dominance)
      else
        PF=ParetoFront_EDSD_class(MLP%NPFS,nX,nY)
      end if  
    end if
    
    ! Initialize the Swarm "SW" 
    SW = New_Swarm(OF,PF,MLP%nSP,nx,ny,MLP%NIM,MLP%Random)

    ! MOPSO_Light_algorithm
    ! Set the number of iterations before updating the leader in the first quarter
    NAT=MLP%NFI1_4      

    ! Set the end of the first and the end of the second quarter of the iterative process      
    iNAT=NAT                
    IqNIT=INT(MLP%NITMOPSO/4)        
    Iq1=1+MLP%NITMOPSO-3*IqNIT  
    Iq2=1+MLP%NITMOPSO-2*IqNIT  

    ! Calculate the constant part of the inertia factor (Omega)
    C_omega = (omegamax-omegamin) / (MLP%NITMOPSO**q)  

    write(*,*)" Number of non-dominated particles in the initial swarm SW  =  ",PF%get_NCSP()
    
    ! Iterative procedure
    do k=1,MLP%NITMOPSO
      ! After a quater of the total cycles, set a new number of iterations before updating the leader in the second quarter of the iterative process
      if (K.eq.Iq1) then  
        NAT=MLP%NFI2_4     
        write (*,*) "first quarter / second quarter cycles =" , k ,"number of particles in Pareto Front = ",PF%get_NCSP()
      end if 
      ! After the second quater of the total cycles, set a new number of iterations before updating the leader in the second half of the iterative process
      if (K.eq.Iq2) then    
        NAT=MLP%NFI2_2
        iNAT=NAT-mod((MLP%NITMOPSO-K),NAT)-1     
        write (*,*) "second quarter / first half cycles =" , k , "number of particles in Pareto Front = ",PF%get_NCSP()
      end if 

      ! Calculate the inertia Factor REF #Chatterjee and Siarry, 2006 [3]       
      omega1=omegamin+C_omega*((MLP%NITMOPSO-K)**q)

      ! Set the global best (leader) "Gb" each 3 to 8 iterations
      if(iNAT.GE.NAT) then
        iNAT=0
        R1=Local_RAND(MLP%Random)
        j=INT(PF%get_NCSP()*R1)  ! Random selection of an initial particle 
        do i=1,MLP%nSP
          SW(i)%Gb=PF%get_X(PF%get_NCSP()-MOD(i+j,PF%get_NCSP()))   ! sequential selection of the global leader
        end do
      end if
      iNat=iNat+1
      
      ! Calculate particles X-velocities
      do i=1,MLP%nSP
        do j=1,Nx
          R1=Local_RAND(MLP%Random)
          R2=Local_RAND(MLP%Random)
          SW(i)%Vx(j)=(omega1*SW(i)%Vx(j)+      &
          R1*C1*(SW(i)%Lb%x(j)-SW(i)%P%x(j))+    &
          R2*C2*(SW(i)%Gb(j)-SW(i)%P%x(j)) )
        end do
      end do

      ! Update swarm position informations      
      do i=1,MLP%nSP

        ! Update particles position     
        SW(i)%P%X=SW(i)%P%X+SW(i)%Vx

        ! Check if the new position is out of bounds and
        !  set a new position inside the bounds for outer particles
        do j=1,Nx
          do while (OUTDOM)
            OUTDOM=.false.
            dpL=-SW(i)%P%X(j)+ OF%liminfX(j)
            if (dpL.GT.0) then                     ! if the particle is less than the lower bound  
              SW(i)%P%X(j)=OF%liminfX(j)+dpL/10    ! The particle dampens ten times. 
              SW(i)%Vx(j)=SW(i)%Vx(j)/10           ! The particle moves back a tenth of the distance it has moved forward out of the bound.
            end if
            dpL=-SW(i)%P%X(j)+OF%limsupX(j)
            if (dpL.LT.0) then                       !  if the particle is greater than the upper bound  
              SW(i)%P%X(j)=OF%limsupX(j)+dpL/10       ! The particle dampens ten times. 
              SW(i)%Vx(j)=SW(i)%Vx(j)/10           ! The particle moves back a tenth of the distance it has moved forward out of the bound.
              OUTDOM=.true.
            end if
          end do
          OUTDOM=.true.
        end do
        
        ! Calculate the function value for particles new position
        call OF%Function_X(SW(i)%P%X,SW(i)%P%Y)

        ! Update particle's "best position"    
        j=1
        do while (j.LE.Ny)
        ! Test if particle is non-dominated by the current best position.
          if (SW(i)%P%Y(j).LT.SW(i)%Lb%Y(j)) then
            j = Ny
            ! Test if the current best position is non-dominated by the particle.
            L = 1
            do while (L.LE.Ny)
              if (SW(i)%Lb%Y(L).LT.SW(i)%P%Y(L)) then
                L=Ny
                ! If particles are not dominated by each other:
                ! Find the dimension with the highest relative variation for the current particle
                maVA=SW(i)%P%Y(1)/SW(i)%Lb%Y(1)-1.0d0
                do m=2,Ny
                  VA=SW(i)%P%Y(m)/SW(i)%Lb%Y(m)-1.0d0
                  if (VA.GT.maVA) maVA=VA
                end do
                m=1
                ! Search for a higher relative variation for the new particle.
                do while (m.LE.Ny)
                  if ((SW(i)%Lb%Y(m)/SW(i)%P%Y(m)-1.0d0).GT.maVa) then   
                    m=Ny                                            
                    EXCHANGE=.true.
                  end if
                  m=m+1
                end do
              end if
              L=L+1
            end do
            if (EXCHANGE)then
              SW(i)%Lb=SW(i)%P
            else
              EXCHANGE=.false.
            end if
          end if
          j=j+1
        end do

        ! Try to insert Particles into the Pareto Front 
        call PF%TrySolution_PFS(SW(i)%P)

      end do

    end do  

    return
  end function MOPSO_Light


  real function Local_RAND(RANDOM)
    implicit none

    logical,intent(IN) :: RANDOM

    if (RANDOM) then
      call RANDOM_NUMBER (Local_RAND)
      return
    end if  
    Local_RAND=RAND()

    return
  end function Local_RAND  

end module MOPSO_Light_mod