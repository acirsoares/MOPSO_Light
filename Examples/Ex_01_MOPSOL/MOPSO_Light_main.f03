! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  teste - Test example for MOPSOLight Optimization Algorithm
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************************
program MOPSO_Light_main
!***********************************************************************
!     ==========================
!     Libraries for MOPSO_Light 
!     ==========================
!      - Precision_defaults_MOPSO_mod
!      - ParetoFront_EDSD_mod  --> (ParetoFront_mod)
!      - Swarm_MOPSO_mod --> (SobolSequence_mod)
!      - MOPSO_Light_mod 
!      - OF_test1_mod  --> (ObjectiveFunction_MOPSO_mod)
!      - File_Manager_ParetoFront_mod
!     ==========================
!     Main - Algorithm 
!     ==========================
!     Instantiate an ObjectiveFunction_MOPSOL object (OF_test1_Class)
!     Read MOPSOL parameters (DimensionsMOPSOL_type) 
!     Call MOPSOLight optimization algorithm 
!     Write Pareto Front to a File.
!
!***********************************************************************
! Modules
  use OF_test1_mod  
  use MOPSO_Light_mod
  use ParetoFront_EDSD_mod
  use File_Manager_ParetoFront_mod
  implicit none

  type(OF_test1_Class), allocatable :: OF !  Objective function object
  type(ParetoFront_EDSD_class), allocatable :: PF
  type(MOPSO_Light_Parameters_type) :: MOPSO_Light_Parameters, Set_MOPSO_Light_Parameters ! Structure that contains MOPSO parameters 
  integer :: i,nX,nY 
  integer :: NRs=10  
  integer :: NTParticles  ! The total particles in the output archive
  real :: T1,T2 
  character (len=10) :: File_id
  character (len=50) :: File_Name 
  character (len=50),allocatable :: All_File_Names (:)
  integer, allocatable :: ALL_nTotal_particles (:)
  call execute_command_line('clear')
 
  allocate(OF)
  OF = OF_test1_Class()

  nx=OF%get_nX() 
  nY=OF%get_nY() 

  MOPSO_Light_Parameters = Set_MOPSO_Light_Parameters(nY)

  ! Initializes the stored solutions counter    
  NTParticles=0

  allocate (PF)
  
  ! * Run the MOPSOL "NRs" times, evaluate the CPU time consuming, 
  ! *  ordenate Pareto Front in respect to the first objective function, 
  ! *  and save results in a file. 

  allocate (All_File_Names(NRs))
  allocate (ALL_nTotal_particles(NRs))

  do i=1,NRs

    call cpu_time(T1)
    PF = MOPSO_Light(OF,MOPSO_Light_Parameters)
    call cpu_time(T2)

    ALL_nTotal_particles(i)=PF%get_NCSP()
    write(*,"(a,i3,a,2i4)")"Number of particles in the",i,"st  Pareto Front = ", PF%get_NCSP()
    write(*,*)"CPU processing time = ", T2-T1

    NTParticles=NTParticles+PF%get_NCSP()

    call PF%Sort_PFS(1)

    if (i.lt.10) then
      write(File_id,"('0',i0)") i
    else
      write(File_id,"(i0)") i
    end if     
    File_Name='MOPSOL_test1_' // trim(adjustl(file_id)) // '.csv'
    All_File_Names(i)=File_Name
    call Save_ParetoFront_in_file(PF,File_Name)

  end do  

  write(*,*) "Particles in all fronts =",NTParticles  
  ! * Merge all Pareto Fronts
  File_Name='MOPSOLight_final.csv'
  PF = ParetoFront_EDSD_class(MOPSO_Light_Parameters%NPFS,nX,nY,MOPSO_Light_Parameters%S_dominance) 
  call Merge_ParetoFront_in_file(PF,NRs,All_File_Names,ALL_nTotal_particles,File_name,NTParticles)
  write(*,*) "Particles in Merged Pareto Front = ",PF%get_NCSP()

  deallocate(ALL_nTotal_particles)
  deallocate(All_File_Names)
  deallocate(PF)  
  deallocate(OF)  

!  call execute_command_line('python3 scatter2d.py')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end program MOPSO_Light_main

function Set_MOPSO_Light_Parameters(nY) result(MP)
  use MOPSO_Light_mod
  
  implicit none

  integer,intent(in):: nY
  type(MOPSO_Light_Parameters_type) :: MP      !     Structure that contains all optimization parameters 
  
  integer i

  open (UNIT=30,FILE='MOPSO_L_Dimension_Data.txt')
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
end function Set_MOPSO_Light_Parameters
