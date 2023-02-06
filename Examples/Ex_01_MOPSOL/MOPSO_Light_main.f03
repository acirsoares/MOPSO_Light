! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def : test - Test the example for MOPSOLight Optimization Algorithm
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************************
program MOPSO_Light_main
!***********************************************************************
!     ==========================
!     Libraries for MOPSO_Light 
!     ==========================
!      - Precision_defaults_MOPSO_mod
!      - ParetoFront_EDSD_mod  --> (ParetoFront_mod)
!      - Swarm_MOPSO_mod  <-- (SobolSequence_mod)
!      - MOPSO_Light_mod 
!      - OF_test1_mod  --> (ObjectiveFunction_MOPSO_mod)
!      - File_Manager_ParetoFront_mod
!     ==========================
!     Main - Algorithm 
!     ==========================
!     Instantiate an ObjectiveFunction_MOPSOL object (OF_test1_Class)
!     Read MOPSOL parameters (DimensionsMOPSOL_type) 
!     Call MOPSO_Light optimization algorithm 
!     Write Pareto Front to a File.
!
!***********************************************************************

! Modules
  use OF_test1_mod  
  use MOPSO_Light_mod
  use ParetoFront_EDSD_mod
  use File_Manager_ParetoFront_mod
  implicit none

! Variables
  type(OF_test1_Class), allocatable :: OF !  Objective function object
  type(ParetoFront_EDSD_class), allocatable :: PF
  type(MOPSO_Light_Parameters_type) :: MOPSO_Light_Parameters, Set_MOPSO_Light_Parameters ! Structure that contains MOPSO parameters 
  integer :: i
  integer :: NRs=10                       ! The number of runs 
  integer :: NTParticles                  ! The total particles in the output archive
  real :: T1,T2                           ! CPU timer variable
  character (len=10) :: File_id
  character (len=16) :: File_base = 'MOPSOL_test1_'
  character (len=16) :: Additional_path = 'results'
  character (len=255) :: File_Name        ! complete file name + path
  character (len=255), allocatable :: All_File_Names (:)    ! A set of file names with resulting Pareto Front
  character (len=200) :: cwd                 
  logical :: dir_e

  ! Set results directory name
  ! get the current working directory path
  call getcwd(cwd) 
  ! check whether there is the additional_path directory
  inquire(file=Additional_path, exist=dir_e)
  if (.not. dir_e) call execute_command_line ('mkdir -p ' // adjustl(trim( Additional_path ) ) )
  ! add a left and a right slash to the additional path name
  Additional_path=adjustl(trim('/'//Additional_path))//'/'


  ! instantiate the objective function object
  allocate(OF)
  OF = OF_test1_Class()

  ! Set parameters for MOPSO_Light algorithm   
  MOPSO_Light_Parameters = Set_MOPSO_Light_Parameters(OF%get_nY())

  ! Initializes the stored solutions counter    
  NTParticles=0

  ! Allocate a Pareto Front 
  allocate (PF)

  ! Allocate a set of file names  
  allocate (All_File_Names(NRs))
  
  ! *  Run "NRs" times 
  do i=1,NRs

    ! Run MOPSO_Light algorithm and evaluate the CPU time consuming
    call cpu_time(T1)
    PF = MOPSO_Light(OF,MOPSO_Light_Parameters)
    call cpu_time(T2)

    write(*,"(a,i3,a,2i4)")"Number of particles in the Pareto Front",i," =", PF%get_NCSP()
    write(*,*)"CPU processing time = ", T2-T1
    NTParticles=NTParticles+PF%get_NCSP()

    ! Ordenate Pareto Front in respect to the first objective function, 
    call PF%Sort_PFS(1)

    ! Set the file name and the writing results directory path
    if (i.lt.10) then
      write(File_id,"('0',i0)") i
    else
      write(File_id,"(i0)") i
    end if     
    File_Name=trim(cwd)//trim(Additional_path)//trim(File_base)//trim(adjustl(file_id))//'.csv'
    All_File_Names(i)=File_Name

    ! *  Save results in a file. 
    call Save_ParetoFront_in_file(PF,File_Name)

  end do  

  write(*,*) "Particles in all fronts =",NTParticles  

  ! Set a Pareto Front Object to merge all Files
  PF = ParetoFront_EDSD_class(MOPSO_Light_Parameters%NPFS,OF%get_nX(),OF%get_nY(),MOPSO_Light_Parameters%S_dominance) 

  ! Set the merger file name
  File_Name=trim(cwd)//trim(Additional_path)//'MOPSOLight_final.csv'

  ! * Merge all Pareto Fronts
  call Merge_ParetoFront_in_file(PF,All_File_Names,File_name,NTParticles)
  write(*,*) "Particles in Merged Pareto Front = ",PF%get_NCSP()

  deallocate(All_File_Names)
  deallocate(PF)  
  deallocate(OF)  

end program MOPSO_Light_main


! ** Set_MOPSO_Light_Parameters
!   Function to set the MOPSO_Light partameters 
function Set_MOPSO_Light_Parameters(nY) result(MP)
  use MOPSO_Light_mod  
  implicit none

  integer,intent(in):: nY
  type(MOPSO_Light_Parameters_type) :: MP      !     Structure that contains all optimization parameters 
  integer i

  open (UNIT=30,FILE='MOPSO_Light_Parameters.txt')
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
