! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  File_Manager_ParetoFront_mod : Write and read Pareto Front from a 
! "ParetoFront_class" to a file  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ***************************************************************************************
! *** File_Manager_ParetoFront_mod.f03 (Version 1.1)  
! ***************************************************************************************
! *** Algorithm Application:
! *     Save_ParetoFront_in_file: recieve a PF and save the Y and X vectors in a file;
! *     Merge_ParetoFront_in_file: Open files containig a PF, merge them all and save in
! *  a new file.
! ***************************************************************************************

module File_Manager_ParetoFront_mod
  use Precision_defaults_MOPSO_mod
  use ParetoFront_mod

  implicit none

  private 

  public :: Save_ParetoFront_in_file, Merge_ParetoFront_in_file

contains

  ! ** Save_ParetoFront_file
  !    save the Pareto Front solutions in a .csv file format      
  subroutine Save_ParetoFront_in_file(PF1,File_Name)
    use ParetoFront_mod
    implicit none
  
    class(ParetoFront_class),intent(in) :: PF1
    character(len=*) :: File_Name
    integer ::nX,nY
    type(SPE_ParetoFront_type):: Paux  ! Auxiliar Particle (nX,nY)
    character(len=28),allocatable :: buffer(:)
    integer i,j
    
    ! Initialize the file to store the results
    nX = PF1%get_nX()
    nY = PF1%get_nY()
    OPEN (UNIT=55,FILE=trim(File_Name),STATUS='REPLACE')
    ! Write real numbers in a .csv file format
    allocate (buffer(nY+nX))
    allocate (Paux%X(nX))
    allocate (Paux%Y(nY))
    do i=1,PF1%get_NCSP()
      Paux%X=PF1%get_X(i)
      Paux%Y=PF1%get_Y(i)
      do j=1,nY
        write( buffer(j), '(E25.18,a)' ) Paux%Y(j)," ,"
      end do  
      do j=1,nX
        write( buffer(nY+j), '(E25.18,a)' ) Paux%X(j)," ,"
      end do  
      buffer = adjustl(buffer)
      write(55,*)buffer
    end do
    close (UNIT=55, STATUS='KEEP')

    return
  end subroutine Save_ParetoFront_in_file

  !** MergeParetoFront
  !   Merge solution vectors stored in different files (.csv format) into a  
  !   single Pareto Front, and save it to a file        
  subroutine Merge_ParetoFront_in_file(PF1,All_File_Names,File_Name,NTPIF)
  Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    class(ParetoFront_class) :: PF1                            
    character(len=*), allocatable :: All_File_Names(:)         ! The array of file names that contains the Pareto Fronts selected to merge
    character(len=*) :: File_Name                              ! Name of the file to save the merged Pareto front
    integer, intent(out) :: NTPIF                              ! Number of Total Points In File with the final merged Pareto Front
    integer :: i,j,k
    integer :: nX,nY
    integer :: error 
    character(len=28),allocatable :: buffer(:)
    type(SPE_ParetoFront_type):: Paux  ! Auxiliar Particle (nX,nY)
    real(kind=rp):: re1,re2,re3,re4
    nX = PF1%get_nX()
    nY = PF1%get_nY()
    allocate (Paux%X(nX))
    allocate (Paux%Y(nY))
    ! Read points in file and try to insert in Pareto Front

    do k=1,size(All_File_Names)
      OPEN (UNIT=55,FILE=trim(All_File_Names(k))) 
      error=0
      do
        read(55,*,iostat=error) Paux%Y,Paux%X
        select case(error)
        case( 0 )
          call PF1%TrySolution_PFS(Paux)
        case(iostat_end)
          exit
        case default
          Write( *, * ) 'Error in reading file'
          Stop
        end select  
      end do
      CLOSE (UNIT=55, STATUS='KEEP')
    end do
    NTPIF=PF1%get_NCSP()
    ! Sort and save the Pareto front to a file
    call PF1%Sort_PFS(1)
    call Save_ParetoFront_in_file(PF1,File_Name)

    return
  end subroutine Merge_ParetoFront_in_file

end module File_Manager_ParetoFront_mod