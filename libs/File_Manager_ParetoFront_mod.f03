! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>
! @def :  Pareto Front - Class (general)
module File_Manager_ParetoFront_mod
  use Precision_defaults_MOPSO_mod
  use ParetoFront_mod

  implicit none

  private 

  public :: Save_ParetoFront_in_file, Merge_ParetoFront_in_file

contains

! Save_ParetoFront_file
! ** adds the Pareto Front solutions to a .csv file      
  subroutine Save_ParetoFront_in_file(PF1,nX,nY,File_Name)
    use ParetoFront_mod
    implicit none
  
    class(ParetoFront_class),intent(in) :: PF1
    integer, intent (in) ::nX,nY
    character(len=*) :: File_Name
    type (SPE_ParetoFront_type) :: P
    integer i,j
    real(kind=rp),allocatable :: PX(:),PY(:)
    character(len=28),allocatable :: buffer(:)
    ! Initializes the file to store the results
    OPEN (UNIT=55,FILE=trim(File_Name),STATUS='REPLACE')
    allocate (buffer(nY+nX))
    allocate (PX(nX))
    allocate (PY(nY))
    do i=1,PF1%get_NCSP()
      PX=PF1%get_X(i)
      PY=PF1%get_Y(i)
      do j=1,nY
        write( buffer(j), '(E25.18,a)' ) PY(j),";"
      end do  
      do j=1,(nX)
        write( buffer(nY+j), '(E25.18,a)' ) PX(j),";"
      end do  
      buffer = adjustl(buffer)
      write(55,*)buffer
    end do
    deallocate (buffer)
    deallocate (PX)
    deallocate (PY)
    close (UNIT=55, STATUS='KEEP')

    return
  end subroutine Save_ParetoFront_in_file

! ===============================================================================
!                    SUBROUTINE         MergeParetoFront
! ===============================================================================
!   Merge solutions vectors stored in file MOPSOLight_complete.csv from all runs 
!   in only one final Pareto Front, and save it to file MOPSOLight_final.csv       
  subroutine Merge_ParetoFront_in_file(PF,nX,nY,n_files,All_File_Names,n_P,File_Name,NTPIF)
    implicit none

    class(ParetoFront_class) :: PF
    integer, intent(in) :: nX,nY,n_files
    character(len=*), allocatable :: All_File_Names(:)
    integer, intent(in), allocatable :: n_P(:)
    character(len=*) :: File_Name
    integer, intent(out) :: NTPIF
    integer i,j,k
    character(len=28),allocatable :: buffer(:)
    type(SPE_ParetoFront_type):: Paux  ! Auxiliar Particle (nX,nY)

    allocate (Paux%X(nX))
    allocate (Paux%Y(nY))
    ! Begin merging processes
    do k=1,n_files
      OPEN (UNIT=55,FILE=trim(All_File_Names(k))) 
      do i=1,n_P(k)
        read(55,*) Paux%Y,Paux%X
        call PF%TrySolution_PFS(Paux)
      end do
      CLOSE (UNIT=55, STATUS='KEEP')
    end do
    NTPIF=PF%get_NCSP()
    call PF%Sort_PFS(1)

    call Save_ParetoFront_in_file(PF,nX,nY,File_Name)

    return
  end subroutine Merge_ParetoFront_in_file

end module File_Manager_ParetoFront_mod

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
