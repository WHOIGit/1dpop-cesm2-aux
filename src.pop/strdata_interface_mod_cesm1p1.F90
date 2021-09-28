module strdata_interface_mod_cesm1p1

!BOP
! !MODULE: strdata_interface_mod_cesm1p1
!
! !DESCRIPTION:
!
!  Provide an interface to shr_strdata to use to read input data sets
!  such as
!  initial data fields for ecosystem tracer fields.
!
!  Initial work done by mlevy (July 2015)

  use shr_strdata_mod, only : shr_strdata_type
  use shr_strdata_mod, only : shr_strdata_create
  use shr_strdata_mod, only : shr_strdata_print
  use shr_strdata_mod, only : shr_strdata_advance
  use shr_pio_mod,     only : shr_pio_getiotype
  use shr_pio_mod,     only : shr_pio_getiosys

  use kinds_mod,        only : char_len
  use kinds_mod,        only : int_kind
  use kinds_mod,        only : r8
  use domain_size,      only : nx_global
  use domain_size,      only : ny_global
  use domain_size,      only : km
  use communicate,      only : my_task
  use communicate,      only : master_task
  use POP_IOUnitsMod,   only : inst_name
  use POP_CommMod,      only : POP_communicator
  use POP_MCT_vars_mod, only : POP_MCT_OCNID
  use POP_MCT_vars_mod, only : POP_MCT_gsMap_o, POP_MCT_gsMap3d_o
  use POP_MCT_vars_mod, only : POP_MCT_dom_o, POP_MCT_dom3d_o

  implicit none
  private

  public :: strdata_input_type
  public :: POP_strdata_init
  public :: POP_strdata_create
  public :: POP_strdata_advance

  type strdata_input_type
    ! Contains arguments for shr_strdata_* that are unique to variable
    ! being
    ! read
    type(shr_strdata_type)  :: sdat
    character(len=char_len) :: field_name
    character(len=char_len) :: short_name
    integer(int_kind)       :: year_first
    integer(int_kind)       :: year_last
    integer(int_kind)       :: year_align
    character(len=char_len) :: file_name
    character(len=char_len) :: field_list_file
    character(len=char_len) :: field_list_model
    integer(int_kind)       :: date
    integer(int_kind)       :: time
    real(r8)                :: dtlimit
    character(len=char_len) :: taxmode   ! time axis mode (cycle, extend, limit)
    character(len=char_len) :: tintalgo  ! time interpolation algorithm
  end type strdata_input_type

!-------------------------------------
contains
!-------------------------------------

  subroutine POP_strdata_init(inputlist)

    type(strdata_input_type), intent(inout) :: inputlist

    inputlist%field_name = 'undefined'
    inputlist%short_name = 'undefined'
    inputlist%year_first = 1
    inputlist%year_last  = 1
    inputlist%year_align = 1
    inputlist%file_name  = 'undefined'
    inputlist%field_list_file = 'undefined'
    inputlist%field_list_model = 'undefined'
    inputlist%date     = 0
    inputlist%time     = 0
    inputlist%dtlimit  = 1.5_r8
    inputlist%taxmode  = 'cycle'
    inputlist%tintalgo = 'linear'

  end subroutine POP_strdata_init

!-------------------------------------

  subroutine POP_strdata_create(inputlist,depthflag,TvarName,XvarName, &
     YvarName,ZvarName,MaskName,AreaName)

    type(strdata_input_type), intent(inout) :: inputlist
    logical, optional, intent(in) :: depthflag    ! true means 3d data expected
    character(len=*), optional, intent(in) :: TvarName
    character(len=*), optional, intent(in) :: XvarName
    character(len=*), optional, intent(in) :: YvarName
    character(len=*), optional, intent(in) :: ZvarName
    character(len=*), optional, intent(in) :: MaskName
    character(len=*), optional, intent(in) :: AreaName

    !--- local ---
    logical :: ldepthflag
    character(len=64) :: domTvarName,domXvarName,domYvarName,domZvarName, &
                         domAreaName,domMaskName

    ldepthflag = .false.
    if (present(depthflag)) then
       ldepthflag = depthflag
    endif

    domTvarName='time'
    domXvarName='TLONG'
    domYvarName='TLAT'
    domZvarName='depth'
    domAreaName='TAREA'
    domMaskName='KMT'
    if (present(TvarName)) then
       domTvarName = trim(TvarName)
    endif
    if (present(XvarName)) then
       domXvarName = trim(XvarName)
    endif
    if (present(YvarName)) then
       domYvarName = trim(YvarName)
    endif
    if (present(ZvarName)) then
       domZvarName = trim(ZvarName)
    endif
    if (present(MaskName)) then
       domMaskName = trim(MaskName)
    endif
    if (present(AreaName)) then
       domAreaName = trim(AreaName)
    endif

    if (ldepthflag) then
       !--- include nzg and domZvarName in call
       call shr_strdata_create(inputlist%sdat,name=trim(inputlist%field_name),&
                            mpicom=POP_communicator, &
                            compid=POP_MCT_OCNID, &
                            gsmap=POP_MCT_gsMap3d_o, ggrid=POP_MCT_dom3d_o,   &
                            nxg=nx_global, nyg=ny_global, nzg=km, &
                            yearFirst=inputlist%year_first, &
                            yearLast=inputlist%year_last, &
                            yearAlign=inputlist%year_align, &
                            offset=0, &
                            domFilePath='', &
                            domFileName=inputlist%file_name, &
                            domTvarName=trim(domTvarName), &
                            domXvarName=trim(domXvarName), &
                            domYvarName=trim(domYvarName), &
                            domZvarName=trim(domZvarName), &
                            domAreaName=trim(domAreaName), &
                            domMaskName=trim(domMaskName), &
                            FilePath='', &
                            FileName=(/trim(inputlist%file_name)/), &
                            fldListFile=inputlist%field_list_file, &
                            fldListModel=inputlist%field_list_model, &
                            pio_subsystem=shr_pio_getiosys(inst_name), &
                            pio_iotype=shr_pio_getiotype(inst_name), &
                            fillalgo='none', mapalgo='none', &
                            taxMode=inputlist%taxMode, &
                            tintalgo=inputlist%tintalgo, &
                            dtlimit=inputlist%dtlimit)
    else
       call shr_strdata_create(inputlist%sdat,name=trim(inputlist%field_name),&
                            mpicom=POP_communicator, &
                            compid=POP_MCT_OCNID, &
                            gsmap=POP_MCT_gsMap_o, ggrid=POP_MCT_dom_o, &
                            nxg=nx_global, nyg=ny_global, &
                            yearFirst=inputlist%year_first, &
                            yearLast=inputlist%year_last, &
                            yearAlign=inputlist%year_align, &
                            offset=0, &
                            domFilePath='', &
                            domFileName=inputlist%file_name, &
                            domTvarName=trim(domTvarName), &
                            domXvarName=trim(domXvarName), &
                            domYvarName=trim(domYvarName), &
                            domAreaName=trim(domAreaName), &
                            domMaskName=trim(domMaskName), &
                            FilePath='', &
                            FileName=(/trim(inputlist%file_name)/), &
                            fldListFile=inputlist%field_list_file, &
                            fldListModel=inputlist%field_list_model, &
                            pio_subsystem=shr_pio_getiosys(inst_name), &
                            pio_iotype=shr_pio_getiotype(inst_name), &
                            fillalgo='none', mapalgo='none', &
                            taxMode=inputlist%taxMode, &
                            tintalgo=inputlist%tintalgo, &
                            dtlimit=inputlist%dtlimit)
    endif

    if (my_task == master_task) then
      call shr_strdata_print(inputlist%sdat)
    endif

  end subroutine POP_strdata_create

!-------------------------------------

  subroutine POP_strdata_advance(inputlist)

    type(strdata_input_type), intent(inout) :: inputlist

    call shr_strdata_advance(inputlist%sdat, &
                             inputlist%date, &
                             inputlist%time, &
                             POP_communicator, &
                             trim(inputlist%short_name))

  end subroutine POP_strdata_advance

!-------------------------------------

end module strdata_interface_mod_cesm1p1
