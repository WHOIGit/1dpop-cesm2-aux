!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module forcing_gterms

! !MODULE: forcing_gterms
!
! Created by Ivan Lima on Mon Mar  8 2021 08:13:59 -0700
! Time-stamp: <2021-04-06 11:05:10 ivan>
!
! !DESCRIPTION:
!  Contains routines and variables used for reading and applying
!  the G-terms forcing for TEMP, SALT UVEL and VVEL.

   use kinds_mod
   use domain
   use constants
   use broadcast
   use io
   ! use forcing_tools
   use time_management
   ! use prognostic
   use grid
   use strdata_interface_mod_cesm1p1
   use tavg
   use forcing_fields, only: FGT, FGS, FGU, FGV
   use timers
   use POP_ErrorMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_HaloMod
   use exit_mod
   use shr_string_mod, only : shr_string_listGetIndexF

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_gterms,      &
             get_gterms_data

! !PUBLIC DATA MEMBERS:

   logical (POP_logical), public :: &
      lsingle_col_force_f1d   ! flag to specify single column f1d external forcing is turned on

   character(char_len), public :: &
      f1d_filename,               &! 1D mode input file for depth fields
      f1d_fieldlistfile,          &! 1D mode file field list for depth fields
      f1d_fieldlistmodel           ! 1D mode model field list for depth fields

!-----------------------------------------------------------------------
!
!  internal module variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_FGT,          &! GT term read in from stream file
      tavg_FGS,          &! GS term read in from stream file
      tavg_FGU,          &! GU term read in from stream file
      tavg_FGV            ! GV term read in from stream file

!-----------------------------------------------------------------------
!
!  1d stream forcing data
!
!-----------------------------------------------------------------------

   type(strdata_input_type) :: &
      F1d_inputlist       ! pop stream datatype for 1d vertical forcing data

   integer (int_kind), public :: &
      F1d_GT_ind,                &! index into the stream av for GT
      F1d_GS_ind,                &! index into the stream av for GS
      F1d_GU_ind,                &! index into the stream av for GU
      F1d_GV_ind                  ! index into the stream av for GV

   integer (int_kind) :: &
      F1d_shr_strdata_advance_timer  ! timer

   logical (log_kind) :: first_call_strdata_create = .true.

!***********************************************************************

 contains

!***********************************************************************

subroutine init_gterms

!   Initilize namelist and set tavg fields

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   namelist /forcing_gterms_nml/f1d_filename, f1d_fieldlistfile, f1d_fieldlistmodel

   integer (POP_i4) :: &
      nml_error           ! namelist i/o error flag

!-----------------------------------------------------------------------
!
!  read input namelist for G-terms forcing options
!
!-----------------------------------------------------------------------

   lsingle_col_force_f1d = .false.
   f1d_filename          = 'unknown_f1d_filename'
   f1d_fieldlistfile     = 'unknown_f1d_fieldlistfile'
   f1d_fieldlistmodel    = 'unknown_f1d_fieldlistmodel'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=forcing_gterms_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR: reading forcing_gterms_nml')
   endif

   call broadcast_scalar(f1d_filename,         master_task)
   call broadcast_scalar(f1d_fieldlistfile,    master_task)
   call broadcast_scalar(f1d_fieldlistmodel,   master_task)

   ! Write namelist to log
   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,*) ' forcing_gterms_nml namelist settings:'
       write(stdout,blank_fmt)
       write(stdout, forcing_gterms_nml)
       write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!
!  tavg fields
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_FGT,'FGT',3,                                   &
                          long_name='Unrepresented large-scale T forcing',    &
                          units='degC/sec', grid_loc='3111',                  &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_FGS,'FGS',3,                                   &
                          long_name='Unrepresented large-scale S forcing',    &
                          units='gram/kilogram/sec', grid_loc='3111',         &
                          scale_factor=1000.0_r8,                             &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_FGU,'FGU',3,                                   &
                          long_name='Unrepresented large-scale U forcing',    &
                          units='centimeter/sec^2', grid_loc='3221',          &
                          coordinates='ULONG ULAT z_t time')

   call define_tavg_field(tavg_FGV,'FGV',3,                                   &
                          long_name='Unrepresented large-scale V forcing',    &
                          units='centimeter/sec^2', grid_loc='3221',          &
                          coordinates='ULONG ULAT z_t time')

!-----------------------------------------------------------------------

   FGT = c0
   FGS = c0
   FGU = c0
   FGV = c0

   call POP_strdata_init(F1d_inputlist)

   F1d_inputlist%field_name = 'F1d data'
   F1d_inputlist%short_name = 'F1d'
   F1d_inputlist%year_first = 1
   F1d_inputlist%year_last  = 1
   F1d_inputlist%year_align = 1
   F1d_inputlist%file_name  = trim(f1d_filename)
   F1d_inputlist%field_list_file  = trim(f1d_fieldlistfile)
   F1d_inputlist%field_list_model = trim(f1d_fieldlistmodel)
   F1d_inputlist%dtlimit    = 1.e30

   ! returns index of 0 if string is not defined
   F1d_GT_ind = shr_string_listGetIndexF(F1d_inputlist%field_list_model,'GT')
   F1d_GS_ind = shr_string_listGetIndexF(F1d_inputlist%field_list_model,'GS')
   F1d_GU_ind = shr_string_listGetIndexF(F1d_inputlist%field_list_model,'GU')
   F1d_GV_ind = shr_string_listGetIndexF(F1d_inputlist%field_list_model,'GV')

   lsingle_col_force_f1d = .false.
   if (F1d_GT_ind > 0 .or. F1d_GS_ind > 0 .or. F1d_GU_ind > 0 .or. F1d_GV_ind > 0) then
      lsingle_col_force_f1d=.true.
   endif

   call get_timer(F1d_shr_strdata_advance_timer, &
                  'F1D_SHR_STRDATA_ADVANCE',1, distrb_clinic%nprocs)
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a,l4)') ' 1D f1d forcing data shr_stream: ',lsingle_col_force_f1d
      write(stdout,'(2a)')   '   f1d_filename        = ',trim(F1d_inputlist%file_name)
      write(stdout,'(2a)')   '   f1d_file_list_file  = ',trim(F1d_inputlist%field_list_file)
      write(stdout,'(2a)')   '   f1d_file_list_model = ',trim(F1d_inputlist%field_list_model)
      write(stdout,'(a,i6)') '   GT index        = ',F1d_GT_ind
      write(stdout,'(a,i6)') '   GS index        = ',F1d_GS_ind
      write(stdout,'(a,i6)') '   GU index        = ',F1d_GU_ind
      write(stdout,'(a,i6)') '   GV index        = ',F1d_GV_ind
   endif

!-----------------------------------------------------------------------

end subroutine init_gterms

!***********************************************************************

subroutine get_gterms_data

!   Read input G-terms forcing data

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,                &! dummy indices for horizontal directions
      n,k,                &! dummy indices for vertical level, tracer
      iblock,             &! counter for block loops
      kp1,km1              ! level index for k+1, k-1 levels

   type (block) ::        &
      this_block           ! block information for current block

!-----------------------------------------------------------------------
!
!  read 1d forcing stream dataset
!
!-----------------------------------------------------------------------

   if (lsingle_col_force_f1d) then
      if (first_call_strdata_create) then
         if (lsingle_col_force_f1d) &
            call POP_strdata_create(F1d_inputlist,depthflag=.true., &
            XvarName='TLONG',YvarName='TLAT',ZvarName='z_t',TvarName='time')
      endif
      first_call_strdata_create = .false.

      call timer_start(F1d_shr_strdata_advance_timer)

      if (lsingle_col_force_f1d) then
         F1d_inputlist%date = iyear*10000 + imonth*100 + iday
         F1d_inputlist%time = isecond + 60 * (iminute + 60 * ihour)
         call POP_strdata_advance(F1d_inputlist)

         n = 0
         do k=1,km
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               if (F1d_GT_ind > 0) &
                  FGT(i,j,k,iblock) = F1d_inputlist%sdat%avs(1)%rAttr(F1d_GT_ind,n)
               if (F1d_GS_ind > 0) &
                  FGS(i,j,k,iblock) = F1d_inputlist%sdat%avs(1)%rAttr(F1d_GS_ind,n)
               if (F1d_GU_ind > 0) &
                  FGU(i,j,k,iblock) = F1d_inputlist%sdat%avs(1)%rAttr(F1d_GU_ind,n)
               if (F1d_GV_ind > 0) &
                  FGV(i,j,k,iblock) = F1d_inputlist%sdat%avs(1)%rAttr(F1d_GV_ind,n)
            enddo
            enddo
            FGT(:,:,k,iblock) = merge(FGT(:,:,k,iblock), c0, FGT(:,:,k,iblock) <= 1.e30)
            FGS(:,:,k,iblock) = merge(FGS(:,:,k,iblock), c0, FGS(:,:,k,iblock) <= 1.e30)
            FGU(:,:,k,iblock) = merge(FGU(:,:,k,iblock), c0, FGU(:,:,k,iblock) <= 1.e30)
            FGV(:,:,k,iblock) = merge(FGV(:,:,k,iblock), c0, FGV(:,:,k,iblock) <= 1.e30)

            call accumulate_tavg_field(FGT(:,:,k,iblock), tavg_FGT, iblock, k)
            call accumulate_tavg_field(FGS(:,:,k,iblock), tavg_FGS, iblock, k)
            call accumulate_tavg_field(FGU(:,:,k,iblock), tavg_FGU, iblock, k)
            call accumulate_tavg_field(FGV(:,:,k,iblock), tavg_FGV, iblock, k)

         enddo
         enddo
      endif

      call timer_stop(F1d_shr_strdata_advance_timer)

   endif

!-----------------------------------------------------------------------

end subroutine get_gterms_data

!***********************************************************************

end module forcing_gterms

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
