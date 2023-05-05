!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up

!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
module m_start_up

    ! Dependencies =============================================================
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

#ifdef MFC_MPI
    use mpi                      !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_check_patches

    use m_helper

    use m_checker
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_start_up_module, &
                     s_read_input_file, &
                     s_check_input_file, &
                     s_read_grid_data_files, &
                     s_read_ic_data_files, &
                     s_read_serial_grid_data_files, &
                     s_read_serial_ic_data_files, &
                     s_read_parallel_grid_data_files, &
                     s_read_parallel_ic_data_files, &
                     s_check_grid_data_files, &
                     s_finalize_start_up_module

    abstract interface ! ===================================================

        subroutine s_read_abstract_grid_data_files()! ----------

        end subroutine s_read_abstract_grid_data_files ! ---------------

        subroutine s_read_abstract_ic_data_files(q_cons_vf) ! -----------

            import :: scalar_field, sys_size

            ! Conservative variables
            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_cons_vf

        end subroutine s_read_abstract_ic_data_files ! -----------------

    end interface ! ========================================================

    character(LEN=path_len + name_len) :: proc_rank_dir !<
    !! Location of the folder associated with the rank of the local processor

    character(LEN=path_len + 2*name_len), private :: t_step_dir !<
    !! Possible location of time-step folder containing preexisting grid and/or
    !! conservative variables data to be used as starting point for pre-process

    procedure(s_read_abstract_grid_data_files), pointer :: s_read_grid_data_files => null()
    procedure(s_read_abstract_ic_data_files), pointer :: s_read_ic_data_files => null()

contains

    !>  Reads the configuration file pre_process.inp, in order to
        !!      populate the parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    subroutine s_read_input_file() ! ---------------------------------------

        character(LEN=name_len) :: file_loc  !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: iostatus
            !! Integer to check iostat of file read

        ! Namelist for all of the parameters to be inputed by the user
        namelist /user_inputs/ case_dir, old_grid, old_ic, &
            t_step_old, t_step_start, m, n, p, x_domain, y_domain, z_domain, &
            stretch_x, stretch_y, stretch_z, a_x, a_y, &
            a_z, x_a, y_a, z_a, x_b, y_b, z_b, &
            model_eqns, num_fluids, &
            adv_alphan, mpp_lim, &
            weno_order, bc_x, bc_y, bc_z, num_patches, &
            hypoelasticity, patch_icpp, fluid_pp, &
            precision, parallel_io, &
            perturb_flow, perturb_flow_fluid, &
            perturb_sph, perturb_sph_fluid, fluid_rho, &
            cyl_coord, loops_x, loops_y, loops_z, &
            rhoref, pref, bubbles, R0ref, nb, &
            polytropic, thermal, Ca, Web, Re_inv, &
            polydisperse, poly_sigma, qbmm, &
            sigR, sigV, dist_type, rhoRV, R0_type

        ! Inquiring the status of the pre_process.inp file
        file_loc = 'pre_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)
            if (iostatus /= 0) then
                call s_mpi_abort('Invalid line in pre_process.inp. It is '// &
                'likely due to a datatype mismatch. Exiting ...')
            end if
            close (1)
            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p
        else
            call s_mpi_abort('File pre_process.inp is missing. Exiting ...')
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !>  Checking that the user inputs make sense, i.e. that the
    !!      individual choices are compatible with the code's options
    !!      and that the combination of these choices results into a
    !!      valid configuration for the pre-process
    subroutine s_check_input_file() ! --------------------------------------

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            print '(A)', 'WARNING: Ensure that compiler flags/choices in Makefiles match your compiler! '
            print '(A)', 'WARNING: Ensure that preprocessor flags are enabled! '
            call s_mpi_abort('Unsupported choice for the value of case_dir.' // &
                'Exiting ...')
        end if

        call s_check_inputs()

        ! Check all the patch properties
        call s_check_patches()

    end subroutine s_check_input_file ! ------------------------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      grid data as well as based on the imported grid, complete
        !!      the necessary global computational domain parameters.
    subroutine s_read_serial_grid_data_files() ! ---

        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Logical variable used to test the existence of folders
        logical :: dir_check

        ! Generic logical used for the purpose of asserting whether a file
        ! is or is not present in the designated location
        logical :: file_check

        ! Setting address of the local processor rank and time-step directory
        write (proc_rank_dir, '(A,I0)') '/p_all/p', proc_rank
        proc_rank_dir = trim(case_dir)//trim(proc_rank_dir)

        write (t_step_dir, '(A,I0)') '/', t_step_start
        t_step_dir = trim(proc_rank_dir)//trim(t_step_dir)

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(t_step_dir)//'/.'
        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the pre-process exits
        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Time-step folder '//trim(t_step_dir)// &
                ' is missing. Exiting ...')
        end if

        ! Reading the Grid Data File for the x-direction ===================

        ! Checking whether x_cb.dat exists
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! If it exists, x_cb.dat is read
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) x_cb(-1:m)
            close (1)
        else 
            call s_mpi_abort('File x_cb.dat is missing in '// &
                trim(t_step_dir)//'. Exiting ...')
        end if

        ! Computing cell-center locations
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0

        ! Computing minimum cell-width
        dx = minval(x_cb(0:m) - x_cb(-1:m - 1))
        if (num_procs > 1) call s_mpi_reduce_min(dx)

        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        ! ==================================================================

        ! Reading the Grid Data File for the y-direction ===================

        if (n > 0) then

            ! Checking whether y_cb.dat exists
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, y_cb.dat is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) y_cb(-1:n)
                close (1)
            else
                call s_mpi_abort('File y_cb.dat is missing in '// &
                    trim(t_step_dir)//'. Exiting ...')
            end if

            ! Computing cell-center locations
            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0

            ! Computing minimum cell-width
            dy = minval(y_cb(0:n) - y_cb(-1:n - 1))
            if (num_procs > 1) call s_mpi_reduce_min(dy)

            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            ! ==================================================================

            ! Reading the Grid Data File for the z-direction ===================

            if (p > 0) then

                ! Checking whether z_cb.dat exists
                file_loc = trim(t_step_dir)//'/z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_check)

                ! If it exists, z_cb.dat is read
                if (file_check) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) z_cb(-1:p)
                    close (1)
                else
                    call s_mpi_abort('File z_cb.dat is missing in '// &
                        trim(t_step_dir)//'. Exiting ...')
                end if

                ! Computing cell-center locations
                z_cc(0:p) = (z_cb(0:p) + z_cb(-1:(p - 1)))/2d0

                ! Computing minimum cell-width
                dz = minval(z_cb(0:p) - z_cb(-1:p - 1))
                if (num_procs > 1) call s_mpi_reduce_min(dz)

                ! Setting locations of domain bounds
                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)

            end if

        end if

        ! ==================================================================

        ! If only the preexisting grid data files are read in and there will
        ! not be any preexisting initial condition data files imported, then
        ! the directory associated with the rank of the local processor may
        ! be cleaned to make room for the new pre-process data. In addition,
        ! the time-step directory that will contain the new grid and initial
        ! condition data are also generated.
        if (old_ic .neqv. .true.) then
            call s_delete_directory(trim(proc_rank_dir)//'/*')
            call s_create_directory(trim(proc_rank_dir)//'/0')
        end if

    end subroutine s_read_serial_grid_data_files ! --------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
    subroutine s_check_grid_data_files() ! -----------------

        ! Cell-boundary Data Consistency Check in x-direction ==============

        if (any(x_cb(0:m) - x_cb(-1:m - 1) <= 0d0)) then
            call s_mpi_abort('x_cb.dat in '//trim(t_step_dir)// &
                ' contains non-positive cell-spacings. Exiting ...')
        end if

        ! ==================================================================

        ! Cell-boundary Data Consistency Check in y-direction ==============

        if (n > 0) then

            if (any(y_cb(0:n) - y_cb(-1:n - 1) <= 0d0)) then
                call s_mpi_abort('y_cb.dat in '//trim(t_step_dir)// &
                    ' contains non-positive cell-spacings. '// &
                    'Exiting ...')
            end if

            ! ==================================================================

            ! Cell-boundary Data Consistency Check in z-direction ==============

            if (p > 0) then

                if (any(z_cb(0:p) - z_cb(-1:p - 1) <= 0d0)) then
                    call s_mpi_abort('z_cb.dat in '//trim(t_step_dir)// &
                        ' contains non-positive cell-spacings'// &
                        ' .Exiting ...')
                end if

            end if

        end if

        ! ==================================================================

    end subroutine s_check_grid_data_files ! -------------------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_serial_ic_data_files(q_cons_vf) ! ---------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc !<
        ! Generic string used to store the address of a particular file

        character(LEN= &
                  int(floor(log10(real(sys_size, kind(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: i !< Generic loop iterator

        ! Reading the Conservative Variables Data Files ====================
        do i = 1, sys_size

            ! Checking whether data file associated with variable position
            ! of the currently manipulated conservative variable exists
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'// &
                       trim(file_num)//'.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! If it exists, the data file is read
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf
                close (1)
            else
                call s_mpi_abort( 'File q_cons_vf'//trim(file_num)// &
                    '.dat is missing in '//trim(t_step_dir)// &
                    '. Exiting ...')
            end if

        end do

        ! ==================================================================

        ! Since the preexisting grid and initial condition data files have
        ! been read in, the directory associated with the rank of the local
        ! process may be cleaned out to make room for new pre-process data.
        ! In addition, the time-step folder that will contain the new grid
        ! and initial condition data are also generated.
        call s_create_directory(trim(proc_rank_dir)//'/*')
        call s_create_directory(trim(proc_rank_dir)//'/0')

    end subroutine s_read_serial_ic_data_files ! ----------------------------------

    !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
    subroutine s_read_parallel_grid_data_files()

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_ALL(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting... ')
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing cell center locations
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2d0
        ! Computing minimum cell width
        dx = minval(x_cb(0:m) - x_cb(-1:(m - 1)))
        if (num_procs > 1) call s_mpi_reduce_min(dx)
        ! Setting locations of domain bounds
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ_ALL(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting... ')
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing cell center locations
            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2d0
            ! Computing minimum cell width
            dy = minval(y_cb(0:n) - y_cb(-1:(n - 1)))
            if (num_procs > 1) call s_mpi_reduce_min(dy)
            ! Setting locations of domain bounds
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ_ALL(ifile, z_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting... ')
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing cell center locations
                z_cc(0:p) = (z_cb(0:p) + z_cb(-1:(p - 1)))/2d0
                ! Computing minimum cell width
                dz = minval(z_cb(0:p) - z_cb(-1:(p - 1)))
                if (num_procs > 1) call s_mpi_reduce_min(dz)
                ! Setting locations of domain bounds
                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)

            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_read_parallel_grid_data_files ! -----------------------

    !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.
        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_ic_data_files(q_cons_vf) ! ------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        integer :: i

        ! Open the file to read
        write (file_loc, '(I0,A)') t_step_start, '.dat'
        file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can read even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            do i = 1, sys_size
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                   MPI_DOUBLE_PRECISION, status, ierr)
            end do

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)

        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting... ')
        end if
        call s_mpi_barrier()

#endif

    end subroutine s_read_parallel_ic_data_files ! -------------------------

    subroutine s_initialize_start_up_module() !-----------------------------

        if (parallel_io .neqv. .true.) then
            s_read_grid_data_files => s_read_serial_grid_data_files
            s_read_ic_data_files => s_read_serial_ic_data_files
        else
            s_read_grid_data_files => s_read_parallel_grid_data_files
            s_read_ic_data_files => s_read_parallel_ic_data_files
        end if

    end subroutine s_initialize_start_up_module ! --------------------------

    subroutine s_finalize_start_up_module() ! ------------------------------

        s_read_grid_data_files => null()
        s_read_ic_data_files => null()

    end subroutine s_finalize_start_up_module ! ----------------------------

end module m_start_up
