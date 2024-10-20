
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_common

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

    implicit none

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: ierr
    !> @}

#ifndef MFC_PRE_PROCESS
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), q_cons_buff_send)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), q_cons_buff_recv)

    !$acc declare link(q_cons_buff_recv, q_cons_buff_send)
#else
    real(kind(0d0)), private, allocatable, dimension(:), target :: q_cons_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! conservative variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(kind(0d0)), private, allocatable, dimension(:), target :: q_cons_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average conservative variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !$acc declare create(q_cons_buff_send, q_cons_buff_recv)
#endif

    integer :: v_size

    !$acc declare create(v_size)
#endif

contains

    !> The subroutine initializes the MPI execution environment
        !!      and queries both the number of processors which will be
        !!      available for the job and the local processor rank.
    subroutine s_mpi_initialize

#ifndef MFC_MPI

        ! Serial run only has 1 processor
        num_procs = 1
        ! Local processor rank is 0
        proc_rank = 0

#else

        ! Initializing the MPI environment
        call MPI_INIT(ierr)

        ! Checking whether the MPI environment has been properly initialized
        if (ierr /= MPI_SUCCESS) then
            print '(A)', 'Unable to initialize MPI environment. Exiting ...'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Querying the number of processors available for the job
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        ! Querying the rank of the local processor
        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)

#endif

    end subroutine s_mpi_initialize

    subroutine s_initialize_mpi_common_module()

#if defined(MFC_MPI) && !defined(MFC_PRE_PROCESS)

        ! Allocating q_cons_buff_send/recv and ib_buff_send/recv. Please note that
        ! for the sake of simplicity, both variables are provided sufficient
        ! storage to hold the largest buffer in the computational domain.

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)))
            end if

            @:ALLOCATE_GLOBAL(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size + 2*nb*4
        else
#endif
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size))
            end if

            @:ALLOCATE_GLOBAL(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size
#if MFC_SIMULATION
        end if
#endif

        !$acc update device(v_size)

#endif
        
    end subroutine s_initialize_mpi_common_module

    subroutine s_finalize_mpi_common_module()

#if defined(MFC_MPI) && !defined(MFC_PRE_PROCESS)
        @:DEALLOCATE_GLOBAL(q_cons_buff_send, q_cons_buff_recv)
#endif

    end subroutine s_finalize_mpi_common_module

    !! @param q_cons_vf Conservative variables
    !! @param ib_markers track if a cell is within the immersed boundary
    subroutine s_initialize_mpi_data(q_cons_vf, ib_markers)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_cons_vf

        type(integer_field), &
            optional, &
            intent(in) :: ib_markers

        integer, dimension(num_dims) :: sizes_glb, sizes_loc
        integer, dimension(1) :: airfoil_glb, airfoil_loc, airfoil_start

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i, j, q, k, l

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m, 0:n, 0:p)
        end do

        !Additional variables pb and mv for non-polytropic qbmm
#ifdef MFC_PRE_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do j = 1, nnode
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j)%sf => pb%sf(0:m, 0:n, 0:p, j, i)
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j + nb*nnode)%sf => mv%sf(0:m, 0:n, 0:p, j, i)
                end do
            end do
        end if
#endif

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do j = 1, nnode
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j)%sf => pb_ts(1)%sf(0:m, 0:n, 0:p, j, i)
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j + nb*nnode)%sf => mv_ts(1)%sf(0:m, 0:n, 0:p, j, i)
                end do
            end do
        end if
#endif
        ! Define global(g) and local(l) sizes for flow variables
        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1
        if (n > 0) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1
            if (p > 0) then
                sizes_glb(3) = p_glb + 1; sizes_loc(3) = p + 1
            end if
        end if

        ! Define the view for each variable
        do i = 1, sys_size
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do

#ifndef MFC_POST_PROCESS
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*4
                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                              MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), ierr)
                call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)

            end do
        end if
#endif

        if (present(ib_markers)) then

#ifdef MFC_PRE_PROCESS
            MPI_IO_IB_DATA%var%sf => ib_markers%sf
#else
            MPI_IO_IB_DATA%var%sf => ib_markers%sf(0:m, 0:n, 0:p)
#endif
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, &
                                          MPI_ORDER_FORTRAN, MPI_INTEGER, MPI_IO_IB_DATA%view, ierr)
            call MPI_TYPE_COMMIT(MPI_IO_IB_DATA%view, ierr)

        end if

#ifndef MFC_POST_PROCESS
        if (present(ib_markers)) then
            do j = 1, num_ibs
                if (patch_ib(j)%c > 0) then

#ifdef MFC_PRE_PROCESS
                    allocate (MPI_IO_airfoil_IB_DATA%var(1:2*Np))
#endif

                    airfoil_glb(1) = 3*Np*num_procs
                    airfoil_loc(1) = 3*Np
                    airfoil_start(1) = 3*proc_rank*Np

#ifdef MFC_PRE_PROCESS
                    do i = 1, Np
                        MPI_IO_airfoil_IB_DATA%var(i)%x = airfoil_grid_l(i)%x
                        MPI_IO_airfoil_IB_DATA%var(i)%y = airfoil_grid_l(i)%y
                    end do
#endif

                    call MPI_TYPE_CREATE_SUBARRAY(1, airfoil_glb, airfoil_loc, airfoil_start, &
                                                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(1), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(1), ierr)

#ifdef MFC_PRE_PROCESS
                    do i = 1, Np
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%x = airfoil_grid_u(i)%x
                        MPI_IO_airfoil_IB_DATA%var(Np + i)%y = airfoil_grid_u(i)%y
                    end do
#endif
                    call MPI_TYPE_CREATE_SUBARRAY(1, airfoil_glb, airfoil_loc, airfoil_start, &
                                                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(2), ierr)
                    call MPI_TYPE_COMMIT(MPI_IO_airfoil_IB_DATA%view(2), ierr)

                end if
            end do

        end if
#endif

#endif

    end subroutine s_initialize_mpi_data

    subroutine mpi_bcast_time_step_values(proc_time, time_avg)

        real(kind(0d0)), dimension(0:num_procs - 1), intent(inout) :: proc_time
        real(kind(0d0)), intent(inout) :: time_avg

#ifdef MFC_MPI

        call MPI_GATHER(time_avg, 1, MPI_DOUBLE_PRECISION, proc_time(0), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine mpi_bcast_time_step_values

    !>  The goal of this subroutine is to determine the global
        !!      extrema of the stability criteria in the computational
        !!      domain. This is performed by sifting through the local
        !!      extrema of each stability criterion. Note that each of
        !!      the local extrema is from a single process, within its
        !!      assigned section of the computational domain. Finally,
        !!      note that the global extrema values are only bookkeept
        !!      on the rank 0 processor.
        !!  @param icfl_max_loc Local maximum ICFL stability criterion
        !!  @param vcfl_max_loc Local maximum VCFL stability criterion
        !!  @param Rc_min_loc Local minimum Rc stability criterion
        !!  @param icfl_max_glb Global maximum ICFL stability criterion
        !!  @param vcfl_max_glb Global maximum VCFL stability criterion
        !!  @param Rc_min_glb Global minimum Rc stability criterion
    subroutine s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                       vcfl_max_loc, &
                                                       ccfl_max_loc, &
                                                       Rc_min_loc, &
                                                       icfl_max_glb, &
                                                       vcfl_max_glb, &
                                                       ccfl_max_glb, &
                                                       Rc_min_glb)

        real(kind(0d0)), intent(in) :: icfl_max_loc
        real(kind(0d0)), intent(in) :: vcfl_max_loc
        real(kind(0d0)), intent(in) :: ccfl_max_loc
        real(kind(0d0)), intent(in) :: Rc_min_loc

        real(kind(0d0)), intent(out) :: icfl_max_glb
        real(kind(0d0)), intent(out) :: vcfl_max_glb
        real(kind(0d0)), intent(out) :: ccfl_max_glb
        real(kind(0d0)), intent(out) :: Rc_min_glb

#ifdef MFC_MPI
#ifdef MFC_SIMULATION

        ! Reducing local extrema of ICFL, VCFL, CCFL and Rc numbers to their
        ! global extrema and bookkeeping the results on the rank 0 processor
        call MPI_REDUCE(icfl_max_loc, icfl_max_glb, 1, &
                        MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                        MPI_COMM_WORLD, ierr)

        if (any(Re_size > 0)) then
            call MPI_REDUCE(vcfl_max_loc, vcfl_max_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                            MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(Rc_min_loc, Rc_min_glb, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MIN, 0, &
                            MPI_COMM_WORLD, ierr)
        end if

#endif
#endif

    end subroutine s_mpi_reduce_stability_criteria_extrema

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the sum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_sum(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_sum

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the minimum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_min(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MIN, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_min

    !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the maximum of all
        !!      values. The reduced variable is recorded back onto the
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
    subroutine s_mpi_allreduce_max(var_loc, var_glb)

        real(kind(0d0)), intent(in) :: var_loc
        real(kind(0d0)), intent(out) :: var_glb

#ifdef MFC_MPI

        ! Performing the reduction procedure
        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_allreduce_max

    !>  The following subroutine takes the inputted variable and
        !!      determines its minimum value on the entire computational
        !!      domain. The result is stored back into inputted variable.
        !!  @param var_loc holds the local value to be reduced among
        !!      all the processors in communicator. On output, the variable holds
        !!      the minimum value, reduced amongst all of the local values.
    subroutine s_mpi_reduce_min(var_loc)

        real(kind(0d0)), intent(inout) :: var_loc

#ifdef MFC_MPI

        ! Temporary storage variable that holds the reduced minimum value
        real(kind(0d0)) :: var_glb

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_MIN, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_min

    !>  The following subroutine takes the first element of the
        !!      2-element inputted variable and determines its maximum
        !!      value on the entire computational domain. The result is
        !!      stored back into the first element of the variable while
        !!      the rank of the processor that is in charge of the sub-
        !!      domain containing the maximum is stored into the second
        !!      element of the variable.
        !!  @param var_loc On input, this variable holds the local value and processor rank,
        !!  which are to be reduced among all the processors in communicator.
        !!  On output, this variable holds the maximum value, reduced amongst
        !!  all of the local values, and the process rank to which the value
        !!  belongs.
    subroutine s_mpi_reduce_maxloc(var_loc)

        real(kind(0d0)), dimension(2), intent(inout) :: var_loc

#ifdef MFC_MPI

        real(kind(0d0)), dimension(2) :: var_glb  !<
            !! Temporary storage variable that holds the reduced maximum value
            !! and the rank of the processor with which the value is associated

        ! Performing reduction procedure and eventually storing its result
        ! into the variable that was initially inputted into the subroutine
        call MPI_REDUCE(var_loc, var_glb, 1, MPI_2DOUBLE_PRECISION, &
                        MPI_MAXLOC, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, MPI_2DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb

#endif

    end subroutine s_mpi_reduce_maxloc

    !> The subroutine terminates the MPI execution environment.
        !! @param prnt error message to be printed
    subroutine s_mpi_abort(prnt)

        character(len=*), intent(in), optional :: prnt

        if (present(prnt)) then
            print *, prnt
            call flush (6)

        end if

#ifndef MFC_MPI

        stop 1

#else

        ! Terminating the MPI environment
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

#endif

    end subroutine s_mpi_abort

    !>Halts all processes until all have reached barrier.
    subroutine s_mpi_barrier

#ifdef MFC_MPI

        ! Calling MPI_BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_barrier

    !> The subroutine finalizes the MPI execution environment.
    subroutine s_mpi_finalize

#ifdef MFC_MPI

        ! Finalizing the MPI environment
        call MPI_FINALIZE(ierr)

#endif

    end subroutine s_mpi_finalize

#ifndef MFC_PRE_PROCESS
    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_variables_buffers(q_cons_vf, &
#ifdef MFC_SIMULATION
                                                pb, mv, &
#endif
                                                mpi_dir, &
                                                pbc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
#ifdef MFC_SIMULATION
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif
        integer, intent(in) :: mpi_dir, pbc_loc

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        type(bc_patch_parameters) :: bc
        
        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        integer :: pack_offset, unpack_offset
        real(kind(0d0)), pointer :: p_send, p_recv
        integer, pointer, dimension(:) :: p_i_send, p_i_recv

        integer :: w, r

#ifdef MFC_MPI

        !$acc update device(v_size)

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            buffer_counts = (/ &
                            buff_size*(sys_size + 2*nb*4)*(n + 1)*(p + 1), &
                            buff_size*(sys_size + 2*nb*4)*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        else
#endif
            buffer_counts = (/ &
                            buff_size*sys_size*(n + 1)*(p + 1), &
                            buff_size*sys_size*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
#ifdef MFC_SIMULATION
        end if
#endif

        do w = 1, num_bc_patches
            bc = patch_bc(w)

            buffer_count = buffer_counts(bc%dir)
            boundary_conditions = (/bc_x, bc_y, bc_z/)
            ! FIXME: this is no longer treu with the patch_bc stuff
            beg_end = (/boundary_conditions(bc%dir)%beg, boundary_conditions(bc%dir)%end/)

            ! Implements:
            ! pbc_loc  bc_x >= 0 -> [send/recv]_tag  [dst/src]_proc
            ! -1 (=0)      0            ->     [1,0]       [0,0]      | 0 0 [1,0] [beg,beg]
            ! -1 (=0)      1            ->     [0,0]       [1,0]      | 0 1 [0,0] [end,beg]
            ! +1 (=1)      0            ->     [0,1]       [1,1]      | 1 0 [0,1] [end,end]
            ! +1 (=1)      1            ->     [1,1]       [0,1]      | 1 1 [1,1] [beg,end]

            send_tag = f_logical_to_int(.not. f_xor(bc%type >= 0, bc%loc == 1))
            recv_tag = f_logical_to_int(bc%loc == 1)

            dst_proc = beg_end(1 + f_logical_to_int(f_xor(bc%loc == 1, bc%type >= 0)))
            src_proc = beg_end(1 + f_logical_to_int(bc%loc == 1))

            #:block IMPLEMENT_BOUNDARY_CONDITION()
                !$acc routine seq
                do i = 1, sys_size
                    q_cons_buff_send(packed_idr) = q_cons_vf(i)%sf(nopack_idx, nopack_idy, nopack_idz)
                end do
            #:endblock

            ! Send/Recv
            #:for rdma_mpi in [False, True]
                if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                    p_send => q_cons_buff_send(0)
                    p_recv => q_cons_buff_recv(0)
                    #:if rdma_mpi
                        !$acc data attach(p_send, p_recv)
                        !$acc host_data use_device(p_send, p_recv)
                    #:else
                        !$acc update host(q_cons_buff_send, ib_buff_send)
                    #:endif

                    call MPI_SENDRECV( &
                        p_send, buffer_count, MPI_DOUBLE_PRECISION, dst_proc, send_tag, &
                        p_recv, buffer_count, MPI_DOUBLE_PRECISION, src_proc, recv_tag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    #:if rdma_mpi
                        !$acc end host_data
                        !$acc end data
                        !$acc wait
                    #:else
                        !$acc update device(q_cons_buff_recv)
                    #:endif
                end if
            #:endfor
        end do

        #:block IMPLEMENT_BOUNDARY_CONDITION()
            !$acc routine seq
            do i = 1, sys_size
                q_cons_vf(i)%sf(j + unpack_offset, k, l) = q_cons_buff_recv(r)
            end do
        #:endblock

#endif

    end subroutine s_mpi_sendrecv_variables_buffers
#endif

    subroutine s_prohibit_abort(condition, message)
        character(len=*), intent(in) :: condition, message

        print *, ""
        print *, "===================================================================================================="
        print *, "                                          CASE FILE ERROR                                           "
        print *, "----------------------------------------------------------------------------------------------------"
        print *, "Prohibited condition: ", trim(condition)
        if (len_trim(message) > 0) then
            print *, "Note: ", trim(message)
        end if
        print *, "===================================================================================================="
        print *, ""
        call s_mpi_abort
    end subroutine s_prohibit_abort

end module m_mpi_common
