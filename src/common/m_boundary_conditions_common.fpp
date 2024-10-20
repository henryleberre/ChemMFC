#:include 'macros.fpp'

module m_boundary_conditions_common

    use m_global_parameters

    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_boundary_conditions_module

#ifndef MFC_PRE_PROCESS
    public :: s_populate_prim_buffers

    @:DEFINE_BOUNDARY_CONDITION_INTERFACE(name=xxx)
#endif

    contains

    subroutine s_initialize_boundary_conditions_module()

        call s_inject_global_user_bc_to_patch_bc()

    end subroutine s_initialize_boundary_conditions_module

    subroutine s_inject_global_user_bc_to_patch_bc()

        integer :: i, j, k, l, shift

        shift = 2*num_dims

        do i = num_bc_patches, 1, -1
            patch_bc(i + shift) = patch_bc(i)
        end do

        do concurrent (i = 1:num_dims, j = 1:2)
            k = (i - 1) * 2 + j

            patch_bc(k)%dir = i
            patch_bc(k)%loc = 2*(j - 1) - 1

            if (num_dims <= 2) then
                patch_bc(k)%geometry = 1
            else
                patch_bc(k)%geometry = 2
            end if
        end do

        patch_bc(1)%vel(1) = bc_x%vb1
        patch_bc(1)%vel(2) = bc_x%vb2
        patch_bc(1)%vel(3) = bc_x%vb3
        patch_bc(2)%vel(1) = bc_x%ve1
        patch_bc(2)%vel(2) = bc_x%ve2
        patch_bc(2)%vel(3) = bc_x%ve3
        patch_bc(3)%vel(1) = bc_y%vb1
        patch_bc(3)%vel(2) = bc_y%vb2
        patch_bc(3)%vel(3) = bc_y%vb3
        patch_bc(4)%vel(1) = bc_y%ve1
        patch_bc(4)%vel(2) = bc_y%ve2
        patch_bc(4)%vel(3) = bc_y%ve3

        patch_bc(1)%type = bc_x%beg; patch_bc(2)%type = bc_x%end
        if (n > 0) then
            patch_bc(3)%type = bc_y%beg; patch_bc(4)%type = bc_y%end

            if (p > 0) then
                patch_bc(5)%type = bc_z%beg; patch_bc(6)%type = bc_z%end
            end if
        end if

        num_bc_patches = num_bc_patches + shift

    end subroutine s_inject_global_user_bc_to_patch_bc

#ifndef MFC_PRE_PROCESS
    !>  The purpose of this procedure is to populate the buffers
    !!      of the primitive variables, depending on the selected
    !!      boundary conditions.
    subroutine s_populate_prim_buffers(q_prim_vf,&
#ifdef MFC_SIMULATION
pb, mv &
#endif
)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_SIMULATION
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif

        procedure(s_apply_xxx_boundary_condition), &
            pointer :: s_apply_boundary_condition => null() 

        integer :: i

        do i = 1, num_bc_patches
            s_apply_boundary_condition => null()

            select case (patch_bc(i)%type)
            case (-13:-3); s_apply_boundary_condition => s_apply_ghost_cell_extrapolation_boundary_condition
            case (-14);    s_apply_boundary_condition => s_apply_axis_boundary_condition
            case (-2);     s_apply_boundary_condition => s_apply_symmetry_boundary_condition
            case (-1);     s_apply_boundary_condition => s_apply_periodic_boundary_condition
            case (-15);    s_apply_boundary_condition => s_apply_slip_wall_boundary_condition
            case (-16);    s_apply_boundary_condition => s_apply_no_slip_wall_boundary_condition
            case default ! Processor BC at beginning
                call s_mpi_sendrecv_variables_buffers( &
                    q_prim_vf, pb, mv, patch_bc(i)%dir, patch_bc(i)%loc)
            end select

            if (associated(s_apply_boundary_condition)) then
                call s_apply_boundary_condition(q_prim_vf, &
#ifdef MFC_SIMULATION
                pb, mv, &
#endif
                patch_bc(i))
            end if

#ifdef MFC_SIMULATION
            if (qbmm .and. .not. polytropic) then
                select case (patch_bc(i)%type)
                case (-13:-3, -15, -16) ! Ghost-cell extrap. BC at beginning
                    call s_apply_qbmm_extrapolation_boundary_condition(q_prim_vf, pb, mv, patch_bc(i))
                end select
            end if
#endif
        end do

    end subroutine s_populate_prim_buffers
#endif

#ifndef MFC_PRE_PROCESS
    #:block DEFINE_BOUNDARY_CONDITION(name="ghost_cell_extrapolation")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
        #:endblock
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="periodic")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(px, py, pz)
        #:endblock

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
                pb(x, y, z, q, i) = pb(px, py, pz, q, i)
                mv(x, y, z, q, i) = mv(px, py, pz, q, i)
            #:endblock
        end if
#endif
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="symmetry")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "momxb + bc%dir - 2")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(sx, sy, sz)
        #:endblock

        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[])
            q_prim_vf(momxb + bc%dir - 1)%sf(x, y, z) = &
                -q_prim_vf(momxb + bc%dir - 1)%sf(sx, sy, sz)
        #:endblock

        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", "momxb + bc%dir", "sys_size")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(sx, sy, sz)
        #:endblock

#ifdef MFC_SIMULATION
        if (qbmm .and. .not. polytropic) then
            #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
                pb(x, y, z, q, i) = pb(sx, sy, sz, q, i)
                mv(x, y, z, q, i) = mv(sx, sy, sz, q, i)
            #:endblock
        end if
#endif
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="axis")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[])
            if (z_cc(k) < pi) then
                !$acc loop seq
                do i = 1, momxb
                    q_prim_vf(i)%sf(x, y, z) = &
                        q_prim_vf(i)%sf(sx, sy, z + ((p + 1)/2))
                end do

                q_prim_vf(momxb + 1)%sf(x, y, z) = &
                    -q_prim_vf(momxb + 1)%sf(sx, sy, z + ((p + 1)/2))

                q_prim_vf(momxe)%sf(x, y, z) = &
                    -q_prim_vf(momxe)%sf(sx, sy, z + ((p + 1)/2))

                !$acc loop seq
                do i = E_idx, sys_size
                    q_prim_vf(i)%sf(x, y, z) = &
                        q_prim_vf(i)%sf(sx, sy, z + ((p + 1)/2))
                end do
            else
                !$acc loop seq
                do i = 1, momxb
                    q_prim_vf(i)%sf(x, y, z) = &
                        q_prim_vf(i)%sf(sx, sy, z - ((p + 1)/2))
                end do

                q_prim_vf(momxb + 1)%sf(x, y, z) = &
                    -q_prim_vf(momxb + 1)%sf(sx, sy, z - ((p + 1)/2))

                q_prim_vf(momxe)%sf(x, y, z) = &
                    -q_prim_vf(momxe)%sf(sx, sy, z - ((p + 1)/2))

                !$acc loop seq
                do i = E_idx, sys_size
                    q_prim_vf(i)%sf(x, y, z) = &
                        q_prim_vf(i)%sf(sx, sy, z - ((p + 1)/2))
                end do
            end if
        #:endblock

#ifdef MFC_SIMULATION
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
            pb(x, y, z, q, i) = pb(sx, sy, z - ((p + 1)/2), q, i)
            mv(x, y, z, q, i) = mv(sx, sy, z - ((p + 1)/2), q, i)
        #:endblock
#endif
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="slip_wall")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            if (i == momxb + bc%dir - 1) then
                q_prim_vf(i)%sf(x, y, z) = &
                    -q_prim_vf(i)%sf(sx, sy, sz) + 2d0 * bc%vel(bc%dir)
            else
                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
            end if
        #:endblock
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="no_slip_wall")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            if (i >= momxb .and. i <= momxe) then
                q_prim_vf(i)%sf(x, y, z) = &
                    -q_prim_vf(i)%sf(sx, sy, sz) + 2d0 * bc%vel(i - momxb + 1)
            else
                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
            end if
        #:endblock
    #:endblock

#ifdef MFC_SIMULATION
    #:block DEFINE_BOUNDARY_CONDITION(name="qbmm_extrapolation")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
            pb(x, y, z, q, i) = pb(ex, ey, ez, q, i)
            mv(x, y, z, q, i) = mv(ex, ey, ez, q, i)
        #:endblock
    #:endblock
#endif
#endif

end module m_boundary_conditions_common