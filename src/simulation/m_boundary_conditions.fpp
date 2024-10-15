!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

#:include 'macros.fpp'

#:def _WRAP_i_LOOP(impl, loops = None)
    !$acc parallel loop collapse(${3 + len(loops or [])}$) gang vector default(present) private(x, y, z, sx, sy, sz, px, py, pz, ex, ey, ez)
    #:for index, lbound, hbound in (loops or [])
        do ${index}$ = ${lbound}$, ${hbound}$
    #:endfor

    $:impl

    #:for i in range(len(loops or []))
        end do
    #:endfor
#:enddef

#:def IMPLEMENT_BOUNDARY_CONDITION(impl, loops = None)
    if (bc%dir == 1) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = loops)
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                             x = -j;           y = k;  z = l ! Regular
                            sx = j - 1;       sy = k; sz = l ! Symmetry
                            px = m - (j - 1); py = k; pz = l ! Periodic
                            ex = 0;           ey = k; ez = l ! Extrapolation

                            $:impl
                        end do
                    end do
                end do
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = loops)
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size

                             x = m + j;        y = k;  z = l ! Regular
                            sx = m - (j - 1); sy = k; sz = l ! Symmetry
                            px = j - 1;       py = k; pz = l ! Periodic
                            ex = m;           ey = k; ez = l ! Extrapolation

                            $:impl

                        end do
                    end do
                end do
            #:endblock

        end if

        !< y-direction =========================================================
    elseif (bc%dir == 2) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = loops)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                             x = l;  y = -j;           z = k ! Regular
                            sx = l; sy = j - 1;       sz = k ! Symmetry
                            px = l; py = n - (j - 1); pz = k ! Periodic
                            ex = l; ey = 0;           ez = k ! Extrapolation

                            $:impl
                        end do
                    end do
                end do
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = loops)
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                             x = l; y  = n + j;        z = k ! Regular
                            sx = l; sy = n - (j - 1); sz = k ! Symmetry
                            px = l; py = j - 1;       pz = k ! Periodic
                            ex = l; ey = n;           ez = k ! Extrapolation

                            $:impl
                        end do
                    end do
                end do
            #:endblock

        end if

        !< z-direction =========================================================
    elseif (bc%dir == 3) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = loops)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                             x = k;  y = l;  z = -j          ! Regular
                            sx = k; sy = l; sz = j - 1       ! Symmetry
                            px = k; py = l; pz = p - (j - 1) ! Periodic
                            ex = k; ey = l; ez = 0           ! Extrapolation

                            $:impl
                        end do
                    end do
                end do
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = loops)
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                             x = k;  y = l;  z = p + j       ! Regular
                            sx = k; sy = l; sz = p - (j - 1) ! Symmetry
                            px = k; py = l; pz = j - 1       ! Periodic
                            ex = k; ey = l; ez = p           ! Extrapolation

                            $:impl
                        end do
                    end do
                end do
            #:endblock

        end if

    end if
#:enddef IMPLEMENT_BOUNDARY_CONDITION

#:def BOUNDARY_CONDITION_INTEGER_DECLARATIONS()
    integer :: i, j, k, l, q

    integer ::  x,  y,  z
    integer :: sx, sy, sz
    integer :: px, py, pz
    integer :: ex, ey, ez
#:enddef

#:def DEFINE_BOUNDARY_CONDITION(impl, name)
subroutine s_${name}$(q_prim_vf, pb, mv, bc)

    type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
    real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
    type(bc_patch_parameters), intent(in) :: bc

    @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

    $:impl

end subroutine s_${name}$
#:enddef

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons
module m_boundary_conditions

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy

    use m_constants

    use m_boundary_conditions_common

    ! ==========================================================================

    implicit none

    private; 
    public :: s_populate_variables_buffers, &
              s_populate_capillary_buffers, &
              s_initialize_boundary_conditions_module

contains

    !>  The purpose of this procedure is to populate the buffers
    !!      of the primitive variables, depending on the selected
    !!      boundary conditions.
    subroutine s_populate_variables_buffers(q_prim_vf, pb, mv)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv

        integer :: i

        do i = 1, num_bc_patches
            select case (patch_bc(i)%type)
            case (-13:-3) ! Ghost-cell extrap. BC at beginning
                call s_ghost_cell_extrapolation(q_prim_vf, pb, mv, patch_bc(i))
            case (-14)    ! Axis BC at beginning
                call s_axis(q_prim_vf, pb, mv, patch_bc(i))
            case (-2)     ! Symmetry BC at beginning
                call s_symmetry(q_prim_vf, pb, mv, patch_bc(i))
            case (-1)     ! Periodic BC at beginning
                call s_periodic(q_prim_vf, pb, mv, patch_bc(i))
            case (-15)    ! Slip wall BC at beginning
                call s_slip_wall(q_prim_vf, pb, mv, patch_bc(i))
            case (-16)    ! No-slip wall BC at beginning
                call s_no_slip_wall(q_prim_vf, pb, mv, patch_bc(i))
            case default ! Processor BC at beginning
                call s_mpi_sendrecv_variables_buffers( &
                    q_prim_vf, pb, mv, patch_bc(i)%dir, patch_bc(i)%loc)
            end select

            if (qbmm .and. .not. polytropic) then
                select case (patch_bc(i)%type)
                case (-13:-3) ! Ghost-cell extrap. BC at beginning
                    call s_qbmm_extrapolation(q_prim_vf, pb, mv, patch_bc(i))
                case (-15)    ! Slip wall BC at beginning
                    call s_qbmm_extrapolation(q_prim_vf, pb, mv, patch_bc(i))
                case (-16)    ! No-slip wall BC at beginning
                    call s_qbmm_extrapolation(q_prim_vf, pb, mv, patch_bc(i))
                end select
            end if
        end do

    end subroutine s_populate_variables_buffers

    #:block DEFINE_BOUNDARY_CONDITION(name="ghost_cell_extrapolation")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
        #:endblock
    #:endblock

    #:block DEFINE_BOUNDARY_CONDITION(name="periodic")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "sys_size")])
            q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(px, py, pz)
        #:endblock

        if (qbmm .and. .not. polytropic) then
            #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
                pb(x, y, z, q, i) = pb(px, py, pz, q, i)
                mv(x, y, z, q, i) = mv(px, py, pz, q, i)
            #:endblock
        end if
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

        if (qbmm .and. .not. polytropic) then
            #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
                pb(x, y, z, q, i) = pb(sx, sy, sz, q, i)
                mv(x, y, z, q, i) = mv(sx, sy, sz, q, i)
            #:endblock
        end if
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

        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
            pb(x, y, z, q, i) = pb(sx, sy, z - ((p + 1)/2), q, i)
            mv(x, y, z, q, i) = mv(sx, sy, z - ((p + 1)/2), q, i)
        #:endblock
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

    #:block DEFINE_BOUNDARY_CONDITION(name="qbmm_extrapolation")
        #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "nb"), ("q", 1, "nnode")])
            pb(x, y, z, q, i) = pb(ex, ey, ez, q, i)
            mv(x, y, z, q, i) = mv(ex, ey, ez, q, i)
        #:endblock
    #:endblock

    subroutine s_populate_capillary_buffers(c_divs)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(bc_patch_parameters) :: bc        

        integer :: r

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        do r = 1, num_bc_patches
            bc = patch_bc(r)

            if (bc%type <= -3) then !< ghost cell extrapolation
                #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "num_dims + 1")])
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(ex, ey, ez)
                #:endblock
            elseif (bc%type == -2) then !< slip wall or reflective
                #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "num_dims + 1")])
                    if (i == 1) then
                        c_divs(i)%sf(x, y, z) = -c_divs(i)%sf(sx, sy, sz)
                    else
                        c_divs(i)%sf(x, y, z) =  c_divs(i)%sf(sx, sy, sz)
                    end if
                #:endblock
            elseif (bc%type == -1) then
                #:block IMPLEMENT_BOUNDARY_CONDITION(loops=[("i", 1, "num_dims + 1")])
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(px, py, pz)
                #:endblock
            else
                call s_mpi_sendrecv_capilary_variables_buffers(c_divs, bc%dir, bc%loc)
            end if
        end do

    end subroutine s_populate_capillary_buffers

end module m_boundary_conditions
