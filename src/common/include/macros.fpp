#:def LOG(expr)
#ifdef MFC_DEBUG
    block
        use iso_fortran_env, only: output_unit

        print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$

        call flush (output_unit)
    end block
#endif
#:enddef

#:def ALLOCATE(*args)
    @:LOG({'@:ALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    allocate (${', '.join(args)}$)
#ifndef CRAY_ACC_WAR
!$acc enter data create(${', '.join(args)}$)
#endif
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    @:LOG({'@:DEALLOCATE(${re.sub(' +', ' ', ', '.join(args))}$)'})
    deallocate (${', '.join(args)}$)
#ifndef CRAY_ACC_WAR
!$acc exit data delete(${', '.join(args)}$)
#endif
#:enddef DEALLOCATE

#:def ALLOCATE_GLOBAL(*args)
    @:LOG({'@:ALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef CRAY_ACC_WAR
    allocate (${', '.join(('p_' + arg.strip() for arg in args))}$)
    #:for arg in args
        ${re.sub('\\(.*\\)','',arg)}$ => ${ 'p_' + re.sub('\\(.*\\)','',arg.strip()) }$
    #:endfor
    !$acc enter data create(${', '.join(('p_' + re.sub('\\(.*\\)','',arg.strip()) for arg in args))}$) &
    !$acc& attach(${', '.join(map(lambda x: re.sub('\\(.*\\)','',x), args))}$)
#else
    allocate (${', '.join(args)}$)
    !$acc enter data create(${', '.join(args)}$)
#endif

#:enddef ALLOCATE_GLOBAL

#:def DEALLOCATE_GLOBAL(*args)
    @:LOG({'@:DEALLOCATE_GLOBAL(${re.sub(' +', ' ', ', '.join(args))}$)'})
#ifdef CRAY_ACC_WAR
    !$acc exit data delete(${', '.join(('p_' + arg.strip() for arg in args))}$) &
    !$acc& detach(${', '.join(args)}$)
    #:for arg in args
        nullify (${arg}$)
    #:endfor
    deallocate (${', '.join(('p_' + arg.strip() for arg in args))}$)
#else
    deallocate (${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#endif

#:enddef DEALLOCATE_GLOBAL

#:def CRAY_DECLARE_GLOBAL(intype, dim, *args)
#ifdef CRAY_ACC_WAR
    ${intype}$, ${dim}$, allocatable, target :: ${', '.join(('p_' + arg.strip() for arg in args))}$
    ${intype}$, ${dim}$, pointer :: ${', '.join(args)}$
#else
    ${intype}$, ${dim}$, allocatable :: ${', '.join(args)}$
#endif
#:enddef CRAY_DECLARE_GLOBAL

#:def CRAY_DECLARE_GLOBAL_SCALAR(intype, *args)
#ifdef CRAY_ACC_WAR
    ${intype}$, target :: ${', '.join(('p_' + arg.strip() for arg in args))}$
    ${intype}$, pointer :: ${', '.join(args)}$
#else
    ${intype}$::${', '.join(args)}$
#endif
#:enddef CRAY_DECLARE_GLOBAL_SCALAR

#:def ACC_SETUP_VFs(*args)
#ifdef CRAY_ACC_WAR
    block
        integer :: macros_setup_vfs_i

        @:LOG({'@:ACC_SETUP_VFs(${', '.join(args)}$)'})

        #:for arg in args
            !$acc enter data copyin(${arg}$)
            !$acc enter data copyin(${arg}$%vf)
            if (allocated(${arg}$%vf)) then
                do macros_setup_vfs_i = lbound(${arg}$%vf, 1), ubound(${arg}$%vf, 1)
                    if (associated(${arg}$%vf(macros_setup_vfs_i)%sf)) then
                        !$acc enter data copyin(${arg}$%vf(macros_setup_vfs_i))
                        !$acc enter data create(${arg}$%vf(macros_setup_vfs_i)%sf)
                    end if
                end do
            end if
        #:endfor
    end block
#endif
#:enddef

#:def ACC_SETUP_SFs(*args)
#ifdef CRAY_ACC_WAR
    block

        @:LOG({'@:ACC_SETUP_SFs(${', '.join(args)}$)'})

        #:for arg in args
            !$acc enter data copyin(${arg}$)
            if (associated(${arg}$%sf)) then
                !$acc enter data create(${arg}$%sf)
            end if
        #:endfor
    end block
#endif
#:enddef

#:def ACC_SETUP_source_spatials(*args)
#ifdef CRAY_ACC_WAR
    block

        @:LOG({'@:ACC_SETUP_source_spatials(${', '.join(args)}$)'})

        #:for arg in args
            !$acc enter data copyin(${arg}$)
            if (allocated(${arg}$%coord)) then
                !$acc enter data create(${arg}$%coord)
            end if
            if (allocated(${arg}$%val)) then
                !$acc enter data create(${arg}$%val)
            end if
            if (allocated(${arg}$%angle)) then
                !$acc enter data create(${arg}$%angle)
            end if
            if (allocated(${arg}$%xyz_to_r_ratios)) then
                !$acc enter data create(${arg}$%xyz_to_r_ratios)
            end if
        #:endfor
    end block
#endif
#:enddef

#:def PROHIBIT(condition, message = None)
    if (${condition}$) then
        call s_prohibit_abort("${condition}$", ${message or '""'}$)
    end if
#:enddef

#define t_vec3   real(kind(0d0)), dimension(1:3)
#define t_mat4x4 real(kind(0d0)), dimension(1:4,1:4)

#:def ASSERT(predicate, message = None)
    if (.not. (${predicate}$)) then
        call s_mpi_abort("${_FILE_.split('/')[-1]}$:${_LINE_}$: "// &
                         "Assertion failed: ${predicate}$. " &
                         //${message or '"No error description."'}$)
    end if
#:enddef


#:def _WRAP_i_LOOP(impl, loops)
    !$acc parallel loop collapse(${len(loops)}$) gang vector default(present) private(x, y, z, sx, sy, sz, px, py, pz, ex, ey, ez, pack_idr, internal_pack_offset, pack_idx, pack_idy, pack_idz, internal_unpack_offset)
    #:for index, lbound, hbound in (loops or [])
        do ${index}$ = ${lbound}$, ${hbound}$
    #:endfor

    $:impl

    #:for i in range(len(loops or []))
        end do
    #:endfor
#:enddef

#:def IMPLEMENT_BOUNDARY_CONDITION(impl, outer_loops = None, inner_loops = None)
    block
        integer, dimension(1:3) :: grid_dims
        grid_dims = (/m, n, p/)

        internal_pack_offset = 0
        if (f_xor(bc%loc == 1, bc%type >= 0)) then
            internal_pack_offset = grid_dims(bc%dir) - buff_size + 1
        end if

        internal_unpack_offset = 0
        if (bc%loc == 1) then
            internal_unpack_offset = grid_dims(bc%dir) + buff_size + 1
        end if
    end block

    if (bc%dir == 1) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("l", 0, "p"), ("k", 0, "n"), ("j", 1, "buff_size")] + (inner_loops or []))
                 x = -j;           y = k;  z = l ! Regular
                sx = j - 1;       sy = k; sz = l ! Symmetry
                px = m - (j - 1); py = k; pz = l ! Periodic
                ex = 0;           ey = k; ez = l ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((j - 1) + buff_size*( k +      (n + 1)*l))
                unpack_idr = (i - 1) + v_size * (-j      + buff_size*((k + 1) + (n + 1)*l))
                  pack_idx =  j - 1 +   internal_pack_offset;   pack_idy = k;   pack_idz = l
                unpack_idx = -j     + internal_unpack_offset; unpack_idy = k; unpack_idz = l

                $:impl
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("l", 0, "p"), ("k", 0, "n"), ("j", 1, "buff_size")] + (inner_loops or []))
                 x = m + j;        y = k;  z = l ! Regular
                sx = m - (j - 1); sy = k; sz = l ! Symmetry
                px = j - 1;       py = k; pz = l ! Periodic
                ex = m;           ey = k; ez = l ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((j - 1) + buff_size*( k +      (n + 1)*l))
                unpack_idr = (i - 1) + v_size * (-j      + buff_size*((k + 1) + (n + 1)*l))
                  pack_idx =  j - 1 +   internal_pack_offset;   pack_idy = k;   pack_idz = l
                unpack_idx = -j     + internal_unpack_offset; unpack_idy = k; unpack_idz = l

                $:impl
            #:endblock

        end if

        !< y-direction =========================================================
    elseif (bc%dir == 2) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("k", 0, "p"), ("j", 1, "buff_size"), ("l", "-buff_size", "m + buff_size")] + (inner_loops or []))
                 x = l;  y = -j;           z = k ! Regular
                sx = l; sy = j - 1;       sz = k ! Symmetry
                px = l; py = n - (j - 1); pz = k ! Periodic
                ex = l; ey = 0;           ez = k ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((l + buff_size) + (m + 2*buff_size + 1) * ((j - 1)          + buff_size*k))
                unpack_idr = (i - 1) + v_size * ((l + buff_size) + (m + 2*buff_size + 1) * ((-j + buff_size) + buff_size*k))
                  pack_idx = l;   pack_idy = (j - 1) + internal_pack_offset;   pack_idz = k
                unpack_idx = l; unpack_idy = -j + internal_unpack_offset;    unpack_idz = k

                $:impl
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("k", 0, "p"), ("j", 1, "buff_size"), ("l", "-buff_size", "m + buff_size")] + (inner_loops or []))
                 x = l; y  = n + j;        z = k ! Regular
                sx = l; sy = n - (j - 1); sz = k ! Symmetry
                px = l; py = j - 1;       pz = k ! Periodic
                ex = l; ey = n;           ez = k ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((l + buff_size) + (m + 2*buff_size + 1) * ((j - 1)          + buff_size*k))
                unpack_idr = (i - 1) + v_size * ((l + buff_size) + (m + 2*buff_size + 1) * ((-j + buff_size) + buff_size*k))
                  pack_idx = l;   pack_idy = (j - 1) + internal_pack_offset;   pack_idz = k
                unpack_idx = l; unpack_idy = -j + internal_unpack_offset;    unpack_idz = k

                $:impl
            #:endblock

        end if

        !< z-direction =========================================================
    elseif (bc%dir == 3) then

        if (bc%loc == -1) then

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("j", 1, "buff_size"), ("l", "-buff_size", "n + buff_size"), ("k", "-buff_size", "m + buff_size")] + (inner_loops or []))
                 x = k;  y = l;  z = -j          ! Regular
                sx = k; sy = l; sz = j - 1       ! Symmetry
                px = k; py = l; pz = p - (j - 1) ! Periodic
                ex = k; ey = l; ez = 0           ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((k + buff_size) + (m + 2*buff_size + 1) * ((l + buff_size) + (n + 2*buff_size + 1)*(j - 1)))
                unpack_idr = (i - 1) + v_size * ((k + buff_size) + (m + 2*buff_size + 1) * ((l + buff_size) + (n + 2*buff_size + 1)*(-j + buff_size)))
                  pack_idx = k;   pack_idy = l;   pack_idz = (j - 1) + internal_pack_offset
                unpack_idx = k; unpack_idy = l; unpack_idz = -j + internal_unpack_offset

                $:impl
            #:endblock

        else

            #:block _WRAP_i_LOOP(loops = (outer_loops or []) + [("j", 1, "buff_size"), ("l", "-buff_size", "n + buff_size"), ("k", "-buff_size", "m + buff_size")] + (inner_loops or []))
                 x = k;  y = l;  z = p + j       ! Regular
                sx = k; sy = l; sz = p - (j - 1) ! Symmetry
                px = k; py = l; pz = j - 1       ! Periodic
                ex = k; ey = l; ez = p           ! Extrapolation

                  pack_idr = (i - 1) + v_size * ((k + buff_size) + (m + 2*buff_size + 1) * ((l + buff_size) + (n + 2*buff_size + 1)*(j - 1)))
                unpack_idr = (i - 1) + v_size * ((k + buff_size) + (m + 2*buff_size + 1) * ((l + buff_size) + (n + 2*buff_size + 1)*(-j + buff_size)))
                  pack_idx = k;   pack_idy = l;   pack_idz = (j - 1) + internal_pack_offset
                unpack_idx = k; unpack_idy = l; unpack_idz = -j + internal_unpack_offset

                $:impl
            #:endblock

        end if

    end if
#:enddef IMPLEMENT_BOUNDARY_CONDITION

#:def BOUNDARY_CONDITION_INTENT_DECLARATIONS()
    type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_SIMULATION
    real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif
    type(bc_patch_parameters), intent(in) :: bc
#:enddef

#:def BOUNDARY_CONDITION_INTEGER_DECLARATIONS()
    integer :: i, j, k, l, q

    integer ::  x,  y,  z
    integer :: sx, sy, sz
    integer :: px, py, pz
    integer :: ex, ey, ez
    integer :: pack_idr, unpack_idr
    integer :: pack_idx, pack_idy, pack_idz
    integer :: unpack_idx, unpack_idy, unpack_idz
    integer :: internal_pack_offset, internal_unpack_offset
#:enddef

#:def DEFINE_BOUNDARY_CONDITION_ROUTINE_STRUCT(impl, name)
subroutine s_apply_${name}$_boundary_condition(q_prim_vf, &
#ifdef MFC_SIMULATION
pb, mv, &
#endif
bc)

    $:impl

end subroutine s_apply_${name}$_boundary_condition
#:enddef DEFINE_BOUNDARY_CONDITION_ROUTINE_STRUCT

#:def DEFINE_BOUNDARY_CONDITION(impl, name)
#:block DEFINE_BOUNDARY_CONDITION_ROUTINE_STRUCT(name=name)

    @:BOUNDARY_CONDITION_INTENT_DECLARATIONS()
    @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

    $:impl

#:endblock
#:enddef DEFINE_BOUNDARY_CONDITION

#:def DEFINE_BOUNDARY_CONDITION_INTERFACE(name)
abstract interface
    #:block DEFINE_BOUNDARY_CONDITION_ROUTINE_STRUCT(name=name)

        import :: scalar_field, sys_size, bc_patch_parameters

#ifdef MFC_SIMULATION
        import :: startx, starty, startz,
#endif

        @:BOUNDARY_CONDITION_INTENT_DECLARATIONS()
        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()
        
    #:endblock
end interface
#:enddef DEFINE_BOUNDARY_CONDITION_INTERFACE
