!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

#:include 'macros.fpp'

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons
module m_boundary_conditions

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Definitions of the MPI proxy

    use m_constants

    use m_boundary_conditions_common

    ! ==========================================================================

    implicit none

    private; 
    public :: s_populate_prim_buffers, &
              s_populate_capillary_buffers, &
              s_initialize_boundary_conditions_module

contains

    subroutine s_populate_capillary_buffers(c_divs)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(bc_patch_parameters) :: bc        

        integer :: r

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        do r = 1, num_bc_patches
            bc = patch_bc(r)

            if (bc%type <= -3) then !< ghost cell extrapolation
                #:block IMPLEMENT_BOUNDARY_CONDITION(outer_loops=[("i", 1, "num_dims + 1")])
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(ex, ey, ez)
                #:endblock
            elseif (bc%type == -2) then !< slip wall or reflective
                #:block IMPLEMENT_BOUNDARY_CONDITION(outer_loops=[("i", 1, "num_dims + 1")])
                    if (i == 1) then
                        c_divs(i)%sf(x, y, z) = -c_divs(i)%sf(sx, sy, sz)
                    else
                        c_divs(i)%sf(x, y, z) =  c_divs(i)%sf(sx, sy, sz)
                    end if
                #:endblock
            elseif (bc%type == -1) then
                #:block IMPLEMENT_BOUNDARY_CONDITION(outer_loops=[("i", 1, "num_dims + 1")])
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(px, py, pz)
                #:endblock
            else
                call s_mpi_sendrecv_capilary_variables_buffers(c_divs, bc%dir, bc%loc)
            end if
        end do

    end subroutine s_populate_capillary_buffers

end module m_boundary_conditions
