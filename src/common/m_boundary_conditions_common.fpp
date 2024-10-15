module m_boundary_conditions_common

    use m_global_parameters

    implicit none

    private; public :: s_initialize_boundary_conditions_module

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

end module m_boundary_conditions_common