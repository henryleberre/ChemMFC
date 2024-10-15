module m_boundary_conditions_common

    use m_global_parameters

    implicit none

    private; public :: s_initialize_boundary_conditions_module

    contains

    subroutine s_initialize_boundary_conditions_module()

        call s_inject_global_user_bc_to_patch_bc()
        call s_populate_grid_extents()

    end subroutine s_initialize_boundary_conditions_module

    subroutine s_inject_global_user_bc_to_patch_bc()

        integer :: i, j, k, l, shift

        shift = 2*num_dims

        do i = num_bc_patches, 1, -1
            patch_bc(i + shift) = patch_bc(i)
        end do

        do i = 1, num_dims
            do j = 1, 2
                k = (i - 1)*2 + j

                patch_bc(k)%loc = i
                patch_bc(k)%dir = 2*(j - 1) - 1

                if (num_dims <= 2) then
                    patch_bc(k)%geometry = 1
                else
                    patch_bc(k)%geometry = 2
                end if

                #:for i, coord in [(1, 'x'), (2, 'y'), (3, 'z')]
                #:for j in [1, 2]
                    if (${i}$ == i .and. ${j}$ == j) then
                        patch_bc(k)%type            = bc_${coord}$%${"end" if j == 2 else "beg"}$
                        patch_bc(k)%centroid(i)     = ${coord}$_domain%${"end" if j == 2 else "beg"}$
                        patch_bc(k)%length(i)       = abs(${coord}$_domain%end - ${coord}$_domain%beg)
                        patch_bc(k)%grid_extents(i) = (/ 0, m /)

                        #:for coord2 in set(['x', 'y', 'z']) - set([coord])
                            l = ${{'x': 1, 'y': 2, 'z': 3}[coord2]}$
                            if (l <= num_dims) then
                                patch_bc(k)%centroid(l) = (${coord2}$_domain%end + ${coord2}$_domain%beg) / 2
                                patch_bc(k)%length  (l) = (${coord2}$_domain%end - ${coord2}$_domain%beg)
                            end if
                        #:endfor
                    end if
                #:endfor
                #:endfor
            end do
        end do

        num_bc_patches = num_bc_patches + shift

        do i = 1, num_bc_patches
            print*, "i", i
            print*, patch_bc(i)%type
            print*, patch_bc(i)%geometry
            print*, patch_bc(i)%centroid(:)
            print*, patch_bc(i)%length(:)
            print*, patch_bc(i)%radius
            print*, patch_bc(i)%dir
            print*, patch_bc(i)%loc
            print*, ""
        end do

    end subroutine s_inject_global_user_bc_to_patch_bc

    function f_binary_search_closest(real_array, value) result(closest_idx)

        real(kind(0d0)), dimension(:), intent(in) :: real_array
        real(kind(0d0)), intent(in) :: value
        integer :: closest_idx, low, high, mid

        low = 1
        high = size(real_array)
        mid = (low + high) / 2

        if (value <= real_array(low)) then
            closest_idx = low
            return
        end if

        if (value >= real_array(high)) then
            closest_idx = high
            return
        end if

        do while (low < high)
            if (value < real_array(mid)) then
                high = mid
            else
                low = mid
            end if

            mid = (low + high) / 2
        end do

        if (abs(value - real_array(low)) < abs(value - real_array(high))) then
            closest_idx = low
        else
            closest_idx = high
        end if

    end function f_binary_search_closest

    subroutine s_populate_grid_extents()

        integer :: i, j, k, l

        real(kind(0d0)), dimension(3) :: corners(2)

        do i = 1, 2*num_dims
            patch_bc(i)%grid_extents 
        end do

        do i = 2*num_dims + 1, num_bc_patches
            corners = (/ &
                patch_bc(i)%centroid - patch_bc(i)%length / 2, &
                patch_bc(i)%centroid + patch_bc(i)%length / 2  &
            /)

            if (patch_bc(i)%geometry == 1) then
                patch_bc(i)%grid_extents(1) = (/ 0, m /)
            else if (patch_bc(i)%geometry == 2) then
                if (patch_bc(i)%loc == 1) then
                    patch_bc(i)%grid_extents(1) = (/ /)
                end if
            end if
        end do

    end subroutine s_populate_grid_extents

end module m_boundary_conditions_common