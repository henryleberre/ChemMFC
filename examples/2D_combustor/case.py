#!/usr/bin/env python3

from mfc.case_utils import *

import json, argparse
import cantera as ct

# Reference:
# + https://doi.org/10.1006/jcph.1996.5622: 2-D Combustor Simulation

parser = argparse.ArgumentParser(
    prog="2D_combustor",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc",     type=str,    default='{}', metavar="DICT")
parser.add_argument("--no-chem", dest='chem', default=True, action="store_false",
                                 help="Disable chemistry.")

args = parser.parse_args()

ctfile     = 'h2o2.yaml'
sol_bg     = ct.Solution(ctfile)
sol_bg.TPX =   700,  36_100, 'O2:2,AR:7'
sol_in     = ct.Solution(ctfile)
sol_in.TPX = 1_166, 121_000, 'H2:4,AR:7'

in_length = 0.4375e-2
in_vely   = 10.0

Nx   = 25
Tend = 1e-4

Lx, Ly = 4e-2, 3e-2
Nx, Ny = 64,   48
dt     = 10e-6/100

NT         = int(Tend / dt)
SAVE_COUNT = 200
NS         = NT // SAVE_COUNT

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0,
    'x_domain%end'                 : Lx,
    'y_domain%beg'                 : 0,
    'y_domain%end'                 : Ly,
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : NS,
    'parallel_io'                  : 'T',

    # Simulation Algorithm Parameters ==========================================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 2,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1E-16,
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 1,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     :-2,
    'bc_x%end'                     :-2,
    'bc_y%beg'                     :-7,
    'bc_y%end'                     :-8,

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    # Motionless Background Patch ==============================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : Lx/2,
    'patch_icpp(1)%y_centroid'     : Ly/2,
    'patch_icpp(1)%length_x'       : Lx,
    'patch_icpp(1)%length_y'       : Ly,
    'patch_icpp(1)%vel(1)'         : 0,
    'patch_icpp(1)%vel(2)'         : 0,
    'patch_icpp(1)%pres'           : sol_bg.P,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : sol_bg.density,
    # ==========================================================================

    # Inlet Patch ==============================================================
    'patch_icpp(2)%geometry'       : 3,
    'patch_icpp(2)%x_centroid'     : Lx / 2,
    'patch_icpp(2)%y_centroid'     : 10 * (Ly / Ny) / 2,
    'patch_icpp(2)%length_x'       : in_length,
    'patch_icpp(2)%length_y'       : 10 * Ly / Ny,
    'patch_icpp(2)%vel(1)'         : 0,
    'patch_icpp(2)%vel(2)'         : in_vely,
    'patch_icpp(2)%pres'           : sol_in.P,
    'patch_icpp(2)%alpha(1)'       : 1,
    'patch_icpp(2)%alpha_rho(1)'   : sol_in.density,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================
}

if args.chem:
    case.update({
        # Chemistry ============================================================
        'chem_wrt_T'            : 'T',
        'chemistry'             : 'T',
        'chem_params%advection' : 'T',
        'chem_params%diffusion' : 'F',
        'chem_params%reactions' : 'T',
        'cantera_file'          : ctfile,
        # ======================================================================
    })

    for i in range(len(sol_bg.Y)):
        case[f'chem_wrt_Y({i+1})'] = 'T'

    for i in range(len(sol_bg.Y)):
        case[f'patch_icpp(1)%Y({i+1})'] = sol_bg.Y[i]
    
    for i in range(len(sol_in.Y)):
        case[f'patch_icpp(2)%Y({i+1})'] = sol_in.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))
