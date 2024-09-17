#!/usr/bin/env python3

# https://www.sciencedirect.com/science/article/pii/S0045793013003976?via=ihub
# 4.6. Perfectly stirred reactor

import re, json, argparse
import cantera as ct

parser = argparse.ArgumentParser(
    prog="nD_Reactor",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
parser.add_argument("--mfc", type=str, metavar="DICT", default='')
parser.add_argument("-nc", "--no-chem", action="store_false",
    help="Disable chemistry.", dest='chem', default=True)
parser.add_argument("-s",  "--scale", type=float, default=1, help="Scale.")
parser.add_argument("-d",  "--ndim", type=int, default=1,
    help="Number of dimensions.")
 
args = parser.parse_args()

ctfile  = 'gri30.yaml'
sol     = ct.Solution(ctfile)

sol.TPX = 1300, ct.one_atm, {'H2': 0.43, 'O2': 1.6*0.43, 'AR': 1 - 0.43 - 1.6*0.43}

Nx   = 20*args.scale
Tend = 1e-6
s    = 1e-2
dt   = 1e-9

NT = int(Tend/dt)
NS = NT//60

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -s/2,
    'x_domain%end'                 : +s/2,
    'y_domain%beg'                 : -s/2,
    'y_domain%end'                 : +s/2,
    'z_domain%beg'                 : -s/2,
    'z_domain%end'                 : +s/2,
    'm'                            : Nx,
    'n'                            : Nx,
    'p'                            : Nx,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : NS,
    'parallel_io'                  : 'T' if ndim > 1 else 'F',

    # Simulation Algorithm Parameters ==========================================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 1,
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
    'bc_x%beg'                     :-1,
    'bc_x%end'                     :-1,
    'bc_y%beg'                     :-1,
    'bc_y%end'                     :-1,
    'bc_z%beg'                     :-1,
    'bc_z%end'                     :-1,

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%z_centroid'     : 0,
    'patch_icpp(1)%length_x'       : s,
    'patch_icpp(1)%length_y'       : s,
    'patch_icpp(1)%length_z'       : s,
    'patch_icpp(1)%vel(1)'         : 0,
    'patch_icpp(1)%vel(2)'         : 0,
    'patch_icpp(1)%vel(3)'         : 0,
    'patch_icpp(1)%pres'           : sol.P,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : sol.density,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================

    # Chemistry ================================================
    'cantera_file'                 : ctfile,
    # ==========================================================
}

if args.chem:
    # Chemistry ================================================================
    case.update({
        'chemistry'             : 'T',
        'chem_params%advection' : 'F',
        'chem_params%diffusion' : 'F',
        'chem_params%reactions' : 'T',
    })

    for i in range(len(sol.Y)):
        case[f'patch_icpp(1)%Y({i+1})'] = sol.Y[i]

rmdims = []
if args.ndim < 3: rmdims += ['z']
if args.ndim < 2: rmdims += ['y']

dkeys = set()
for key in case.keys():
    for dim in rmdims:
        dirid = 3 if dim == 'z' else 2

        if (  re.match(f'.+_{dim}', key) or re.match(f'{dim}_.+', key)
           or re.match(f'%{dim}', key)   or f'%vel({dirid})' in key):
            dkeys.add(key)
            break

case = {k: v for k, v in case.items() if k not in dkeys}

if __name__ == '__main__':
    print(json.dumps(case))
