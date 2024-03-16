#!/usr/bin/env python3

import re
import json
import argparse


N=200
s=1e-2
v=100000*s
dt=10e-9
NT=int(((2*s)/v)/dt)

parser = argparse.ArgumentParser(
    prog="nD_chem_binary",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("dict", type=str, metavar="DICT", help=argparse.SUPPRESS)
parser.add_argument("ndim", type=int, metavar="NDIM", default=1)

ARGS = vars(parser.parse_args())
ndim = ARGS['ndim']

y1=f'(1d0+sin(2d0*3.141592d0*x/(2d0*{s})))/2d0'
y2=f'1.0d0-{y1}'

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -s,
    'x_domain%end'                 : +s,
    'y_domain%beg'                 : -s if ndim > 1 else 0,
    'y_domain%end'                 : +s if ndim > 1 else 0,
    'z_domain%beg'                 : -s if ndim > 2 else 0,
    'z_domain%end'                 : +s if ndim > 2 else 0,
    'm'                            : N,
    'n'                            : N if ndim > 1 else 0,
    'p'                            : N if ndim > 2 else 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : int(NT),
    't_step_save'                  : int(NT/(5*30)),
    't_step_print'                 : int(NT/(5*30)),

    # Simulation Algorithm Parameters ==========================================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 1,
    'adv_alphan'                   : 'F',
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
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    'patch_icpp(1)%geometry'       : 3**(ndim - 1),
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%z_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 2*s,
    'patch_icpp(1)%length_y'       : 2*s,
    'patch_icpp(1)%length_z'       : 2*s,
    'patch_icpp(1)%vel(1)'         : v,
    'patch_icpp(1)%vel(2)'         : v,
    'patch_icpp(1)%vel(3)'         : v,
    'patch_icpp(1)%pres'           : 10000,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : 1000,
    'patch_icpp(1)%Y(1)'           : y1,
    'patch_icpp(1)%Y(2)'           : y2,

    # Chemistry ================================================================
    'chemistry'                    : 'T',
    'chem_params%advection'        : 'T',
    'chem_params%diffusion'        : 'T',
    'chem_params%reactions'        : 'F',
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.0E+08/(4.4E+00-1.E+00),
    # ==========================================================================

    # Chemistry ================================================
    'cantera_file'                 : 'h2o2.yaml',
    # ==========================================================
}

rmdims = []
if ndim < 3:
    rmdims += ['z']
if ndim < 2:
    rmdims += ['y']

dkeys = set()
for key in case.keys():
    for dim in rmdims:
        dirid = 3 if dim == 'z' else 2

        rm = False
        rm = rm or re.match(f'.+_{dim}', key)
        rm = rm or re.match(f'{dim}_.+', key)
        rm = rm or re.match(f'%{dim}', key)
        rm = rm or f'%vel({dirid})' in key

        if rm:
            dkeys.add(key)
            break

case = {k: v for k, v in case.items() if k not in dkeys}

print(json.dumps(case))
