#!/usr/bin/env python3

import re
import json
import argparse


parser = argparse.ArgumentParser(
    prog="nD_chem_binary",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("dict", type=str,   metavar="DICT", help=argparse.SUPPRESS)
parser.add_argument("--ndim", type=int,   metavar="NDIM", default=1)
parser.add_argument("--ts",   type=float, metavar="TS", default=0.001)
parser.add_argument("--te",   type=float, metavar="TE", default=0.010)

ARGS = vars(parser.parse_args())
ndim = ARGS['ndim']

ts, te = ARGS['ts'], ARGS['te']
pa=997
pb=789
D=0.0003
# The domain is 10 mixing layers wide to avoid boundary effects.
s=10*(2*D*te)**0.5
# Have at least 100 cells across the mixing layer.
N=100*int(s/((2*D*ts)**0.5))
v=0
dt=1e-9
NT=int((te-ts)/dt)

r = '((x**2d0'
if ndim > 1:
    r += '+y**2d0'
if ndim > 2:
    r += '+z**2d0'
r += ')**0.5d0)'

# Left - Right binary
r='x'
erf=f'(1.0d0 - erf({r}/(2.0d0*({(D*ts):.15f})**(0.5d0))))'
y1=f'{pa/2.0}*{erf}/({pb}+({pa - pb})/2.0d0*{erf})'
a1=f'{pb} + ({pa - pb})/2.0d0*{erf}'

# 2D polar
#a=s/8
#erfp=f'erf(({r}+{a})/(2.0d0*({(D*ts):.15f})**(0.5d0)))'
#erfm=f'erf(({r}-{a})/(2.0d0*({(D*ts):.15f})**(0.5d0)))'
#a1=f'({pb} + ({pa - pb})/2.0d0*({erfp} - {erfm}))'
#y1=f'({pa}*({pb} - {a1}))/({a1}*({pb} - {pa}))'

y2=f'1.0d0-{y1}'

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -s/2,
    'x_domain%end'                 : +s/2,
    'y_domain%beg'                 : -s/2 if ndim > 1 else 0,
    'y_domain%end'                 : +s/2 if ndim > 1 else 0,
    'z_domain%beg'                 : -s/2 if ndim > 2 else 0,
    'z_domain%end'                 : +s/2 if ndim > 2 else 0,
    'm'                            : N,
    'n'                            : N if ndim > 1 else 0,
    'p'                            : N if ndim > 2 else 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : int(NT),
    't_step_save'                  : int(NT/(30*5)),
    't_step_print'                 : int(NT/(30*5)),
    'parallel_io'                  : 'T' if ndim != 1 and json.loads(ARGS['dict'])['mpi'] else 'F',

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
    'bc_x%beg'                     :-3,
    'bc_x%end'                     :-3,
    'bc_y%beg'                     :-3,
    'bc_y%end'                     :-3,
    'bc_z%beg'                     :-3,
    'bc_z%end'                     :-3,

    # Chemistry ================================================================
    'chemistry'                    : 'T',
    'chem_params%advection'        : 'F',
    'chem_params%diffusion'        : 'T',
    'chem_params%reactions'        : 'F',
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
    'patch_icpp(1)%alpha_rho(1)'   : a1,
    'patch_icpp(1)%Y(1)'           : y1,
    'patch_icpp(1)%Y(2)'           : y2,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
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
