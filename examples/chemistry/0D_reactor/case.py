#!/usr/bin/env python3

# https://www.sciencedirect.com/science/article/pii/S0045793013003976?via=ihub
# 4.6. Perfectly stirred reactor

import json
import cantera as ct

ctfile  = 'h2o2.yaml'
sol     = ct.Solution(ctfile)

sol.TPY = 1200, ct.one_atm, {'H2': 0.1, 'O2': 0.2, 'N2': 0.7}
dt = 1e-8
Tend = 1e-4

NT=int(Tend/dt)
SAVE_COUNT=60
NS=NT//SAVE_COUNT
Nx=25
s=1e-2

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -s/2,
    'x_domain%end'                 : +s/2,
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : NS,
    'parallel_io'                  : 'F',

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

    # Chemistry ================================================================
    'chemistry'                    : 'T',
    'chem_params%advection'        : 'F',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'T',
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%length_x'       : s,
    'patch_icpp(1)%vel(1)'         : 0,
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

for i in range(len(sol.Y)):
    case[f'patch_icpp(1)%Y({i+1})'] = sol.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))
