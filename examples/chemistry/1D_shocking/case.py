#!/usr/bin/env python3

# https://www.sciencedirect.com/science/article/pii/S0045793013003976?via=ihub
# 4.6. Perfectly stirred reactor

import json
import cantera as ct

ctfile  = 'h2o2.yaml'
sol     = ct.Solution(ctfile)
sol.TPX = 1400, ct.one_atm, 'H2:2,O2:1,N2:7'

rho_l = sol.density
u_l   = 0
p_l   = sol.P

# Speeds of sound:
# N_2 = 968 m/s
#  Ar = 630 m/s
# 
# Test low temperature 500 k
# Track number of get_temperature iterations
# Build up

rho_r = sol.density
u_r   = -487.34
p_r   = sol.P
L  = 0.12
Nx = 1000
dx = L/Nx
dt = dx/abs(u_r)*0.1
Tend=230e-6

NT=int(Tend/dt)
SAVE_COUNT=100
NS=NT//SAVE_COUNT

chemistry = True

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -L/2,
    'x_domain%end'                 : +L/2,
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS, # 1592
    't_step_print'                 : NS, # 1592
    'parallel_io'                  : 'F',

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
    'bc_x%end'                     :-3,

    # Chemistry ================================================================
    'chemistry'                    : 'F' if not chemistry else 'T',
    'chem_params%advection'        : 'T',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'T',
    'cantera_file'                 : ctfile,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : -L/4,
    'patch_icpp(1)%length_x'       : L/2,
    'patch_icpp(1)%vel(1)'         : u_l, #f'{u_l} + {u_r-u_l}/2d0 + ({u_r-u_l}/2d0)*(1+tanh(x/0.005))',
    'patch_icpp(1)%pres'           : p_l, #f'{p_l} + {p_r-p_l}/2d0 + ({p_r-p_l}/2d0)*(1+tanh(x/0.005))',
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : rho_l, #f'{rho_l} + {rho_r-rho_l}/2d0 + ({rho_r-rho_l}/2d0)*(1+tanh(x/0.005))',
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : L/4,
    'patch_icpp(2)%length_x'       : L/2,
    'patch_icpp(2)%vel(1)'         : u_r, #f'{u_l} + {u_r-u_l}/2d0 + ({u_r-u_l}/2d0)*(1+tanh(x/0.005))',
    'patch_icpp(2)%pres'           : p_r, #f'{p_l} + {p_r-p_l}/2d0 + ({p_r-p_l}/2d0)*(1+tanh(x/0.005))',
    'patch_icpp(2)%alpha(1)'       : 1,
    'patch_icpp(2)%alpha_rho(1)'   : rho_r, #f'{rho_l} + {rho_r-rho_l}/2d0 + ({rho_r-rho_l}/2d0)*(1+tanh(x/0.005))',
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================

    # Chemistry ================================================
    # ==========================================================
}

if chemistry:
    for i in range(len(sol.Y)):
        case[f'patch_icpp(1)%Y({i+1})'] = sol.Y[i]
        case[f'patch_icpp(2)%Y({i+1})'] = sol.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))
