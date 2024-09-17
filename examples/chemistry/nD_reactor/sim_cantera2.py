import copy
import numpy as np
import cantera as ct
import pyrometheus as pyro
import matplotlib.pyplot as plt

sol = ct.Solution('h2o2.yaml')

pyro_cls = pyro.codegen.python.get_thermochem_class(sol)
pyro_gas = pyro_cls(np)

pressure    = pyro_gas.one_atm
temperature = 1000

# Cantera
sol.TPY = temperature, pressure, 'H2:2,O2:1,N2:7'

# Constant density (using pyro)
density = pyro_gas.get_density(pressure, temperature, sol.Y)

print(f"Density: {density} v. {sol.density}")

def rhs(y, density, temperature, use_pyro: bool):
    if not use_pyro:
        sol.TDY = temperature, density, y
        wdot = sol.net_production_rates
        qdot = -np.sum(wdot * sol.partial_molar_enthalpies) / (sol.density * sol.cp_mass)
        return wdot * sol.molecular_weights / sol.density, qdot

    return pyro_gas.get_net_production_rates(
        density, temperature, y
    ) * pyro_gas.molecular_weights / density

def runge_kutta(dt, y, temperature, use_pyro: bool):
    global density
    k_1, q_1 = rhs(y, density, temperature, use_pyro)
    k_2, q_2 = rhs(y + 0.5 * dt * k_1, density, temperature + 0.5 * dt * q_1, use_pyro)
    k_3, q_3 = rhs(y + 0.5 * dt * k_2, density, temperature + 0.5 * dt * q_2, use_pyro)
    k_4, q_4 = rhs(y + dt * k_3, density, temperature + dt * q_3, use_pyro)
    y_new = y + dt * (k_1 + 2 * (k_2 + k_3) + k_4) / 6
    temperature_new = temperature + dt * (q_1 + 2 * (q_2 + q_3) + q_4) / 6
    return y_new, temperature_new

def time_march(dt, y_init, t_final, use_pyro: bool):
    global density, pressure
    Ys, rhos, temperatures, pressures = [], [], [], []
    t = 0
    step = 0
    temp = temperature
    temperatures.append(temp)
    y = copy.deepcopy(y_init)
    Ys.append(copy.deepcopy(y))
    rhos.append(density)
    pressures.append(pressure)
    while t < t_final:
        y, temp = runge_kutta(dt, y, temp, use_pyro)
        #y = np.maximum(y, 0)
        #y = np.minimum(y, 1)
        #y /= np.sum(y)
        if not use_pyro:
            density = sol.density
            sol.TDY = temp, density, y
            pressure = sol.P
        else:
            pass
        #print(f"Step {step}: {y}")
        t += dt
        step = step + 1
        Ys.append(copy.deepcopy(y))
        rhos.append(density)
        pressures.append(pressure)
        temperatures.append(temp)
    return Ys, rhos, temperatures, pressures

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(15, 20))
# make the entire image a square
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

y_init  = copy.deepcopy(sol.Y)
t_final = 0.12e-3
dt      = 1e-8

Ys_ct, rhos_ct, temps_ct, pressures_ct = time_march(dt, y_init, t_final, False)
#Ys_py, rhos_py, temps_py = time_march(dt, y_init, t_final, True)

axes[0,0].plot(np.linspace(0, t_final, len(Ys_ct)), Ys_ct, marker='o', markersize=1)
axes[0,0].set_title('Cantera')
axes[0,0].legend(sol.species_names)
axes[0,0].set_xlabel('Time (s)')
axes[0,0].set_ylabel('Mass Fraction')

every_nth = 2
for n, label in enumerate(axes[0,0].xaxis.get_ticklabels()):
    if n % every_nth != 0:
        label.set_visible(False)

#axes[0,1].plot(np.linspace(0, t_final, len(Ys_py)), Ys_py)
#axes[0,1].set_title('Pyro')
#axes[0,1].legend(sol.species_names)
#axes[0,1].set_xlabel('Time (s)')
#axes[0,1].set_ylabel('Mass Fraction')

axes[1,0].plot(np.linspace(0, t_final, len(rhos_ct)), rhos_ct)
axes[1,0].set_title('Cantera')
axes[1,0].set_xlabel('Time (s)')
axes[1,0].set_ylabel('Density')

#axes[1,1].plot(np.linspace(0, t_final, len(rhos_py)), rhos_py)
#axes[1,1].set_title('Pyro')
#axes[1,1].set_xlabel('Time (s)')
#axes[1,1].set_ylabel('Density')

axes[2,0].plot(np.linspace(0, t_final, len(temps_ct)), temps_ct)
axes[2,0].set_title('Cantera')
axes[2,0].set_xlabel('Time (s)')
axes[2,0].set_ylabel('Temperature')

#axes[2,1].plot(np.linspace(0, t_final, len(temps_py)), temps_py)
#axes[2,1].set_title('Pyro')
#axes[2,1].set_xlabel('Time (s)')
#axes[2,1].set_ylabel('Temperature')

axes[3,0].plot(np.linspace(0, t_final, len(pressures_ct)), pressures_ct)
axes[3,0].set_title('Cantera')
axes[3,0].set_xlabel('Time (s)')
axes[3,0].set_ylabel('Pressure')

plt.savefig('py_Ys.png')
