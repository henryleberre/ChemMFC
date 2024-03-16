import glob
import pathlib
import matplotlib.pyplot as plt

import cantera as ct
import numpy as np

# Constants

from case import dt, Tend, SAVE_COUNT, sol

BG_COLOR = '#1a1a1a'
TX_COLOR = '#FFFFFF'

import seaborn as sns

sns.set_style('dark', {
    'axes.facecolor':    '#121212',
    'axes.edgecolor':    BG_COLOR,
    'axes.labelcolor':   TX_COLOR,
    'text.color':        TX_COLOR,
    'xtick.color':       TX_COLOR,
    'ytick.color':       TX_COLOR,
    'grid.color':        BG_COLOR,
    'figure.facecolor':  BG_COLOR,
    'figure.edgecolor':  BG_COLOR,
    'savefig.facecolor': BG_COLOR,
    'savefig.edgecolor': BG_COLOR,
})

# Initialize lists
timesteps = []

casedir = pathlib.Path(__file__).parent.resolve()

# Get timesteps from file names
for filename in glob.glob(f"{casedir}/D/cons.5.*.dat"):
    timesteps.append(int(filename.split(".")[-2]))

# Sort timesteps
timesteps.sort()

# Calculate corresponding times
times = [timestep * dt for timestep in timesteps]

# Initialize DATA list
DATA = []

# Read data for each variable and timestep
for timestep in timesteps:
    values = []
    for var in range(5, 14 + 1):
        filepath = f"{casedir}/D/prim.{var}.00.{timestep:06d}.dat"
        with open(filepath) as f:
            values.append(float(f.readline().split()[1]))
    DATA.append(values)

# Transpose DATA to separate each variable's data
MFC_DATA = list(map(list, zip(*DATA)))

# Legends dictionary
legends = {
    5: 'H_2',
    6: 'H',
    7: 'O',
    8: 'O_2',
    9: 'OH',
    10: 'H_2O',
    11: 'HO_2',
    12: 'H_2O_2',
    13: 'AR',
    #14: 'N2',
}

# Create a reactor and add the gas mixture to it
reactor = ct.IdealGasReactor(sol)

# Create a reactor network
reactor_network = ct.ReactorNet([reactor])

# Time parameters
time      = 0.0             # Initial time
time_step = Tend/SAVE_COUNT # Time step (seconds)

# Lists to store time and temperature for plotting
time_history = []
temperature_history = []
pressure = []
density = []
Ys = []
energys = []

# Integrate the reactor network over time
while time < Tend:
    #time = reactor_network.step()
    time_history.append(time)
    Ys.append(reactor.thermo.Y)
    density.append(reactor.thermo.density)
    pressure.append(reactor.thermo.P)
    temperature_history.append(reactor.T)
    energys.append(reactor.thermo.int_energy_mass)

    reactor_network.advance(time + time_step)
    time += time_step

time_history.append(time)
Ys.append(reactor.thermo.Y)
density.append(reactor.thermo.density)
pressure.append(reactor.thermo.P)
temperature_history.append(reactor.T)
energys.append(reactor.thermo.int_energy_mass)

# Plot each variable over time
for i, var in enumerate(range(5, 14 + 1)):
    if var in legends:
        plt.plot(times, MFC_DATA[i], label=f"${{{legends[var]}}}$")

# Now plot cantera results in dotted lines, use same colors and dashed lines.
for i, var in enumerate(range(5, 14 + 1)):
    if var in legends:
        plt.plot(time_history, [y[var-5] for y in Ys], linestyle='dashed', color=plt.gca().lines[i+1].get_color())

# Y
plt.xlabel("Time (s)")
plt.ylabel("Mass Fractions")
plt.legend()
plt.title("Mass Fractions (MFC v. Cantera)")
plt.savefig(f"{casedir}/compare_Y.png", dpi=700)
plt.close()

pressure_data = []
energy_data = []

# Read pressure data for each timestep
for timestep in timesteps:
    filepath = f"{casedir}/D/prim.{3}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        pressure_data.append(float(f.readline().split()[1]))
    filepath = f"{casedir}/D/cons.{3}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        energy_data.append(float(f.readline().split()[1]))

# Pressure
plt.plot(times, pressure_data)
plt.plot(time_history, pressure, linestyle='dashed')
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.legend(["MFC", "Cantera"])
plt.title("Pressure (MFC v. Cantera)")
plt.savefig(f"{casedir}/compare_p.png", dpi=700)
plt.close()

# Energy
plt.plot(times, energy_data)
plt.plot(time_history, energys, linestyle='dashed')
plt.xlabel("Time (s)")
plt.ylabel("Internal Energy")
plt.legend(["MFC", "Cantera"])
plt.title("Internal Energy (MFC v. Cantera)")
plt.savefig(f"{casedir}/compare_e_t.png", dpi=700)
plt.close()

# Temperature
temperatures = []

# Read pressure data for each timestep
for timestep in timesteps:
    filepath = f"{casedir}/D/cons.{15}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        temperatures.append(float(f.readline().split()[1]))

print(len(times), len(temperatures))
print(len(time_history), len(temperature_history))

plt.plot(times, temperatures)
plt.plot(time_history, temperature_history, linestyle='dashed')
plt.xlabel("Time (s)")
plt.ylabel("Temperature")
plt.legend(["MFC", "Cantera"])
plt.title("Temperature (MFC v. Cantera)")
plt.ticklabel_format(style='plain')
plt.ticklabel_format(useOffset=False)
plt.savefig(f"{casedir}/compare_T.png", dpi=700)
plt.close()
