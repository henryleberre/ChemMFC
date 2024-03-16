import cantera as ct
import numpy as np
import pathlib
import matplotlib.pyplot as plt

casedir = pathlib.Path(__file__).parent.resolve()

from case import sol as gas
from case import Tend as end_time

# Create a reactor and add the gas mixture to it
reactor = ct.IdealGasReactor(gas)

# Create a reactor network
reactor_network = ct.ReactorNet([reactor])

# Time parameters
time      = 0.0  # Initial time
time_step = end_time/100 # Time step (seconds)

# Lists to store time and temperature for plotting
time_history = []
temperature_history = []
pressure = []
density = []
Ys = []

# Integrate the reactor network over time
while time < end_time:
    #time = reactor_network.step()
    time_history.append(time)
    temperature_history.append(reactor.T)
    Ys.append(reactor.thermo.Y)
    density.append(reactor.thermo.density)
    pressure.append(reactor.thermo.P)

    reactor_network.advance(time + time_step)
    time += time_step

    print(time, reactor.thermo.Y)

# Plot the results
plt.figure()
plt.plot(time_history, temperature_history)
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Evolution Over Time')
plt.savefig(f"{casedir}/ct_T.png")

# Plot the results
plt.figure()
plt.plot(time_history, Ys)
plt.xlabel('Time (s)')
plt.ylabel('Mass Fractions')
plt.title('Mass Fractions Evolution Over Time')
plt.legend(gas.species_names)
# only show half the ticks
plt.xticks([x*100000 for x in time_history[::int(8*end_time/time_step)]])
plt.savefig(f"{casedir}/ct_Y.png")

# Plot the results
plt.figure()
plt.plot(time_history, pressure)
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure Evolution Over Time')
plt.savefig(f"{casedir}/ct_p.png")

# Plot the results
plt.figure()
plt.plot(time_history, density)
plt.xlabel('Time (s)')
plt.ylabel('Density (kg/m^3)')
plt.title('Density Evolution Over Time')
plt.savefig(f"{casedir}/ct_rho.png")

# 7.999999999999997e-05