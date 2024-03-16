import glob
import pathlib
import matplotlib.pyplot as plt

# Constants

dt=7e-10
Tend=0.12e-3
NT=int(Tend/dt)
NS=NT//60

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
DATA = list(map(list, zip(*DATA)))

# Legends dictionary
legends = {
    5: 'H2',
    6: 'H',
    7: 'O',
    8: 'O2',
    9: 'OH',
    10: 'H2O',
    11: 'HO2',
    12: 'H2O2',
    13: 'AR',
    14: 'N2',
}

# Plot each variable over time
for i, var in enumerate(range(5, 14 + 1)):
    plt.plot(times, DATA[i], label=f"prim.{var} : Y_{legends[var]}")

plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.legend()
plt.savefig(f"{casedir}/mfc_Y.png")
plt.close()

# Initialize DATA list for pressure
pressure_data = []

# Read pressure data for each timestep
for timestep in timesteps:
    filepath = f"{casedir}/D/prim.{3}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        pressure_data.append(float(f.readline().split()[1]))

# Plot pressure over time
plt.plot(times, pressure_data)
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.legend(["prim.3 : Pressure"])
plt.savefig(f"{casedir}/mfc_p.png")
plt.close()

# Initialize DATA list for pressure
temperatures = []

# Read pressure data for each timestep
for timestep in timesteps:
    filepath = f"{casedir}/D/cons.{15}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        temperatures.append(float(f.readline().split()[1]))

# Plot pressure over time
plt.plot(times, temperatures)
plt.xlabel("Time (s)")
plt.ylabel("Temperature")
plt.legend(["cons.15 : Temperature"])
plt.ticklabel_format(style='plain')
plt.ticklabel_format(useOffset=False)
plt.savefig(f"{casedir}/mfc_T.png")
plt.close()

# Initialize DATA list for pressure
densities = []

# Read pressure data for each timestep
for timestep in timesteps:
    filepath = f"{casedir}/D/cons.{1}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        densities.append(float(f.readline().split()[1]))

# Plot pressure over time
plt.plot(times, densities)
plt.xlabel("Time (s)")
plt.ylabel("Density")
plt.legend(["cons.1 : Density"])
plt.ticklabel_format(useOffset=False)
plt.savefig(f"{casedir}/mfc_rho.png")
plt.close()

velocity = []

for timestep in timesteps:
    filepath = f"{casedir}/D/cons.{2}.00.{timestep:06d}.dat"
    with open(filepath) as f:
        velocity.append(float(f.readline().split()[1]))

plt.plot(times, velocity)
plt.xlabel("Time (s)")
plt.ylabel("Velocity")
plt.legend(["cons.2 : Velocity"])
plt.ticklabel_format(useOffset=False)
plt.savefig(f"{casedir}/mfc_u.png")
