import os
import glob
import typing
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BG_COLOR = '#1a1a1a'
TX_COLOR = '#FFFFFF'

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

class Case:
    def __init__(self, dirpath: str, is_0d = False):
        self._dirpath   = dirpath
        self._data      = {}
        self._procs     = set()
        self._timesteps = set()
        self._ndims     = 0
        self._coords    = [set(), set(), set()]

        for f in glob.glob(os.path.join(self._dirpath, 'D', f'cons.1.*.*.dat')):
            self._procs.add(int(f.split('.')[-3]))
            self._timesteps.add(int(f.split('.')[-2]))

        for i, t_step in enumerate(self._timesteps):
            df = pd.DataFrame()
            for proc in self._procs:
                df = pd.concat([
                    df,
                    pd.read_csv(
                        os.path.join(self._dirpath, 'D', f'cons.1.{proc:02d}.{t_step:06d}.dat'),
                        sep=r'\s+', header=None, names=['x'], usecols=['x']
                    )
                ])

            self._data[t_step] = df

            if i == 0:
                self._coords[0] = set(df['x'])
                self._ndims     = 1 + min(len(self._coords[1]), 1) + min(len(self._coords[2]), 1)

    def get_ndims(self)     -> int:  return self._ndims
    def get_coords(self)    -> set:  return self._coords
    def get_timesteps(self) -> set:  return self._timesteps
    def get_procs(self)     -> set:  return self._procs
    def get_data(self)      -> dict: return self._data

    def define_variable(self, name: str, func: typing.Callable):
        for t_step in self._data.keys():
            xs = self._data[t_step]['x']

            self._data[t_step] = pd.merge(
                self._data[t_step],
                pd.DataFrame({
                    'x':  xs,
                    name: [ func(x, t_step, self._data[t_step]) for x in xs ]
                })
            )

    def load_variable(self, name: str, path: str):
        for t_step in self._timesteps:
            dfs = []
            for proc in self._procs:
                dfs.append(
                    pd.read_csv(
                        os.path.join(self._dirpath, 'D', f'{path}.{proc:02d}.{t_step:06d}.dat'),
                        sep=r'\s+', header=None, names=['x', name]
                    )
                )

            self._data[t_step] = pd.merge(self._data[t_step], pd.concat(dfs))

    def plot(self, t_step: int, varname: str):
        sns.lineplot(
            data=self._data[t_step], x='x', y=varname,
            color="white", linewidth=1.5, alpha=1,
            ax=None, dashes=False, label=varname)


from case import sol

case = Case(".")

for i, name in enumerate(sol.species_names):
    case.load_variable(name, f"prim.{5+i}")

case.load_variable("rho", f"cons.1")

for step in case.get_timesteps():
    for i, name in enumerate(sol.species_names):
        case.plot(step, name)

    plt.show()
    plt.cla()

    case.get_data()[step].to_csv(f"step-{step}.csv", sep='\t')

case.define_variable("test", lambda x, t, df: x*x*x*x*x*x*x*x)

print(case._data[0].to_csv('test.csv', sep='\t'))
