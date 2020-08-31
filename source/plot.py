import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

filename ='results/{}_{}.csv'
column_names = ['x', 't', 'n']

def plot_pt_comp(mol):
    no_pt = pd.read_csv(filename.format(mol, 'False'), names=column_names)
    no_pt['t_cs'] = no_pt['t'].cumsum()
    no_pt['n_cs'] = no_pt['n'].cumsum()
    pt = pd.read_csv(filename.format(mol, 'True'), names=column_names)
    pt['t_cs'] = pt['t'].cumsum()
    pt['n_cs'] = pt['n'].cumsum()

    data = no_pt.join(pt.set_index('x'), on='x', rsuffix='_pt')


    print(data)

    fig, axes = plt.subplots(2, 2, figsize=(13,10))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)

    #data.plot(x='x', subplots=True, title='H2 VQE PES time, Parameter Transfer (PT)')
    data.plot(x='x', y=['t', 't_pt'], ax=axes[0,0], ylabel='VQE Time [s]', xlabel='Seperation [Å]', title=f'{mol} VQE PES time')
    data.plot(x='x', y=['n', 'n_pt'], ax=axes[0,1], ylabel='Optimization iterations', xlabel='Seperation [Å]', title=f'{mol} VQE PES iterations')
    data.plot(x='x', y=['t_cs', 't_cs_pt'], ax=axes[1,0], ylabel='VQE Time [s]', xlabel='Seperation [Å]', title=f'{mol} VQE PES cumulative time')
    data.plot(x='x', y=['n_cs', 'n_cs_pt'], ax=axes[1,1], ylabel='Optimization iterations', xlabel='Seperation [Å]', title=f'{mol} VQE PES cumulative iterations')

    for ax in axes.flatten(): 
        ax.grid(True)
        ax.legend(['W/o PT', 'With PT'])


plot_pt_comp('H2')
plot_pt_comp('LiH')
plot_pt_comp('BeH2')

plt.show()
