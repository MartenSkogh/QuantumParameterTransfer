import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import sys

filename ='results/{}_{}.csv'
column_names = ['x', 't', 'n', 'E']

def plot_pt_comp(mol):
    no_pt = pd.read_csv(filename.format(mol, 'False'), names=column_names)
    no_pt['t_cs'] = no_pt['t'].cumsum()
    no_pt['n_cs'] = no_pt['n'].cumsum()
    pt = pd.read_csv(filename.format(mol, 'True'), names=column_names)
    pt['t_cs'] = pt['t'].cumsum()
    pt['n_cs'] = pt['n'].cumsum()

    data = no_pt.join(pt.set_index('x'), on='x', rsuffix='_pt')


    print(data)

    fig = plt.figure(figsize=(13,10))

    gs = gridspec.GridSpec(2, 3)
    axes = []


    #fig, axes = plt.subplots(3, 2, figsize=(13,10))

    #data.plot(x='x', subplots=True, title='H2 VQE PES time, Parameter Transfer (PT)')
    axes.append(fig.add_subplot(gs[0,0]))
    data.plot(x='x', y=['t', 't_pt'], ax=axes[-1], ylabel='VQE Time [s]', xlabel='Seperation [Å]', title=f'{mol} VQE PES time')

    axes.append(fig.add_subplot(gs[1,0]))
    data.plot(x='x', y=['n', 'n_pt'], ax=axes[-1], ylabel='Optimization iterations', xlabel='Seperation [Å]', title=f'{mol} VQE PES iterations')
    
    axes.append(fig.add_subplot(gs[0,1]))
    data.plot(x='x', y=['t_cs', 't_cs_pt'], ax=axes[-1], ylabel='VQE Time [s]', xlabel='Seperation [Å]', title=f'{mol} VQE PES cumulative time')
    
    axes.append(fig.add_subplot(gs[1,1]))
    data.plot(x='x', y=['n_cs', 'n_cs_pt'], ax=axes[-1], ylabel='Optimization iterations', xlabel='Seperation [Å]', title=f'{mol} VQE PES cumulative iterations')
    
    axes.append(fig.add_subplot(gs[:,2]))
    data.plot(x='x', y=['E', 'E_pt'], ax=axes[-1], ylabel='VQE Time [s]', xlabel='Seperation [Å]', title=f'{mol} VQE PES cumulative time')

    #plt.subplots_adjust(wspace=0.3, hspace=0.3)
    #plt.tight_layout()

    for ax in axes: 
        ax.grid(True)
        ax.legend(['W/o PT', 'With PT'])


    gs.tight_layout(fig)

plot_pt_comp('H2')
plot_pt_comp('LiH')
#plot_pt_comp('BeH2')

plt.show()
