import os
import numpy as np
import numpy.ma as ma
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.basemap import Basemap
import string
import speedy.grads as grads
import speedy.noleap as nl

plt.switch_backend('agg')

# define color bar
mycm = plt.cm.RdBu_r
mycm.set_bad('silver',1.0)

ft_start_list = np.arange('1870-01','1875-01',1,dtype='datetime64[M]')
n = grads.grads('/data15/okazaki/speedy_ver41.5/work/nature_icsea-2/out_monthly/attm000.ctl')
Xm, Ym = np.meshgrid(n.lon, n.lat)

# parameters
ft_max = 1
var_name_list = ['temp','sst']
res_list = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
exp_list = ['ctrl','long']
dirref = {'Serr' : '../../error_growth_Serr', 'Lerr' : '../../error_growth_Lerr', 'Eerr' : '../../error_growth_Eerr', 'noint':'../../error_growth_noint/'}
err_list = ['Serr','Lerr','Eerr','noint']
ls = {'Serr':'dotted', 'Lerr':'dashed', 'Eerr':'solid', 'noint':'solid'}
lw = {'Serr':2, 'Lerr':2, 'Eerr':2, 'noint':4}
col = {'ctrl':'', 'long':''}
col['ctrl'] = {'Serr':'k', 'Lerr':'k', 'Eerr':'k', 'noint':'lightgray'}
col['long'] = {'Serr':'r', 'Lerr':'r', 'Eerr':'r', 'noint':'lightpink'}
title = {'temp':'SAT', 'sst':'SST'}

nt = 60
ne = 29
nx = 96
ny = 48
nz = 8
df = nt-2
recl = 4 *(((nx * ny + 2) * nz * 9) + ((nx * ny + 2) * 32 ))

cor = np.empty((ny,nx),dtype=float)
cor_glb = np.empty((len(err_list),len(exp_list),len(res_list),ft_max),dtype=float)

fig = plt.figure(figsize=(12,4))

for vvv,var_name in enumerate(var_name_list) :

    npyFileName = 'cor_glb_' + var_name + '.npy'

    if os.path.isfile( npyFileName ) :
        cor_glb = np.load( npyFileName )
    else :
        print('file not found --skip')
        print(npyFileName)

    r_thresh = 0.2143825406137366


    ### global mean
    ax = fig.add_subplot(1,3,vvv+1)

    for i,err in enumerate(err_list) :

        for j,exp in enumerate(exp_list) :

            #--- plot
            if exp == 'long' and i < 2 :
                continue 
            plt.plot([0,1,2,3,4,5,6],cor_glb[i,j,:,0], color=col[exp][err], linestyle=ls[err], linewidth=2)
            plt.plot([0,1,2,3,4,5,6],ma.masked_less(cor_glb[i,j,:,0],r_thresh), color='w', linewidth=0, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col[exp][err], alpha=1.0)
            plt.ylim(0,1)
            ax.set_xticklabels(res_list)
            ax.set_xticks([0,1,2,3,4,5,6])
            plt.yticks(color='k')
            plt.ylabel('Correlation')
            plt.xlabel('Time averaging')
            
            #--- title
            plt.title('(' + string.ascii_lowercase[vvv] + ')', loc='left')
            plt.title(title[var_name], loc='center')
    

### legend
ax = fig.add_subplot(1,3,3)
d1, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['ctrl']['Serr'], ls=ls['Serr'], linewidth=2, label=r'intrinsic($\tau_{ocn}$=90d)'  )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['ctrl']['Serr'], markerfacecolor='w')  
d2, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['ctrl']['Lerr'], ls=ls['Lerr'], linewidth=2, label=r'potential($\tau_{ocn}$=90d)'  )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['ctrl']['Lerr'], markerfacecolor='w')
d3, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['ctrl']['Eerr'], ls=ls['Eerr'], linewidth=2, label=r'practical($\tau_{ocn}$=90d)'  )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['ctrl']['Eerr'], markerfacecolor='w')
d4, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['ctrl']['noint'],ls=ls['noint'], linewidth=2,label=r'noint($\tau_{ocn}$=90d)'      )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['ctrl']['noint'], markerfacecolor='w')
d5, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['long']['Eerr'], ls=ls['Eerr'], linewidth=2, label=r'practical($\tau_{ocn}$=360d)' )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['long']['Eerr'], markerfacecolor='w')  
d6, = plt.plot(np.arange(7),cor_glb[0,0,:,0], color=col['long']['noint'],ls=ls['noint'], linewidth=2,label=r'noint($\tau_{ocn}$=360d)'     )#, marker='o', markersize=5, markeredgewidth=2, markeredgecolor=col['long']['noint'], markerfacecolor='w')
ax.patch.set_alpha(1.0)
ax.axis('off')

legend = plt.legend(bbox_to_anchor=(-0.1,0.5), loc='center left', fontsize=12, frameon=True, fancybox=False, ncol=1, borderaxespad=0.5)
legend.get_frame().set_linewidth(1.)
legend.get_frame().set_edgecolor('k')
d1.set_visible(False)
d2.set_visible(False)
d3.set_visible(False)
d4.set_visible(False)
d5.set_visible(False)
d6.set_visible(False)
            

plt.tight_layout()
plt.savefig('figure_3.pdf',dpi=300)
plt.close()



    
