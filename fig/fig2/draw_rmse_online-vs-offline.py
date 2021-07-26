import numpy as np
import matplotlib.pyplot as plt
import string
import os
from scipy import stats
#plt.switch_backend('agg')


#
# params
#
exp = 'ctrl'
tmax = 99
da_list = ['offline','online','aa','oa','ao','oo','noda']
res_list = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
var_list = ['SAT','SST']
loc_list = ['L1000','L2000','L3000','L4000']
inf_list_off = ['none0.9']
inf_list_on  = ['RTPS0.0','RTPS0.6','RTPS0.9']
inf_list = {'noda':inf_list_off,'offline':inf_list_off, 'online':inf_list_on, 'aa':inf_list_on, 'ao':inf_list_on, 'oa':inf_list_on, 'oo':inf_list_on} 
label = {'offline':'OFFLINE','online':'ONLINE','aa':'ONLINE_AA','ao':'ONLINE_AO','oa':'ONLINE_OA','oo':'ONLINE_OO'} 
col  = {'offline':'royalblue', 'online':'lightcoral', 'aa':'lightgreen', 'ao':'lightgreen', 'oo':'lightgreen', 'oa':'lightgreen'}
hatch = {'offline':'None', 'online':'None', 'aa':'--', 'ao':'xxx', 'oo':'///', 'oa':'ooo'}

#
# vars
#
# vars
rmse_off_all  = np.empty((len(res_list),len(loc_list)*len(inf_list_off),len(var_list),tmax),dtype=float) 
rmse_on_all   = np.empty((len(res_list),len(loc_list)*len(inf_list_on) ,len(var_list),tmax),dtype=float) 
rmse_nd_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax),dtype=float) # var,loc*inf,time
rmse_off_save = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax),dtype=float) # var,loc*inf,time
rmse_on_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_on),tmax),dtype=float) # var,loc*inf,time
rmse_aa_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_on),tmax),dtype=float) # var,loc*inf,time
rmse_ao_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_on),tmax),dtype=float) # var,loc*inf,time
rmse_oa_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_on),tmax),dtype=float) # var,loc*inf,time
rmse_oo_save  = np.empty((len(var_list),len(loc_list)*len(inf_list_on),tmax),dtype=float) # var,loc*inf,time
rmse_nd  = np.empty((len(res_list),len(loc_list)*len(inf_list_off),len(var_list)),dtype=float) # tres,loc*inf,var
rmse_off = np.empty((len(res_list),len(loc_list)*len(inf_list_off),len(var_list)),dtype=float) # tres,loc*inf,var
rmse_on  = np.empty((len(res_list),len(loc_list)*len(inf_list_on),len(var_list)),dtype=float) # res,loc*inf,var
rmse_aa  = np.empty((len(res_list),len(loc_list)*len(inf_list_on),len(var_list)),dtype=float) # res,loc*inf,var
rmse_ao  = np.empty((len(res_list),len(loc_list)*len(inf_list_on),len(var_list)),dtype=float) # res,loc*inf,var
rmse_oa  = np.empty((len(res_list),len(loc_list)*len(inf_list_on),len(var_list)),dtype=float) # res,loc*inf,var
rmse_oo  = np.empty((len(res_list),len(loc_list)*len(inf_list_on),len(var_list)),dtype=float) # res,loc*inf,var
rmse_save = {'noda':rmse_nd_save, 'offline':rmse_off_save, 'online':rmse_on_save, 'aa':rmse_aa_save, 'ao':rmse_ao_save, 'oa' : rmse_oa_save, 'oo' : rmse_oo_save}
rmse      = {'noda':rmse_nd,      'offline':rmse_off,      'online':rmse_on     , 'aa':rmse_aa     , 'ao':rmse_ao     , 'oa' : rmse_oa     , 'oo' : rmse_oo     }

#
# read data
#
for k,da in enumerate(da_list) :

    for j,res in enumerate(res_list) :

        for l1,roi in enumerate(loc_list) :
            for l2,inf in enumerate(inf_list[da]) : 

                l = len(inf_list[da]) * l1 + l2

                for i,var in enumerate(['temp','sst']) :

                    npz_file_name = './data/'+exp+'_'+da+'_'+res+'_'+roi+'.'+inf+'_'+var+'.npz'
                    npz = np.load( npz_file_name )
                    rmse_save[da][i,l,:] = npz['arr_0']
                    rmse[da][j,l,i] = np.nanmean(rmse_save[da][i,l,:], axis=0)
                    # for statistical test
                    if da == 'online' :
                        rmse_on_all[j,l,i,:] = rmse_save[da][i,l,:]
                    elif da == 'offline' :
                        rmse_off_all[j,l,i,:] = rmse_save[da][i,l,:]

#
# plot
#

fig = plt.figure(figsize=(9,8),dpi=200)

for i,var in enumerate(var_list) :

    ax = fig.add_subplot(len(var_list),1,i+1)

    for j,res in enumerate(res_list) :

        for k,da in enumerate(da_list) :

            if da == 'noda' :
                continue

            left = 8 * j + k + 1
            if da == 'online' or da == 'offline' :
                l_am = np.argmin(rmse[da][j,:,i],axis=0)
            if da == 'offline' :
                l_am_save = l_am
            data = rmse[da][j,l_am,i] / rmse['noda'][j,0,i] - 1.
            #print(var,res,da,l_am,data)
            if da == 'online' or da == 'offline' or da == 'aa' :
                plt.bar(left, data, color=col[da], edgecolor='k', label=label[da])
            else :
                plt.bar(left, data, color=col[da], edgecolor='k', hatch=hatch[da], label=label[da])

            # statistical test
            if da == 'online' :
                t, p = stats.ttest_ind( rmse_on_all[j,l_am,i,5:], rmse_off_all[j,l_am_save,i,5:])
                if p < 0.05 :
                    ax.text(left, 0.25, '*', ha='center', va='center',size=9)#, style='italic')#transform=ax.transAxes)

            plt.ylim(-0.8,0.3)
            plt.xlim(0,56)

        if j == 0 :
            plt.legend(frameon=True,fancybox=False,bbox_to_anchor=(1,0),loc='lower right',borderaxespad=0.5,ncol=3)

        if i == 0 and j == 1 :
            ax.annotate('OFFLINE',
                        xy=(9,-0.33),xycoords='data',
                        xytext=(-50,-40),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=.1"))
            ax.annotate('ONLINE',
                        xy=(10,-0.39),xycoords='data',
                        xytext=(-30,-50),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=.1"))
            ax.annotate('AA',
                        xy=(11,0.02),xycoords='data',
                        xytext=(-30,20),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=-.1"))
            ax.annotate('OA',
                        xy=(12,0.02),xycoords='data',
                        xytext=(-30,30),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=-.1"))
            ax.annotate('AO',
                        xy=(13,0.02),xycoords='data',
                        xytext=(20,30),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=.1"))
            ax.annotate('OO',
                        xy=(14,0.02),xycoords='data',
                        xytext=(30,20),textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                        connectionstyle="arc3,rad=.1"))
            ax.annotate('impact of\nanalyses',
                        xy=(11.5,-0.42),xycoords='data',
                        xytext=(20,-54),textcoords='offset points',
                        arrowprops=dict(arrowstyle="-[,widthB=0.7,lengthB=0.5",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10"))
            ax.annotate('impact of\nforecasts',
                        xy=(13.5,-0.05),xycoords='data',
                        xytext=(20,-89),textcoords='offset points',
                        arrowprops=dict(arrowstyle="-[,widthB=0.7,lengthB=0.5",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10"))

        if i == 1 and j == 1 :
            ax.annotate('impact of\nforecasts',
                        xy=(11.5,-0.18),xycoords='data',
                        xytext=(20,-70),textcoords='offset points',
                        arrowprops=dict(arrowstyle="-[,widthB=0.7,lengthB=0.5",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10"))
            ax.annotate('impact of\nanalyses',
                        xy=(13.5,-0.08),xycoords='data',
                        xytext=(20,-50),textcoords='offset points',
                        arrowprops=dict(arrowstyle="-[,widthB=0.7,lengthB=0.5",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10"))

    ax.set_xticklabels(res_list)
    ax.set_xticks([4,12,20,28,36,44,52])
    plt.title('('+string.ascii_lowercase[i]+') '+var, loc='center')
    plt.ylabel(r'Normalized $\Delta$RMSE')
    plt.hlines(0,0,56,linestyles='solid')

#plt.tight_layout()
#plt.show()
plt.savefig('figure_2.pdf',dpi=300)
#plt.savefig('tmp.png')
plt.close()
