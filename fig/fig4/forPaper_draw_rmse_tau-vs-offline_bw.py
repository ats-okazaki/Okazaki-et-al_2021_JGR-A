import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import speedy.grads as grads
import os
from scipy import stats
import string


# params
method_list = ['offline','tau']
exp_list = ['ctrl','long']
res_list = ["1dy","3dy","10dy","1mo","3mo","6mo","12mo"]
var_list = ['temp','sst']
loc_list = ['L1000','L2000','L3000','L4000']
inf_list_off = ['none0.9']
inf_list_tau = ['RTPS0.0','RTPS0.6','RTPS0.9']
inf_list = {'offline' : inf_list_off, 'tau' : inf_list_tau} 
xlabels = {'ctrl':r'$\tau_{ocn}$=90d', 'long':r'$\tau_{ocn}=360d'}
dir_list = {'ctrl':{'truth':'','offline':'','tau':''},'long':{'truth':'','offline':'','tau':''}}
dir_list['ctrl']['truth'] = '/data/user/speedy_ver41.5/work/nature_icsea-2/obs/ctl'
dir_list['ctrl']['offline'] = '/data/user/speedy_ver41.5/work/offline_all_snr2/'
dir_list['ctrl']['tau'] = '/data/user/speedy_ver41.5/work/tau_all_snr2/'
dir_list['long']['truth'] = '/data/user/speedy_ver41.5/work/nature_icsea-2-long/obs/ctl'
dir_list['long']['offline'] = '/data/user/speedy_ver41.5/work/offline_all_long_snr2'
dir_list['long']['tau'] = '/data/user/speedy_ver41.5/work/tau_all_long_snr2/'
color = {'offline' : {'none0.9' : 'k'}, 'tau' : {'RTPS0.0' : 'r', 'RTPS0.6' : 'g', 'RTPS0.9' : 'b'}}
var_title = {'temp':'SAT', 'sst':'SST'}
col = {'ctrl':'lightcoral', 'long':'limegreen'}
tmax = 99

# init
rmse_off = np.empty((len(res_list),len(loc_list)*len(inf_list_off),len(var_list)),dtype=float) # tres,loc,var
rmse_tau = np.empty((len(res_list),len(loc_list)*len(inf_list_tau),len(var_list)),dtype=float) # tres,loc,var
rmse_off_all = np.empty((len(res_list),len(loc_list)*len(inf_list_off),len(var_list),tmax),dtype=float) # tres,loc,var,time
rmse_tau_all = np.empty((len(res_list),len(loc_list)*len(inf_list_tau),len(var_list),tmax),dtype=float) # tres,loc,var,time
rmse_off_save = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax),dtype=float) # var,time
rmse_tau_save = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax),dtype=float) # var,time
sprd_off_save = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax),dtype=float) # var,time
sprd_tau_save = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax),dtype=float) # var,time
rmse      = {'offline' : rmse_off,      'tau' : rmse_tau}
rmse_all  = {'offline' : rmse_off_all,  'tau' : rmse_tau_all}
rmse_save = {'offline' : rmse_off_save, 'tau' : rmse_tau_save}
sprd_save = {'offline' : sprd_off_save, 'tau' : sprd_tau_save}


fig = plt.figure(dpi=200)

for exp in exp_list :

    for i,res in enumerate(res_list) :

        if res[-2:] == 'mo' :
            dtunit = 'M' 
        elif res[-2:] == 'dy' :
            dtunit = 'D'
        dt = res[:-2]

        ctlfile_truth = dir_list[exp]['truth'] +'/attm' + res + "_snr10.ctl"
        t = grads.grads(ctlfile_truth)
        X, Y = np.meshgrid(t.lon, t.lat)
        weights=abs(np.cos(Y*np.pi/180.)).reshape(t.nx*t.ny)
        time_start = t.date[0]
        time_end = t.date[99]
        dates = np.arange(time_start,time_end,int(dt),dtype='datetime64['+dtunit+']')

        #--- read data
        for method in method_list :


            npz_file_name = './data/'+ exp +'_'+ method +'_'+ res +'.npz'

            if os.path.isfile( npz_file_name ) :

                npz = np.load( npz_file_name )
                rmse_save[method] = npz['arr_0']
                sprd_save[method] = npz['arr_1']

            else :

                for j1,roi in enumerate(loc_list) :

                    for j2,inf in enumerate(inf_list[method]) : 

                        subexp = roi +'.'+ inf
                        j = len(inf_list[method]) * j1 + j2

                        # read
                        ctlfile_mean = dir_list[exp][method] +'/'+ res +'/'+ subexp +'/ctl/anal_mean.ctl'
                        ctlfile_sprd = dir_list[exp][method] +'/'+ res +'/'+ subexp +'/ctl/anal_sprd.ctl'
                        
                        a = grads.grads(ctlfile_mean)
                        sa = grads.grads(ctlfile_sprd)
    
                        tmpa = np.empty((np.shape(dates)[0],a.ny,a.nx),dtype=float)
                        tmpt = np.empty((np.shape(dates)[0],a.ny,a.nx),dtype=float)
                        tmpsa = np.empty((np.shape(dates)[0],a.ny,a.nx),dtype=float)

                        
                        for k,var in enumerate(var_list) :
                        
                            if var in t.var3d_mylist :
                                for it, date in enumerate(dates) :
                                    tmpt[it,:,:] = t.read(var, date)[0,:,:]
                                    tmpa[it,:,:] = a.read(var, date)[0,:,:]
                                    tmpsa[it,:,:] = sa.read(var, date)[0,:,:]
                            elif var in t.var2d_mylist : 
                                for it, date in enumerate(dates) :
                                    tmpt[it,:,:] = t.read(var, date)[:,:]
                                    tmpa[it,:,:] = a.read(var, date)[:,:]
                                    tmpsa[it,:,:] = sa.read(var, date)[:,:]
    
                            weights=abs(np.cos(Y*np.pi/180.)).reshape(t.nx*t.ny)
                            if var == 'sst' or var == 'lst' :
                                flag = tmpt[0,:,:].reshape(a.nx*a.ny)
                                weights = np.where((flag < 200) | (flag > 400), 0, weights)
                                if var == 'sst' :
                                    weights = np.where(abs(Y.reshape(a.nx*a.ny))>60, 0, weights)
                            
                            #
                            # rmse
                            #
                            tmptmpa = tmpa[:,:,:].reshape( np.shape(tmpa)[0], a.nx * a.ny )
                            tmptmpt = tmpt[:,:,:].reshape( np.shape(tmpa)[0], a.nx * a.ny )
                            rmse_save[method][k,j,:] = np.sqrt(np.average(np.square(tmptmpa-tmptmpt),weights=weights, axis=1))

                            #
                            # sprd
                            #
                            tmptmpa = tmpsa[:,:,:].reshape( np.shape(tmpa)[0], a.nx * a.ny )
                            sprd_save[method][k,j,:] = np.sqrt(np.average(np.square(tmptmpa),weights=weights, axis=1))

                #
                # save data
                #
                np.savez( npz_file_name, rmse_save[method], sprd_save[method] )


            #--- plot rmse & sprd timeseries
            for j1,roi in enumerate(loc_list) :

                for j2,inf in enumerate(inf_list[method]) : 

                    j = len(inf_list[method]) * j1 + j2

                    for k,var in enumerate(var_list) :
        
                        # temporal mean
                        rmse[method][i,j,k] = np.nanmean(rmse_save[method][k,j,:], axis=0)
                        rmse_all[method][i,j,k,:] = rmse_save[method][k,j,:]

    print('----------------')


    #
    # plot data
    # 
    
    for k,var in enumerate(var_list) :

        ax = fig.add_subplot( len(var_list), 1, k+1 )
        ax.set_title( '('+ string.ascii_lowercase[k] +') '+ var_title[var])


        for i,res in enumerate(res_list) :

            data_tau = rmse_tau_all[i,np.argmin(rmse_tau[i,:,k],axis=0),k,:]
            data_off = rmse_off_all[i,np.argmin(rmse_off[i,:,k],axis=0),k,:]

            if exp == 'ctrl' :
                xpos = [i*3+1.1]
            else :
                xpos = [i*3+1.9]

            data = [ data_tau[5:]/np.mean(data_off[5:]) ] # discard spinup
            bp = ax.boxplot(data, 
                            positions=xpos,
                            patch_artist=True, # fill with color
                            showfliers=False,
                            widths=0.6)

            # colors
            plt.setp(bp['medians'],color='k')
            bp['boxes'][0].set_facecolor(col[exp])

            plt.tick_params(labelbottom=True,labelleft=True,labelright=False,labeltop=False)
            ax.set_ylabel('Normalized RMSE')

            ax.set_ylim(0.5,1.3)
            ax.set_xlim(0,21)


            # statistical test
            t, p = stats.ttest_ind(data_tau[5:], data_off[5:]) 

            if p < 0.001 :
                str = '***'
            elif p < 0.01 :
                str = '**'
            elif p < 0.05 : 
                str = '*'
            else :
                str =''

            if exp == 'ctrl' :
                xpos = i*3+1.1
            else :
                xpos = i*3+1.9
            ax.text(xpos, 0.55, str, ha='center', va='center',size=7,fontweight='bold')#, style='italic')#transform=ax.transAxes)
            
        ax.set_xticklabels(res_list)
        ax.set_xticks([1.5,4.5,7.5,10.5,13.5,16.5,19.5])
        plt.hlines(1., 0, 21, 'k', linestyles='dashed', linewidth=2)

        # legend
        if exp == 'ctrl' and k == 0 :
            d1 = plt.bar(100,1,color=col['ctrl'],edgecolor='k',label=r'$\tau_{ocn}$=90d') 
            d2 = plt.bar(200,1,color=col['long'],edgecolor='k',label=r'$\tau_{ocn}$=360d') 
            legend = plt.legend(frameon=True,fancybox=False,bbox_to_anchor=(1,0),loc='lower right',borderaxespad=0.5)
            legend.get_frame().set_linewidth(1.)
            legend.get_frame().set_edgecolor('k')

plt.subplots_adjust(hspace=0.4)
plt.savefig('bw_rmse_tau-offline.pdf',dpi=300)
#plt.savefig('bw_rmse_tau-offline.png')
plt.close()

