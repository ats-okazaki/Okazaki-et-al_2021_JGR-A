import os
import numpy as np
from scipy import stats
import speedy.grads as grads
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import string


colors = plt.cm.rainbow(np.linspace(0, 1, 20))
# color bar
cmap = plt.cm.PuOr_r
#cmap = plt.cm.coolwarm
cmap.set_bad('gainsboro',1.0)

# init
exp_list = ['ctrl']
res_list = ['1dy','3dy','10dy','1mo']
method_list = ['tau','offline','noda','aa','oa','ao','oo']
var_list = ['temp','sst']
loc_list = ['L1000','L2000','L3000','L4000']
inf_list_off = ['none0.9']
inf_list_tau = ['RTPS0.0','RTPS0.6','RTPS0.9']
inf_list = {'offline':inf_list_off,'tau':inf_list_tau,'noda':inf_list_off,'aa':inf_list_tau,'oa':inf_list_tau,'ao':inf_list_tau,'oo':inf_list_tau }
dir_list = {'ctrl':'','long':''}
dir_list['ctrl'] = {'truth':'','offline':'','tau':'','noda':'','aa':'','oa':'','ao':'','oo':''}
dir_list['long'] = {'truth':'','offline':'','tau':'','noda':'','aa':'','oa':'','ao':'','oo':''}
dir_list['ctrl']['truth'] = '/data/user/speedy_ver41.5/work/nature_icsea-2/obs/ctl'
dir_list['ctrl']['offline'] = '/data/user/speedy_ver41.5/work/offline_all_snr2/'
dir_list['ctrl']['noda'] = '/data/user/speedy_ver41.5/work/offline_all_snr2/'
dir_list['ctrl']['tau'] = '/data/user/speedy_ver41.5/work/tau_all_snr2/'
dir_list['ctrl']['aa'] = '/data/user/speedy_ver41.5/work/tau_all_snr2_aa/'
dir_list['ctrl']['oa'] = '/data/user/speedy_ver41.5/work/tau_all_snr2_oa/'
dir_list['ctrl']['ao'] = '/data/user/speedy_ver41.5/work/tau_all_snr2_ao/'
dir_list['ctrl']['oo'] = '/data/user/speedy_ver41.5/work/tau_all_snr2_oo/'
var_title = {'temp':'SAT','sst':'SST'}
method_title = {'offline':'OFFLINE','tau':'ONLINE','aa':'ONLINE_AA','ao':'ONLINE_AO','oa':'ONLINE_OA','oo':'ONLINE_OO'}
tmax = 99


for exp in exp_list :

    fig = plt.figure(figsize=(10,10))

    for i,res in enumerate(res_list) :

        if res[-2:] == 'mo' :
            dtunit = 'M' 
            dtunit_title = 'mon'
        elif res[-2:] == 'dy' :
            dtunit = 'D'
            dtunit_title = 'day'
        dt = res[:-2]


        #
        # read data
        #
        rmse_tau = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax,48,96),dtype=float)
        rmse_off = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax,48,96),dtype=float)
        rmse_nod = np.empty((len(var_list),len(loc_list)*len(inf_list_off),tmax,48,96),dtype=float)
        rmse_aa = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax,48,96),dtype=float)
        rmse_oa = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax,48,96),dtype=float)
        rmse_ao = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax,48,96),dtype=float)
        rmse_oo = np.empty((len(var_list),len(loc_list)*len(inf_list_tau),tmax,48,96),dtype=float)
        rmse_tau_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_tau)),dtype=float)
        rmse_off_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_off)),dtype=float)
        rmse_nod_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_off)),dtype=float)
        rmse_aa_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_tau)),dtype=float)
        rmse_oa_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_tau)),dtype=float)
        rmse_ao_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_tau)),dtype=float)
        rmse_oo_glb = np.empty((len(var_list),len(loc_list)*len(inf_list_tau)),dtype=float)
        rmse = {'offline':rmse_off,'tau':rmse_tau,'noda':rmse_nod,'aa':rmse_aa,'oa':rmse_oa,'ao':rmse_ao,'oo':rmse_oo}
        rmse_glb = {'offline':rmse_off_glb,'tau':rmse_tau_glb,'noda':rmse_nod_glb,'aa':rmse_aa_glb,'oa':rmse_oa_glb,'ao':rmse_ao_glb,'oo':rmse_oo_glb}


        ctlfile_truth = dir_list[exp]['truth'] +'/attm' + res + "_snr10.ctl"
        t = grads.grads(ctlfile_truth)
        X, Y = np.meshgrid(t.lon,t.lat)
        weights=abs(np.cos(Y*np.pi/180.)).reshape(t.nx*t.ny)
        time_start = t.date[0]
        time_end = t.date[99]
        dates = np.arange(time_start,time_end,int(dt),dtype='datetime64['+dtunit+']')


        for method in method_list :

            npz_file_name = './data/'+ exp +'_'+ method +'_'+ res +'_map.npz'

            if os.path.isfile( npz_file_name ) :

                npz = np.load( npz_file_name )
                rmse[method] = npz['arr_0']
                rmse_glb[method] = npz['arr_1']

            else :


                for j1,roi in enumerate(loc_list) :

                    for j2,inf in enumerate(inf_list[method]) : 

                        subexp = roi +'.'+ inf
                        j = len(inf_list[method]) * j1 + j2

                        for k,var in enumerate(var_list) :

                            if method == 'noda' :
                                ctlfile_mean = dir_list[exp][method] +'/'+ res +'/'+ subexp +'/ctl/gues_mean.ctl'
                            else :
                                ctlfile_mean = dir_list[exp][method] +'/'+ res +'/'+ subexp +'/ctl/anal_mean.ctl'
                            
                            a = grads.grads(ctlfile_mean)
    
                            tmpa = np.empty((np.shape(dates)[0],a.ny,a.nx),dtype=float)
                            tmpt = np.empty((np.shape(dates)[0],a.ny,a.nx),dtype=float)
                        
                        
                            if var in t.var3d_mylist :
                                for it, date in enumerate(dates) :
                                    tmpt[it,:,:] = t.read(var, date)[0,:,:]
                                    tmpa[it,:,:] = a.read(var, date)[0,:,:]
                            elif var in t.var2d_mylist : 
                                for it, date in enumerate(dates) :
                                    tmpt[it,:,:] = t.read(var, date)[:,:]
                                    tmpa[it,:,:] = a.read(var, date)[:,:]
    
                            weights=abs(np.cos(Y*np.pi/180.)).reshape(t.nx*t.ny)
                            if var == 'sst' or var == 'lst' :
                                flag = tmpt[0,:,:].reshape(a.nx*a.ny)
                                weights = np.where((flag < 200) | (flag > 400), 0, weights)
                                if var == 'sst' :
                                    weights = np.where(abs(Y.reshape(a.nx*a.ny))>60., 0, weights)
                            
                            # rmse
                            rmse[method][k,j,:,:,:] = np.sqrt(np.square(tmpa-tmpt))
                            rmse_glb[method][k,j] = np.sqrt(np.average(np.square(tmpa-tmpt).reshape(np.shape(tmpa)[0],a.nx*a.ny),weights=weights, axis=1)).mean()

                # save data
                np.savez( npz_file_name, rmse[method], rmse_glb[method] )

        #
        # draw
        #
        var = 'sst'
        k = 1
        #var = 'temp'
        #k = 0

        # ice mask
        icefile = grads.grads('../../../data/bc/t30/clim//seaice_7908clim.t30.sea.ctl')
        if var == 'sst' :
            ice = icefile.read('sice',t=0) / 12
            for it in np.arange(1,12) :
                ice = ice + icefile.read('sice', t=it) / 12

        for mm,method in enumerate(['offline','tau','aa','oa','ao','oo']) :

            ax = fig.add_subplot(len(['offline','tau','aa','oa','ao','oo']),len(res_list),i+mm*len(res_list)+1) 

            if method == 'offline' or method == 'tau' :
                j_best = np.nanargmin(rmse_glb[method][k,:])

            drmse_mean = ( rmse[method][k,j_best,:,:,:].mean(axis=0) / rmse['noda'][k,0,:,:,:].mean(axis=0) - 1. ) * 100
            # plot drmse
            m = Basemap(projection = 'cyl', lon_0 = 180, resolution = 'c', ax = ax)
            m.drawcoastlines( color = 'k', linewidth = 0.5)
            if var == 'sst' :
                drmse_mean = np.where( drmse_mean == 0., np.nan, drmse_mean )
                drmse_mean = ma.masked_where(np.isnan(drmse_mean),drmse_mean)
            clevs = [-80,-40,-20,-10,-5,0,5,10,20,40,80]
            im = ax.contourf( X, Y, drmse_mean, clevs, latlon=True, extend='both', cmap=cmap)

            # observation point
            ctlfile = "/home/okazaki/data15/speedy_ver41.5/work/obsnet/obsnet_flag.ctl"
            o = grads.grads(ctlfile)
            dat = o.read('temp',o.date[0])
            m.scatter( X[dat[0,:,:]==1.], Y[dat[0,:,:]==1.], s=5, marker='.', c='w', alpha=0.9,linewidths=0.5, edgecolors='lawngreen')
            dat = o.read('sst',o.date[0])
            m.scatter( X[dat[:,:]==1.], Y[dat[:,:]==1.], s=5, marker='.', c='w', alpha=0.9,linewidths=0.5, edgecolors='lawngreen')

            # sea ice mask
            if var == 'sst' :
                ax.pcolor( X, Y, ma.masked_less(ice,0.3), hatch='///', alpha=0.0, lw=0.1)


            # title
            plt.title(str(dt) +'-'+ dtunit_title, loc='center', fontsize=12)
            plt.title('('+ string.ascii_lowercase[i+mm*len(res_list)] +')', loc='left', fontsize=12)
            if i == 0 :
                plt.ylabel(method_title[method], fontsize=12)


    # colorbar
    divider = make_axes_locatable(ax)
    fig.subplots_adjust(bottom=0.15, hspace=0.07, wspace=0.05)
    cax = fig.add_axes([0.15, 0.12, 0.7, 0.02])
    clb = fig.colorbar(im, orientation='horizontal', cax=cax)
    clb.ax.set_title('(%)', fontsize=12, y=-2.3)
    plt.savefig('figure_5v2.pdf', dpi=300)
    plt.close()
    plt.clf()

