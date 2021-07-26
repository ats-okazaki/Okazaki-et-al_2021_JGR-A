import os
import numpy as np
import numpy.ma as ma
import scipy.stats as st
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import speedy.grads as grads
import speedy.geomet as geo
import speedy.noleap as nl
import string
import random
import matplotlib.cm as cm

# output
pngFileName = 'cov_sat-sst_dist.png'

# parameters
var_list = ['temp','sst']
res_list = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
exp_list = ['ctrl']
nx = 96
ny = 48
nz = 8
year_str = 1900
year_end = 1950
recl = 4 *(((nx * ny + 2) * nz * 9) + ((nx * ny + 2) * 32 ))
nSample = 100
length_list = np.arange(0,5500,500)
lengthM_list = np.arange(250,5250,500)

# data
dirref = {'ctrl': {'D':'', 'M':''}, 'long':''}
dirref['ctrl']['D'] = '/data/user/speedy_ver41.5/work/nature_icsea-2/out_daily'
dirref['ctrl']['M'] = '/data/user/speedy_ver41.5/work/nature_icsea-2/out_monthly'
vrec_temp = 4 * ((nx * ny + 2) * nz )
vrec_sst  = 4 * ((nx * ny + 2) * nz * 9 + (nx * ny + 2) * 15 )
vrec = {'temp' : vrec_temp, 'sst': vrec_sst }
nstp_clim = {'D': 365, 'M':12}
nstp = {'D': 365*9, 'M':1800}

# defile colors
colors = cm.rainbow(np.linspace(0, 1, len(res_list)))

# distance
distFileName = 'dist_speedy.npy'
if os.path.isfile(distFileName) :
    dist = np.load( distFileName )
else :
    n = grads.grads('/data/user/speedy_ver41.5/work/nature_icsea-2/out_monthly/attm000.ctl')
    X, Y = np.meshgrid(n.lon, n.lat)
    Xf, Yf = X.flatten(), Y.flatten()
    dist = np.empty((n.ny*n.nx,n.ny*n.nx),dtype=float)
    for ii in range(len(Xf)) :
        ii_x = ii % nx
        if ii_x == 0 :
            ii_base = ii
            for jj in range(len(Xf)) :
                dist[ii,jj] = geo.distance( Xf[ii], Yf[ii], Xf[jj], Yf[jj] )
        else :
            tmp = np.roll( dist[ii_base,:].reshape(ny,nx), ii_x, axis=1 )
            dist[ii,:] = tmp.reshape(ny*nx)
    np.save( distFileName, dist )
    exit()

# dfine vars
per50 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per25 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per75 = np.empty((len(res_list),len(length_list)-1),dtype=float)


# covariance
for i,exp in enumerate(exp_list) :

    #fig = plt.figure(dpi=200)

    for j, res in enumerate(res_list) :

        #
        # time resolution
        #
        tres = int(res[:-2])
        tres_unit = res[-2].upper()

        #
        # read data
        #
        npy_file_name = 'cov_sat-sst_'+ exp +'_'+ res +'.npy'

        if os.path.isfile( npy_file_name ) :

            cov = np.load( npy_file_name )

        else :

            # data
            dat = np.empty((2,nstp[tres_unit],ny*nx),dtype=float)
            sat = np.empty((nstp[tres_unit],ny*nx),dtype=float)
            sst = np.empty((nstp[tres_unit],ny*nx),dtype=float)
            a = grads.grads(dirref[exp][tres_unit] + '/attm000.ctl')
            X,Y = np.meshgrid(a.lon,a.lat)
            time_start = a.date[0]
            #time_end = a.date[tres*nstp[tres_unit]]
            time_end = a.date[nstp[tres_unit]]
            dates = nl.arange(time_start,time_end,1,tres_unit)
            for it, date in enumerate(dates) :
                dat[0,it,:] = a.read('temp',t=it)[0,:,:].reshape(ny*nx)
                dat[1,it,:] = a.read('sst', t=it)[:,:].reshape(ny*nx)
            dat = np.ma.masked_where( abs(dat) > 500., dat )


            icefile = grads.grads('../../../data/bc/t30/clim//seaice_7908clim.t30.sea.ctl')
            ice = icefile.read('sice',t=0) / 12
            for it in np.arange(1,12) :
                ice = ice + icefile.read('sice', t=it) / 12
            for it, date in enumerate(dates) :
                dat[1,it,:] = np.ma.masked_where( ice.flatten() > 0., dat[1,it,:])


            # climatology
            clim = np.zeros((2,nstp_clim[tres_unit],ny*nx),dtype=float)
            climSmth = np.zeros((2,nstp_clim[tres_unit],ny*nx),dtype=float)

            for year in np.arange(year_str,year_end) :

                datafile = dirref[exp][tres_unit] + '/attm000_' +str(year)+ '.grd'

                with open(datafile,'rb') as f :

                    for k,var in enumerate(var_list) :

                        for it in np.arange(nstp_clim[tres_unit]) :
                            rec = recl * it + vrec[var] + 4
                            f.seek( rec )
                            clim[k,it,:] += np.fromfile( f, '>f', nx*ny) / float(year_end - year_str)



            if tres_unit == 'D' :
                #weight = np.ones(11)/11.
                #climSmth = np.convolve( clim, weight, mode='same' )
                dt = 10
                for it in np.arange(dt) :
                    climSmth[:,it,:] = (np.sum(clim[:,-dt+it:,:],axis=1) + np.sum(clim[:,0:dt+1+it,:],axis=1))/(2*dt+1)
                for it in np.arange(dt,365-dt+1) :
                    climSmth[:,it,:] = clim[:,it-dt:it+dt,:].mean(axis=1)
                for it in np.arange(365-dt+1,365) :
                    climSmth[:,it,:] = (np.sum(clim[:,it-dt:,:],axis=1) + np.sum(clim[:,0:it+dt+1-365,:],axis=1))/(2*dt+1)
            elif tres_unit == 'M' :
                climSmth = clim


            # anomaly
            for it,date in enumerate(dates) :
                dat[:,it,:] = dat[:,it,:] - climSmth[:,(it%nstp_clim[tres_unit]),:]


            # average for res
            nt = int(nstp[tres_unit]/tres)
            print(nt)
            dat2 = np.empty((2,nt,ny*nx),dtype=float)
            for it in np.arange(nt) :
                dat2[:,it,:] = np.mean(dat[:,it*tres:(it+1)*tres,:],axis=1)


            # random sampling
            index_all = range(nt)
            index = random.sample(index_all,nSample)
            samples = np.empty((2,nSample,ny*nx),dtype=float)
            for jj,ii in enumerate(index) :
                if j >= nSample :
                    break
                samples[:,jj,:] = dat2[:,ii,:]

            '''
            plt.plot(np.arange(100),samples[0,:,3060])
            plt.plot(np.arange(100),samples[1,:,3060])
            plt.show()
            '''



            # correlation
            ave_sat = np.mean( samples[0,:,:], axis=0 )
            std_sat = np.std ( samples[0,:,:], axis=0 )
            ave_sst = np.mean( samples[1,:,:], axis=0 )
            std_sst = np.std ( samples[1,:,:], axis=0 )
            '''
            plt.imshow(std_sat.reshape(ny,nx))
            plt.colorbar()
            plt.show()
            plt.clf()
            plt.imshow(ave_sat.reshape(ny,nx))
            plt.colorbar()
            plt.show()
            plt.clf()
            '''
            #for it in np.arange(nstp[tres_unit]) :
            for it in np.arange(nSample) :
                sat[it,:] = (samples[0,it,:] - ave_sat[:]) / std_sat[:]
                sst[it,:] = (samples[1,it,:] - ave_sst[:]) / std_sst[:]

            '''
            plt.plot(np.arange(nstp[tres_unit]),sst[:,96*10+48])
            plt.plot(np.arange(nstp[tres_unit]),sst[:,96*20+48])
            plt.plot(np.arange(nstp[tres_unit]),sst[:,96*30+48])
            plt.plot(np.arange(nstp[tres_unit]),sst[:,96*40+48])
            plt.show()
            '''

            #cov = np.dot(sat[:,:].T,sat[:,:]) / float(nSample)
            cov = np.dot(sat[:,:].T,sst[:,:]) / float(nSample)

            '''
            ij =    3040
            for ij in np.arange(6*96,48*96,96+12) :
                plt.clf()
                plt.imshow(cov[ij,:].reshape(ny,nx)[:,:], cmap=cm,vmin=-1,vmax=1)
                plt.plot(ij%96,ij//96,'o',c='k')
                plt.colorbar()
                plt.show()
            '''

            np.save( npy_file_name, cov )


        #----------------------------------------------------------------

        
        #
        # distance
        #
        '''
        plt.imshow( np.ma.masked_where(dist[3030,:]>3000, cov[3030,:]).reshape(ny,nx) )
        plt.colorbar()
        plt.show()
        '''
        print(lengthM_list)
        #if j > 3 :
        if j == 1 or j == 2 or j == 4 :
            continue
        for l,length in enumerate(length_list) :
            if l == len(length_list)-1 :
                break
            tmpSample = cov[ (dist >= length_list[l]) & (dist <= length_list[l+1]) ]
            per50[j,l] = np.nanmedian(tmpSample)
            per75[j,l] = np.nanpercentile(tmpSample,75)
            per25[j,l] = np.nanpercentile(tmpSample,25)
            print(j,l,per50[j,l])

        plt.fill_between( lengthM_list[:], per25[j,:], per75[j,:], facecolor=colors[j], alpha=0.3 )
        plt.plot( lengthM_list[:], per50[j,:], color=colors[j], linewidth=3, label=res )

    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation Coefficient')
    plt.title('SPEEDY-SOM')
    plt.xlim(0,5000)
    plt.ylim(-0.2,1.0)
    legend = plt.legend(frameon=True, fancybox=False)
    #legend = plt.legend(bbox_to_anchor=(0.5,-1.0), loc='center', fontsize=18, frameon=True, fancybox=False)
    legend.get_frame().set_linewidth(1.)
    legend.get_frame().set_edgecolor('k')

    plt.savefig(pngFileName)
