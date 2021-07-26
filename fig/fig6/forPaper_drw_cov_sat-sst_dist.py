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

# output
pngFileName = 'cov_sat-sst_dist.png'
pdfFileName = 'cov_sat-sst_dist.pdf'

# parameters
var_list = ['temp','sst']
res_list = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
exp_list = ['ctrl']
length_list = np.arange(0,5500,500)
lengthM_list = np.arange(250,5250,500)
mdl_list = ['speedy','miroc5']

#---------------- SPEEDY ---------------------

#---------------- MIROC5 ------------------------
exp = 'piControl'
model = 'MIROC5'
ensemble = 'r1i1p1'


# dfine vars
per50 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per25 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per75 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per99 = np.empty((len(res_list),len(length_list)-1),dtype=float)
per01 = np.empty((len(res_list),len(length_list)-1),dtype=float)

# stats
nSample = 100
df = nSample - 2
t_thresh = st.t.ppf(0.95,df)
r_thresh = np.sqrt(t_thresh * t_thresh / (df + t_thresh * t_thresh))

# covariance
for i,exp in enumerate(exp_list) :

    fig = plt.figure(figsize=(13,5))

    for Nmdl, mdl in enumerate(mdl_list) :

        for j, res in enumerate(res_list) :

            #
            # time resolution
            #
            tres = int(res[:-2])
            tres_unit = res[-2].upper()
            if tres_unit == 'D' :
                tres_unit_title = 'day'
            elif tres_unit == 'M' :
                tres_unit_title = 'mon'

            #
            # read data
            #
            if mdl == 'speedy' :
                npyFileName = 'cov_sat-sst_'+ exp +'_'+ res +'.npy'
                distFileName = 'dist_speedy.npy'

            elif mdl == 'miroc5' :
                npyFileName = 'cov_sat-sst_MIROC5_piControl_r1i1p1_'+ res +'.npy'
                distFileName = 'dist_miroc5.npy'
                print(npyFileName)
                print(distFileName)

            if os.path.isfile( npyFileName ) and os.path.isfile(distFileName) :
                cov = np.load( npyFileName )
                dist = np.load( distFileName )
            else :
                print('data file does not exist!')
                exit()


            #----------------------------------------------------------------
            
            #
            # plot
            #
            '''
            plt.imshow( np.ma.masked_where(dist[3030,:]>3000, cov[3030,:]).reshape(ny,nx) )
            plt.colorbar()
            plt.show()
            '''
            nPanel = Nmdl * len(res_list) + j + 1
            ax = fig.add_subplot(2,len(res_list),nPanel)
            #if j > 3 :
            for l,length in enumerate(length_list) :
                if l == len(length_list)-1 :
                    break
                tmpSample = cov[ (dist >= length_list[l]) & (dist <= length_list[l+1]) ]
                per50[j,l] = np.nanmedian(tmpSample)
                per99[j,l] = np.nanpercentile(tmpSample,99)
                per75[j,l] = np.nanpercentile(tmpSample,75)
                per25[j,l] = np.nanpercentile(tmpSample,25)
                per01[j,l] = np.nanpercentile(tmpSample, 1)

            plt.fill_between( lengthM_list[:], per01[j,:], per99[j,:], facecolor='gainsboro')#, alpha=0.3 )
            plt.fill_between( lengthM_list[:], per25[j,:], per75[j,:], facecolor='darkgray')#, alpha=0.3 )
            plt.plot( lengthM_list[:], per50[j,:], color='k', linewidth=3, label=res )
            plt.hlines(0,0,5000,linestyles='solid',linewidths=1)
            plt.hlines(r_thresh,0,5000,linestyles='dashed',linewidths=2)
            plt.xlim(0,5000)
            plt.ylim(-0.2,1.0)
            if j == 0 :
                plt.ylabel('Correlation')
            plt.xlabel('Distance [km]')
            plt.title('('+ string.ascii_lowercase[nPanel-1] +') '+ str(tres) +'-'+ tres_unit_title)
            plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

    plt.subplots_adjust(hspace=0.5,wspace=0.33)
    plt.savefig(pdfFileName, dpi=300)
