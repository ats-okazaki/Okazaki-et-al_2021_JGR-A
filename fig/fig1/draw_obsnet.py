import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
import speedy.grads as grads
import string
from matplotlib.colors import LogNorm
plt.switch_backend('agg')


resList = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
vmin=0.01
vmax=5
alpha=0.5


fig = plt.figure(figsize=(11,6.5),dpi=200)
gs = gridspec.GridSpec(3,6)

#
# observation points
#
plt.subplot(gs[0, 1:3])
# observation point
ctlfile = "/data/user/speedy_ver41.5/work/obsnet/obsnet_flag.ctl"
o = grads.grads(ctlfile)
X, Y = np.meshgrid(o.lon, o.lat)
X = np.where( X>180., X-360., X)
# plot basemap
m = Basemap(projection = 'cyl', lon_0 = 0, resolution = 'c')
m.drawcoastlines()
# plot obs point
obspointsat = o.read('temp',o.date[0])
m.scatter( X[obspointsat[0,:,:]==1.], Y[obspointsat[0,:,:]==1.], s=20, marker='x', c='r', alpha=1.0, label='SAT',linewidths=0.3, facecolors='None')
obspointsst = o.read('sst',o.date[0])
m.scatter( X[obspointsst[:,:]==1.],   Y[obspointsst[:,:]==1.],   s=20, marker='x', c='b', alpha=1.0, label='SST',linewidths=0.3, facecolors='None')
# legend
plt.legend(loc='center right',bbox_to_anchor=(-0.01,0.5),fontsize=15,frameon=True,fancybox=False)
# title
plt.title('('+string.ascii_lowercase[0]+')',loc='left',fontsize=15)
plt.title('Observation',loc='center',fontsize=15)


#
# error size
#
for r,res in enumerate(resList) :
    
    row = (r+2) // 3
    if r == 0 :
        col = 3
    else :
        col = (r-1)*2%6


    ax = plt.subplot(gs[row,col:col+2])
    ctlfile = '../../nature_icsea-2/obs/ctl/std' +res+ '_snr2.ctl'
    s = grads.grads(ctlfile)
    # plot basemap
    m = Basemap(projection = 'cyl', lon_0 = 0, resolution = 'c', ax = ax)
    m.drawcoastlines()
    # sst
    dat = s.read('sst',s.date[0])
    im = ax.scatter( X[obspointsst[:,:]==1.], Y[obspointsst[:,:]==1.],  c=dat[obspointsst[:,:]==1.],  s=20, marker='o', cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax), alpha=alpha)#, linewidths=1, edgecolors='k')
    # sat
    dat = s.read('temp',s.date[0])[0,:,:]
    im = ax.scatter( X[obspointsat[0,:,:]==1.], Y[obspointsat[0,:,:]==1.],  c=dat[obspointsat[0,:,:]==1.],  s=20, marker='o', cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax), alpha=alpha)#, linewidths=1, edgecolors='k')
     # title
    plt.title('('+string.ascii_lowercase[1+r]+')',loc='left',fontsize=15)
    plt.title(r'$\sigma^{o}$ $\tau$='+res,loc='center',fontsize=15)


fig.subplots_adjust(right=0.9,wspace=0.1,hspace=0.2)
cax = fig.add_axes([0.91, 0.25, 0.020, 0.5])
clb = fig.colorbar(im, cax=cax)
clb.ax.set_title('(K)',fontsize=12)

#plt.show()
plt.savefig('obsnet.pdf',dpi=300)
plt.close()

