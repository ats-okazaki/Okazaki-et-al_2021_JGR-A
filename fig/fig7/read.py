import numpy as np
import glob,os
import matplotlib.pyplot as plt
import speedy.grads as grads


# files & exp. params
expName = "tau_all_snr2"
dataDir = "/data/user/speedy_ver41.5/work/tau_all_snr2/"
rmseDir = "/data/user/speedy_ver41.5/work/fig/online-vs-offline/data"
obsCtlFile = "/data/user/speedy_ver41.5/work/obsnet/obsnet_flag.ctl"
resList = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
locList = ['L1000','L2000','L3000','L4000']
infList = ['RTPS0.0','RTPS0.6','RTPS0.9']

# speedy params
nx = 96
ny = 48 
nv =  2 # temp, sst
nobs = 447 #   1-413 : temp
           # 414-447 : sst

# plot params (map)
vmax = [0.3,0.5,0.5,1.0,1.0,1.0,1.0]
#cmap = 'viridis'
cmap = 'cividis'
c_om = 'lawngreen'
c_ol = 'lawngreen'
s_om = 10

# plot params (bar)
bar_width = 0.3

# define vars
so_mean = np.empty((len(resList),nobs,nv,ny,nx),dtype=float)
idx_argmin = np.empty((len(resList),nv,2),dtype=int)
impact_of_sst = np.empty((len(resList),nv,ny,nx),dtype=float)
impact_of_sat = np.empty((len(resList),nv,ny,nx),dtype=float)

# rmse
for i,res in enumerate(resList) :
    npz_file_name = rmseDir +'/ctrl_tau_'+ res +'_map.npz'
    if os.path.isfile( npz_file_name ) :
        npz = np.load( npz_file_name )
        rmse_glb = npz['arr_1']
        for j in range(nv) :
            idx_argmin[i,j,:] = np.unravel_index(np.nanargmin(rmse_glb[j,:]),(len(locList),len(infList)))

# analysis sensitivy (2d)
npyFileName = './as2d.npy'
if os.path.isfile( npyFileName ) :

    so_mean = np.load(npyFileName)

else :
    for i,res in enumerate(resList) :

        # for sst (j=1)
        loc = locList[idx_argmin[i,0,0]]
        inf = infList[idx_argmin[i,0,1]]

        inFileList = glob.glob(dataDir +'/'+ res +'/'+ loc +'.'+ inf +'/anal/mean/as_*')
        so = np.empty((len(inFileList),nobs,nv,ny,nx),dtype=float)

        for l,inFile in enumerate(inFileList) :

            with open(inFile,'rb') as f :

                so[l,:,:,:,:] = np.fromfile(f, '>f', nobs*ny*nx*nv).reshape(nobs,nv,ny,nx)

        # temporal mean
        so_mean[i,:,:,:,:] = np.absolute(so).mean(axis=0)
        # save
        np.save(npyFileName,so_mean)

# obs points
o = grads.grads(obsCtlFile)
obs_temp = o.read('temp',o.date[0])
obs_sst  = o.read('sst',o.date[0])

#----------------------------------------------------------------------------------------------------
# plot
import cartopy.crs as ccrs
import cartopy.feature as cftr
from matplotlib import gridspec
import string

# x,y,weights
a = grads.grads('../../tau_all_snr2/1dy/L1000.RTPS0.0/ctl/anal_mean.ctl')
flag = a.read('sst',t=0).reshape(ny*nx)

f = grads.grads('../../../data/bc/t30/clim//seaice_7908clim.t30.sea.ctl')
X, Y = np.meshgrid(f.lon,f.lat)
ice = f.read('sice',t=0) / 12
for it in np.arange(1,12) :
    ice = ice + f.read('sice', t=it) / 12
weights = abs(np.cos(Y*np.pi/180.)).reshape(ny*nx)
weights = np.where(flag > 400, 0, weights)
weights = np.where(ice.reshape(ny*nx)>0, 0, weights)

# obs impact
impact_of_sat = so_mean[:,:414,:,:,:].sum(axis=1)
impact_of_sst = so_mean[:,414:,:,:,:].sum(axis=1)

# divide figure
spec = gridspec.GridSpec(ncols=3, nrows=2, width_ratios=[5,5,0.2], height_ratios=[1,1])
fig = plt.figure(figsize=(10,5))

for i,res in enumerate(resList) :
    if i == 0 :
        pn = 0
        pn2 = 0
    elif i == 6 :
        pn = 3
        pn2 = 2
    else :
        continue
    # impact of temp on sst
    ax = fig.add_subplot(spec[pn], projection=ccrs.Robinson(central_longitude=180.))
    ax.add_feature(cftr.LAND, facecolor='gainsboro',zorder=100)
    ax.coastlines(zorder=101)
    im = ax.pcolormesh(X,Y, np.ma.masked_where(so_mean[i,1,1,:,:]>10000, impact_of_sat[i,1,:,:]), cmap=cmap, vmax=vmax[i], transform=ccrs.PlateCarree())
    im2 = ax.pcolor(X,Y, np.ma.masked_less(ice,0.3), hatch='///', alpha=0.0, lw=0.1, transform=ccrs.PlateCarree())
    ax.scatter( X[obs_temp[0,:,:]==1.], Y[obs_temp[0,:,:]==1.], s=s_om, marker='.', c=c_om, alpha=1,linewidths=0.5, edgecolors=c_ol, transform=ccrs.PlateCarree(), zorder=102)
    ax.set_title('('+ string.ascii_lowercase[pn2] +')', loc='left')
    ax.set_title('SAT', loc='center')
    ax.text(-0.15, 0.5, res, ha='center', va='center', rotation='vertical', fontsize=12,bbox={'boxstyle':'square','fc':'w'}, transform=ax.transAxes)
    
    
    # impact of sst on sst
    ax = fig.add_subplot(spec[pn+1], projection=ccrs.Robinson(central_longitude=180.))
    ax.add_feature(cftr.LAND, facecolor='gainsboro',zorder=100)
    ax.coastlines(zorder=101)
    im = ax.pcolormesh(X,Y, np.ma.masked_where(so_mean[i,1,1,:,:]>10000, impact_of_sst[i,1,:,:]), cmap=cmap, vmax=vmax[i], transform=ccrs.PlateCarree())
    im2 = ax.pcolor(X,Y, np.ma.masked_less(ice,0.3), hatch='///', alpha=0.0, lw=0.1, transform=ccrs.PlateCarree())
    ax.scatter( X[obs_sst[:,:]==1.], Y[obs_sst[:,:]==1.], s=s_om, marker='.', c=c_om, alpha=1,linewidths=0.5, edgecolors=c_ol, transform=ccrs.PlateCarree(), zorder=102)
    ax.set_title('('+ string.ascii_lowercase[pn2+1] +')', loc='left')
    ax.set_title('SST', loc='center')


    # colorbar
    ax = fig.add_subplot(spec[pn+2])
    clb = fig.colorbar(im, cax=ax, extend='max')
#
plt.savefig('anlsens_map_' +expName+ '_sst.pdf', dpi=300)
#plt.savefig('anlsens_map_' +expName+ '_sst.png', dpi=200)
plt.close()

#----------------------------------------------------------------------------------------------------
# plot 2
x = np.arange(len(resList))
y_sat = np.average(impact_of_sat.reshape(len(resList),nv,ny*nx), weights=weights, axis=2)
y_sst = np.average(impact_of_sst.reshape(len(resList),nv,ny*nx), weights=weights, axis=2)

fig = plt.figure(figsize=(10,4))
ax = fig.add_subplot(1,2,1)
ax.bar( x - bar_width/2., y_sat[:,1], width=bar_width, lw=1, fc='w', ec='k', label='SAT')
ax.bar( x + bar_width/2., y_sst[:,1], width=bar_width, lw=1, fc='gray', ec='k', label='SST')
ax.set_ylabel('Analysis Sensitivity to Obs')
ax.set_xticklabels(resList)
ax.set_xticks(x)
ax.set_title('('+ string.ascii_lowercase[0] +')', loc='left')
ax.set_title('Total impact', loc='center')
legend = ax.legend(frameon=True, fancybox=False)
legend.get_frame().set_linewidth(1.)
legend.get_frame().set_edgecolor('k')


#
x = np.arange(len(resList))
y_sat = np.average(impact_of_sat.reshape(len(resList),nv,ny*nx), weights=weights, axis=2) / 413
y_sst = np.average(impact_of_sst.reshape(len(resList),nv,ny*nx), weights=weights, axis=2) / 37

ax = fig.add_subplot(1,2,2)
ax.bar( x - bar_width/2., y_sat[:,1], width=bar_width, lw=1, fc='w', ec='k', label='SAT')
ax.bar( x + bar_width/2., y_sst[:,1], width=bar_width, lw=1, fc='gray', ec='k', label='SST')
#ax.set_ylabel('Analysis Sensitivity to Obs')
ax.set_xticklabels(resList)
ax.set_xticks(x)
ax.set_title('('+ string.ascii_lowercase[1] +')', loc='left')
ax.set_title('Per-obs impact', loc='center')
#plt.show()
plt.savefig('anlsens_bar_' +expName+ '_sst.pdf', dpi=300)
#plt.savefig('anlsens_bar_' +expName+ '_sst.png', dpi=200)


