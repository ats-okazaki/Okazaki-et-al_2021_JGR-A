'''
note : momory size is too large to be ran in head node. use computational node.
note : anaconda3/2019.10 or higher version should be used
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import speedy.noleap as nl
import speedy.geomet as geo
from netCDF4 import Dataset
import random
from mpl_toolkits.basemap import Basemap



def getMon(iday) :
    daystot = [31 ,59 ,90 ,120 ,151 ,181 ,212 ,243 ,273 ,304 ,334 ,365]
    for mon in np.arange(12) :
        if iday < daystot[mon] :
            break
    return mon
#------------------------------------------------------------------------
#
# file, parameter, variable
#
tresList = ['1dy','3dy','10dy','1mo','3mo','6mo','12mo']
lengthList = np.arange(0,5500,500)
lengthListM = np.arange(250,5250,500)
per50 = np.empty((len(tresList),len(lengthListM)),dtype=float)
per25 = np.empty((len(tresList),len(lengthListM)),dtype=float)
per75 = np.empty((len(tresList),len(lengthListM)),dtype=float)
colors = cm.rainbow(np.linspace(0, 1, len(tresList)))

exp = 'piControl'
model = 'MIROC5'
ensemble = 'r1i1p1'
timeFreq = 'day'
realm = {'tas':'atmos', 'tos':'ocean'}
varList = ['tas','tos']
xIndx = {'tas':'lon', 'tos':'rlon'}
yIndx = {'tas':'lat', 'tos':'rlat'}

dirBase = '/data/user/CMIP5/'+ exp +'/'+ timeFreq +'/'
distFilePath = 'dist_miroc5.npy'
yearStart = 2400
yearEnd = 2499
dyearFile = 10
nYear = yearEnd - yearStart + 1
nFiles = int(nYear / dyearFile)
nSample = 100

#------------------------------------------------------------------------
#
# read data (get indices)
#
# tas
sDateStart = str(yearStart + dyearFile * 0) + '0101'
sDateEnd = str(yearStart + dyearFile * (0 + 1) - 1) + '1231'
var = 'tas'
filePath = dirBase + realm[var] +'/'+ var +'/'+ model +'/'+ ensemble +'/'
filePath += var +'_'+ timeFreq +'_'+ model +'_'+ exp +'_'+ ensemble +'_'+ sDateStart +'-'+ sDateEnd +'.nc'
if not os.path.exists(filePath) :
    print('File does not exist! ',filePath)
nc = Dataset(filePath,'r')
nt = len(nc.dimensions['time'])
nxA = len(nc.dimensions[xIndx[var]])
nyA = len(nc.dimensions[yIndx[var]])
lonA = nc.variables[xIndx[var]][:]
latA = nc.variables[yIndx[var]][:]
nc.close()
# tos
var = 'tos'
filePath = dirBase + realm[var] +'/'+ var +'/'+ model +'/'+ ensemble +'/'
filePath += var +'_'+ timeFreq +'_'+ model +'_'+ exp +'_'+ ensemble +'_'+ sDateStart +'-'+ sDateEnd +'.nc'
if not os.path.exists(filePath) :
    print('File does not exist! ',filePath)
nc = Dataset(filePath,'r')
nt = len(nc.dimensions['time'])
nxO = len(nc.dimensions[xIndx[var]])
nyO = len(nc.dimensions[yIndx[var]])
print(filePath)
lonO = nc.variables['lon'][:]
latO = nc.variables['lat'][:]
nc.close()
datA = np.empty((nt*nFiles,nyA,nxA),dtype=float)
datO = np.empty((nt*nFiles,nyO,nxO),dtype=float)
    
#
# file check
#
existNpyFiles = True
for tres in tresList :
    
    npyFileName = 'cov_sat-sst_'+ model +'_'+ exp +'_'+ ensemble +'_'+ tres +'.npy'

    if not os.path.exists(npyFileName) :
        existNpyFiles = False

#
# read data
#
if not existNpyFiles :
    print('file does not exist')

    
    #
    # read data (data)
    #
    for nf in np.arange(nFiles) :
        sDateStart = str(yearStart + dyearFile * nf) + '0101'
        sDateEnd = str(yearStart + dyearFile * (nf + 1) - 1) + '1231'
    
        # tas
        var = 'tas'
        filePath = dirBase + realm[var] +'/'+ var +'/'+ model +'/'+ ensemble +'/'
        filePath += var +'_'+ timeFreq +'_'+ model +'_'+ exp +'_'+ ensemble +'_'+ sDateStart +'-'+ sDateEnd +'.nc'
        if not os.path.exists(filePath) :
            print('File does not exist! ',filePath)
        nc = Dataset(filePath,'r')
        datA[nf*nt:(nf+1)*nt,:,:] = nc.variables[var][:,:,:]
        nc.close()
    
        # tos
        var = 'tos'
        filePath = dirBase + realm[var] +'/'+ var +'/'+ model +'/'+ ensemble +'/'
        filePath += var +'_'+ timeFreq +'_'+ model +'_'+ exp +'_'+ ensemble +'_'+ sDateStart +'-'+ sDateEnd +'.nc'
        if not os.path.exists(filePath) :
            print('File does not exist! ',filePath)
        nc = Dataset(filePath,'r')
        datO[nf*nt:(nf+1)*nt,:,:] = nc.variables[var][:,:,:]
        nc.close()
    
        '''
        print(filePath)
        plt.imshow(np.ma.masked_where(datO[0,:,:]<200,datO[0,:,:]))
        plt.colorbar()
        plt.show()
        '''
    datO = np.ma.masked_where( datO < 200., datO )
    
    
    #
    # climatology (daily)
    #
    climDayA = np.zeros((365,nyA,nxA),dtype=float)
    climDayO = np.zeros((365,nyO,nxO),dtype=float)
    climSmthDayA = np.zeros((365,nyA,nxA),dtype=float)
    climSmthDayO = np.zeros((365,nyO,nxO),dtype=float)
    
    for it in np.arange(365 * nYear) :
        climDayA[it%365,:,:] += datA[it,:,:] / nYear 
        climDayO[it%365,:,:] += datO[it,:,:] / nYear 
    
    dt = 10
    for it in np.arange(dt) :
        climSmthDayA[it,:,:] = (np.sum(climDayA[-dt+it:,:,:],axis=0) + np.sum(climDayA[0:dt+1+it,:,:],axis=0))/(2*dt+1)
        climSmthDayO[it,:,:] = (np.sum(climDayO[-dt+it:,:,:],axis=0) + np.sum(climDayO[0:dt+1+it,:,:],axis=0))/(2*dt+1)
    for it in np.arange(dt,365-dt+1) :
        climSmthDayA[it,:,:] = climDayA[it-dt:it+dt,:,:].mean(axis=0)
        climSmthDayO[it,:,:] = climDayO[it-dt:it+dt,:,:].mean(axis=0)
    for it in np.arange(365-dt+1,365) :
        climSmthDayA[it,:,:] = (np.sum(climDayA[it-dt:,:,:],axis=0) + np.sum(climDayA[0:it+dt+1-365,:,:],axis=0))/(2*dt+1)
        climSmthDayO[it,:,:] = (np.sum(climDayO[it-dt:,:,:],axis=0) + np.sum(climDayO[0:it+dt+1-365,:,:],axis=0))/(2*dt+1)
        
    '''
    plt.imshow(np.mean(climDayO,axis=0)[-1::-1,:])
    plt.colorbar()
    plt.show()
    plt.plot(np.arange(365),climDayA[:,80,80])
    plt.show()
    '''
    
    #
    # climatology (monthly)
    #
    daysMon = [31 ,28 ,31 ,30 ,31 ,30 ,31 ,31 ,30 ,31 ,30 ,31]
    datMonA = np.zeros((12*nYear,nyA,nxA),dtype=float)
    datMonO = np.zeros((12*nYear,nyO,nxO),dtype=float)
    climMonA = np.zeros((12,nyA,nxA),dtype=float)
    climMonO = np.zeros((12,nyO,nxO),dtype=float)
    # data
    for it in np.arange(365 * nYear) :
        mon = getMon(it%365) 
        datMonA[ it//365*12 + mon,:,: ] += datA[it,:,:] 
        datMonO[ it//365*12 + mon,:,: ] += datO[it,:,:] 
    for mon in np.arange(12 * nYear) :
        datMonA[mon,:,:] /= daysMon[mon%12] 
        datMonO[mon,:,:] /= daysMon[mon%12] 
    # clim
    for mon in np.arange(12 * nYear) :
        climMonA[mon%12,:,:] += datMonA[mon,:,:] / nYear 
        climMonO[mon%12,:,:] += datMonO[mon,:,:] / nYear 
    
    '''
    plt.plot(np.arange(365),climSmthDayA[:,80,80])
    plt.show()
    plt.plot(np.arange(12),climMonA[:,80,80])
    plt.show()
    '''
    
    #
    # anomaly from climatology (daily)
    #
    for it in np.arange(365 * nYear) :
        datA[it,:,:] -= climSmthDayA[it%365,:,:]
        datO[it,:,:] -= climSmthDayO[it%365,:,:]
    for it in np.arange(12 * nYear) :
        datMonA[it,:,:] -= climMonA[it%12,:,:]
        datMonO[it,:,:] -= climMonO[it%12,:,:]
    
    '''
    plt.plot(np.arange(365*nYear),datO[:,80,80])
    plt.show()
    plt.plot(np.arange(12),climMonO[:,80,80])
    plt.show()
    print(' ')
    print('climSmthDayA')
    print(climSmthDayA[:,0,0])
    print(' ')
    print('climMonA')
    print(climMonA[:,0,0])
    print(' ')
    print('datA')
    print(datA[:,0,0])
    print(' ')
    print('datMonA')    
    print(datMonA[:,0,0])
    '''
    
    #
    # time-average
    #
    for i,tres in enumerate(tresList) :
        
        res = int(tres[:-2])
        resUnit = tres[-2].upper()
    
        if resUnit == 'D' :
            nt = int(365*nYear/res)
        elif resUnit == 'M' :
            nt = int(12*nYear/res)
    
        datAveA = np.zeros((nt,nyA,nxA),dtype=float)
        datAveO = np.zeros((nt,nyO,nxO),dtype=float)
        
        if resUnit == 'D' :
            for it in np.arange(nt) :
                datAveA[it,:,:] = np.mean(datA[it*res:(it+1)*res,:,:],axis=0)
                datAveO[it,:,:] = np.mean(datO[it*res:(it+1)*res,:,:],axis=0)
        elif resUnit == 'M' :
            for it in np.arange(nt) :
                datAveA[it,:,:] = np.mean(datMonA[it*res:(it+1)*res,:,:],axis=0)
                datAveO[it,:,:] = np.mean(datMonO[it*res:(it+1)*res,:,:],axis=0)
            
    
    #
    # sampling
    #        
        sampleA = np.empty((nSample,nyA,nxA),dtype=float)
        sampleO = np.empty((nSample,nyO,nxO),dtype=float)
        index = random.sample(range(nt),nSample)
        for jj,ii in enumerate(index) :
            sampleA[jj,:,:] = datAveA[ii,:,:]
            sampleO[jj,:,:] = datAveO[ii,:,:]
        '''
        print('datAveA')
        print(datAveA[:,0,0])
        print('sampleA')
        print(sampleA[:,0,0])
        '''
    
    #
    # covariance
    #
        sampleAveA = np.mean(sampleA,axis=0)
        sampleAveO = np.mean(sampleO,axis=0)
        sampleStdA = np.std(sampleA,axis=0)
        sampleStdO = np.std(sampleO,axis=0)
        sampleNormA = np.empty((nSample,nyA,nxA),dtype=float)
        sampleNormO = np.empty((nSample,nyO,nxO),dtype=float)
        '''
        plt.imshow(sampleAveA.reshape(nyA,nxA))
        plt.colorbar()
        plt.show()
        plt.imshow(sampleAveA.reshape(nyA,nxA))
        plt.colorbar()
        plt.show()
        plt.imshow(sampleAveO.reshape(nyO,nxO))
        plt.colorbar()
        plt.show()
        plt.imshow(sampleAveO.reshape(nyO,nxO))
        plt.colorbar()
        plt.show()
        exit()
        '''
        for it in np.arange(nSample) :
            sampleNormA[it,:,:] = (sampleA[it,:,:] - sampleAveA[:,:]) / sampleStdA[:,:]
            sampleNormO[it,:,:] = (sampleO[it,:,:] - sampleAveO[:,:]) / sampleStdO[:,:]
        cov = np.dot( sampleNormA.reshape(nSample,nyA*nxA).T, sampleNormO.reshape(nSample,nyO*nxO) ) / float(nSample)
        '''
        print('sampleNormA')
        print(sampleNormA[:,0,0])
        print(np.shape(cov))
        plt.imshow(cov)
        plt.colorbar()
        plt.show()
        '''
        npyFileName = 'cov_sat-sst_'+ model +'_'+ exp +'_'+ ensemble +'_'+ tres +'.npy'
        np.save( npyFileName, cov )

print('--------')

#
#------------------- data plot --------------------
#
#
# distance 
#
#'''
XA,YA = np.meshgrid(lonA,latA)
XA,YA = XA.flatten(),YA.flatten()
XO,YO = lonO, latO
XO,YO = XO.flatten(),YO.flatten()
if os.path.exists( distFilePath ) :
    dist = np.load( distFilePath )
else :
    dist = np.empty((nyA*nxA,nyO*nxO),dtype=float)
    for ii in range(len(XA)) :
        for jj in range(len(XO)) :
            dist[ii,jj] = geo.distance( XA[ii], YA[ii], XO[jj], YO[jj] )
            #print(XA[ii],YA[ii], XO[jj], YO[jj], dist[ii,jj])
    np.save( distFilePath, dist )

'''
plt.imshow(dist)
plt.colorbar()
plt.show()
'''
exit()
#'''
#
# read data (*.npy)
#
#for j,tres in enumerate(['1mo']) :
for j,tres in enumerate(tresList) :
    npyFileName = 'cov_sat-sst_'+ model +'_'+ exp +'_'+ ensemble +'_'+ tres +'.npy'
    cov = np.load( npyFileName )
    '''
    plt.imshow(cov)
    #plt.imshow(cov[nxO*100+150,:].reshape(nyO,nxO))
    plt.colorbar()
    plt.show()
    '''
    '''
    iyA = 100
    ixA = 130
    ijA = iyA * nxA + ixA
    ijO = (np.abs(YO-YA[ijA])+np.abs(XO-XA[ijA])).argmin()
    print(XA[ijA],YA[ijA])
    print(XO[ijO],YO[ijO])
    m = Basemap(projection='cyl', lon_0=180, resolution='c')
    m.drawcoastlines(color='k',linewidth=0.5)
    #m.scatter(XO,YO,c=dist[ijA,:], latlon=True, cmap='jet')
    m.scatter(XO,YO,c=np.ma.masked_where( dist[ijA,:]>1000,cov[ijA,:]), latlon=True, cmap='jet')
    plt.plot(XO[ijO],YO[ijO],'o',c='k')
    plt.colorbar()
    plt.show()
    '''
    #if j > 3 :
    if j == 1 or j == 2 or j == 4 or j == 6 :
        continue
    for l,length in enumerate(lengthList) :
        if l == len(lengthList)-1 :
            break
        tmpSample = cov[ (dist >= lengthList[l]) & (dist <= lengthList[l+1]) ]
        per50[j,l] = np.nanmedian(tmpSample)
        per75[j,l] = np.nanpercentile(tmpSample,75)
        per25[j,l] = np.nanpercentile(tmpSample,25)
        print(j,l,per50[j,l], per75[j,l], per25[j,l])

    plt.fill_between( lengthListM[:], per25[j,:], per75[j,:], facecolor=colors[j], alpha=0.3 )
    plt.plot( lengthListM[:], per50[j,:], color=colors[j], linewidth=3, label=tres )

plt.xlabel('Distance [km]')
plt.ylabel('Correlation Coefficient')
plt.title('MIROC5')
plt.xlim(0,5000)
plt.ylim(-0.2,1.0)
legend = plt.legend(frameon=True, fancybox=False)
#legend = plt.legend(bbox_to_anchor=(0.5,-1.0), loc='center', fontsize=18,frameon=True,fancybox=False)
legend.get_frame().set_linewidth(1.)
legend.get_frame().set_edgecolor('k')
plt.savefig('cov_miroc5_sat-sst_dist.png')

