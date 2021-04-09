import numpy as np
import matplotlib.pyplot as plt
import colorcet as cc # use this for choosing cmap
from astropy.table import Table

####################################
## 23 Mar 2021
## edited v1 with updates to
## 1. color map
## 2. log bin for the 2D histo

####################################
## July 2, 2020:
## This is copy of ZItools.py
## edited to incorporate my changes
###################################

# Incorporating an additional plotting scheme
# using 2D histo, instead of a scatter plot
#
# CHANGE HISTORY
# 1. New plot function which does 2D histo
# 2. New median function which returns xEdge as well
# 3. New kw arguments, nBinX, nBinY
###################

# robust standard deviation
def sigG(arr):
    return 0.741*(np.quantile(arr, 0.75)-np.quantile(arr, 0.25))



################################
### PLOTS  
################################

# wrapper around plotdelMag below to use only two vectors instead of astropy Table 
def plotdelMagArr(xVec, yVec, kw):
    t = Table()
    t[kw['Xstr']] = xVec
    t[kw['Ystr']] = yVec
    plotdelMag(t, kw)
    return 

################################
# THIS DOES THE SCATTER PLOT
################################
def plotdelMag(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    # binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBin'], 1)
    # plotting
    if (kw['offset'] >= 0):
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c='red', linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c='red', linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c='red', linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='yellow', alpha=0.3)
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='yellow', linewidth=3)
    ax.plot(xBinM, TwoSigM, c='yellow', linewidth=3)
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='cyan', linewidth=3)
    ax.plot(xBinM, rmsM, c='cyan', linewidth=3)
    # 
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return
 
################################
# THIS DOES THE 2D HISTO
################################
def plotdelMag_KT(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    # get binned medians and std devs
    # get med, sig, as well as x-binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBinX'], 1)
    # get x-binning
    xMin,xMax,nBinX = kw['Xmin'], kw['Xmax'], kw['nBinX']
    # xEdge = np.linspace(yMin, yMax, (nBinY+1))    
    # get y-binning
    yMin,yMax,nBinY = kw['Ymin'], kw['Ymax'], kw['nBinY']
    # form the x,y edges vectors for 2d histo
    xedges = np.linspace(xMin,xMax,nBinX+1,endpoint=False,dtype=float)
    yedges = np.linspace(yMin,yMax,nBinY+1,endpoint=False,dtype=float)

    # create the 2D histo
    xvect,yvect = d[kw['Xstr']], d[kw['Ystr']]
    Histo2D, xedges, yedges = np.histogram2d(xvect,yvect, bins=(xedges, yedges))
    Histo2D = Histo2D.T
    
    # plotting
    # get the colormap from kw
    cmap = kw['cmap']
    fig, ax = plt.subplots(figsize=(12, 8))
    # ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    X, Y = np.meshgrid(xedges, yedges)
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap='Greys')
    cs = ax.pcolormesh(X,Y, Histo2D, cmap=cmap)
    ## Rearrange the plotting order in order to make the lines stand out
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='black', linewidth=2, linestyle='dashdot')
    ax.plot(xBinM, rmsM, c='black', linewidth=2, linestyle='dashdot')
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='blue', linewidth=3)
    ax.plot(xBinM, TwoSigM, c='blue', linewidth=3)
    # 
    if (kw['offset'] >= 0):
        offsetColor0 = 'magenta'
        offsetColor1 = 'black'
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c=offsetColor0, linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c=offsetColor1, linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c=offsetColor1, linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='yellow', alpha=0.3)
    #
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return


################################
# THIS DOES THE LOG 2D HISTO
# THIS ALSO PLOTS A COLORBAR
################################
def plotdelMag2log_KT(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    # get binned medians and std devs
    # get med, sig, as well as x-binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBinX'], 1)
    # get x-binning
    xMin,xMax,nBinX = kw['Xmin'], kw['Xmax'], kw['nBinX']
    # xEdge = np.linspace(yMin, yMax, (nBinY+1))    
    # get y-binning
    yMin,yMax,nBinY = kw['Ymin'], kw['Ymax'], kw['nBinY']
    # form the x,y edges vectors for 2d histo
    xedges = np.linspace(xMin,xMax,nBinX+1,endpoint=True,dtype=float)
    yedges = np.linspace(yMin,yMax,nBinY+1,endpoint=True,dtype=float)

    # create the 2D histo
    xvect,yvect = d[kw['Xstr']], d[kw['Ystr']]
    Histo2D, xedges, yedges = np.histogram2d(xvect,yvect, bins=(xedges, yedges))
    Histo2D = Histo2D.T
    logHisto2D = Histo2D # initial assignment
    pidx = Histo2D > 0 # values > 0, can take log10
    nidx = Histo2D <= 0 # values <= 0, cannot take log10, leave a min value
    logHisto2D[pidx] = np.log10(Histo2D[pidx])
    minval = np.min(logHisto2D[pidx]) # find the min val of transformed pix
    logHisto2D[nidx] = minval # set these pix to min val
    
    
    # plotting
    # get the colormap from kw
    cmap = kw['cmap']
    fig, ax = plt.subplots(figsize=(12, 8))
    # ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    X, Y = np.meshgrid(xedges, yedges)
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap='Greys')
    cs = ax.pcolormesh(X,Y, logHisto2D, cmap=cmap)
    # fig.colorbar(cs,location='top')
    cbar = fig.colorbar(cs, spacing='proportional',
                        shrink=0.95, ax=ax)
    cbar.set_label(label='Log$_{10}$(counts)',size=20)
    cbar_tix = 16 # Adjust as appropriate.
    cbar.ax.tick_params(labelsize=cbar_tix)
    
    ## Rearrange the plotting order in order to make the lines stand out
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='black', linewidth=2, linestyle='dotted')
    ax.plot(xBinM, rmsM, c='black', linewidth=2, linestyle='dotted')
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='black', linewidth=3, linestyle='dotted')
    ax.plot(xBinM, TwoSigM, c='black', linewidth=3, linestyle='dotted')
    # 
    if (kw['offset'] >= 0):
        offsetColor0 = 'black'
        offsetColor1 = 'black'
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c=offsetColor0, linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c=offsetColor1, linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c=offsetColor1, linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='yellow', alpha=0.3)
    #
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return

####################################
# THIS DOES THE LOG 2D HISTO
# AND CAN PLOT CONTOURS IF KW['contur'] = True
####################################
def plotdelMag2logDC_KT(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    # get binned medians and std devs
    # get med, sig, as well as x-binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBinX'], 1)
    # get x-binning
    xMin,xMax,nBinX = kw['Xmin'], kw['Xmax'], kw['nBinX']
    # xEdge = np.linspace(yMin, yMax, (nBinY+1))    
    # get y-binning
    yMin,yMax,nBinY = kw['Ymin'], kw['Ymax'], kw['nBinY']
    # form the x,y edges vectors for 2d histo
    xedges = np.linspace(xMin,xMax,nBinX+1,endpoint=True,dtype=float)
    yedges = np.linspace(yMin,yMax,nBinY+1,endpoint=True,dtype=float)

    # create the 2D histo
    xvect,yvect = d[kw['Xstr']], d[kw['Ystr']]
    Histo2D, xedges, yedges = np.histogram2d(xvect,yvect, bins=(xedges, yedges))
    Histo2D = Histo2D.T
    logHisto2D = Histo2D # initial assignment
    pidx = Histo2D > 0 # values > 0, can take log10
    nidx = Histo2D <= 0 # values <= 0, cannot take log10, leave a min value
    logHisto2D[pidx] = np.log10(Histo2D[pidx])
    minval = np.min(logHisto2D[pidx]) # find the min val of transformed pix
    logHisto2D[nidx] = minval # set these pix to min val
    
    
    # plotting
    # get the colormap from kw
    cmap = kw['cmap']
    fig, ax = plt.subplots(figsize=(12, 8))
    # ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    X, Y = np.meshgrid(xedges, yedges)
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap='Greys')
    cs = ax.pcolormesh(X,Y, logHisto2D, cmap=cmap)
    # fig.colorbar(cs,location='top')
    cbar = fig.colorbar(cs, spacing='proportional',
                        shrink=0.95, ax=ax)
    cbar.set_label(label='Counts [dex]',size=20)
    cbar_tix_font = 16 # Adjust as appropriate.
    cbar.ax.tick_params(labelsize=cbar_tix_font)
    
    # plot contours if kw['contur'] = True
    if kw['contur']:
        # check if contour levels have been set
        # else use from color bar
        cbar_tix = cbar.get_ticks()
        cbar_tix = cbar_tix[1:] # skip the min level for clarity
        if kw.get('cbar_tix'): # if kw set, then use that value
            cbar_tix = kw['cbar_tix']
        # cbar.location('top')
        CS = plt.contour(X[1:,1:],Y[1:,1:],logHisto2D,cbar_tix,
                    colors='brown',linestyles='dashdot')
        ax.clabel(CS, CS.levels, inline=True, fmt='%d', fontsize=20,colors='black')
    
    ## Rearrange the plotting order in order to make the lines stand out
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='black', linewidth=2, linestyle='dotted')
    ax.plot(xBinM, rmsM, c='black', linewidth=2, linestyle='dotted')
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='black', linewidth=3, linestyle='dotted')
    ax.plot(xBinM, TwoSigM, c='black', linewidth=3, linestyle='dotted')
    # 
    if (kw['offset'] >= 0):
        offsetColor0 = 'black'
        offsetColor1 = 'black'
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c=offsetColor0, linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c=offsetColor1, linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c=offsetColor1, linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='yellow', alpha=0.3)
    #
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return


################################
# THIS DOES THE 2D HISTO, in Black & White
# CHANGING ALL LINES TO BLACK, WITH GREYS CMAP
################################
def plotdelMagBW_KT(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    # get binned medians and std devs
    # get med, sig, as well as x-binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBinX'], 1)
    # get x-binning
    xMin,xMax,nBinX = kw['Xmin'], kw['Xmax'], kw['nBinX']
    # xEdge = np.linspace(yMin, yMax, (nBinY+1))    
    # get y-binning
    yMin,yMax,nBinY = kw['Ymin'], kw['Ymax'], kw['nBinY']
    # form the x,y edges vectors for 2d histo
    xedges = np.linspace(xMin,xMax,nBinX+1,endpoint=False,dtype=float)
    yedges = np.linspace(yMin,yMax,nBinY+1,endpoint=False,dtype=float)

    # create the 2D histo
    xvect,yvect = d[kw['Xstr']], d[kw['Ystr']]
    Histo2D, xedges, yedges = np.histogram2d(xvect,yvect, bins=(xedges, yedges))
    Histo2D = Histo2D.T
    
    # plotting
    # get the colormap from kw
    cmap = kw['cmap']
    fig, ax = plt.subplots(figsize=(12, 8))
    # ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    X, Y = np.meshgrid(xedges, yedges)
    cs = ax.pcolormesh(X,Y, Histo2D, cmap='Greys')
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap=cmap)
    ## Rearrange the plotting order in order to make the lines stand out
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='black', linewidth=2, linestyle='dashdot')
    ax.plot(xBinM, rmsM, c='black', linewidth=2, linestyle='dashdot')
    # 
    if (kw['offset'] >= 0):
        offsetColor0 = 'black'
        offsetColor1 = 'black'
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c=offsetColor0, linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c=offsetColor1, linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c=offsetColor1, linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='black', alpha=0.3)
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='silver', linewidth=3)
    ax.plot(xBinM, TwoSigM, c='silver', linewidth=3)
    #
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return

################################
# THIS DOES THE LOG 2D HISTO, in Black & White
# CHANGING ALL LINES TO BLACK, WITH GREYS CMAP
# THIS ALSO PLOTS A COLORBAR
################################
def plotdelMag2logBW_KT(d, kw):
 
    print('medianAll:', np.median(d[kw['Ystr']]), 'std.dev.All:', sigG(d[kw['Ystr']]))
    print('N=', np.size(d[kw['Ystr']]), 'min=', np.min(d[kw['Ystr']]), 'max=', np.max(d[kw['Ystr']]))

    # get binned medians and std devs
    # get med, sig, as well as x-binning
    xBinM, nPtsM, medianBinM, sigGbinM = fitMedians(d[kw['Xstr']], \
                                         d[kw['Ystr']], kw['XminBin'], kw['XmaxBin'], kw['nBinX'], 1)
    # get x-binning
    xMin,xMax,nBinX = kw['Xmin'], kw['Xmax'], kw['nBinX']
    # xEdge = np.linspace(yMin, yMax, (nBinY+1))    
    # get y-binning
    yMin,yMax,nBinY = kw['Ymin'], kw['Ymax'], kw['nBinY']
    # form the x,y edges vectors for 2d histo
    xedges = np.linspace(xMin,xMax,nBinX+1,endpoint=True,dtype=float)
    yedges = np.linspace(yMin,yMax,nBinY+1,endpoint=True,dtype=float)

    # create the 2D histo
    xvect,yvect = d[kw['Xstr']], d[kw['Ystr']]
    Histo2D, xedges, yedges = np.histogram2d(xvect,yvect, bins=(xedges, yedges))
    Histo2D = Histo2D.T
    logHisto2D = Histo2D # initial assignment
    pidx = Histo2D > 0 # values > 0, can take log10
    nidx = Histo2D <= 0 # values <= 0, cannot take log10, leave a min value
    logHisto2D[pidx] = np.log10(Histo2D[pidx])
    minval = np.min(logHisto2D[pidx]) # find the min val of transformed pix
    logHisto2D[nidx] = minval # set these pix to min val
    
    # plotting
    # get the colormap from kw
    cmap = kw['cmap']
    fig, ax = plt.subplots(figsize=(12, 8))
    # ax.scatter(d[kw['Xstr']], d[kw['Ystr']], s=kw['symbSize'], c='black') 
    X, Y = np.meshgrid(xedges, yedges)
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap='Greys')
    # cs = ax.pcolormesh(X,Y, Histo2D, cmap=cmap)
    cs = ax.pcolormesh(X,Y, logHisto2D, cmap='Greys')
    # fig.colorbar(cs,location='top')
    cbar = fig.colorbar(cs, spacing='proportional',
                        shrink=0.95, ax=ax)
    cbar.set_label(label='Log$_{10}$(counts)',size=20)
    cbar_tix = 16 # Adjust as appropriate.
    cbar.ax.tick_params(labelsize=cbar_tix)

    ## Rearrange the plotting order in order to make the lines stand out
    #
    rmsBin = np.sqrt(nPtsM) / np.sqrt(np.pi/2) * sigGbinM
    rmsP = medianBinM + rmsBin
    rmsM = medianBinM - rmsBin
    ax.plot(xBinM, rmsP, c='black', linewidth=2, linestyle='dashdot')
    ax.plot(xBinM, rmsM, c='black', linewidth=2, linestyle='dashdot')
    #
    TwoSigP = medianBinM + kw['Nsigma']*sigGbinM
    TwoSigM = medianBinM - kw['Nsigma']*sigGbinM 
    ax.plot(xBinM, TwoSigP, c='silver', linewidth=3, linestyle='dotted')
    ax.plot(xBinM, TwoSigM, c='silver', linewidth=3, linestyle='dotted')
    # 
    if (kw['offset'] >= 0):
        offsetColor0 = 'black'
        offsetColor1 = 'black'
        xL = np.linspace(kw['XminBin'], kw['XmaxBin'])
        ax.plot(xL, 0*xL, '-', c=offsetColor0, linewidth=3)
        ax.plot(xL, 0*xL+kw['offset'], '--', c=offsetColor1, linewidth=3)
        ax.plot(xL, 0*xL-kw['offset'], '--', c=offsetColor1, linewidth=3)
    if (0):
        ax.scatter(xBinM, medianBinM, s=30.0, c='black', alpha=0.8)
        ax.scatter(xBinM, medianBinM, s=15.0, c='black', alpha=0.3)
    #
    ax.set_xlabel(kw['Xlabel'], fontsize=22)
    ax.set_ylabel(kw['Ylabel'], fontsize=22)
    ax.set_xlim(kw['Xmin'], kw['Xmax'])
    ax.set_ylim(kw['Ymin'], kw['Ymax'])
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(kw['plotName'], dpi=600)
    print('saved plot as:', kw['plotName']) 
    plt.show()
    return


################################
# given vectors x and y, fit medians in bins from xMin to xMax, with Nbin steps,
# and return xBin, medianBin, medianErrBin 
def fitMedians(x, y, xMin, xMax, Nbin, verbose=1): 

    # first generate bins
    xEdge = np.linspace(xMin, xMax, (Nbin+1)) 
    xBin = np.linspace(0, 1, Nbin)
    nPts = 0*np.linspace(0, 1, Nbin)
    medianBin = 0*np.linspace(0, 1, Nbin)
    sigGbin = -1+0*np.linspace(0, 1, Nbin) 
    for i in range(0, Nbin): 
        xBin[i] = 0.5*(xEdge[i]+xEdge[i+1]) 
        yAux = y[(x>xEdge[i])&(x<=xEdge[i+1])]
        if (yAux.size > 0):
            nPts[i] = yAux.size
            medianBin[i] = np.median(yAux)
            # robust estimate of standard deviation: 0.741*(q75-q25)
            sigmaG = 0.741*(np.percentile(yAux,75)-np.percentile(yAux,25))
            # uncertainty of the median: sqrt(pi/2)*st.dev/sqrt(N)
            sigGbin[i] = np.sqrt(np.pi/2)*sigmaG/np.sqrt(nPts[i])
        else:
            nPts[i] = 0
            medianBin[i] = 0
            sigGbin[i] = 0
        
    if (verbose):
        print('median:', np.median(medianBin[nPts>0]), 'std.dev:', np.std(medianBin[nPts>0]))

    return xBin, nPts, medianBin, sigGbin

################################
# THIS ALSO RETURNS XEDGE
# given vectors x and y, fit medians in bins from xMin to xMax, with Nbin steps,
# and return xBin, medianBin, medianErrBin 
def fitMedians_KT(x, y, xMin, xMax, Nbin, verbose=1): 

    # first generate bins
    xEdge = np.linspace(xMin, xMax, (Nbin+1)) 
    xBin = np.linspace(0, 1, Nbin)
    nPts = 0*np.linspace(0, 1, Nbin)
    medianBin = 0*np.linspace(0, 1, Nbin)
    sigGbin = -1+0*np.linspace(0, 1, Nbin) 
    for i in range(0, Nbin): 
        xBin[i] = 0.5*(xEdge[i]+xEdge[i+1]) 
        yAux = y[(x>xEdge[i])&(x<=xEdge[i+1])]
        if (yAux.size > 0):
            nPts[i] = yAux.size
            medianBin[i] = np.median(yAux)
            # robust estimate of standard deviation: 0.741*(q75-q25)
            sigmaG = 0.741*(np.percentile(yAux,75)-np.percentile(yAux,25))
            # uncertainty of the median: sqrt(pi/2)*st.dev/sqrt(N)
            sigGbin[i] = np.sqrt(np.pi/2)*sigmaG/np.sqrt(nPts[i])
        else:
            nPts[i] = 0
            medianBin[i] = 0
            sigGbin[i] = 0
        
    if (verbose):
        print('median:', np.median(medianBin[nPts>0]), 'std.dev:', np.std(medianBin[nPts>0]))

    return xEdge,xBin, nPts, medianBin, sigGbin


