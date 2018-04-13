# Given a file with directory and LC file information,
# loop over LCs and do some computations. 

# Libraries
import os, sys, math
import numpy as np
from astropy.table import Table
from astropy.io import ascii 
from matplotlib import pyplot as plt
from astroML.plotting import scatter_contour
from astroML.plotting import hist as histML

def processLC(LCfile): 
	# read into astropy Table 
	data = ascii.read(LCfile) 
        rmag = data['col31']
        rerr = data['col36']
        tai = data['col51'] + data['col56']
	Ndata = rmag.size 
        ### compute and return
	# dT = taiMax - taiMin
        # min, max, rms, and chi2 for the r band
	dT = np.max(tai) - np.min(tai) 
	rmed = np.median(rmag)
	# robust standard deviation
	rms = 0.7413*(np.percentile(rmag, 75)-np.percentile(rmag, 25))
	# robust chi2dof
	chi = (rmag-rmed)/rerr
	chi2 = 0.7413*(np.percentile(chi, 75)-np.percentile(chi, 25))
        ### eventualy, run Lomb-Scargle here
	# bestPeriods = LombScargle(tai, rmag) 

        return np.min(tai), dT, Ndata, rmed, rms, chi2 

def getLCdata(LCfile): 
	# read into astropy Table 
	data = ascii.read(LCfile) 
        rmag = data['col31']
        rerr = data['col36']
        tai = data['col51'] + data['col56']
        return rmag, rerr


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
            nPts[i], medianBin[i], sigGBin[i] = 0
        
    if (verbose):
        print('median:', np.median(medianBin[Npts>0]))

    return xBin, nPts, medianBin, sigGbin


def analyzeSTATSfile(statsFile, outdir): 
    """Make a 4-panel plot illustrating basic stats"""
     # Ndata histogram 
     # rchi2 histogram 
     # rrms vs. rmed 
     # rchi2 vs. rmed 

    # <TAImin     dT   Ndata   rmed   rrms   rchi2>
    dataT = np.loadtxt(statsFile, skiprows=1, usecols = (0,1,2,3,4,5))
    data = dataT.transpose(1,0)
    TAImin = data[0]
    dT = data[1]
    Ndata = data[2]
    rmag = data[3]
    rrms = data[4]
    rchi2 = data[5]

    ### PLOTTING ###
    plot_kwargs = dict(color='k', linestyle='none', marker='.', markersize=1)
    plt.subplots_adjust(bottom=0.10, top=0.92, left=0.11, right=0.9, wspace=0.41, hspace=0.32)

    ax1 = plt.subplot(2,2,1)
    hist, bins = np.histogram(Ndata, bins=60)
    center = (bins[:-1]+bins[1:])/2
    ax1.plot(center, hist, drawstyle='steps')   
    #histML(Ndata, bins='knuth', ax=ax1, histtype='stepfilled', ec='k', fc='#AAAAAA')
    ax1.set_xlim(0, 60)
    ax1.set_ylim(0, 1.1*np.max(hist))
    ax1.set_xlabel(r'$\mathrm{Ndata}$')
    ax1.set_ylabel(r'$\mathrm{dN/dNdata}$') 
    ax1.set_title('Standard stars from SDSS Stripe 82')

    rchi2ok = rchi2[rchi2<3.0]
    hist, bins = np.histogram(rchi2ok, bins=50)
    center = (bins[:-1]+bins[1:])/2
    ax2 = plt.subplot(2,2,2)
    ax2.plot(center, hist, drawstyle='steps')   
    ax2.set_xlim(0, 3.0)
    ax2.set_ylim(0, 1.1*np.max(hist))
    ax2.set_xlabel(r'$\mathrm{\chi^2_{dof}(r)}$')
    ax2.set_ylabel(r'$\mathrm{dN/d\chi^2_{dof}}$')

    ax3 = plt.subplot(2,2,3)
    ax3.set_xlim(13, 23)
    ax3.set_ylim(0, 0.1)
    ax3.set_xlabel(r'$\mathrm{rMag}$')
    ax3.set_ylabel(r'$\mathrm{rRms}$')
    if (0):
        ax3.scatter(rmag, rrms, s=3, alpha=0.5)
    else: 
        col0 = rmag
	col1 = rrms 
        # 2D-histogram 
	im3 = ax3.hexbin(col0, col1, bins='log', cmap=plt.cm.viridis,
               mincnt=1, extent=(13, 22.0, 0, 0.1))
        # color bar
        # cb = plt.colorbar(im,label='log(N)')

    ax4 = plt.subplot(2,2,4)
    ax4.set_xlim(13, 22)
    ax4.set_ylim(0, 3.0)
    ax4.set_xlabel(r'$\mathrm{rMag}$')
    ax4.set_ylabel(r'$\mathrm{\chi^2_{dof}(r)}$')
    if (0):
        ax4.scatter(rmag, rchi2, s=6, alpha=0.5)
    else: 
        col0 = rmag[rchi2<3.0]
	col1 = rchi2ok
        # 2D-histogram 
	im4 = ax4.hexbin(col0, col1, bins='log', cmap=plt.cm.viridis,
               mincnt=1, extent=(13, 22.0, 0, 3.0))
        # color bar
        # cb = plt.colorbar(im,label='log(N)')

    outfile = outdir + "/statsFile.png" 	
    plt.savefig(outfile)
    # plt.show() 

    return



def processLCs(LCdir):
    # Read Karun's file with LC data 
    filein = LCdir + "/LC_dirs_fils.lst" 
    fileout = LCdir + "/LC_dirs_fils.stats" 
    Nlines = 0 
    NlineMin = 2
    NlineMax = 2000001
    EpochSum = 0 
    mssgStep = 1000
    fout = open(fileout, "w")
    fout.write("   TAImin     dT   Ndata   rmed   rrms   rchi2 LCfile\n")
    with open(filein, "r") as f:
        for line in f:
            Nlines += 1
	    if (Nlines == Nlines/mssgStep*mssgStep):
		print Nlines
	    lineS = line.split()
	    if ((Nlines >= NlineMin) and (Nlines <= NlineMax)):
		RAdir = lineS[3]
		DECdir = lineS[4]
		LCfile = lineS[5]
		Nepoch = int(lineS[6])
		EpochSum += Nepoch
		LCpath = RAdir + "/" + DECdir + "/" + LCfile
		LCfullpath = LCdir + "/" + LCpath
		# test that it's not a file without data (i.e. only 2 header lines)
		lcf = open(LCfullpath, "r")
		lines = lcf.readlines()
		if (len(lines)>2):
                    r1, r2, r3, r4, r5, r6 = processLC(LCfullpath)
		    s = str("%.4f " % r1) + str("%.2f  " % r2) + str("%3.0f  " % r3) + str("%.3f  " % r4)
		    s = s + str("%.3f  " % r5) + str("%6.2f  " % r6) + LCpath + "\n"
                    fout.write(s)
		else:
	            print 'EMPTY FILE: ', LCpath 		
    fout.close() 
    print 'EpochSum = ', EpochSum 
    return 

if(0):
    LCdir = "/Users/ivezic/Work/Science/CalibrationV2/Data/testLCdir"
    processLCs(LCdir)
    statsfile = LCdir + "/LC_dirs_fils.stats" 
    analyzeSTATSfile(statsfile, LCdir)
    # plot in statsFile_r.png 
    # 481,000 with 15<r<20 and Ndata>10  (162,000 with Ndata>20) 



