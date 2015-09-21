#!/usr/bin/python

### 8 Dec 2014
### First attempt to match Zeljko's 2007 SS catalog
### with each of the S82 field catalogs 


############ Libraries #########################

# Import libraries
import os,sys,math,numpy,urllib,pyfits,csv,time,math,string
from numpy import *
import string
import sys
import copy

sys.stdout.flush()

############ Function: Read in main SS catalog  #########################
def readZI_SSC(indir,SSC,SSCra,SSCdec,nSS):
        success = 0L
        try:
                f1= open(indir+SSC,"r")

                ## define running counters
                nlin = 0L
                nSS = 0L
                ## parse file line by line
                for line in f1:
                        ## check for comment lines
                        ## li=line.strip()
                        if not line.startswith("#"):
                                data = line.split()
                                SSCra[nSS] = float(data[1])
                                SSCdec[nSS]  = float(data[2])
                                nSS += 1L
                        nlin += 1L
                f1.close()
                print 'Number of lines read: %d' % (nlin)
                print 'Number of SS read: %d' % (nSS)
                # trim input vectors
                SSCra = SSCra[0:nSS-1]
                SSCdec = SSCdec[0:nSS-1]
                success = nSS
                return success
        except:
                print 'Problem opening input SSC: %s/%s ' % (indir,SSC)
                return success

###################################################################################
## http://www.stsci.edu/~ferguson/software/pygoodsdist/pygoods/angsep.py
# angsep.py
# Program to calculate the angular separation between two points
# whose coordinates are given in RA and Dec
# From angsep.py Written by Enno Middelberg 2001
##################################################################################
def angsep(ra1deg,dec1deg,ra2deg,dec2deg):
    """ Determine separation in degrees between two celestial objects 
        arguments are RA and Dec in decimal degrees. 
    """
    ra1rad=ra1deg*pi/180
    dec1rad=dec1deg*pi/180
    ra2rad=ra2deg*pi/180
    dec2rad=dec2deg*pi/180
    
    # calculate scalar product for determination
    # of angular separation
    
    x=cos(ra1rad)*cos(dec1rad)*cos(ra2rad)*cos(dec2rad)
    y=sin(ra1rad)*cos(dec1rad)*sin(ra2rad)*cos(dec2rad)
    z=sin(dec1rad)*sin(dec2rad)
    
    rad=arccos(x+y+z) # Sometimes gives warnings when coords match
    
    # use Pythargoras approximation if rad < 1 arcsec
    sep = choose( rad<0.000004848 , (
        sqrt((cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2),rad))
        
    # Angular separation
    sep=sep*180/pi

    return sep

##################################################################################
### python match pros taken from:
### http://www.stsci.edu/~ferguson/software/pygoodsdist/pygoods/match.py
##################################################################################

def matchsorted(ra,dec,ra1,dec1,tol):
    """ Find closest ra,dec within tol to a target in an ra-sorted list of ra,dec.
        Arguments:
          ra - Right Ascension decimal degrees (numpy sorted in ascending order)
          dec - Declination decimal degrees (numpy array)
          ra1 - RA to match (scalar, decimal degrees)
          ra1 - Dec to match (scalar, decimal degrees)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - index of the best match within tol; -1 if no match within tol
          sep - separation (defaults to tol if no match within tol)
    """
    i1 = searchsorted(ra,ra1-tol)-1
    i2 = searchsorted(ra,ra1+tol)+1
    if i1 < 0:
        i1 = 0
    sep = angsep.angsep(ra[i1:i2],dec[i1:i2],ra1,dec1)
    # print "tolerance ",tol
    indices=argsort(sep)
    # print sep
    if sep[indices[0]] > tol:
        return -1, tol
    ibest = indices[0] + i1
    return ibest, sep[indices[0]]
    
##################################################################################
def matchpos(ra1,dec1,ra2,dec2,tol):
    """ Match two sets of ra,dec within a tolerance.
        Longer catalog should go first
        Arguments:
          ra1 - Right Ascension decimal degrees (numpy array)
          dec1 - Declination decimal degrees (numpy array)
          ra2 - Right Ascension decimal degrees (numpy array)
          dec2 - Declination decimal degrees (numpy array)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - indices of the best matches within tol; -1 if no match within tol
          sep - separations (defaults to tol if no match within tol)
    """
    indices = argsort(ra1)
    rasorted = ra1[indices]
    decsorted = dec1[indices]
    ibest = []
    sep = []
    for i in range(len(ra2)):
        j,s = matchsorted(rasorted,decsorted,ra2[i],dec2[i],tol)
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
        sep += [s]
    return array(ibest),array(sep)

##################################################################################
def matchjoin(si1,si2,ibest,sep=[],dict1={},dict2={}):
    """ Keep only elements that match in both catalogs. 
        Arguments:
          si1 -- object with data arrays as attributes (e.g. si1.mag_auto, si1.id)
          si2 -- object with data arrays as attributes 
          ibest -- indices of si1 that match si2 (in order of si2)
        Keyword Arguments:
          sep -- array of separations, will be added as an attribute to second object
          dict1 -- dictionary of attributes for si1 (e.g. {'mag_auto':0,'id':1} )
          dict2 -- dictionary of attributes for si2 
        The objects si1 and si2 would normally be sextractor objects returned
        by sextutils.sextractor(). In that case, the column names are stored in
        si1._d and si2._d. If the objects are not sextractor objects, the user can
        provide separate dictionaries whose keys are the object attributes 
	that correspond to the numpy data arrays. 
       
        Returns:
          s1, s2 -- object arrays that include only the matching objects contained
              in both s1 and s2.
    """
    if len(dict1) == 0:
        dict1 = si1._d
    if len(dict2) == 0:
        dict2 = si2._d
    indices = compress(ibest > 0,ibest)
    # print indices
    s1 = copy.deepcopy(si1)
    s2 = copy.deepcopy(si2)
    for k in dict1.keys():
        if type(s1.__dict__[k]) == type([]):
             s1.__dict__[k] = char.array(s1.__dict__[k])
        s1.__dict__[k] = s1.__dict__[k][indices]
    flag = ibest > 0
    for k in dict2.keys():
        if type(s2.__dict__[k]) == type([]):
             s2.__dict__[k] = char.array(s2.__dict__[k])
        s2.__dict__[k] = compress(flag,s2.__dict__[k])
    if len(sep) > 0:
        s2.separation = compress(flag,sep)
    return s1,s2

def matchids(id1,id2):
    """ Match two sets of ids. 
        Returns: 
          ibest -- array of indices of i1 that match i2; -1 if no match
    """
    indices = argsort(id1)
    idsorted = id1[indices]
    ibest = []
    for i in range(len(id2)):
        j = matchidsorted(idsorted,id2[i])
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
    return array(ibest)

##################################################################################
def matchidsorted(ids,targetid):
    """ Find id matches, return index in i1 that matches targetid; -1 if no match. """
    i1 = searchsorted(ids,targetid)
    if targetid == ids[i1]:
        ibest = i1
    else:
        ibest = -1 
    return ibest


##################################################################################
############ Main Pro #########################
# Read in Zeljko's SS catalog
indir = '/villa2/karun/ALTAIR/SDSScalib/'
SSC = 'stripe82calibStars_v2.6.dat'
# define input arrays
SSCra = zeros(1500000)
SSCdec = zeros(1500000)
nSS = 0L
nSS = readZI_SSC(indir,SSC,SSCra,SSCdec,nSS)
print 'Number of SS read: %9d' % (nSS)
if (nSS > 0):
        print 'SS 0 ra, dec: %14.9f %14.9f ' % (SSCra[0],SSCdec[0])
        print 'SS nSS ra,dec: %14.9f %14.9f ' % (SSCra[nSS-1],SSCdec[nSS-1])

############### all done ################
print "Done",time.strftime("%d %b %Y %H:%M:%S")
sys.stdout.flush()
sys.exit()


