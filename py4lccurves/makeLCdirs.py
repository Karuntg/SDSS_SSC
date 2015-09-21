# Make two-levels deep directories of the form RAm12/Dp6 where the top level 
# runs over 1 deg step in RA and the bottom level over 0.1 deg steps in Dec 
# (RA/Dec limits are set below).
# The code, e.g. RAm12, means all wrapped RA (RAw = RA - 360 when RA > 180) 
# is in the range -13 < RAw < -12 (deg) and Dp6 means +0.6 < Dec < +0.7 (deg)
# Directories are made relative to directory specified below by rootName 
# The code below also includes LCfilename(RA,Dec) for returning the names of
# the two diretories and also light curve file name for provided (RA, Dec) 
# of a given object. E.g. an object at RA=342.31214581 and Dec = -0.0441024
# would have its light curve stored in file with path 
#    RAname, Dname, LCname = LCfilename(342.31214581, -0.0441024)
#    print RAname, Dname, LCname
#    RAm17 Dm0 LC_RA342.312146_Dm0.044102.dat
#  that is, SDSSdata/RAm17/Dm0/LC_RA342.312146_Dm0.044102.dat

# Libraries
import os, sys, math, numpy, urllib

def LCfilename(RA, Dec): 
        # given (RA, Dec), return RA/Dec directory names and LC file name 
	RAdir = getRAname(RA)
        DecDir = getDname(Dec)

        # for filename itself, use 6 decimal places for RA and Dec
        sRA = '%0.6f'%(RA)        
        sDec = '%0.6f'%(abs(Dec))
        sgnString = "p" 
	if Dec < 0:
           sgnString = "m"
	LCfilename = "LC_RA"+sRA+"_D"+sgnString+sDec+".dat"

        return RAdir, DecDir, LCfilename 


def getRAname(RA):         
	if RA > 180:
           RAname = "RAm" + str(int(abs(RA-360)))
	else:
           RAname = "RAp" + str(int(RA))
	return RAname

def getDname(Dec): 
	if Dec < 0:
           Dname = "Dm" + str(int(abs(10*Dec)))
	else:
           Dname = "Dp" + str(int(10*Dec))
	return Dname


# RAw (RA-360 when RA>180) and Dec boundaries for Stripe 82 standard star catalog
RAwMin = -52
RAwMax = 60
DecMin = -1.3
DecMax = 1.3 
rootName = "SDSSdata"
numberDirs = 0 

for RAw in range(RAwMin, RAwMax+1):
    if RAw < 0:
       RA = RAw + 360.5
    else:
       RA = RAw
    RAname = getRAname(RA)
    print RAname
    # make directory rootName/RAname
    newdir = "%s/%s" % (rootName,RAname)
    try:
        os.stat(newdir)
    except:
        os.makedirs(newdir)
    for Dec10 in range(int(round(DecMin*10)), int(round(DecMax*10))):
        if Dec10 < 0: 
           Dec10 += 0.05
        Dec = Dec10/10.0
	Dname = getDname(Dec) 
        # print Dname
        # now make directory rootName/RAname/Dname
        newdir2 = "%s/%s" % (newdir,Dname)
        try:
            os.stat(newdir2)
        except:
            os.makedirs(newdir2)
        numberDirs += 1
  
print("Made %s directories" % numberDirs) 
