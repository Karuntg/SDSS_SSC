# Read a list of (run, minField, maxField) from a text file and
# download photoObj*fits files in directories
# workingDir/rootName/run/camcol/
# (and optionally jpg images, but they are downloaded to the working directory) 

        # the final volume of light curves is about 20 GB, but all fits files are about 400 GB
        # it's about 4 GB per run worth of fits files, thus download a run, extract photometry
        # for standard stars and then purge fits files 
        # 3 mins for 50 files: for ~200,000: 8 days 

# Libraries
import os, sys, math, numpy, urllib



def dumpSDSSfile(outdir, run, camcol, field, downloadJPG):
        success = 0 
        try: 
            outfile = "photoObj-%06d-%d-%04d.fits" % (run,camcol,field)
            filepath = "%s/%s" % (outdir,outfile)
            if not os.path.exists(filepath):
                # this file doesn't exist - try to fetch it
                success = fetchSDSSfile(filepath,run,camcol,field,downloadJPG)
        except: 
            print "some problem with run=%06d camCol=%d field=%04d" % (run,camcol,field) 
            print "outdir=", outdir
            print "filepath=", filepath
            os.remove(outfile)
            pass 
        return success


def fetchSDSSfile(outfile,run, camcol, field, downloadJPG):
        try: 
            infile = "http://data.sdss3.org/sas/dr9/boss/photoObj/301/%d/%d/photoObj-%06d-%d-%04d.fits" % (run, camcol, run, camcol, field)
            # Download the fits file
            urllib.urlretrieve (infile,outfile)

            if downloadJPG:  
                outfile = "frame-irg-%06d-%d-%04d.jpg" % (run,camcol,field)
                infile = "http://data.sdss3.org/sas/dr9/boss/photoObj/frames/301/%d/%d/frame-irg-%06d-%d-%04d.jpg" % (run, camcol, run, camcol, field)
                # Download the jpg file
                urllib.urlretrieve (infile,outfile)
                print "downloaded photoObj-%06d-%d-%04d.fitsframe-irg-%06d-%d-%04d.jpg" % (run,camcol,field)
        except: 
            print "some problem with run=%06d camCol=%d field=%04d" % (run,camcol,field) 
            os.remove(outfile)
            pass 
        # if the file doesn't exist, urllib still makes an (almost) empty file, remove it...
        statinfo = os.stat(outfile)
        filesize = statinfo.st_size
        if filesize < 300: 
            os.remove(outfile)
            return 0 
        else:
            # print "downloaded photoObj-%06d-%d-%04d.fits" % (run,camcol,field)
            return 1


# Read in file with SDSS run data
## objlist = numpy.loadtxt('Stripe82RunList.dat')
objlist = numpy.loadtxt('S82runs_test.txt')
downloadJPG = 0
print objlist

numberFields = 0 
numberFieldsOK = 0 
rootName = "SDSSdata"
for line in objlist:
    print line
    run = int(line[0])
    outdir = "%s/%d" % (rootName,run)
    try:
        os.stat(outdir)
    except:
        os.makedirs(outdir)

    minField = int(line[1])
    maxField = int(line[2])
    for camcol in range(1,7):
        outdir = "SDSSdata/%d/%d" % (run,camcol)
        try:
            os.stat(outdir)
        except:
            os.mkdir(outdir)
        for iField in range(minField,maxField+1):
            print run, camcol, iField 
            numberFields += 1                
            numberFieldsOK += dumpSDSSfile(outdir, run, camcol, iField, downloadJPG)
            # processSDSSphotoObj(run, camcol, iField, downloadJPG, stdStars) 
    print "***** done with run", run

print "tried numberFields:", numberFields
print "downloaded:", numberFieldsOK


