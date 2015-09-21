#!/usr/bin/python

### v1: 21 Aug 2015
### Modifying v0 to run for just one run
### input pars are /SDSSdata, <run #>, minfield, maxfield
### this may be easier for batch processing
### no input file to vcp from vos, also easier to check/manage the output
###

### 12 SEPT 2014
### Modifying DownloadSDSScsv.py to run on
### a Canfar VMOD in batch mode
### needs input/output dir and input run list to be
### given as input parameters
###
### 6 Aug 2014
### Extending Zeljko's version of the SDSS photoobj fits file downloader
### to also process the fits files and write out a csv file with the reqd cols only


## Zeljko's code
# Read a list of (run, minField, maxField) from a text file and
# download photoObj*fits files in directories
# workingDir/rootName/run/camcol/
# (and optionally jpg images, but they are downloaded to the working directory) 
## KT's code
# python script to read in Stripe 82 PhotoObj fits tables,
## then parse for required cols, output to ASCII csv table
## >>>>> cols reqd: 'RA', 'DEC','MJD','OBJC_ROWC','OBJC_COLC','NCHILD','PSFMAG','PSFMAGERR','TAI','AIRMASS','PSF_FWHM','SKYFLUX'


######### Notes ######################
        # 1. The final volume of light curves is about 20 GB, but all fits files are about 400 GB
        # it's about 4 GB per run worth of fits files, thus download a run, extract photometry
        # for standard stars and then purge fits files 
        # Estimated proc time = 3 mins for 50 files: for ~200,000: 8 days
	# 2. Processing the fits tables in situ and writing to csv (then deleting the fits)
	# may reduce storage requirements, but will increase processing time
	# 3. Establish a VM on Canfar for processing and storage
#####################################
	

# Import libraries
import os,sys,math,numpy,urllib,pyfits,csv,time,math,string
sys.stdout.flush()

############ Function Defs #########################
def dumpSDSSfile(outdir, run, camcol, field, downloadJPG):
        success = 0
	## print 'dumpSDSSfile: ',success
        try: 
            outfile = "photoObj-%06d-%d-%04d.fits" % (run,camcol,field)
            filepath = "%s/%s" % (outdir,outfile)
	    success = fetchSDSSfile(filepath,run,camcol,field,downloadJPG)
	    print "dumpSDSSfile done: %s" % (filepath)
        except: 
            print "some problem with run=%06d camCol=%d field=%04d" % (run,camcol,field) 
            print "outdir=", outdir
            print "filepath=", filepath
	    if os.path.exists(filepath):
		    os.remove(filepath)
	    success = 0
            pass

        print 'dumpSDSSfile: ',success
        return success


############ Function Defs #########################
def fetchSDSSfile(outfile,run, camcol, field, downloadJPG):
        success = 0
        try: 
            infile = "http://data.sdss3.org/sas/dr9/boss/photoObj/301/%d/%d/photoObj-%06d-%d-%04d.fits" % (run, camcol, run, camcol, field)
            # Download the fits file
            urllib.urlretrieve (infile,outfile)
            print "downloaded photoObj-%06d-%d-%04d.fits" % (run,camcol,field)
            success = 1
        except: 
            print "some problem with run=%06d camCol=%d field=%04d" % (run,camcol,field) 
            if os.path.exists(outfile):
                    os.remove(outfile)
            success = 0
            pass 
        # if the file doesn't exist, urllib still makes an (almost) empty file, remove it...
        if os.path.exists(outfile):
		statinfo = os.stat(outfile)
		filesize = statinfo.st_size
                print "filesize %06d" % (filesize)
		if filesize < 300: 
			os.remove(outfile)
                        success = 0
        return success

############ Function Defs #########################
def chk4FitsTable(indir,infits):
        	success = 0 
		fitshd = pyfits.getheader(indir+infits,1)
		## fitshd['XTENSION'].lower()
		if fitshd['XTENSION'].lower() == 'bintable':
			success = 1
		## print success
		return success

############ Function Defs #########################
def convFits2CSV(indir,infits):
        	success = 0 
        ## try: 
		## create the out cat name from the fits name
		fitsnam = infits.split('.')
		outcat = fitsnam[0]+'.cat'
		## print "%s %s %s" % (indir,infits,outcat)
		## open the fits table, get info
		hdulist = pyfits.open(indir+infits)

		## read in table data to array
		tbdata = hdulist[1].data

		## get col names for later ref
		cols = hdulist[1].columns

		## assign all required fields to vectors/arrays
		ra = tbdata['RA']
		dec = tbdata['DEC']
		mjd = tbdata['MJD']
		obj_row = tbdata['OBJC_ROWC']
		obj_col = tbdata['OBJC_COLC']
		nchild = tbdata['NCHILD']
		psfmag = tbdata['PSFMAG']
		psferr = tbdata['PSFMAGERR']
		tai = tbdata['TAI']
		amass = tbdata['AIRMASS']
		psffwhm = tbdata['PSF_FWHM']
		skyflux = tbdata['SKYFLUX']
		nobj = (ra.shape)[0]
		nmag = (psfmag.shape)[1]
		## print "Num obj/mag: %d %d" % (nobj,nmag)
		## close fits table
		hdulist.close()

		## open output ASCII catalog
		f1= open(indir+outcat,"w")
		f1.write("## RA,Dec,MJD,objc_colc,objc_rowc,nchild,psfmag[5],psfmagerr[5],tai[5],airmass[5],psf_fwhm[5],skyflux[5]\n")

		j = 0
		while j <= nobj-1:
## while j <= 1:
			outstr = str(j+1)+','+str(ra[j])+','+str(dec[j])+','+str(mjd[j])+','+str(obj_col[j])+','+str(obj_row[j])+','+str(nchild[j])

			## make mag string
			magstr = ''
			for tmpmag in psfmag[j][:]:
				tmpstr = ','+str(tmpmag)
				magstr += tmpstr

			## make err string
			errstr = ''
			for tmperr in psferr[j][:]:
				tmpstr = ','+str(tmperr)
				errstr += tmpstr
		
			## make tai string
			taistr = ''
			for tmptai in tai[j][:]:
				tmpstr = ','+str(tmptai)
				taistr += tmpstr

			## make airmass string
			amstr = ''
			for tmpam in amass[j][:]:
				tmpstr = ','+str(tmpam)
				amstr += tmpstr

			## make fwhm string
			fwstr = ''
			for tmpfw in psffwhm[j][:]:
				tmpstr = ','+str(tmpfw)
				fwstr += tmpstr

			## make skyflux string
			sfstr = ''
			for tmpsf in skyflux[j][:]:
				tmpstr = ','+str(tmpsf)
				sfstr += tmpstr


			## output string to outfile
			outstr += magstr+errstr+taistr+amstr+fwstr+sfstr
			f1.write(outstr+'\n')
			## increment while counter
			j += 1


		## all done, close output catalog and quit
		f1.close()
		success = 1

        ## except:
		## print "Problem with Fits to CSV conv in %s/%s" % (indir,infits)
		## success = 0
		return success

############ Main Pro #########################

## output root dir name: SDSSdata
## input run number and min/max fields given as input args
starttim = time.time()
indir = sys.argv[1]
run = int(sys.argv[2])
minField = int(sys.argv[3])
maxField = int(sys.argv[4])
downloadJPG = 0

numberFields = 0 
numberFieldsOK = 0 
success = 0
rootName = indir
## Main loop for run
outdir = "%s/%d" % (rootName,run)
try:
        os.stat(outdir)
except:
        os.makedirs(outdir)
## Loop for camcol
for camcol in range(1,7):
        outdir = "%s/%d/%d/" % (rootName,run,camcol)
        try:
            os.stat(outdir)
        except:
            os.mkdir(outdir)
	outglog = outdir+'Success.log'
	outblog = outdir+'Failure.log'
	flog1 = open(outglog,"w")
	flog2 = open(outblog,"w")
	## Loop for min -> max fields
        for iField in range(minField,maxField+1):
            ## print run, camcol, iField 
            numberFields += 1
	    fits_success = 0
            fits_success = dumpSDSSfile(outdir, run, camcol, iField, downloadJPG)
	    print 'fits success: ',fits_success
            ## success += fits_success
            if fits_success:
                    indir = "%s/%d/%d/" % (rootName,run,camcol)
                    infits = "photoObj-%06d-%d-%04d.fits" % (run, camcol, iField)
                    ## check if fits table, if so call Fits2CSV, else log error
                    if chk4FitsTable(indir,infits):
                            ## print "Fits2CSV %s %s" % (indir,infits)
                            numberFieldsOK += convFits2CSV(indir,infits)
                            flog1.write(infits+'\n')
                            ## delete the fits table
                            inpath = '%s/%s' % (indir,infits)
                            print 'Deleting converted fits table: %s' % (inpath)
                            os.remove(inpath)
                    else:
                            print "Bad fits: %s" % (infits)
                            flog2.write(infits+'\n')
            else:
                     infits = "photoObj-%06d-%d-%04d.fits" % (run, camcol, iField)
                     print "Bad fits: %s" % (infits)
                     flog2.write(infits+'\n')

                     
    	flog1.close()
    	flog2.close()
        # processSDSSphotoObj(run, camcol, iField, downloadJPG, stdStars) 
    	## print "***** done with run", run

## print "tried numberFields:", numberFields
print "downloaded:", numberFieldsOK
endtim = time.time()
print "Total time: %12.2f" % ((endtim - starttim)/3600.) 

