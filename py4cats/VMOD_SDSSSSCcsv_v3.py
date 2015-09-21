#!/usr/bin/python

### v3: 4 Sept 2015
### (Adding cuts suggested by Zeljko, see his email Sept 1:
### 64 < objc_rowc <= 1425
### nchild = 0 
### taiFrac =  tai - 24 * 3600 * (mjd + 3.0471e-07 * objc_rowc) 
 
### v2: 28 Aug 2015
### (Following telecon with Zeljko/ Douglas today)
### Modifying v1 to return additional columns for aux science
### cols needed:  OBJC_TYPE, ROWC[5], ROWCERR[5], COLC[5], COLCERR[5], MODELMAG[5], MODELMAGERR[5]
###

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
import numpy as np
import os,sys,math,urllib,pyfits,csv,time,math,string
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
def getSelIndx(objc_rowc,nchild):
### selection function for the foll two conditions
### returns a vector of selected indices
### 64 < objc_rowc <= 1425
### nchild = 0
	selindx = (np.logical_and(np.logical_and((objc_rowc > 64),(objc_rowc <= 1425)),(nchild == 0))).nonzero()
	## convert selindx from tuple to ndarray, extract only row 0
	selindx = np.asarray(selindx)
	selindx = selindx[0,:]
	return selindx

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
		objc_type = tbdata['OBJC_TYPE']
		objc_rowc = tbdata['OBJC_ROWC']
		objc_colc = tbdata['OBJC_COLC']
		nchild = tbdata['NCHILD']
		rowc = tbdata['ROWC']
		rowcerr = tbdata['ROWCERR']
		colc = tbdata['COLC']
		colcerr = tbdata['COLCERR']
		psfmag = tbdata['PSFMAG']
		psferr = tbdata['PSFMAGERR']
		modelmag = tbdata['MODELMAG']
		modelmagerr = tbdata['MODELMAGERR']
		tai = tbdata['TAI']
		amass = tbdata['AIRMASS']
		psffwhm = tbdata['PSF_FWHM']
		skyflux = tbdata['SKYFLUX']
		nobj = (ra.shape)[0]
		nmag = (psfmag.shape)[1]
		## print "Num obj/mag: %d %d" % (nobj,nmag)
		## close fits table
		hdulist.close()

		## get indices of sel objs
		selindx = getSelIndx(objc_rowc,nchild)

		## calc tai_frac
		taiFrac = tai * 0.
		## print tai.shape
		## print taiFrac.shape
		corr_fact = mjd + 3.0471e-07 * objc_rowc
		## print corr_fact.shape
		for colnum in xrange(0,4):
			tmptai = tai[:,colnum]
			tmptaiFrac = tmptai - 24. * 3600. * corr_fact
			taiFrac[:,colnum] =  tmptaiFrac
		
		## open output ASCII catalog
		f1= open(indir+outcat,"w")
		f1.write("## obj_ctr,RA,Dec,MJD,objc_type,objc_rowc,objc_colc,nchild,rowc[5],rowcerr[5],colc[5],colcerr[5],psfmag[5],psfmagerr[5],modelmag[5],modelmagerr[5],tai[5],taifrac[5],airmass[5],psf_fwhm[5],skyflux[5]\n")

		j = 0
		## write out the list of sel obj
		nsel = len(selindx)
		print "num of obj sel: %d" % nsel
		while j <= nsel-1:
## while j <= 1:
			tmpindx = selindx[j]
			outstr = str(j+1)+','+str(ra[tmpindx])+','+str(dec[tmpindx])+','+str(mjd[tmpindx])+','+str(objc_type[tmpindx])+','+str(objc_rowc[tmpindx])+','+str(objc_colc[tmpindx])+','+str(nchild[tmpindx])

			## make rowc string
			rowcstr = ''
			for tmprowc in rowc[tmpindx][:]:
				tmpstr = ','+str(tmprowc)
				rowcstr += tmpstr

			## make rowcerr string
			rowcerrstr = ''
			for tmprowcerr in rowcerr[tmpindx][:]:
				tmpstr = ','+str(tmprowcerr)
				rowcerrstr += tmpstr

			## make colc string
			colcstr = ''
			for tmpcolc in colc[tmpindx][:]:
				tmpstr = ','+str(tmpcolc)
				colcstr += tmpstr

			## make colcerr string
			colcerrstr = ''
			for tmpcolcerr in colcerr[tmpindx][:]:
				tmpstr = ','+str(tmpcolcerr)
				colcerrstr += tmpstr

			## make psf mag string
			magstr = ''
			for tmpmag in psfmag[tmpindx][:]:
				tmpstr = ','+str(tmpmag)
				magstr += tmpstr

			## make psf mag err string
			errstr = ''
			for tmperr in psferr[tmpindx][:]:
				tmpstr = ','+str(tmperr)
				errstr += tmpstr
		
			## make model mag string
			modstr = ''
			for tmpmod in modelmag[tmpindx][:]:
				tmpstr = ','+str(tmpmod)
				modstr += tmpstr

			## make psf mag err string
			moderrstr = ''
			for tmpmoderr in modelmagerr[tmpindx][:]:
				tmpstr = ','+str(tmpmoderr)
				moderrstr += tmpstr

			## make tai string
			taistr = ''
			for tmptai in tai[tmpindx][:]:
				tmpstr = ','+str(tmptai)
				taistr += tmpstr


			## make taifrac string
			taifracstr = ''
			for tmptaiFrac in taiFrac[tmpindx][:]:
				tmpstr = ','+str(tmptaiFrac)
				taifracstr += tmpstr

			## make airmass string
			amstr = ''
			for tmpam in amass[tmpindx][:]:
				tmpstr = ','+str(tmpam)
				amstr += tmpstr

			## make fwhm string
			fwstr = ''
			for tmpfw in psffwhm[tmpindx][:]:
				tmpstr = ','+str(tmpfw)
				fwstr += tmpstr

			## make skyflux string
			sfstr = ''
			for tmpsf in skyflux[tmpindx][:]:
				tmpstr = ','+str(tmpsf)
				sfstr += tmpstr


			## output string to outfile
			outstr += rowcstr+rowcerrstr+colcstr+colcerrstr+magstr+errstr+modstr+moderrstr+taistr+taifracstr+amstr+fwstr+sfstr
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

