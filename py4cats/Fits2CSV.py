#!/usr/bin/python
##
## python script to read in Stripe 82 PhotObjAll fits tables,
## then parse for required cols, output to ASCII csv table
## >>>>> cols reqd: 'RA', 'DEC','MJD','OBJC_ROWC','OBJC_COLC','NCHILD','PSFMAG','PSFMAGERR','TAI','AIRMASS','PSF_FWHM','SKYFLUX'

import os,sys,numpy,pyfits,csv,time,math,string
sys.stdout.flush()

## define input and output dirs
indir = "/Users/Karun/Science/Altair/SDSScalib/SDSSdata/94/1/"
infits = "photoObj-000094-1-0012.fits"
outcat = "photoObj-000094-1-0012.cat"

## get input dir and input fits name from command line
indir = sys.argv[1]
infits = sys.argv[2]

## create the out cat name from the fits name
fitsnam = infits.split('.')
outcat = fitsnam[0]+'.cat'

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
## print 'num obj: ',nobj
## print 'num mag: ',nmag

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

print "Done",time.strftime("%d %b %Y %H:%M:%S")
sys.stdout.flush()
sys.exit()
