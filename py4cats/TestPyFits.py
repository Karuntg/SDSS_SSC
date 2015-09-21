#!/usr/bin/python
##
## python script to read in Stripe 82 PhotObjAll fits tables,
## then parse for required cols
## output to ASCII csv table
 
import os,sys,numpy,pyfits,csv,time,math,string
sys.stdout.flush()

## define input and output dirs
indir = "/Users/Karun/Science/Altair/SDSScalib/SDSSdata/94/1/"
infits = "photoObj-000094-1-0012.fits"
outcat = "photoObj-000094-1-0012.cat"
outcsv = "photoObj-000094-1-0012.csv"

## open the fits table, get info
hdulist = pyfits.open(indir+infits)

## Get info regarding table
## hdulist.info()
## No.    Name         Type      Cards   Dimensions   Format
## 0    PRIMARY     PrimaryHDU      39   ()           int16   
## 1                BinTableHDU    


## read in table data to array
tbdata = hdulist[1].data

## get col names for later ref
cols = hdulist[1].columns
## print cols.names
## ['OBJID', 'PARENTID', 'FIELDID', 'SKYVERSION', 'MODE', 'CLEAN', 'RUN', 'RERUN', 'CAMCOL', 'FIELD', 'ID', 'PARENT', 'NCHILD', 'OBJC_TYPE', 'OBJC_PROB_PSF', 'OBJC_FLAGS', 'OBJC_FLAGS2', 'OBJC_ROWC', 'OBJC_ROWCERR', 'OBJC_COLC', 'OBJC_COLCERR', 'ROWVDEG', 'ROWVDEGERR', 'COLVDEG', 'COLVDEGERR', 'ROWC', 'ROWCERR', 'COLC', 'COLCERR', 'PETROTHETA', 'PETROTHETAERR', 'PETROTH50', 'PETROTH50ERR', 'PETROTH90', 'PETROTH90ERR', 'Q', 'QERR', 'U', 'UERR', 'M_E1', 'M_E2', 'M_E1E1ERR', 'M_E1E2ERR', 'M_E2E2ERR', 'M_RR_CC', 'M_RR_CCERR', 'M_CR4', 'M_E1_PSF', 'M_E2_PSF', 'M_RR_CC_PSF', 'M_CR4_PSF', 'THETA_DEV', 'THETA_DEVERR', 'AB_DEV', 'AB_DEVERR', 'THETA_EXP', 'THETA_EXPERR', 'AB_EXP', 'AB_EXPERR', 'FRACDEV', 'FLAGS', 'FLAGS2', 'TYPE', 'PROB_PSF', 'NPROF', 'PROFMEAN_NMGY', 'PROFERR_NMGY', 'STAR_LNL', 'EXP_LNL', 'DEV_LNL', 'PSP_STATUS', 'PIXSCALE', 'RA', 'DEC', 'CX', 'CY', 'CZ', 'RAERR', 'DECERR', 'L', 'B', 'OFFSETRA', 'OFFSETDEC', 'PSF_FWHM', 'MJD', 'AIRMASS', 'PHI_OFFSET', 'PHI_DEV_DEG', 'PHI_EXP_DEG', 'EXTINCTION', 'SKYFLUX', 'SKYFLUX_IVAR', 'PSFFLUX', 'PSFFLUX_IVAR', 'PSFMAG', 'PSFMAGERR', 'FIBERFLUX', 'FIBERFLUX_IVAR', 'FIBERMAG', 'FIBERMAGERR', 'FIBER2FLUX', 'FIBER2FLUX_IVAR', 'FIBER2MAG', 'FIBER2MAGERR', 'CMODELFLUX', 'CMODELFLUX_IVAR', 'CMODELMAG', 'CMODELMAGERR', 'MODELFLUX', 'MODELFLUX_IVAR', 'MODELMAG', 'MODELMAGERR', 'PETROFLUX', 'PETROFLUX_IVAR', 'PETROMAG', 'PETROMAGERR', 'DEVFLUX', 'DEVFLUX_IVAR', 'DEVMAG', 'DEVMAGERR', 'EXPFLUX', 'EXPFLUX_IVAR', 'EXPMAG', 'EXPMAGERR', 'APERFLUX', 'APERFLUX_IVAR', 'CLOUDCAM', 'CALIB_STATUS', 'NMGYPERCOUNT', 'NMGYPERCOUNT_IVAR', 'TAI', 'RESOLVE_STATUS', 'THING_ID', 'IFIELD', 'BALKAN_ID', 'NOBSERVE', 'NDETECT', 'NEDGE', 'SCORE']

## >>>>> cols reqd: 'RA', 'DEC','MJD','OBJC_ROWC','OBJC_COLC','NCHILD','PSFMAG','PSFMAGERR','TAI','AIRMASS','PSF_FWHM','SKYFLUX'

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
print 'num obj: ',nobj
nmag = (psfmag.shape)[1]
print 'num mag: ',nmag

## print psfmag.shape
## print psfmag[0:2][:]

## close fits table
hdulist.close()

## write to csv output file
with open(indir+outcsv, 'wb') as csvfile:
    datwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    datwriter.writerow(['Spam'] * 5 + ['Baked Beans'])

csvfile.close()

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
