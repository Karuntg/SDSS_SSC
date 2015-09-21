#!/usr/bin/python
##
## python script to download fits stamp images from SDSS DR7 
## centered on ra, dec of SNe Host galaxies
 
import os,sys,numpy,sqlcl_dr9,sqlcl_dr7,time,math,string

print "Querying for SDSS SNe Host Images ...",time.strftime("%d %b %Y %H:%M:%S")
sys.stdout.flush()

## define input and output dirs
indir = "/copiapo1/karun/SNe_yan/SDSSSNeHosts/Images/Ver1/"
incat = "HostGalImCat.dat"
outdir = "/copiapo1/karun/SNe_yan/SDSSSNeHosts/Images/Ver1/"

# Read SNe host galaxy cat
## index, SDSS ID, ra, dec, petro_r, R25, SNeHostdist_asec, SNeHostdist_pix 
f1= open(indir+incat,"r")

# input arrays
indx =[0]*2000
HostID =[0]*2000
Hostra =[0]*2000
Hostdec =[0]*2000
Hostmag =[0]*2000
HostR25 =[0]*2000
Distasec =[0]*2000
Distpix =[0]*2000
nHost= 0
nlin = 0

for line in f1:
	data = line.split(' | ')
	indx[nHost] = int(data[0])
	HostID[nHost] = int(data[1])
	Hostra[nHost]  = float(data[2])
	Hostdec[nHost]  = float(data[3])
	Hostmag[nHost]  = float(data[4])
	HostR25[nHost]  = float(data[5])
	Distasec[nHost]  = float(data[6])
	Distpix[nHost]  = float(data[7])
	nHost += 1
f1.close()
print "Number of Host read in: ",nHost

# trim input vectors
HostID = HostID[0:nHost]
Hostra = Hostra[0:nHost]
Hostdec = Hostdec[0:nHost]
Hostmag = Hostmag[0:nHost]
HostR25 = HostR25[0:nHost]
Distasec = Distasec[0:nHost]
Distpix = Distpix[0:nHost]

## for each host ra, dec 
## get stamp image
## keep track of processed hosts, out fits names, etc
outlog = outdir+'SNeHostsIm.log'
f1= open(outlog,"w")
f1.write("## ctr,Host ID, fits image\n")
j = 0
while j <= nHost-1:
## while j <= 1:
	objid_str = str(HostID[j])
	print objid_str
	lines = sqlcl_dr7.query("select objID,ra,dec,run,rerun,camcol,field,obj from BESTDR7..PhotoObjAll where objID="+objid_str).readlines()
	print "number of dr7 lines read: ",len(lines)
	print lines

	hostname = os.uname()[1]
	for line in lines[1:]:
    		data = line.split(',')
    		objID = data[0]
    		ra = float(data[1])
    		dec = float(data[2])
    		run = int(float(data[3]))
    		rerun = int(float(data[4]))
    		camcol = int(float(data[5]))
    		field = int(float(data[6]))
    		ID = int(float(data[7]))

    		try:
                        		pix_pos = sqlcl_dr7.query("select rowc_g,colc_g from BESTDR7..PhotoPrimary where objID = "+objID).readlines()
    		except:
                        		os.system("sleep 902")
    		while len(pix_pos) != 2:
        			os.system("sleep 902")
        			try:
                            			pix_pos = sqlcl_dr7.query("select rowc_g,colc_g from BESTDR7..PhotoPrimary where objID = "+objID).readlines()
        			except:
                            			continue

		## make wget command
    		for line2 in pix_pos[1:]:
                        		rowc_g = string.split(line2,",")[0]
                        		colc_g = string.split(line2,",")[1][:-1]
    		print rowc_g,colc_g
    		tmp_run = str(run)
    		while len(tmp_run) < 6: tmp_run = "0"+tmp_run 
    		tmp_field = str(field)
    		while len(tmp_field) < 4: tmp_field = "0"+tmp_field
    		sdss_path_c = "imaging/"+str(run)+"/"+str(rerun)+"/corr/"+str(camcol)+"/"

		## get g-band image
    		## sdss_corr_name_g = "fpC-"+tmp_run+"-g"+str(camcol)+"-"+tmp_field+".fit"
    		## wgetcmd = 'wget http://das.sdss.org/'+sdss_path_c+sdss_corr_name_g+'.gz'
		## get r-band image
    		sdss_corr_name_r = "fpC-"+tmp_run+"-r"+str(camcol)+"-"+tmp_field+".fit"
    		wgetcmd = 'wget http://das.sdss.org/'+sdss_path_c+sdss_corr_name_r+'.gz'
		print wgetcmd

		## run wget, gunzip downloaded *.fit.gz file
    		os.system(wgetcmd)
    		## os.system("gunzip "+sdss_corr_name_g+'.gz')
    		os.system("gunzip "+sdss_corr_name_r+'.gz')
		jplus1_str = str(j+1)
		## f1.write(jplus1_str+"  "+objid_str+"  "+sdss_corr_name_g+"\n")
		f1.write(jplus1_str+"  "+objid_str+"  "+sdss_corr_name_r+"\n")

		## increment while counter
		j += 1

## all done, close and quit
f1.close()
print "Done",time.strftime("%d %b %Y %H:%M:%S")
sys.stdout.flush()
sys.exit()
