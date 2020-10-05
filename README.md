# SDSS_SSC
Standard star catalogs with SDSS data
This repository contains all the Python code used to generate the light curves for standard stars observared as part the SDSS Stripe 82 photometric data. 

# Analysis_2020
## Create a new dir to hold the notebooks from the 2020 analysis
Uploading the following notebooks
1. Match_NvsO_SSCv2.ipynb
This is the code to match the 2007 and 2020 catalogs. See the selection cuts applied to the 2020 catalog

2. sdss_gaia_matching_IZv2.ipynb
This is ZI's original code to match 2007 SSC and the Gaia catalog. I have reorganized the code somewhat

3. NEWsdss_gaia_matching_v0.ipynb
This is my version of ZI's code for Gaia matching


### How-to produce paper plots ### 

A. Note that SDSS_SSC/Data directory is intentionally empty in GitHub repo because data files 
   are rather large (~1 GB). Instead of storing them on GitHub, they need to be downloaded
   to this directory from the following links:

http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/stripe82calibStars_v2.6.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/N2020_stripe82calibStars.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/stripe82calibStars_v3.1_noheader_final.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/stripe82calibStars_v3.1.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/stripe82calibStars_v3.4.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/Stripe82_GaiaDR2.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/Stripe82_GaiaDR2_BPRP.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/cfis_stripe82.dat
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/match_stripe82calibStars_v3.2_des_dr1_cleaned.csv
http://faculty.washington.edu/ivezic/sdss/calib82/dataV2/nSSC2PS1_matched.csv
 

B. Execute notebooks stored in Analysis_2020 directory from PaperPlots directory:

1. recalibration_gray.ipynb 
      astroVSpm_RA_pm.png
      astroVSpm_RA_r.png
      GrVSgi.png
      GmagCorrection_RA_Hess.png  
      GmagCorrection_Dec_Hess.png	
      GmagCorrection_Gmag_Hess.png
  
2. recalibration_colors.ipynb  
      colorResidGaiaColorsB_ri_Dec_Hess.png
      colorResidGaiaColors_gr_RA_Hess.png
      colorResidGaiaColors_gr_Dec_Hess.png

3. test26vs33.ipynb 
      testV26vsV33_r_dr_r_mag_Hess.png 
      testV26vsV33_r_dr_Dec_Hess.png 
      testV26vsV33_r_w_old_RA_Hess.png
      testV26vsV33_r_w_new_RA_Hess.png
      testV26vsV33_r_w_old_Dec_Hess.png
      testV26vsV33_r_w_new_Dec_Hess.png
      testV26vsV33_snew_u_s_new_RA_Hess.png
      testV26vsV33_snew_u_s_new_Dec_Hess.png  
      testV26vsV33_ynew_z_y_new_RA_Hess.png 
      testV26vsV33_ynew_z_y_new_Dec_Hess.png 

4. compareDES.ipynb  
      colorResidDES2bright_dr_RA_Hess.png  
      colorResidDES2bright_di_RA_Hess.png 
      colorResidDES2bright_dz_RA_Hess.png 
      colorResidDES2bright_dr_Dec_Hess.png  
      colorResidDES2bright_di_Dec_Hess.png  
      colorResidDES2bright_dz_Dec_Hess.png 
      colorResidDES2_dr_rmag_Hess.png  
      colorResidDES2_di_rmag_Hess.png  
 
5. comparePanSTARRS.ipynb
      colorResidPSbright_dr_RA_Hess.png
      colorResidPSbright_di_RA_Hess.png
      colorResidPSbright_dz_RA_Hess.png
      colorResidPSbright_dr_Dec_Hess.png
      colorResidPSbright_di_Dec_Hess.png
      colorResidPSbright_dz_Dec_Hess.png
   
6. compareCFIS.ipynb
      colorResidCFISug_Dec_Hess.png 


=== Notebooks for figures selected for the paper (same order as in tex file) ===  
astroVSpm_RA_pm.png: recalibration_gray.ipynb
astroVSpm_RA_r.png: recalibration_gray.ipynb 
GrVSgi.png: recalibration_gray.ipynb
GmagCorrectionTest_Gmag_Hess.png: recalibration_gray.ipynb      
GmagCorrection_RA_Hess.png: recalibration_gray.ipynb 
GmagCorrection_Dec_Hess.png: recalibration_gray.ipynb 
colorResidGaiaColorsB_ri_Dec_Hess.png: recalibration_colors.ipynb  
testV26vsV33_r_dr_r_mag_Hess.png: test26vs33.ipynb  
testV26vsV33_r_dr_Dec_Hess.png: test26vs33.ipynb   
testV26vsV33_r_w_old_RA_Hess.png: test26vs33.ipynb  
testV26vsV33_r_w_new_RA_Hess.png: test26vs33.ipynb  
testV26vsV33_r_w_old_Dec_Hess.png: test26vs33.ipynb  
testV26vsV33_r_w_new_Dec_Hess.png: test26vs33.ipynb  
testV26vsV33_snew_u_s_new_RA_Hess.png: test26vs33.ipynb  
testV26vsV33_snew_u_s_new_Dec_Hess.png: test26vs33.ipynb    
testV26vsV33_ynew_z_y_new_RA_Hess.png: test26vs33.ipynb   
testV26vsV33_ynew_z_y_new_Dec_Hess.png: test26vs33.ipynb   
colorResidGaiaColors_gr_RA_Hess.png: recalibration_colors.ipynb  
colorResidGaiaColors_gr_Dec_Hess.png: recalibration_colors.ipynb
colorResidDES2bright_dr_RA_Hess.png: compareDES.ipynb 
colorResidPSbright_dr_RA_Hess.png: comparePanSTARRS.ipynb
colorResidDES2bright_di_RA_Hess.png: compareDES.ipynb 
colorResidPSbright_di_RA_Hess.png: comparePanSTARRS.ipynb
colorResidDES2bright_dz_RA_Hess.png: compareDES.ipynb 
colorResidPSbright_dz_RA_Hess.png: comparePanSTARRS.ipynb
colorResidDES2bright_dr_Dec_Hess.png: compareDES.ipynb 
colorResidPSbright_dr_Dec_Hess.png: comparePanSTARRS.ipynb
colorResidDES2bright_di_Dec_Hess.png: compareDES.ipynb  
colorResidPSbright_di_Dec_Hess.png: comparePanSTARRS.ipynb
colorResidDES2bright_dz_Dec_Hess.png: compareDES.ipynb 
colorResidPSbright_dz_Dec_Hess.png: comparePanSTARRS.ipynb
colorResidDES2_dr_rmag_Hess.png: compareDES.ipynb 
colorResidDES2_di_rmag_Hess.png: compareDES.ipynb 
colorResidCFISug_Dec_Hess.png: compareCFIS.ipynb 


C. All *png files are removed by hand from Analysis_2020 directory and are NOT uploaded to GitHub 


