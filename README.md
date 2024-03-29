# BOLD-CSF coupling

## The theory
The BOLD-CSF coupling method was proposed by [Fultz and colleagues (2019)](https://www.science.org/doi/full/10.1126/science.aax5440). The coupling monitors a cerebral blood volume-driven CSF inflow at the fourth ventricle, linked to fluctuations of neuronal activity. An anti-correlation is observed between cortical grey matter activity (BOLD signal) and CSF waves, in awake participants, but even stronger during sleep. The coupling between the two signals is explained by the following mechanism: activity fluctuations in the cortex give rise to cerebral blood volume changes, which in turn push CSF out, or allow for CSF inflow, due to the restricted space inside the skull. The BOLD-CSF coupling method has been applied in healthy people during sleep (Fultz et al., 2019), Alzheimer’s Disease during wake (Han et al., 2021) and Cerebral Amyloid Angiopathy during wake (Hirschler et al., in prep).

## Overview of general workflow
To retrieve the BOLD-CSF coupling, separate processing of the BOLD and CSF signal is required:

BOLD signal
- T1w-image processing: segmentation of cortical grey matter
- functional image processing: slice-time and motion correction
- register both images to the same space to retrieve the BOLD signal in the mask

CSF signal
- functional image processing: NO corrections!
- select ROI location suitable for your dataset and draw manual ROIs on the bottom slice
- retrieve the CSF signal in the manual ROIs

Both signals are normalised, detrended, and low-pass filtered. The cross-correlation is calculated by taking the peak closest to 0. Calculate the BOLD-CSF correlation and dBOLD-CSF correlation amplitude and lag for statistical analysis.

<img src="https://github.com/evavanheese/BOLD-CSF/assets/82935174/66465b47-5b9f-4e4b-97b8-77f38e4b384e" width="500">

## Required software and dependencies
We recommend the use of standardised software pipelines for the image processing steps. In this study, T1w images were processed with FreeSurfer to segment the cortex. The fMRI data were processed using HALFpipe (Waller et al., 2022). Registration was done with ANTs. 
  
## Retrieving the BOLD signal
### HALFpipe
We used [HALFpipe](https://github.com/HALFpipe/HALFpipe) to pre-process the fMRI signal and perform a slice-time and motion correction. HALFpipe can be accessed as a container through Docker/Apptainer on a server. The apptainer command to open the HALFpipe GUI our slurm server looked like this:
```
apptainer run /path/to/container/halfpipe-latest.sif --keep none --verbose
```

A GUI will open that allows you to specify settings for preprocessing. You can type as usual, but only use ‘delete’ (not ‘backspace’) to remove text. Use ‘enter’ to confirm an option and the spacebar to select or deselect options in a multiple-choice menu. Go back to the previous question using ctrl + ‘C’. We used the following preprocessing settings:
```
Specify working directory:
/path/to/workdir/

Is the data available in BIDS format?
yes

Specify data directory:
/path/to/bidsdir/
```
The GUI now prints how many T1w and BOLD images have been found in that data directory. Confirm the numbers or go back to your data.
```
Do slice timing?
yes

Missing slice acquisition direction values. Specify slice acquisition direction
Inferior to superior

Check slice timing values. Proceed with these values?
yes

Remove initial volumes from scans? 
[0]

Specify first-level features? 
No

Output a preprocessed image?
Yes

Specify image name. 
[preproc]
```
This name will be used to save your output files, if you choose to run different preprocessing settings at the same time to compare, specify the difference in this name.
```
Apply smoothing? 
No

Do grand mean scaling? 
No

Apply a temporal filter? 
No

Remove confounds? 
Motion parameters (only)

Output another preprocessed image? 
No
```
Preprocessing takes approximately 1-3 hours per subject depending on the available memory and parallel processing options.

In the specified output folder, data is saved in a ‘reports’, ‘rawdata’, and ‘derivatives’ folder. The output we use for the next steps in this analysis are saved in derivatives/halfpipe/sub-XX/func/sub-XX_task-rest_setting_{your_label}_bold.nii.gz.

### Cortical grey matter segmentation
We used [FreeSurfer recon-all](https://freesurfer.net/fswiki/recon-all) for whole-brain segmentation and took the voxels classified as grey matter in the cortical ribbon (ribbon.mgz in {subject}/mri). 

### Image registration
As the motion- and slice-time corrected functional images generated by HALFpipe are in MNI standard space ([ICBM-152 atlases](https://nist.mni.mcgill.ca/icbm-152-nonlinear-atlases-2009/)), we will register the grey matter mask (currently in T1w image space) to the same MNI space using [ANTs](https://github.com/ANTsX/ANTs)'s antsREgistrationSyN.sh function. 

An example command to register the T1w image of one participant to MNI space:
```
antsRegistrationSyN.sh -d 3 -f {MNI_template}.nii.gz -m {sub_T1}.nii.gz -o {T1_MNI_} -n 2 -t s -y 1 -e 2
```

Consequently, to apply the transformations to the grey matter mask:
```
antsApplyTransforms \
  -d 3 \
  -i {sub_GM_mask}.nii.gz \
  -r {MNI_template}.nii.gz \
  -t {T1_MNI_}1Warp.nii.gz \
  -t {T1_MNI_}0GenericAffine.mat \
  -n GenericLabel \
  -o {sub_GM_mask_MNI}.nii.gz
```

### Signal extraction
Before extracting the BOLD signal from the grey matter mask, ensure that:
1) both images (grey matter mask and functional) are of the same dimensions and orientations. We resampled the grey matter mask to the functional image using AFNI:
```
3dresample -master {sub_func_preprocessed}.nii.gz -input {sub_GM_mask_MNI}.nii.gz -prefix {sub_GM_mask_MNI_resampled}.nii.gz
```
2) you remove voxels from the mask that are outside the functional image (as they will drag the mean signal down). We multiplied the grey matter mask with the functional image using fslmaths and binarised the mask again:
```
fslmaths {sub_GM_mask_MNI_resampled}.nii.gz -mul {sub_func_preprocessed}.nii.gz {sub_GM_mask_MNI_resampled_nonzero}.nii.gz

fslmaths {sub_GM_mask_MNI_resampled_nonzero}.nii.gz -thr 0.5 -bin {sub_GM_mask_MNI_resampled_nonzero_bin}.nii.gz
```

We extracted the BOLD signal using fslmeants (which outputs a .txt file of the timeseries):
```
fslmeants -i {sub_func_preprocessed}.nii.gz -m {sub_GM_mask_MNI_resampled_nonzero_bin}.nii.gz -o [sub_output_BOLD_timeseries].txt 
```

All previously described steps to retrieve the BOLD signal (except for running HALFpipe) are performed by the [get_BOLD_timeseries.sh script](https://github.com/evavanheese/BOLD-CSF/blob/main/get_BOLD_timeseries_feb24.sh)

## Retrieving the CSF signal
The aim of this step is to draw an ROI on a location that shows sufficient CSF signal. In case of an existing dataset: depending on how consistently the field-of-view is planned during fMRI acquisition, this could vary little to greatly between participants. Preferably, the CSF ROI is NOT placed in a larger area, i.e. where CSF moves more isotropically (for example in the middle of the fourth ventricle) but rather at the bottom or top of the fourth ventricle or in the cerebral aqueduct/central canal. It is imperative that the ROI is drawn at the bottom slice (see inflow effect explained in Fultz et al. 2019, [supplementary material](https://www.science.org/doi/full/10.1126/science.aax5440#supplementary-materials)) and that the image is not processsed. 

When drawing the regions of interest:
- Select the bottom slice on the functional image, keep an anatomical scan for reference.
- Scroll through the volumes (i.e. time) and look at the CSF signal in the transverse view.
- Select at least 3 voxels that are generally included in the CSF area, preferably more.
- In case of much motion, select a somewhat larger area to make sure there is always CSF included in the mask over time.

We extracted the CSF signal using fslmeants (which outputs a .txt file of the timeseries):
```
fslmeants -i {sub_func_raw}.nii.gz -m {sub_CSF_manual_mask}.nii.gz -o [sub_output_CSF_timeseries].txt
```

## Signal processing and cross-correlation
We used Matlab functions for:
- Signal normalisation
- Detrending of the signal (using `spline_detrend` from chronux)
- Low-pass filtering (<0.1 Hz)
- Cross-correlation

Scripts to perform these steps and make BOLD-CSF output figures with the same layout and formatting as published in the original paper are available upon request.

## Recommended checks to perform along the way
- Open T1w images and grey matter masks to confirm accurate segmentation by FreeSurfer (remove mask voxels close to or inside CSF)
- Open T1w and functional images to confirm accurate registration to MNI space
- Plot raw / normalised / detrended / low-pass filtered signal to confirm these processing steps worked well
- Plot TF plots to confirm filtering
- Plot both processed signals and see whether they show any kind of coupling
- Plot cross-correlation (lags over x axis) for some subjects to assess lag (does it make sense)

## Sources
Han, F., Chen, J., Belkin-Rosen, A., Gu, Y., Luo, L., Buxton, O. M., Liu, X., & the Alzheimer’s Disease Neuroimaging Initiative. (2021). Reduced coupling between cerebrospinal fluid flow and global brain activity is linked to Alzheimer disease–related pathology. PLoS Biology, 19(6), e3001233. https://doi.org/10.1371/journal.pbio.3001233

Fultz, N. E., Bonmassar, G., Setsompop, K., Stickgold, R. A., Rosen, B. R., Polimeni, J. R., & Lewis, L. D. (2019). Coupled electrophysiological, hemodynamic, and cerebrospinal fluid oscillations in human sleep. Science, 366(6465), 628–631. https://doi.org/10.1126/science.aax5440

Waller, L., Erk, S., Pozzi, E., & Toenders, Y. J. (2022). ENIGMA HALFpipe: Interactive, reproducible, and efficient analysis for resting‐state and task‐based fMRI data. Human Brain Mapping. https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.25829
