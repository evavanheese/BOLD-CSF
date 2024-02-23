# BOLD-CSF coupling

## The theory
The BOLD-CSF coupling method was proposed by Fultz and colleagues (2019). The coupling monitors a cerebral blood volume-driven CSF inflow at the fourth ventricle, linked to fluctuations of neuronal activity. An anti-correlation is observed between cortical grey matter activity (BOLD signal) and CSF waves, in awake participants, but even stronger during sleep. The coupling between the two signals is explained by the following mechanism: activity fluctuations in the cortex give rise to cerebral blood volume changes, which in turn push CSF out, or allow for CSF inflow, due to the restricted space inside the skull. The BOLD-CSF coupling method has been applied in healthy people during sleep (Fultz et al., 2019), Alzheimer’s Disease during wake (Han et al., 2021) and Cerebral Amyloid Angiopathy during wake (Hirschler et al., in prep).

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

## Required software and dependencies
We recommend the use of standardised software pipelines for the image processing steps. In this study, T1w images were processed with FreeSurfer to segment the cortex. The fMRI data were processed using HALFpipe (Waller et al., 2022). Registration was done with ANTs. 
  
## Retrieving the BOLD signal
### HALFpipe
We used HALFpipe to pre-process the fMRI signal and perform a slice-time and motion correction. HALFpipe can be accessed as a container through Docker/Apptainer on a server. The apptainer command to open the HALFpipe GUI our slurm server looked like this:
```
apptainer run /path/to/container/halfpipe-latest.sif --keep none --verbose
```

A GUI will open that allows you to specify settings for preprocessing. You can type as usual, but only use ‘delete’ (not ‘backspace’) to remove text. Use ‘enter’ to confirm an option and the spacebar to select or deselect options in a multiple-choice menu. Go back to the previous question using ctrl + ‘C’. We used the following preprocessing settings:
```
Specify working directory: /path/to/workdir/
Is the data available in BIDS format?
yes
Specify data directory: /path/to/bidsdir/
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

## Retrieving the CSF signal

## Signal processing

## Sources
Han, F., Chen, J., Belkin-Rosen, A., Gu, Y., Luo, L., Buxton, O. M., Liu, X., & the Alzheimer’s Disease Neuroimaging Initiative. (2021). Reduced coupling between cerebrospinal fluid flow and global brain activity is linked to Alzheimer disease–related pathology. PLoS Biology, 19(6), e3001233. https://doi.org/10.1371/journal.pbio.3001233

Fultz, N. E., Bonmassar, G., Setsompop, K., Stickgold, R. A., Rosen, B. R., Polimeni, J. R., & Lewis, L. D. (2019). Coupled electrophysiological, hemodynamic, and cerebrospinal fluid oscillations in human sleep. Science, 366(6465), 628–631. https://doi.org/10.1126/science.aax5440

Waller, L., Erk, S., Pozzi, E., & Toenders, Y. J. (2022). ENIGMA HALFpipe: Interactive, reproducible, and efficient analysis for resting‐state and task‐based fMRI data. Human Brain Mapping. https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.25829
