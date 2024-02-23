#!/bin/bash
##script to register fMRI to grey matter mask & extract fMRI timeseries based on mask
##Eva van Heese
##Feb 2024, version 2.0
##see steps described on Github page: https://evavanheese.github.io/BOLD-CSF/

##input requirements
#T1w image (.nii.gz) in T1 space
#Grey matter mask of the cortex in T1 space (.nii.gz)
#processed fMRI image in MNI space (.nii.gz)
#MNI template (.nii.gz) best resolution possible (here we use 1mm)
#list of subjects to be processed subj_list.txt (01, 02, 03) in workdir

##required adjustments to this script
#change paths under define paths
#change original name t1w image in line 43 to match
#change original name GM mask in line 46 to match
#change original name functional in line 47 to match
#skip step 1 if not relevant
#adjust file names in step 2 and 3 if step 1 was skipped

##expected output for each subject
#in output folder: 
#T1w image in MNI space
#final GM mask in MNI space
#text file with timeseries extracted from GM mask

#save terminal output to get_BOLD_timeseries_log.txt
#xxx

#required software and packages
module load ANTs

##define paths
workdir="/home/anw/evanheese/my-scratch/brain_clear_narco/bold-csf/BOLD_signal/test123"
t1dir="/home/anw/evanheese/my-scratch/brain_clear_narco/bold-csf/t1_freesurfer_output" #location T1 scans
GMmaskdir="/home/anw/evanheese/my-scratch/brain_clear_narco/bold-csf/BOLD_signal/GM_masks"
fmridir="/home/anw/evanheese/my-scratch/brain_clear_narco/bold-csf/BOLD_signal/processed_rsfMRI_data"
MNI_template="/home/anw/evanheese/my-scratch/brain_clear_narco/bold-csf/BOLD_signal/mni_icbm152_t1_tal_nlin_asym_09c.nii"
module load FreeSurfer

mkdir ${workdir}/output

##read subjects from subj_list.txt and register T1 to MNI standard space
while read -r sub; do

	echo "__________working on subject ${sub}__________"
	mkdir ${workdir}/output/sub-${sub}
	
	#get all required files in the temporary dir
	mkdir ${workdir}/temp
	cp ${t1dir}/00${sub}_smri_t1/mri/T1.mgz ${workdir}/temp/sub-${sub}_t1.mgz
	
	mri_convert ${workdir}/temp/sub-${sub}_t1.mgz ${workdir}/temp/sub-${sub}_t1.nii.gz
	cp ${GMmaskdir}/00${sub}_smri_t1_bilat_GM_mask.nii.gz ${workdir}/temp/sub-${sub}_GM_mask_t1space.nii.gz
	cp ${fmridir}/sub-${sub}_task-rest_setting-preproc_bold.nii.gz ${workdir}/temp/
	cp ${MNI_template} ${workdir}/temp/

	cd ${workdir}/temp
	
	#step 1 [AFNI]: de-oblique t1 and GM mask (specific to our dataset which had oblique images - skip if not relevant)
	#echo "***deoblique t1w image and GM mask***"
	#apptainer run --writable /mnt/scratch-01/anw/share-np/AFNIr/ 3dWarp -deoblique -prefix ${workdir}/temp/sub-${sub}_t1_deobl.nii.gz ${workdir}/temp/sub-${sub}_t1.nii.gz
	#apptainer run --writable /mnt/scratch-01/anw/share-np/AFNIr/ 3dWarp -deoblique -prefix ${workdir}/temp/sub-${sub}_GM_mask_t1space_deobl.nii.gz ${workdir}/temp/sub-${sub}_GM_mask_t1space.nii.gz

	#registration t1 to MNI space
	echo "***starting registration to MNI space***"

	#step 2 [ANTs]: 3-stage registration of T1w image to MNI 2009c template	
	antsRegistrationSyN.sh -d 3 -f ${MNI_template} -m sub-${sub}_t1.nii.gz -o sub-${sub}_t1_MNI_ -n 2 -t s -y 1 -e 2

	#step 3 [ANTs]: apply warp from step 2 to GM mask
	antsApplyTransforms \
  	-d 3 \
  	-i sub-${sub}_GM_mask_t1space.nii.gz \
  	-r ${MNI_template} \
  	-t sub-${sub}_t1_MNI_1Warp.nii.gz \
  	-t sub-${sub}_t1_MNI_0GenericAffine.mat \
  	-n GenericLabel \
  	-o sub-${sub}_GM_mask_MNIspace.nii.gz

	#step 4 [AFNI]: resample to fMRI (2mm slices)
	echo "--resampling image to fmri"
	apptainer run --writable /mnt/scratch-01/anw/share-np/AFNIr/ 3dresample -master sub-${sub}_task-rest_setting-preproc_bold.nii.gz -input sub-${sub}_GM_mask_MNIspace.nii.gz -prefix sub-${sub}_GM_mask_MNIspace_resampled2fmri.nii.gz

	#step 5 [FSL]: get rid of GM mask voxels outside fmri FOV -> values are 0 and drag mean signal down (multiply mask with functional image)
	echo "----excluding 0 values outside field of view"
	fslmaths sub-${sub}_GM_mask_MNIspace_resampled2fmri.nii.gz -mul sub-${sub}_task-rest_setting-preproc_bold.nii.gz sub-${sub}_GM_mask_MNIspace_resampled2fmri_nonzero.nii.gz

	#step 9 [FSL]: binarise again
	fslmaths sub-${sub}_GM_mask_MNIspace_resampled2fmri_nonzero.nii.gz -thr 0.5 -bin sub-${sub}_GM_mask_MNIspace_resampled2fmri_nonzero_bin.nii.gz

	#step 10 [FSL]: extract timeseries
	echo "------extracting BOLD timeseries"
	fslmeants -i ${workdir}/temp/sub-${sub}_task-rest_setting-preproc_bold.nii.gz -m sub-${sub}_GM_mask_MNIspace_resampled2fmri_nonzero_bin.nii.gz -o sub-${sub}_cort_BOLD_signal_in_mask.txt

	#copy useful output to output folder
	cp sub-${sub}_t1_MNI_Warped.nii.gz sub-${sub}_GM_mask_MNIspace_resampled2fmri_nonzero_bin.nii.gz sub-${sub}_cort_BOLD_signal_in_mask.txt ${workdir}/output/sub-${sub}

	#comment the next line to avoid removing the temp folder containing output from intermediate steps
	#rm -r ${workdir}/temp
	
	echo "--------done with subject ${sub}"
done < ${workdir}/subj_list.txt
