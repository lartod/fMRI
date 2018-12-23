Preprocessing pipeline:

1. transfer anatomical and field maps from .IMA 2 .nii>>> raw2nii_anat_job.m + raw2nii_fieldmap.m

2. transfer functionals from .IMA 2 .nii + combine within each block 
>>> echo_combination.m (C*.nii)
<first 30 .nii deleted for 1st task only>

3. fieldmap correction, part A >>> [NB. 1 image = phase; 2 images = magnitude with long and shirt TE; 4.30; 6.76ms;
https://www.mccauslandcenter.sc.edu/crnl/tools/fieldmap
Total EPI read out = echo spacing * # of phase encoding lines / acceleration factor
.vdm == presubstracted ]
fmap_estimation (VDM calculation) >>> fieldmap_vdmEstimation.m (vdms per session)

4. realign all images between blocks (to the 1st image of the first block, and all within a block to the 1st) 
AND unwarp (apply wdm here) >>> realign_and_unwarp_job.m (uC*.nii)
<SNR check>

5. coregistration (anat/to the mean image) >>> coregistration.m (cors*.nii) 
6. make a grey matter mask (segment coregistered structural) >>> segmentation_job.m (c1cors*.nii)

7. normalization >>> normalization.m (wuC*.nii)
8. smoothing >>> smoothing.m (swuC*.nii)

THE END