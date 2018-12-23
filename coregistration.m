%this is coregistration script that takes mean image after realignment and
%coregisters anatomical to it

clear all;close all; clc; 
addpath('/home/common/matlab/spm12')
study_path = '/project/3014036.01/preprocessed_data'; cd(study_path);
folder_names = dir; subject_ids = folder_names(3:end); 

for subj = 2 %:size(subject_ids,1)% for each subject
    
    cur_subj = fullfile(study_path, subject_ids(subj).name);
    content_cur_subject = dir(cur_subj); content_cur_subject = content_cur_subject(3:end);
        
        for cont = 1:size(content_cur_subject,1)
            task1(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task1')); 
            task2(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task2')); %1 if not empty
            task3(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task3')); %1 if not empty
        end      
            raw_folders_index = task1 - task2 -task3;
            all_content_folders = {content_cur_subject(:).name};
            all_nii_blocks = all_content_folders(logical(raw_folders_index)); 

            %go to the ref path (only first block) and read out meanImage
            dirREF = fullfile(cur_subj, all_nii_blocks(1));    
            file_REF = dir(char((fullfile(dirREF,'mean*.nii'))));   
            ref_path = fullfile(dirREF, file_REF(1).name);
            
            %go to the source image (struct T1)
            dirSOURCE = fullfile(cur_subj, '/anat/');  %/T1/  
            file_SOURCE  = dir(char((fullfile(dirSOURCE ,'s*.nii'))));   
            struct_path = {fullfile(dirSOURCE , file_SOURCE (1).name)}; 
            
            %go to C1
            %dirSOURCE = fullfile(cur_subj, '/anat/');   
            %file_C1 = dir(char((fullfile(dirSOURCE ,'c1*.nii'))));   
            %C1_path = {fullfile(dirSOURCE , file_C1(1).name)}; 

        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = ref_path;
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = struct_path;
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'cor';

        
         spm('defaults', 'fmri');
         spm_jobman('initcfg');
         spm_jobman('run', matlabbatch);
         
        disp([' Coregistration is done for SUBJECT', num2str(subj)]);


end