%-----------------------------------------------------------------------
% Normalize: realigned + coregistered structural
%-----------------------------------------------------------------------
clear all;close all; clc; 
addpath('/home/common/matlab/spm12')
study_path =   '/project/3014036.01/preprocessed_data'; cd(study_path);
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
            
            %for all the data in the existing blocks (1-2 or 1-3)
            for bl = 1:size(all_nii_blocks,2)
                each_subj_bl = fullfile(cur_subj,all_nii_blocks{bl});cd(each_subj_bl);
                cur_block_COMBO = dir(char(fullfile(each_subj_bl, 'uC*.nii')));
                out =  fullfile(each_subj_bl, {cur_block_COMBO(:).name}); CC = out'; 
                tasks{bl} = strcat(CC,',1');
            end
             
             dir_str = fullfile(cur_subj, '/anat/'); 
             files_str = dir(char((fullfile(dir_str,'cors*.nii')))); %pick the rs
             struc_cor = {strcat(fullfile(dir_str, files_str(1).name), ',1')};              
                
      %combine them all
      combined  = [tasks,{struc_cor}];
      ALLCOMBO = cat(1, combined{:});
      
      matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = struc_cor;  %{'/home/predatt/lartod/fMRI/FaceGender/SUBJECT02/STRUCTURAL/T1/rs220151205.nii,1'};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = ALLCOMBO;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/home/common/matlab/spm12/tpm/TPM.nii'};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                     78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

 spm('defaults', 'fmri');
 spm_jobman('initcfg');
 spm_jobman('run', matlabbatch);
        
         %delete all realigned images
          for bl = 1:size(all_nii_blocks,2)
              each_subj_bl = fullfile(cur_subj,all_nii_blocks{bl});cd(each_subj_bl);
              delete uC*.nii
          end
          
 disp([' Normalization is done for SUBJECT', num2str(subj)]);
       
end