%-----------------------------------------------------------------------
% Smoothing 
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
                cur_block_COMBO = dir(char(fullfile(each_subj_bl, 'wuC*.nii')));
                out =  fullfile(each_subj_bl, {cur_block_COMBO(:).name}); CC = out'; 
                tasks{bl} = strcat(CC,',1');
            end
            
    
     CC_smooth = cat(1, tasks{:});
      
    matlabbatch{1}.spm.spatial.smooth.data = CC_smooth; 
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6]; %change degree of smoothng
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);
    
    %delete all normalized
          for bl = 1:size(all_nii_blocks,2)
              each_subj_bl = fullfile(cur_subj,all_nii_blocks{bl});cd(each_subj_bl);
              delete wuC*.nii
          end

  disp([' Smoothing is done for SUBJECT', num2str(subj)]);

 end