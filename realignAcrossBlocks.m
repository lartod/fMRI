%-----------------------------------------------------------------------
% Realigh the sessions (2 or 3 depending on the Subject) to the first .nii
% of the first block
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
                cur_block_COMBO = dir(char(fullfile(each_subj_bl, 'C*.nii')));
                out =  fullfile(each_subj_bl, {cur_block_COMBO(:).name}); CC = out'; 
                tasks{bl} = CC;
            end
                
                matlabbatch{1}.spm.spatial.realign.estwrite.data = tasks;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % to the first of the first session
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                
                %start the batch
                spm('defaults', 'fmri');
                spm_jobman('initcfg');
                spm_jobman('run', matlabbatch);
                
         %here you get the MEAN image (it is only in the first session)
         
                    %delete all Combined images
                    for bl = 1:size(all_nii_blocks,2)
                            each_subj_bl = fullfile(cur_subj,all_nii_blocks{bl});cd(each_subj_bl);
                            %cur_block_COMBO = dir(char(fullfile(each_subj_bl, 'C*.nii')));
                            delete C*.nii
                    end

               disp([' RealAcrossBlocks is done for SUBJECT', num2str(subj)]);

end