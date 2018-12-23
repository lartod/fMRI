%-----------------------------------------------------------------------
% FMap: VDM estimaton: blip direction -1 [from A>>P]
%-----------------------------------------------------------------------
clear all;close all; clc; 
addpath('/home/common/matlab/spm12')
study_path =   '/project/3014036.01/preprocessed_data'; cd(study_path);
folder_names = dir; subject_ids = folder_names(3:end); 

for subj = 1:size(subject_ids,1)% for each subject
    
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
                all_blocks{bl} = strcat(CC,',1');
            end
             
            
            %pick up the MAG/PHASE
             dir_fmap1 = fullfile(cur_subj, '/fmap1/'); 
             files_MAG = dir(char((fullfile(dir_fmap1,'s*.nii')))); 
             MAG = {strcat(fullfile(dir_fmap1, files_MAG(1).name), ',1')};  
             dir_fmap2 = fullfile(cur_subj, '/fmap2/'); 
             files_PHASE = dir(char((fullfile(dir_fmap2,'s*.nii')))); 
             PHASE= {strcat(fullfile( dir_fmap2, files_PHASE(1).name), ',1')};
             
             %pick up anatomical
             dir_str = fullfile(cur_subj, '/anat/'); 
             anatomical_info = dir(char((fullfile(dir_str,'s*.nii')))); 
             anatomical = {strcat(fullfile(dir_str, anatomical_info(1).name), ',1')};
            
            %sometimes 2, sometimes 3 sessions:
            if size(all_blocks,2) == 3
                block1 = all_blocks{1}(1);
                block2 = all_blocks{2}(1);
                block3 = all_blocks{3}(1);
             
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = PHASE;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = MAG;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.3 6.76];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1; %1
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 27.72;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/home/common/matlab/spm12_r6470_20150506/toolbox/FieldMap/T1.nii'};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = block1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(2).epi = block2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(3).epi = block3;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = anatomical;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
                
            elseif size(all_blocks,2) == 2
                block1 = all_blocks{1}(1);
                block2 = all_blocks{2}(1);
                
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = PHASE;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = MAG;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.3 6.76];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1; %1
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 27.72;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/home/common/matlab/spm12_r6470_20150506/toolbox/FieldMap/T1.nii'};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = block1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(2).epi = block2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = anatomical;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
                
            elseif size(all_blocks,2) == 4
                block1 = all_blocks{1}(1);
                block2 = all_blocks{2}(1);
                block3 = all_blocks{3}(1);
                block4 = all_blocks{4}(1);
                
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = PHASE;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = MAG;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.3 6.76];
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1; %1
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 27.72;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/home/common/matlab/spm12_r6470_20150506/toolbox/FieldMap/T1.nii'};
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = block1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(2).epi = block2;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(3).epi = block3;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(4).epi = block4;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = anatomical;
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
                
            end
                 spm('defaults', 'fmri');
                 spm_jobman('initcfg');
                 spm_jobman('run', matlabbatch);
                 
          clear content_cur_subject cur_subject all_blocks
          disp([' Fmap is done for SUBJECT', num2str(subj), ' subjname:', subject_ids(subj).name]);

end




