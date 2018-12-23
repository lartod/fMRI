%-----------------------------------------------------------------------
% This script transsfers raw .IMA to .nii for anatomical scan per each
% subject
%-----------------------------------------------------------------------
%%
clear all; close all; clc;

addpath('/home/common/matlab/spm12'); 
study_path =  '/project/3014036.01/raw';cd(study_path);
target_out = '/project/3014036.01/preprocessed_data';

%get subject folder
folder_names = dir; subject_ids = folder_names(3:end); 

for subj= 26 %:size(subject_ids,1) % 31 subjects
   
    
    data_path=fullfile(study_path, subject_ids(subj).name, 'ses-mri01');
    subj_content = dir(data_path);subj_content = subj_content(3:end);
     
    disp(strcat('Subj_counter:', '  ', num2str(subj), '   Subj_name:', '  ', num2str(subject_ids(subj).name)))
    
        %find mprage
        for cont = 1:size(subj_content,1)
            out_cont(cont) = ~isempty(strfind(subj_content(cont).name, 'mprage')); %1 if not empty
        end        
        all_content_folders = {subj_content(:).name};
        cur_anatomical = all_content_folders(logical(out_cont));        
        anatomical_files = dir(char(fullfile(data_path, cur_anatomical{1}, '*.IMA')));
        
        str2append = fullfile(data_path, cur_anatomical{1});                 
        out =  fullfile(str2append, {anatomical_files(:).name}); CC = out'; 

        target_folder = fullfile(target_out,subject_ids(subj).name, 'anat');
        mkdir(target_folder)
        clear subj_content out_cont  all_content_folders cur_anatomical anatomical_files out
        
matlabbatch{1}.spm.util.import.dicom.data = CC;
%%
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {target_folder};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

%start the batch
spm('defaults', 'fmri');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end
