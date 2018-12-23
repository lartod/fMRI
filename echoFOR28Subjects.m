% Echo Combination
% Parts copied from Rasim's ME_Combine_GUI.m
% 2019.11.21 LT // no realignment before weighting commbination

%get all folders from preprocessed_data folder+ DELETE
%{
pathPREP = '/project/3014036.01/preprocessed_data';cd(pathPREP)
ACTIVEfolder_names = dir; PR_names = ACTIVEfolder_names(3:end); 
  for subj = 1:size(PR_names,1)% for each subject
    cur_subject = fullfile(pathPREP, PR_names(subj).name);
    content_cur_subject = dir(cur_subject); content_cur_subject = content_cur_subject(3:end);
    for cont = 1:size(content_cur_subject,1)
        out_cont(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task')); %1 if not empty
    end      
        all_content_folders = {content_cur_subject(:).name};
        all_raw_blocks = all_content_folders(logical(out_cont)); 
        
        for i = 1:size(all_raw_blocks,2)
            CUR = fullfile(cur_subject, all_raw_blocks{i});
            rmdir(CUR, 's')
        end
  end
  %}
%% 
clear all;close all; clc; 
addpath('/home/common/matlab/spm12')
%% fixed parameters
nechos = 3;               % number of echos
numOnlyWeightScan = 0;    % number of volumes to be discarded for SPM12 
numberWeightScan  = 30;   % number of volumes used for calculating the weighting
startVolume = 1;
KernelSize = 3;                 %default kernel size for smoothing
WeightVolumes = numberWeightScan;

base_path = '/project/3014036.01/raw';
root_path = '/project/3014036.01';
prepr_subj = '/project/3014036.01/preprocessed_data';cd(prepr_subj)
folder_names = dir;folder_names = folder_names(3:end); %get number of subjects

%seriesFunc_start = 31; 

for subj = 1:size(folder_names,1)% for each subject
   
    cur_subject = fullfile(base_path, folder_names(subj).name, 'ses-mri01');
    content_cur_subject = dir(cur_subject); content_cur_subject = content_cur_subject(3:end);
    %find the raw folders data: contains TaskPart, dont contain SBRef
    for cont = 1:size(content_cur_subject,1)
        out_cont(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'TaskPart')); %1 if not empty
        ref_out_cont(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'SBRef')); %1 if not empty
    end      
        raw_folders_index = out_cont - ref_out_cont;
        all_content_folders = {content_cur_subject(:).name};
        all_raw_blocks = all_content_folders(logical(raw_folders_index)); 
        

       %% 1. DICOM to .NII part (functional runs)
        disp(strcat('converting DICOM to NIFTI started for functional scans','  for Subject:  ', num2str(subj)))
        
    for bl = 1:size(all_raw_blocks,2)        
        parsed_block= strsplit(all_raw_blocks{bl}, 'Part'); 
        seq_parser = strsplit(parsed_block{1}, '-');
        block_out = strcat(seq_parser{1},'_task',parsed_block{2}); %part of .nii directory
        dir_data = fullfile(cur_subject,  all_raw_blocks{bl});
       %transfer those to .nii
        raw_files = dir(char(fullfile(dir_data, '*.IMA')));         
        out =  fullfile(dir_data, {raw_files(:).name}); 
        CC = out'; 
        
        %discard first 30 volumes: ONLY FOR first BLOCK!!
        %if bl == 1
        %   CC2 = CC(31:end,:);
        %else
         %  CC2 = CC;
        %end
        
        for i = 1:3 %check the TE values
            hdr = spm_dicom_headers(CC{i,:});
            TE(i) = hdr{1}.EchoTime;
            clear hdr
        end
        
        nii_folder = fullfile(root_path,'preprocessed_data',folder_names(subj).name, block_out); % save .nii per block depending on the # of blocks
        mkdir(nii_folder)
        %clear subj_content out_cont  all_content_folders cur_anatomical anatomical_files out
        
        matlabbatch{1}.spm.util.import.dicom.data = CC;
        matlabbatch{1}.spm.util.import.dicom.root = 'flat';
        matlabbatch{1}.spm.util.import.dicom.outdir = {nii_folder};
        matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
        
        %start the batch
        spm('defaults', 'fmri');
        spm_jobman('initcfg');
        spm_jobman('run', matlabbatch);
        
        clear nii_folder CC CC2 out parsed_block block_out dir_data raw_files
    end
        
        disp(strcat('end of DICOM conversion','  for Subject:  ', num2str(subj)))
        
        %% 2. echo volumes are realigned to the first volume of first echo
        %%and then the rest of the echo is set to the same as the first echo
        %%and all volumes are resliced        
        
        %disp('Realignment for functional volumes started')
        %fetch .niis for each run for each subject
        nii_cur_subj = fullfile(root_path,'preprocessed_data', folder_names(subj).name);
        content_cur_subject = dir(nii_cur_subj); content_cur_subject = content_cur_subject(3:end);
        
        %find the raw folders data: contains TaskPart, dont contain SBRef
        for cont = 1:size(content_cur_subject,1)
            task1(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task1')); 
            task2(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task2')); %1 if not empty
            task3(cont) = ~isempty(strfind(content_cur_subject(cont).name, 'task3')); %1 if not empty
        end      
            raw_folders_index = task1 - task2 -task3;
            all_content_folders = {content_cur_subject(:).name};
            all_nii_blocks = all_content_folders(logical(raw_folders_index)); 
            
        %FOR EACH BLOCK: 1. realign first echo to the first volume
        for bl = 1:size(all_nii_blocks,2) % for each block
           each_subj_bl = fullfile(nii_cur_subj,all_nii_blocks{bl});cd(each_subj_bl);
           filesTemp = dir('f*01.nii');%here you indicate all vols for the 1st echo
            files = char(zeros(length(filesTemp),length(filesTemp(1).name)+2,nechos));
            for i=startVolume:size(files,1)
                files(i,1:length(filesTemp(i).name),1) = filesTemp(i).name;
            end
            
            clear filesTemp
          %  spm_jobman('initcfg');
          % spm_realign(files(:,:,1)); 
            %rp file is created here (after realignment) == transformation
            %matrix
            
            %get echos 01,02, 03 in 3D : 1183 1st
            for j=2:nechos
                filesTemp = dir(['*' num2str(j) '.nii']); %% assuming number of echoes is less than 10!
               for i=startVolume:size(files(:,:,j),1)
                    files(i,1:length(filesTemp(i).name),j) = filesTemp(i).name;
                end
            end
            
            % Transformation matrices of all volumes of all echoes 
            % (except first echo) are changed to the matrix of first echo,
            % thus, realigned.
            %for i=1:size(files,1)
            %    V{1} = spm_get_space(files(i,:,1)); % matrix of the first echo
            %    for j=2:nechos
            %       spm_get_space(files(i,:,j),V{1});
            %   end
            %end

            %resliceFiles = dir('*.nii'); %% reslicing of all volumes
            %resliceFiles = char(resliceFiles.name);
            %spm_reslice(resliceFiles);
            
            %files prefixed with 'r'== transformation matrix applied to
            %echoes 02, 03;
            %mean file created with 'mean' before original files and first file
            %name
            
             %FOR EACH BLOCK: STEP 2 weight calculation, 
            % (a) extract volumes (optional, set e.g. N=30), 
            % (b) apply smoothing (optional)
            % (c) calculate weights
            % (d) apply weights to all volumes
            
            %% Smoothing of weight calculation volumes %%
            smoothingPrefix = 's';cd(each_subj_bl)
            %cd([targetPath '/converted_Volumes']);
            disp('Smoothing is applied to weight calculation volumes')
            
            %for first 30 vols of each echo
            for j=1:nechos  
                    if size(files,1) < 30
                        WeightVolumes = 15;
                    else
                        WeightVolumes = 30;
                    end
                    
                    for i=startVolume:startVolume+WeightVolumes-1
                        spm_smooth([files(i,:,j)],['s' files(i,:,j)],KernelSize);
                    end
            end  
            %files created with 's' at the beginning (not equal to total number of
            %files)
            
            %% Weight Calculation %%
            %cd([targetPath '/converted_Volumes']);
            disp('Weight calculation started...')

            dimVolume = spm_vol(files(1,:,1));
            dim = dimVolume.dim;

            clear volume4D;
            for i=1:nechos
                volume4D(:,:,:,:,i) = zeros(dim(1),dim(2),dim(3),WeightVolumes);
            end
            clear V;

            for i=startVolume:startVolume + WeightVolumes-1
                for j=1:nechos
                    V{j} = spm_vol([smoothingPrefix files(i,:,j)]);
                    volume4D(:,:,:,i-(startVolume-1),j) = spm_read_vols(V{j});       
                end
            end

            for j=1:nechos
                 tSNR(:,:,:,j) = mean(volume4D(:,:,:,:,j),4)./std(volume4D(:,:,:,:,j),0,4);
                 CNR(:,:,:,j) = tSNR(:,:,:,j) * TE(1,j); %% assuming all runs have the same TEs!!
            end

            CNRTotal = sum(CNR,4);

            for i=1:nechos
                weight(:,:,:,i) = CNR(:,:,:,i) ./ CNRTotal;
            end

            % apply weighting and combining the functional volumes
            for i=numOnlyWeightScan+1:size(files,1)
                clear V
                for j=1:nechos
                    V{j} = spm_vol([files(i,:,j)]);
                end    

                newVolume = V{1};
                newVolume.fname = ['CombinedEchos_', all_nii_blocks{bl}, '_', sprintf('%04d',i),'.nii'];           

                I_weighted = zeros(newVolume.dim);
                for j=1:nechos
                    I(:,:,:,j) = spm_read_vols(V{j});
                    I_weighted = I_weighted + I(:,:,:,j).*weight(:,:,:,j); 
                end        

                spm_create_vol(newVolume);
                spm_write_vol(newVolume,I_weighted);
                vols(:,:,:,i) = I_weighted;  %for to create a mean image of the EPI
            end

            % Create a mean EPI image
            %mean_vol = mean(vols,4);
            %newVolume.fname = ['mean_image_',all_nii_blocks{bl},'.nii'];  
            %cd(each_subj_bl);
            %spm_write_vol(newVolume,mean_vol);  % save meanEPI in the functional directory

            disp('Volumes are combined!')

            %% copy functional data and rp-file with the proper length
            % SPM realignment output file
            %cd (dir_temp)
            %FileSPMOutput = 'spm_*';
            %filename = dir(FileSPMOutput);
            %dir_data_info = fullfile(each_subj_bl, 'info');
            %if ~exist(dir_data_info,'dir'); mkdir(dir_data_info); end   % make directory for output if it does not exist
            %if isempty(filename)==0
            %    destFile=fullfile(dir_data_info,filename.name);
            %    movefile(filename.name,destFile);
            %    disp(['SPM output.ps copied for ' subj]);
            %else
            %    disp('SPM output.ps was not found. Did not start with SPM open?');
            %end


            % Realignment Parameter file (rp_ file)
            %cd (dir_temp)
            %FileRP = ['rp_*'];
            %filename = dir(FileRP);
            %[x,y,z,t1,t2,t3] = textread(filename.name, '%n %n %n %n %n %n ','delimiter','\t');
            %rp = [x,y,z,t1,t2,t3];
            %R = rp([numOnlyWeightScan+1]:length(rp),:);
            %rp_filename = ['rp_' subj all_nii_blocks{bl}];
            %save (fullfile(dir_data_info,rp_filename), 'R');
            %copyfile(filename.name, fullfile(dir_data_info,filename.name));
            %disp([ rp_filename '.mat is copied to ' dir_data_info]);

            % EPI files
            %cd (dir_temp)
            %disp(['copying MEAN functional files to ' dir_data_info]);
            %cd(each_subj_bl)    
            %filesToBeMoved = dir(['mean', '*']);
            %for n=1:size(filesToBeMoved,1);
            %    sourceFile= filesToBeMoved(n).name;
            %    destFile=fullfile(dir_data_info);
            %    copyfile(sourceFile, destFile);
            %end;

            %disp([num2str(n) ' functional files copied for subject' num2str(subj), ' ',all_nii_blocks{bl}]);
            if bl == 1
                each_subj_bl = fullfile(nii_cur_subj,all_nii_blocks{bl});
                prepscanDIR = fullfile(each_subj_bl, 'prepscans');
                mkdir(prepscanDIR)
                temp_files = dir(char(fullfile(each_subj_bl, '*.nii')));         
                allfiles =  fullfile(each_subj_bl, {temp_files(:).name}); files2move = allfiles(1:30); 
                for ii = 1:30
                    movefile(files2move{ii}, prepscanDIR)
                end
            end
            
            delete sf*.nii
            delete f*.nii
            
            clear files each_subj_bl
            disp([' EchoCombination is done for SUBJECT', num2str(subj), ' ',all_nii_blocks{bl}]);

        end % for each block
        
end 