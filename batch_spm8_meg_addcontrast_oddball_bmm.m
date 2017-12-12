function batch_spm8_meg_addcontrast_oddball_bmm(list,ep)
% this batch will take directory listing of source analyzed MEG
% files, add a contrast and produce a new set of smoothed
% images for the contrast.

% options to change as necessary
invnum  =  1; % can change to add source analyses, but do not rerun this script to do so!
cwoi    = [0 500]; % time window of interest for contrast/image in ms
ctype   = 'evoked'; % can be evoked or induced if adding to epoch file, only evoked for average
cfboi   = [0]; % time-frequency to analyze, 0 for ERP all frequencies, [start stop] for tf window
csmooth = 8; % fwhm in mm for smoothing contrast into volume


% Select working directories

if exist('list')
    pth_subjdirs = list;
else
    
cwd             = spm_select(1,'dir','Select root directory for studies',...
                  '',pwd);
cd(cwd);
pth_subjdirs    = spm_select([1,Inf],'dir','Select subject directories to process',...
                  '',pwd);
end
nsub = size(pth_subjdirs,1);
fprintf('The following %d subject(s) will be preprocessed:\n',nsub);
disp(pth_subjdirs);

if ~exist('ep')
    ep='ep10'; %default
end
filestr = ['_spm_' ep '.mat']; % change if you use a different file naming convention than id_spm.mat

% loop through directories
for sub=1:nsub
    fprintf('Processing file %d of %d\n', sub, nsub);
      
    % change working directory
    cd(pth_subjdirs(sub,:))
    
    %get mat name from dat name to avoid other mat files
    if strcmp(pth_subjdirs(sub,end),filesep)
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end-1)); % get subject id from path
    else
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end));
    end
    megfile             = dir([id filestr]);
    
    %load D structure from mat file
    D = spm_eeg_load(megfile.name);
    
    if isfield(D,'inv')
        %create new contrast here (assumes at least one contrast exists)
        D.inv{D.val}.contrast.woi   = cwoi; % time window to use
        D.inv{D.val}.contrast.fboi  = cfboi; % 0 = no tft, [5 50] = tft over 5-50 Hz
        D.inv{D.val}.contrast.type  = ctype;
        D                           = spm_eeg_inv_results(D);
    
        %create images and save results
        D.inv{D.val}.contrast.smooth  = csmooth;
        D.inv{D.val}.contrast.display = 0;
        D.inv{D.val}.contrast.space   = 1;
        D                             = spm_eeg_inv_Mesh2Voxels(D);
        
        %save file
        save([id, ep,'_spm_1_t0_500_f_1.mat'], 'D'); %BMM changed nii to mat, b/c saving the D structure to a structure doesn't feed forward to 2nd level analysis 
        %movefile([id,'_spm_1_t0_500_f_1.nii'],[id,ep,'_spm_1_t0_500_f_1.nii']);%BMM...
        %to save the "1st"-level contrast. Moving function no longer needed, though, b/c
        %the ep is in the id.
        
    else
        disp('No inverse structure in file!');
    end
    
    clear D;
end