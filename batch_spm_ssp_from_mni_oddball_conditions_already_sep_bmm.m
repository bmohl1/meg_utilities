function batch_spm_ssp_from_mni_oddball_condition_already_sep_bmm(list,ep,coords,labels)
% batch to do ssp from spm8 results - uses spm8 routines but saves result
% in format consistent with megtools ssp.m output. Could be changed in
% future to save result as spm8 formatted file with chantype 'LFP'

% options
%Kristina's coordinates
%mni_coords = [52 -10 4;
%              -56 -12 4;
%              40 18 20;
%              -40 18 20;
%              40 8 36;
%              -40 8 36;
%              56 -44 32;
%              -56 -44 32;
%              28 -62 60;
%             -28 -62 60;

%              54 -4 0;
%                -48 -2 -14;
%                46 -56 -8;
%                54 -2 0;
%                -54 2 10];

%labels     = {'R_Aud',...
%              'L_Aud',...
%              'R_Finf',...
%              'L_Finf'...
%              'R_Fmid'...
%              'L_Fmid'...
%              'R_smg'...
%              'L_smg'...
%              'R_Spar'...
%              'L_Spar'...
%              'R_Audb',...
%              'L_Audb',...
%              'R_infTL',...
%              'R_Aud20',...
%              'L_Aud20'};

%Brianne's coordinates

mni_coords = [coords];
labels = labels{1,1};
ind        = '1'; % number of inverse

% defaults
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/spm_meg/');
spm('defaults','eeg');
spmDir = [spm('dir') filesep];

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


% loop through directories
for sub=1:nsub

    % change working directory
    cd(pth_subjdirs(sub,:));

    % get file name
    if strcmp(pth_subjdirs(sub,end),filesep)
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end-1)); % get subject id from path
    else
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end));
    end


    mstr = [id '_spm_' ep '.mat'];
    megfile = dir(mstr);
    fprintf('Working on subject: %s\n',megfile.name);

    for i=1:size(mni_coords,1)
        saveFile = [pth filesep id filesep id '_' char(labels{i}) '_' ep '_ssp.mat'];
        if ~exist (saveFile,'file');
            % read in MEG data and get some information
            D           = spm_eeg_load(megfile.name);
            nsamp       = D.nsamples;
            ntrials     = D.ntrials;
            megchn      = find(ismember(D.chantype,'MEG'));
            ssp.Q       = zeros(ntrials,nsamp);
            ssp.time    = D.time;
            ssp.epdur   = abs(D.time(1))+abs(D.time(end));
            vert        = D.inv{1}.forward.mesh.vert;
            face        = D.inv{1}.forward.mesh.face;
            norms       = spm_mesh_normals(struct('faces',face,'vertices',vert),true);
            start       = 1;
            [t0 stop]   = min(abs(ssp.time) - 0);

            % NOTE: SPMgainmatrix scaling is 1e-6 or so due to scaling of data in
            % fT and mm? Multiplying by 1e6 produces good result. See spm_eeg_lgainmat.m
            % and spm_eeg_dipole_waveforms.m

            % get vertices closest to coordinates and do ssp for each one


            fprintf('Projecting to %s\n',char(labels{i}));
            [dist meshind] = meg_spm_mni2vert(D, mni_coords(i,:));
            load(['SPMgainmatrix_' id '_spm_' ep '_' ind '.mat']);
            f              = G(:,meshind)*1e6;           % adjust leadfield scaling
            rf             = repmat(f,1,nsamp);
            ssp.W          = f;
            ssp.dipole     = mni_coords(i,:);
            for j=1:ntrials
                data           = squeeze(D(megchn,:,j)); % SPM data scale is fT
                trialproj      = dot(data,rf);
                base           = mean(trialproj(start:stop)); % baseline for subtraction
                ssp.Q(j,:)     = trialproj-base; % scale is nAm when data are fT
            end
            ssp.Q = ssp.Q/1e9; % change scale to Am from nAm for consistency with megtools ssp.m routine
            save([id '_' char(labels{i}) '_' ep '_ssp.mat'], 'ssp'); %change hard code
        end
    end
end
end
