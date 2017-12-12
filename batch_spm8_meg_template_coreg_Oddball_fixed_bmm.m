function batch_spm8_meg_template_coreg_Oddball_fixed_bmm(list,ep)
% batch file to process meg data in SPM8. The 4D files are read in, converted to
% SPM format and coregistered with the SPM8 standard MRI template. A forward model
% is also selected and stored for each subject. If continuous files are input, they
% are epoched and the epoch files are prepared. This batch can be used to prepare
% a large number of MEG files for a group MEG inversion in SPM8. Minor
% modifications could be added to allow averaging, filtering and baseline
% correction.
global correct
set_spm 8
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/spm_meg');
addpath(genpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/'));
rmpath(genpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/ft_meg'));
% options
invname =  'ODD'; % name for your inverse solution
invnum  =  1; % can change to add source analyses, but do not rerun this script to do so!
%prestim =  200; % pre and post-stim only used if input files are continuous
%poststim = 800;
usehs    = 1; % set to 1 to use headshape for coreg, 0 to use only fiducials

% defaults
spm('defaults','eeg');
%ft_defaults;
spmDir=which('spm');
spmDir=spm('dir');

%child_check = dir([spmDir, filesep, 'canonical/CHILD_TEMPLATES.txt']);
%if ~isempty(fieldnames(child_check()))
%    disp('Do you know you are using CHILD TEMPLATES for coregistration?')
%end

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

if strcmp(correct,'yes')
    ep20ver = 'correp20';
else
    ep20ver = 'ep20';
end

if ~exist('ep')
    ep='ep10'; %default
end
filestr = strcat('_', ep, '.mat');

%epoched         = spm_input('Type of Input File?',1,'b',...
%                  'Epoched|Continuous',[1 2], 1);

epoched  = 1;
              nsub = size(pth_subjdirs,1);
fprintf('The following %d subject(s) will be preprocessed:\n',nsub);
disp(pth_subjdirs);

% loop through directories
for sub=1:nsub

    fprintf('Working on subject %d of %d\n', sub, nsub);

    % change working directory
    cd(pth_subjdirs(sub,:))

    % get file name
    if strcmp(pth_subjdirs(sub,end),filesep)
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end-1)); % get subject id from path
    else
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end));
    end

    megfile          = dir ([id filestr]);

    % load filtered epoched file created using batch_OddballKM_convert_ica.m
    load ([id filestr])

    %convert to spm format
    if strcmp(ep,'ep10')
            D = meg2spm(ep10,[id, '_spm_ep10']); %Error here? did you add the path?
                clear ep10;
    else
            D = meg2spm(ep20,[id, '_spm_',ep20ver]); %BMM 4.17 added "_ep?0" so that the steps didn't have to keep being re-estimated.
                clear ep20;
    end


    % prep template MRI, not individual MRI
    D.inv                   = {struct('mesh', [])};
    D.inv{invnum}.date      = strvcat(date,datestr(now,15));
    D.inv{invnum}.comment 	= {invname}; % name of your inverse solution
    Msize = 3;  % 1 = coarse, 2 = normal 3 = fine
    sMRI  = []; % empty if using template
    D.inv{invnum}.mesh      = spm_eeg_inv_mesh(sMRI, Msize);

    % coregister template with MEG fids - note: if you convert to spm8 directly
    % from raw data, not using MEG meg2spm.m function, then you will have 5
    % fiducials, most likely, because of the 2 extra coils on head used in our
    % laboratory. Next code sets options and input fields.
    meegfid                  = D.fiducials;
    mrifid                   = D.inv{invnum}.mesh.fid;
    S                        = [];
    S.sourcefid              = meegfid;
    S.sourcefid.fid.pnt      = meegfid.fid.pnt(1:3,:); %use 1st 3 coils only
    S.sourcefid.fid.label    = meegfid.fid.label(1:3,:);
    S.targetfid              = mrifid;
    S.targetfid.fid.pnt      = [];
    S.targetfid.fid.label    = {};
    S.targetfid.fid.pnt      = mrifid.fid.pnt(1:3,:);
    S.targetfid.fid.label    = S.sourcefid.fid.label;
    S.targetfid.fod.label    = S.targetfid.fid.label(1:3, :);
    S.useheadshape           = usehs; % try setting to 0 if poor result or bad/no hs info
    S.template               = 2;
    M1                       = spm_eeg_inv_datareg(S); % co-register

    % set values of D struct from coregistration
    ind                                    = 1;
    D.inv{invnum}.datareg                  = struct([]);
    D.inv{invnum}.datareg(ind).sensors     = D.sensors('MEG');
    D.inv{invnum}.datareg(ind).fid_eeg     = S.sourcefid;
    D.inv{invnum}.datareg(ind).fid_mri     = ft_transform_headshape(inv(M1), S.targetfid);
    D.inv{invnum}.datareg(ind).toMNI       = D.inv{invnum}.mesh.Affine*M1;
    D.inv{invnum}.datareg(ind).fromMNI     = inv(D.inv{invnum}.datareg(ind).toMNI);
    D.inv{invnum}.datareg(ind).modality    = 'MEG';

    % forward model type here
    D.inv{invnum}.forward = struct([]);
    D.inv{invnum}.forward(invnum).voltype ='Single Shell'; % also 'Single Sphere'
    D=spm_eeg_inv_forward(D);
    spm_eeg_inv_checkforward(D, invnum);

    % save and go back to root directory
    D.save;
%    cd(cwd);
end
