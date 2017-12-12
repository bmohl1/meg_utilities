function indvd_meg_inversion(list,ep)

  set_spm 8 %Normalization and new segmentation methods done with SPM 8, could be upgraded
  spm_home = fileparts(which('spm'));
  addpath(fullfile(spm('dir'), 'EEGtemplates'));
  addpath(genpath(fullfile(spm('dir'), 'external','fieldtrip'))); %cannot batch without

  if nargin < 1;
    cwd = spm_select(1,'dir','Select root directory for studies',...
    '',pwd);
    cd(cwd);
    pth_subjdirs    = cellstr(spm_select([1,Inf],'dir','Select subject directories to process',...
    '',pwd));

  else
    pth_subjdirs = cellstr(list);
  end

  if ~exist('ep')
    ep='ep20'; %default
  end

  %% Start setting up the individual's script
  for n = 1:length(pth_subjdirs)
    subj_pth = pth_subjdirs{n}
    if strcmp(subj_pth(end),filesep)
      [subj_dir subj ~]     = fileparts(subj_pth(1:end-1)); % get subject id from path
    else
      [subj_dir subj ~]     = fileparts(subj_pth(1:end));
    end

    input_img = dir(strcat(subj_pth,filesep,subj,'_spm_',ep,'.mat'));
    input_img = {strcat(subj_pth,filesep,input_img.name)};

    clear matlabbatch
    spm_jobman('initcfg');
    matlabbatch{1}.spm.meeg.source.headmodel.D = input_img;
    %%
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 3;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    matlabbatch{2}.spm.meeg.source.invert.D(1) = cfg_dep;
    matlabbatch{2}.spm.meeg.source.invert.D(1).tname = 'M/EEG datasets';
    matlabbatch{2}.spm.meeg.source.invert.D(1).tgt_spec{1}.name = 'filter';
    matlabbatch{2}.spm.meeg.source.invert.D(1).tgt_spec{1}.value = 'mat';
    matlabbatch{2}.spm.meeg.source.invert.D(1).sname = 'M/EEG head model specification: M/EEG dataset(s) with a forward model';
    matlabbatch{2}.spm.meeg.source.invert.D(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.meeg.source.invert.D(1).src_output = substruct('.','D');
    matlabbatch{2}.spm.meeg.source.invert.val = 1;
    matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1;
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = [0 500];
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 0;
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.locs = [];
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
    matlabbatch{2}.spm.meeg.source.invert.modality = {'MEG'};
    savefile = strcat(subj_pth,filesep,'meg_inv_',subj);

    save(savefile, 'matlabbatch');
    spm_jobman('run',matlabbatch)
  end
