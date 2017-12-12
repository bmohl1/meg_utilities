%1.	batch_OddballKM_convert_ica_fixed.m
%2.	Oddball_batch_epoch_extractor.m


%% If you want to separate out correct trials, set correct to 'yes'
clear -global correct
%global correct
%correct = 'yes';

%%
if exist('correct','var')
    ep20ver = 'correp20';
else
    ep20ver = 'ep20';
end

%% Toggle switches for variations on analysis themes
%Break up which part of the script runs
preproc = ('full'); % options: full, first, next, stats, depending on which parts of the pipeline you want to (re)run.
types = { 'ep10' ep20ver}; % set names for variables later in script
grpPrefix = 'child'; %what did you put in front of the coords.txt and label.txt files?

%% Select the files to analyze
cwd     = spm_select(1,'dir','Select root directory for studies',...
    '',pwd);
cd(cwd);
list    = spm_select([1,Inf],'dir','Select subject directories to process',...
    '',pwd);

%% Begin pre-processing
if strcmpi (preproc, 'full') || strcmpi (preproc, 'first')
    disp('Checking coreg and estimation')
    for s = 1:size(list,1)
        subj = list(s, 1:end);
        for iType = 1:length(types)
            ep = types{iType};

            touch = [subj,filesep, 'touch_',ep,'_estCorrTrials.txt'];
            if ~exist (touch, 'file')
                fprintf('Crunching numbers for %s\n', ep);
                rescue_headfile(subj,ep); %checks that the fiducials are available
                batch_spm8_meg_template_coreg_Oddball_fixed_bmm (subj,ep);
                indvd_meg_inversion(subj,ep);
                batch_spm8_meg_addcontrast_oddball_bmm (subj,ep);
                fclose(fopen(touch, 'w+'))
                close all
            end
        end
    end
end

%% Define where to get coords and labels, rather than hard coding in multiple scripts
if ~strcmpi (preproc, 'first')
  fidPath  = 'path_to_your_MEG_dir/results/';
    labelFid = fopen([fidPath, grpPrefix,'_labels.txt'],'r');
    labelFormatSpec = '%s';
    labels = textscan(labelFid,labelFormatSpec);
    fclose(labelFid);
end

%% SSP, if completing 'full' or 'next' processing pipelines
if strcmpi (preproc, 'stats') || strcmpi (preproc, 'first')
    disp('Moving to last section')
else
    coordsFid = fopen([fidPath, grpPrefix,'_coords.txt'],'r');
    coordsFormatSpec = '%d';
    sizeC = [3 Inf];
    coords = fscanf(coordsFid,coordsFormatSpec,sizeC)';
    fclose(coordsFid);

    %% Continue processing

    for iType = 1:length(types)
        ep = types{iType};
        batch_spm_ssp_from_mni_oddball_conditions_already_sep_bmm(list,ep,coords,labels);
        Oddball_batch_tft_bmm(list,coords,labels, ep);
    end
    script_convert_tf2fieldtrip_oddball_batched_bmm (list,labels)

end

%% Generate the statistical tests on the TF plots.
if strcmpi(preproc, 'first') || strcmpi(preproc, 'next')
    disp('"Last" section.');
    disp('Presuming you wanted to run alternate stats or to check something.')
    disp('Exiting without running stats.')
else
    %% Final stats
    script_clusteranalysis_ft_oddball_condCompNoCorr_bmm_grps(labels)
end
