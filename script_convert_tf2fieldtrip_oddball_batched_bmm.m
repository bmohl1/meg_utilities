function script_convert_tf2fieldtrip_oddball_batched_bmm(list,labels)
  % convert all tf structs in directory to Fieldtrip tfr
  %edited to batch

  %Necessary directory
  addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/ft_meg')
  global correct

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

  nsub = size(pth_subjdirs,1);
  fprintf('The following %d subject(s) will be preprocessed:\n',nsub);
  disp(pth_subjdirs);

  if ~exist('labels','var')
    labelFid = fopen('/path_to_your_MEG_data_dir/results/labels.txt','r');
    labelFormatSpec = '%s';
    labels = textscan(labelFid,labelFormatSpec);
    fclose(labelFid);
  end

  rois = labels{1,1};

  msrs = {'nepower' 'mplf'};
  conds = {'ep10' ep20ver};

  %rois = {'rAud',...
  %              'lAud',...
  %              'rFinf',...
  %              'lFinf'...
  %              'rFmid'...
  %              'lFmid'...
  %              'rsmg'...
  %              'lsmg'...
  %              'rSpL'...
  %              'lSpL'...
  %              'raudb'...
  %    'laudb',...
  %    'rInfTL',...
  %    'raud20',...
  %    'laud20'};

  % loop through directories
  for sub=1:nsub
    fprintf('working on subject %d of %d\n',sub,nsub);

    % change working directory
    cd(pth_subjdirs(sub,:));

    % get subject id from path
    if strcmp(pth_subjdirs(sub,end),filesep)
      [pth id ext]     = fileparts(pth_subjdirs(sub,1:end-1)); % get subject id from path
    else
      [pth id ext]     = fileparts(pth_subjdirs(sub,1:end));
    end

    % all for nepower first
    % covert options


    for h = 1: length(msrs)
      msr = msrs{h};
      cfg.measure = msr;


      for i = 1:length(conds)
        cond = conds{i};

        fprintf('converting to fieldtrip...\n');
        for j = 1:length(rois)
          roi=rois{j} ;

          cfg.label = roi;
          cfg.baselinetype = 'percent';

          % load file
          file = ([id '_' cond '_' roi '_tft.mat']);
          load(file);
          tfr = meg2ft_tfr(cfg,tf);
          save([id cond '_' cfg.label '_' cfg.measure '_tfr.mat'],'tfr');
          clear tf;
          clear cfg.label;

        end
      end
    end
  end
