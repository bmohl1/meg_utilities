%batch tf stats on already done tf analysis (see QTF.m and batch_tf.m)
%change dir to tf directory and run after setting options below

%% OPTIONS TO SET %%

% NOTE: these woi points need to be in the range of the tf file output (i.e., you can't specify a
% frequency range of 20 - 50 if you only have 10-50 in your tf files!
fwoi     = [30 50]; % frequency window of interest for statistics
% e.g., [40 40] for single frequency
twoi     = [.2 .5]; % time window of interest for statistics -
% e.g., [50 50] for single timepoint

% scaling option - note that tft program outputs data in raw amplitude units (e.g., fT, uV or nA-m)
% if you want a different scale for stats, enter it here
%scale = '.^2)';
scale    = '*1e9)'; % scale to nA-m from A-m, for example: 1e9
%scale    = '*1e15)';
%scale     = '*1)'; % uncomment for no scaling
cd ('/path_to_your_MEG_data_dir/tft_files');
% END OF OPTIONS TO SET %

%% Load the available tft files
%rois = {'rAud',...
%    'lAud',...
%    'rFinf',...
%    'lFinf'...
%    'rFmid'...
%    'lFmid'...
%    'rsmg'...
%    'lsmg'...
%    'rSpL'...
%    'lSpL'...
%    'raudb'...
%    'laudb',...
%    'rInfTL',...
%    'raud20',...
%    'laud20'};


labelFid = fopen('/path_to_your_MEG_data_dir/results/child_labels.txt','r');
labelFormatSpec = '%s';
labels = textscan(labelFid,labelFormatSpec);
fclose(labelFid);
rois = labels{1,1};

%names from the Oddball_tft script
eps = {'ep10', 'ep20'};
for iEp = 1:length(eps)
  ep = eps{iEp}
  for iRoi = 1:length(rois)
    roi = rois{iRoi};
    prefix =  [ep, '_', roi]
    tftNames = strcat('*_',prefix,'_tft.mat');
    files    = dir(tftNames); % change this if you altered file naming convention
    output   = [prefix,'_200_500.txt']; %filename for output

    %% BEGIN SCRIPT %%
    fid = fopen(output,'w');
    fprintf(fid,...
    'ID\tmPLF\tmtpower\tmepower\tmipower\tmntpower\tmnepower\tmnipower\tmbaseline\tPeakHz\n');

    stats = zeros(length(files), 12); %array for group results
    for i=1:length(files)
      %read the file
      if length(files) == 1
        cur = files.name;
      else
        cur = files(i).name;
      end
      fprintf('Processing %s\n',cur);

      %load file
      load(cur);

      %get nearest indices of requested time frequency windows
      if fwoi(1) < min(tf.freq) || fwoi(2) > max(tf.freq)
        disp('You must specify a frequency range within your actual data range!');
        return;
      end
      if twoi(1) < min(tf.time) || twoi(2) > max(tf.time)
        disp('You must specify a time range within your actual data range!');
        return;
      end
      [diff fstart] = min(abs(tf.freq - fwoi(1)));
      [diff fstop]  = min(abs(tf.freq - fwoi(2)));
      [diff tstart] = min(abs(tf.time - twoi(1)));
      [diff tstop]  = min(abs(tf.time - twoi(2)));
      [diff onset]  = min(abs(tf.time));

      fvals            = tf.freq(fstart:fstop);
      tvals            = tf.time(tstart:tstop);

      %compute mean stats for window on all measures
      mmplf   = mean(mean(tf.mplf(fstart:fstop,tstart:tstop)));
      mtpower = eval(['mean(mean(tf.tpower(fstart:fstop,tstart:tstop)' scale ')']);
      mepower = eval(['mean(mean(tf.epower(fstart:fstop,tstart:tstop)' scale ')']);
      mipower = eval(['mean(mean(tf.ipower(fstart:fstop,tstart:tstop)' scale ')']);
      mntpower = mean(mean(tf.ntpower(fstart:fstop,tstart:tstop)));
      mnepower = mean(mean(tf.nepower(fstart:fstop,tstart:tstop)));
      mnipower = mean(mean(tf.nipower(fstart:fstop,tstart:tstop)));
      mbaseline = eval(['mean(mean(tf.tpower(1:onset))' scale]);
      [val ind]  = max(max(tf.nepower(fstart:fstop,tstart:tstop),...
      [],2));
      peakf      = fvals(ind);

      fprintf(fid,'%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n',...
      strtrim(cur),mmplf,mtpower,mepower,mipower,mntpower,mnepower,mnipower,...
      mbaseline,peakf);
    end
    fclose('all');
  end
end
