function batch_OddballKM_convert_ica_fixed_bmm(list,ep)
% 1st batch run on Oddball data
% Does following:
% 1. Data conversion
% 2. Automatic bad channel detection and removal via FFT
% 3. ICA to remove eye blink and eye movement artifacts
% 4. Saves output of ica, original and trimmed datasets as well as quality
%    control picture and text of removed components

addpath('path_to_your_MEG_tools_dir/tools/fieldtrip-20130302/')

addpath('path_to_your_MEG_tools_dir/tools/FastICA_25/')
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/')
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/ft_meg/')
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/megtools/')

from_file   = 0; % use on restart after crash or force quit
%subdirs     = {'left' 'right'};
file        = 'c,rfhp0.1Hz';
dn_freq     = 55;
%run_no      = '_Run3';
ps_file     = ['Oddball_ICA_' '.ps'];
algorithm   = 'fastica'; % could be runica, fastica or binica

% defaults
ft_defaults;

% check if eeglab and fastica are properly represented on path
cfg_ica.method = algorithm;
switch algorithm
    case 'fastica'
        if ~ft_hastoolbox('fastica')
            error('You need fastica on the path to run script!');
        end
        cfg_ica.fastica.lasteig  = 25;
        cfg_ica.fastica.approach = 'symm';
    case 'runica'
        if ~ft_hastoolbox('eeglab')
            if which('eeglab')
                eeglab;
            else
                error('You need eeglab on the path to run script!');
            end
        end
        cfg_ica.runica.pca  = 25;
    case 'binica'
    if ~ft_hastoolbox('eeglab')
        if which('eeglab')
            eeglab;
        else
            error('You need eeglab on the path to run script!');
        end
    end
    cfg_ica.binica.pca  = 25;
end

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

if ~exist('ep','var')
    ep='ep10'; %default
end

h = figure('name','ICA result','color','w');
p = get(h,'position');

% loop through directories
for sub=1:nsub
    fprintf('working on subject %d of %d\n',sub,nsub);

    % change working directory
    cd(pth_subjdirs(sub,:));

    % get subject id from path
    [pth id ext] = fileparts(pth_subjdirs(sub,1:end-1));

    cnt  = get4D_bmm(file,'reference','yes');
    rescue_headfile(id,ep);
%% BAD Trigger files addition
if length(cnt.events) < 300
    fprintf('Need to replace the trial type info.\n Correcting now.\n')
    load(strcat(cwd, 'cnt_events_triggers.mat'))
    cnt.events=[cnt.events,cnt_tmp] %JUST A TEMP FIX FOR THE BAD TRIGGER FILES
end
%%
    % detect bad channels and delete them
   [bad fftrat] = fft_detect_bad_chn(cnt,2);
   if ~isempty(bad)
       cnt_orig  = deleter(cnt,bad);
       %cnt_orig = deleter (cnt_orig,{'A1','A10', 'A33','A86','A115','A172','A175','A178','A179','A204','A209'});
       %fprintf('deleting additional bad channels');
   else
       cnt_orig  = cnt; clear cnt;
       %cnt_orig = deleter (cnt,{'A1','A10', 'A33','A86','A115','A172','A175','A178','A179','A204','A209'});
   end


        % convert and denoise using highpassed ref channels so that denoising
        % operation is gentle and does not compromise lower frequency signal
        % strength
        ft                  = meg2ft(cnt_orig);
        refchans            = ft_channelselection('MEGREF',ft.label);
        cfg_ref             = [];
        cfg_ref.channel     = refchans;
        cfg_ref.hpfreq      = dn_freq;
        cfg_ref.hpfiltord   = 2;
        cfg_ref.hpfilter    = 'yes';
        ft_ref              = ft_preprocessing(cfg_ref,ft);
        cfg_denoise         = [];
        ft                  = ft_denoise_pca(cfg_denoise,ft,ft_ref);

        % do ica in Fieldtrip
        if strcmp(algorithm,'fastica')
            ft.trial = {ft.trial{1}*1e15}; % scale to ft
        end
        cfg_ica.channel         = 'MEG';
        ic_data                 = ft_componentanalysis(cfg_ica,ft);

        % plot component topography
        cfg_comp             = [];
        cfg_comp.component   = 1:25;
        cfg_comp.layout      = '4D248.lay';
        cfg_comp.comment     = 'no';
        ft_topoplotIC(cfg_comp,ic_data);

        % print to postscript file
        set(h,'position',[p(1) p(2) round(p(3)*1.5), p(4)*2],'PaperPositionMode','auto');
        axes('position',[0 0 1 1],'xlim',[0 1],'ylim',[0,1],'Box','off',...
            'Visible','off','Units','normalized','clipping','off');
        text(0.5,1,['\bf',id],'HorizontalAlignment','center',...
            'verticalalignment','top');
        print(h, '-append', '-dpsc2', fullfile(cwd,ps_file));

        % prompt at command line in MATLAB for components to remove
        prompt = {'Components to remove (comma separated list):'};
        dlg_title = 'ICA Component removal';
        num_lines = 1;
        def = {'1'};
        win = inputdlg(prompt,dlg_title,num_lines,def);
        noise = str2num(char(win(1)));
        cfg_rem             = [];
        cfg_rem.component   = noise;
        ft_rej              = ft_rejectcomponent(cfg_rem,ic_data);

        % convert back, epoch and average to view final
        if strcmp(algorithm,'fastica')
            ft_rej.trial = {ft_rej.trial{1}/1e15}; % scale to T
        end
        ft_rej.hdr  = ft.hdr;
        cnt         = ft2meg(ft_rej);

        if exist(fullfile(pth,'fiducials.mat'),'file')
            load(fullfile(pth,'fiducials.mat'));
            cnt.fiducials = fiducials;
        end

        % save original data and ica corrected data
        save([id '_orig_cnt.mat'],'cnt_orig');
        save([id '_cnt.mat'],'cnt');
        save([id '_icadat.mat'],'ic_data');

        % save numbers of components removed
        csvwrite([id '_icarej.txt'],noise);

        % delete intermediate files from binica if present
        tmpfiles=[dir('fastica*');dir('bias*');dir('temp*')];
        if ~isempty(tmpfiles)
            for tfile=1:length(tmpfiles)
                delete(tmpfiles(tfile).name);
            end
        end

        % change working directory
        cd(pth_subjdirs(sub,:));

        % clean up workspace
        %clear cnt_orig ft ft_ref ft_rej ic_data;

    % offset, epoch and visualize data, then filter
    cnt = offset (cnt);
    eps = epocher(cnt,'trigger',200,800,'OFFSET',20.5,'THRESHOLD',2500);
    %avg = offset(averager(eps));
    %figure;meg_dataplot('data',avg);
    %figure;plot_hs_sens(eps);
    %feps = filterer (eps,'band',[1 50]);
    %favg = offset(averager(feps));
    %save([id '_feps.mat'],'feps');
    save([id '_neweps.mat'],'eps');
    clear cnt;
    clear eps;
    clear cnt_orig ft ft_ref ft_rej ic_data id;


    % prompt to process another subject
    %another = menu('Continue with next subject?','Yes','No');
    %if ~another
        %cd(cwd);
        %pth_subjdirs(1:sub,:)=[]; % remove subjects already processed
        %save(['alleps_'], 'cwd','pth_subjdirs');
        %return;
    %else
        %continue;
    %end

end
