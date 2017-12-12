function script_clusteranalysis_ft_oddball_condCompNoCorr_bmm_grps (labels)
% script to test cluster analysis approach on converted tfr data
% Both ep10 and ep20 together; mplf and nepower separate
addpath('path_to_your_MEG_tools_dir/tools/fieldtrip-20130302');
set_spm 8
covary = 'yes';
covaryLabel = '_covaryAge';%_covaryAge
cwd = '/path_to_your_MEG_data_dir';
cd (cwd);
ft_defaults;

global correct

if ~exist('labels','var');
    labelFid = fopen('/path_to_your_MEG_data_dir/child_labels.txt','r');
    labelFormatSpec = '%s';
    labels = textscan(labelFid,labelFormatSpec);
    fclose(labelFid);
end
rois = labels{1,1};

[~,~, grp_chars] = xlsread('/path_to_your_MEG_data_dir/subj_lists.xlsx');
idIx = find(strcmp('ID',grp_chars(1,:)));
dxIx = find(strcmp('StudyDiagnosis',grp_chars(1,:)));
ageIx = find(strcmp('Age',grp_chars(1,:)));
%dxIx = find(strcmp('relatives',grp_chars(1,:))); % Goes with second
%goi_pairs definition to get the groups collasped across age

goi_pairs = [ {'Autism';'TDC'} {'Autism';'Sibling'} {'Sibling';'TDC'} ];%  {'ASD_Parent'; 'Adult_Control'} ];
%goi_pairs = [{'relative';'control'}];
cond1 = 'ep10';
if strcmp(correct,'yes')
    cond2 = 'correp20';
else
    cond2 = 'ep20'
end
%rois = {'lSTG','rSTG'};
bounds = 'yes';
for nGP = 1:size(goi_pairs,2)
    fprintf('Found %i comparisons', size(goi_pairs,2));
    clear goi y grp1members grp2members
    goi = goi_pairs(:,nGP); % "Group of interest"
    disp ('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    fprintf('Comparison = %s\n',goi{:});
    grp1Name = char(goi(1));
    grp2Name = char(goi(2));
    for k = 1:length(goi)
        y = arrayfun(@(x) ismember(x,goi{k,1}),grp_chars(:,dxIx),'UniformOutput',false);  %09/17 Can throw error, if NaN's present
        members = cellfun(@(a,b) (a*b),  grp_chars(:,idIx),y,'UniformOutput',false);
        possMembers = cellfun(@(x) (x(x ~= 0)),members,'UniformOutput',false);
        grpMembers = members(~cellfun(@isempty,possMembers));
        ages = cellfun(@(a,b) (a*b),  grp_chars(:,ageIx),y,'UniformOutput',false);
        possAges = cellfun(@(x) (x(x ~= 0)),ages,'UniformOutput',false);
        grpAges = ages(~cellfun(@isempty,possAges));

        if k == 1
            grp1members = grpMembers;
            grp1_ages = grpAges;
        else
            grp2members = grpMembers;
            grp2_ages = grpAges;
        end
    end

    msrs = {  'nepower' 'mplf' }; %  %currently not resetting something between msrs and causing the stats to not run...
    tests = { 't-test' 'diffs' 'GrandAvg'};% 'F-test' 'Correlation'};%     }; %
    clusters = {''}; % 'mcc'


    for i = 1: length(rois)
        roi = rois{i}
        for j = 1:length(msrs)
            msr = msrs{j};
            clear *_cond1 *_cond2 *diff *data nGRP*;
            if exist('bounds','var')
                zlimited = [];
            if strcmp(msr,'mplf')
                zlimited.max = .2;
                zlimited.min = -.05;
            else
                zlimited.max = 2;
                zlimited.min = -0.5;
            end
            end
            fprintf('%s had %d participants\n', grp1Name, length(grp1members));
            for g = 1:length(grp1members)
                grp1_files_cond1{g} = [cwd,filesep,num2str(grp1members{g}),filesep,num2str(grp1members{g}), cond1,'_',roi, '_', msr ,'_tfr.mat'];
                grp1_files_cond2{g} = [cwd,filesep,num2str(grp1members{g}),filesep,num2str(grp1members{g}), cond2, '_',roi, '_', msr ,'_tfr.mat'];
                %y = arrayfun(@(x) ismember(x,goi{1,1}),grp_chars(:,dxIx),'UniformOutput',false);
                %grp1_ages{g} = cellfun(@(a,b) (a*b),  grp_chars(:,ageIx),y,'UniformOutput',false);
            end

            %% Controls files
            fprintf('%s had %d participants\n', grp2Name, length(grp2members));
            for g2 = 1:length(grp2members)
                grp2_files_cond1{g2} = [cwd,filesep,num2str(grp2members{g2}),filesep,num2str(grp2members{g2}), cond1, '_',roi, '_', msr ,'_tfr.mat'];
                grp2_files_cond2{g2} = [cwd,filesep,num2str(grp2members{g2}),filesep,num2str(grp2members{g2}), cond2, '_',roi, '_', msr ,'_tfr.mat'];
            end

            %% number of files
            nGRP1_cond1 = length(grp1_files_cond1);
            nGRP1_cond2 = length(grp1_files_cond2);
            nGRP1_diff = length(grp1_files_cond2);
            nGRP2_cond1  = length(grp2_files_cond1);
            nGRP2_cond2  = length(grp2_files_cond2);
            nGRP2_diff  = length(grp2_files_cond2);



            %% Load data
            for ii=1:nGRP1_cond1
                load(grp1_files_cond1{ii});
                grp1_cond1{ii}=tfr;
                grp1_cond1data{ii}=tfr.powspctrm;
            end

            for ii=1:nGRP1_cond2
                load(grp1_files_cond2{ii});
                grp1_cond2{ii}=tfr;
                grp1_cond2data{ii}=tfr.powspctrm;
                grp1_diffdata{ii}=grp1_cond2data{ii} - grp1_cond1data{ii};
                grp1_diff{ii}=grp1_cond2{ii};
                grp1_diff{ii}.powspctrm = grp1_diffdata{ii};
                if mean2(~isnan(grp1_cond2{ii}.powspctrm)) == 0
                    fprintf('Not enough values: %s\n', grp1_files_cond2{ii} )
                end

                if mean2(~isnan(grp1_cond1{ii}.powspctrm)) == 0
                    fprintf('Not enough values: %s\n', grp1_files_cond1{ii} )
                end
            end

            for ii=1:nGRP2_cond1
                load(grp2_files_cond1{ii});
                grp2_cond1{ii}=tfr;
                grp2_cond1data{ii}=tfr.powspctrm;
            end

            for ii=1:nGRP2_cond2
                load(grp2_files_cond2{ii});
                grp2_cond2{ii}=tfr;
                grp2_cond2data{ii}=tfr.powspctrm;
                grp2_diffdata{ii}=grp2_cond2data{ii} - grp2_cond1data{ii};
                grp2_diff{ii}=grp2_cond2{ii};
                grp2_diff{ii}.powspctrm = grp2_diffdata{ii};
                if mean2(~isnan(grp2_cond2{ii}.powspctrm)) == 0
                    fprintf('Not enough values: %s\n', grp2_files_cond2{ii} )
                end
                if mean2(~isnan(grp2_cond1{ii}.powspctrm)) == 0
                    fprintf('Not enough values: %s\n', grp2_files_cond1{ii} )
                end
            end


            grp1a_cond1 = ft_freqgrandaverage([],grp1_cond1{:});
            grp1a_cond2 = ft_freqgrandaverage([],grp1_cond2{:});
            grp2a_cond1 = ft_freqgrandaverage([],grp2_cond1{:});
            grp2a_cond2 = ft_freqgrandaverage([],grp2_cond2{:});

            %differences
            grp1a_diff = ft_freqgrandaverage([],grp1_diff{:});
            grp2a_diff = ft_freqgrandaverage([],grp2_diff{:});

            single_groups = {'grp1';'grp2'};

            %% side-by-side
            resultsDir = ['path_to_your_MEG_tools_dir/asd_meg_kl/results/'];
            if ~exist (resultsDir,'dir')
                mkdir (resultsDir)
                mkdir ([resultsDir 'fTest'])
                %mkdir ([cwd '/results/2nd_level/statFiles'])
                %mkdir ([cwd '/results/2nd_level/fTest'])
            end

            for cl = 1:length(clusters)
                cluster = clusters{cl};
                for r = 1: length(tests)
                    test = tests{r};

                    grp1 = char(single_groups(1));
                    grp2 = char(single_groups(2));
                    grp1Name = char(goi (1))
                    grp2Name = char(goi (2))

                    switch test
                        %% Difference plots
                        case 'diffs'
                            disp ('')
                            disp (' > Portraying group differences <')
                            pairs = [single_groups];
                            v=1;
                            while v < length(goi)
                                if strcmp(cluster,'mcc')
                                    fid = [resultsDir 'mcc_' grp1Name 'V' grp2Name '_respDiff_' roi '_' msr covaryLabel '.eps'];
                                    fid2 = [resultsDir 'statFiles/' 'mcc_'  grp1Name 'V' grp2Name '_respDiff_' roi '_' msr covaryLabel '_stat.mat'];
                                else
                                    fid = [resultsDir grp1Name 'V' grp2Name '_respDiff_' roi '_' cond2 '_' msr covaryLabel '.eps'];
                                    fid2 = [resultsDir 'statFiles/' grp1Name 'V' grp2Name '_respDiff_' roi '_' cond2 '_' msr covaryLabel '_stat.mat'];
                                end
                                check_done = rdir(fid);
                                if isempty(check_done)

                                    [cfg_stat,cfg_plot,h] = load_cfg(cluster);
                                    clear cfg_stat.parameter; %Just for this series
                                    grp1 = char(pairs (v));
                                    grp2 = char(pairs (v+1));
                                    grp1Name = char(goi (v));
                                    grp2Name = char(goi (v+1));
                                    diff1    = eval(strcat(grp1, '_diff'));
                                    diff2    = eval(strcat(grp2, '_diff'));
                                    diff1Avg = eval(strcat(grp1, 'a_diff'));
                                    diff2Avg = eval(strcat(grp2, 'a_diff'));
                                    n1       = eval(strcat('n', upper(grp1), '_diff'));
                                    n2       = eval(strcat('n', upper(grp2), '_diff'));

                                    cfg_stat.design = [%1:n1+n2;
                                        ones(1,n1) ones(1,n2)*2];
                                    if strcmpi(covary,'yes')
                                        cfg_stat.design(2,:) = cell2mat ([grp1_ages; grp2_ages])';
                                        cfg_stat.cvar = 2;
                                    end
                                    cfg_stat.numrandomization = 500;

                                    stat = ft_freqstatistics(cfg_stat, diff1{:}, diff2{:});
                                    ind = find(stat.prob<.05);
                                    mask = stat.prob;
                                    mask(:)=0;
                                    mask(ind)=1;
                                    zmax = max(diff1Avg.powspctrm(:));
                                    zmin = min(diff1Avg.powspctrm(:));
                                    if ~exist('zlimited','var');
                                        cfg_plot.zlim = [zmin zmax];
                                    else
                                        cfg_plot.zlim = [zlimited.min zlimited.max];
                                    end

                                    if zmax > 0
                                        subplot(3,1,1); ft_singleplotTFR(cfg_plot,diff1Avg);
                                        title([(grp1Name) ' Diff (' cond1 '-' cond2 ') Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                        subplot(3,1,2); ft_singleplotTFR(cfg_plot,diff2Avg);
                                        title([(grp2Name) ' Diff (' cond1 '-' cond2 ') Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');

                                        %cfg_stat               = [];
                                        cfg_stat.parameter     = 'stat';
                                        subplot(3,1,3); ft_singleplotTFR(cfg_stat,stat);
                                        title(['T-statistic map: ' msr ' ' cluster]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                        hold on;
                                        %% MCC
                                        if strcmp(cluster,'mcc')
                                            if exist('stat.posclusters')
                                                pos_pvals=[stat.posclusters(:).prob];
                                                pind = find(pos_pvals<stat.cfg.alpha);
                                                pos = squeeze(ismember(stat.posclusterslabelmat, pind));
                                                %over plot a line with the cluster corrected p < .05 result
                                                contour(stat.time,stat.freq,squeeze(pos),[1 1],'w-','linewidth',3);
                                            end
                                            if exist('stat.negclusters')
                                                neg_pvals=[stat.negclusters(:).prob];
                                                nind = find(neg_pvals<stat.cfg.alpha);
                                                neg = squeeze(ismember(stat.negclusterslabelmat, nind));
                                                %over plot a line with the cluster corrected p < .05 result
                                                contour(stat.time,stat.freq,squeeze(neg),[1 1],'k-','linewidth',3);
                                            end
                                        else

                                            contour(stat.time,stat.freq,squeeze(mask),[1 1],'w-','linewidth',2);
                                        end

                                        %% Save output
                                        print(h,fid,'-depsc');
                                        save(fid2,'stat');
                                    end
                                    close (h);
                                else
                                    fprintf('Already plotted %s. Skipping \n', test)
                                end
                                v = v+1;
                            end

                            %% Grand Average plot
                        case 'GrandAvg'
                            if ~strcmp(cluster,'mcc')
                                fid = [cwd '/results/grandAvg_' roi '_' cond2 '_' msr covaryLabel '.eps'];
                                check_done = rdir(fid);
                                if isempty(check_done)
                                    [cfg_stat,cfg_plot,h] = load_cfg(cluster);
                                    zmax = max(grp2a_cond1.powspctrm(:));
                                    zmin = min(grp2a_cond1.powspctrm(:));

                                    if ~exist('zlimited','var');
                                        cfg_plot.zlim = [zmin zmax];
                                    else
                                        cfg_plot.zlim = [zlimited.min zlimited.max];
                                    end

                                    subplot(2,2,1); ft_singleplotTFR(cfg_plot,grp1a_cond1);
                                    title([grp1Name cond1 ' Grand Average: ' msr]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    subplot(2,2,3); ft_singleplotTFR(cfg_plot,grp2a_cond1);
                                    title([grp2Name cond1 ' Grand Average: ' msr]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    subplot(2,2,2); ft_singleplotTFR(cfg_plot,grp1a_cond2);
                                    title([grp1Name cond2 ' Grand Average: ' msr]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    subplot(2,2,4); ft_singleplotTFR(cfg_plot,grp2a_cond2);
                                    title([grp2Name cond2 ' Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    %subplot(3,2,5); ft_singleplotTFR(cfg_plot,grp2a_cond2);
                                    %title(['HC cond2 Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    %subplot(3,2,6); ft_singleplotTFR(cfg_plot,grp2a_cond2);
                                    %title(['HC cond2 Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                    print(h,fid,'-depsc');
                                    close (h);
                                else
                                    fprintf ('already plotted %s. skipping \n', test)
                                end
                            end
                        otherwise
                            trials = {cond1 cond2 'diff'};
                            for t = 1: length(trials);
                                trial = trials{t};
                                grpComp = [grp1Name 'V' grp2Name];
                                fDir = [resultsDir 'tTest/' grpComp];
                                sDir = [resultsDir 'statFiles/' grpComp];
                                if ~exist (fDir)
                                    mkdir (fDir)
                                    mkdir (sDir)
                                end

                                if strcmp(cluster,'mcc')
                                    fid = [fDir '/mcc_' test '_' trial '_' roi '_' cond2 '_' msr covaryLabel '.eps'];
                                    fid2 = [sDir '/mcc_' test '_' trial '_' roi '_' cond2 '_' msr covaryLabel '_stat.mat'];

                                else
                                    fid = [fDir '/' test '_' trial  '_' roi '_' cond2 '_' msr covaryLabel '.eps'];
                                    fid2 = [sDir '/' test '_' trial  '_' roi '_' cond2 '_' msr covaryLabel '_stat.mat'];
                                end
                                check_done = rdir(fid);

                                if isempty(check_done)


                                    [cfg_stat,cfg_plot,h] = load_cfg(cluster);

                                    if strcmp('diff', trial)
                                        cndtn = strcat('_',trial);
                                    else
                                        cndtn = ['_cond',int2str(t)];
                                    end
                                    n1 = eval(strcat('n', upper(grp1),cndtn));
                                    n2 = eval(strcat('n',upper(grp2),cndtn));
                                    files1 = eval(strcat(grp1,cndtn));
                                    files2 = eval(strcat(grp2,cndtn));
                                    inputAvg1 = eval(strcat(grp1,'a', cndtn));
                                    inputAvg2 = eval(strcat(grp2,'a', cndtn));
                                    if any(strcmp(test ,'t-test') || strcmp(test, 'F-test'))
                                        cfg_stat.clustercritval = .05; %uncomment for F-test
                                        cfg_stat.design = [%1:n1+n2+n3;
                                            ones(1,n1) ...
                                            ones(1,n2)*2];
                                        cfg_stat.cvar = '';
                                        if strcmpi(covary,'yes')
                                            cfg_stat.design(2,:) = cell2mat ([ grp1_ages; grp2_ages])';
                                            cfg_stat.cvar = 2;
                                        end

                                        cfg_stat.tail = 1; %BMM addition to make the stat script happy
                                        if strcmp(test,'F-test')
                                            cfg_stat.statistic = 'indepsamplesF'; %Swap the default
                                        end
                                        cfg_stat.correctm   = 'fdr'; % Monte Carlo not appropriate for F-test, per Don (3/17)... I think.
                                    elseif test == 'Correlation'
                                        disp('Running correlation')
                                        cfg_stat.design = cell2mat ([ grp1_ages; grp2_ages] )'; %uncenters so the variance isn't the sqrt of a negative
                                        cfg_stat.statistic = 'indepsamplesregrT';
                                        cfg_stat.cvar = '';
                                    end

                                    fprintf('\n\nRunning %s stats || %s || %s || %s || %s \n',test, trial, roi, msr, grpComp);
                                    stat = ft_freqstatistics(cfg_stat,files1{:},files2{:});
                                    ind = find(stat.prob<.05);
                                    mask = stat.prob;
                                    mask(:)=0;
                                    mask(ind)=1;


                                    save(fid2,'stat'); %save before cfg_stat is overwritten

                                    result_check = find(~isnan(stat.stat));
                                    if isreal(stat.stat) == 1 & result_check > 0 ;
                                        zmax = max(inputAvg1.powspctrm(:));
                                        zmin = min(inputAvg1.powspctrm(:));
                                        if ~exist('zlimited','var');
                                            cfg_plot.zlim = [zmin zmax];
                                        else
                                            cfg_plot.zlim = [zlimited.min zlimited.max];
                                        end

                                        set(h,'Position', [300 300 600 600]); %Example of sizing!!
                                        subplot(4,1,1); ft_singleplotTFR(cfg_plot,inputAvg1);
                                        title([grp1Name ' ' trial ': ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                        subplot(4,1,2); ft_singleplotTFR(cfg_plot,inputAvg2);
                                        title([grp2Name ' ' trial ': ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                        %cfg_stat               = [];
                                        cfg_stat.parameter     = 'stat';
                                        subplot(4,1,3); ft_singleplotTFR(cfg_stat,stat);
                                        title([test ' map: ' msr cluster]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                        hold on;
                                        if strcmp(cluster,'mcc');
                                            if exist('stat.posclusters')
                                                pos_pvals=[stat.posclusters(:).prob];
                                                pind = find(pos_pvals<stat.cfg.alpha);
                                                pos = squeeze(ismember(stat.posclusterslabelmat, pind));
                                                %over plot a line with the cluster corrected p < .05 result
                                                contour(stat.time,stat.freq,squeeze(pos),[1 1],'w-','linewidth',3);
                                            end
                                            if exist('stat.negclusters')
                                                neg_pvals=[stat.negclusters(:).prob];
                                                nind = find(neg_pvals<stat.cfg.alpha);
                                                neg = squeeze(ismember(stat.negclusterslabelmat, nind));
                                                %over plot a line with the cluster corrected p < .05 result
                                                contour(stat.time,stat.freq,squeeze(neg),[1 1],'k-','linewidth',3);
                                            end

                                        else
                                            contour(stat.time,stat.freq,squeeze(mask),[1 1],'w-','linewidth',2);

                                        end
                                        print(h,fid,'-depsc');

                                        fprintf('End of %s for: %s \nCorrection was: %s \n', test, trial, cluster);
                                    else
                                        fprintf('%s analysis did not yield significant results', test);
                                    end
                                    close (h);
                                else
                                    fprintf('Already plotted %s. Skipping\n', test);
                                end
                            end

                    end

                end




                %Cluster-correction
                %% Single group comparisons

                for k = 1: length(single_groups)
                    grp = single_groups{k};
                    grpName = goi{k};
                    grpUpper = upper(grp);

                    for cl = 1:length(clusters)
                        cluster = clusters{cl};
                        if strcmp(cluster,'mcc')
                            fid = [resultsDir grpName '_' roi '_' cond2 '_' msr covaryLabel '_corr.eps'];
                            fid2 = [resultsDir 'statFiles/' grpName '_' roi '_' cond2 '_' msr covaryLabel '_corrStat.mat'];
                        else
                            fid = [resultsDir grpName '_' roi '_' cond2 '_' msr covaryLabel '.eps'];
                            fid2 = [resultsDir 'statFiles/' grpName '_' roi '_' cond2 '_' msr covaryLabel '_stat.mat'];
                        end
                        check_done = rdir(fid);

                        if isempty(check_done)
                            %Defaults
                            [cfg_stat,cfg_plot,h] = load_cfg(cluster);
                            if strcmp (grp, 'grp1')
                                cfg_stat.design = [%1:npar_cond1+npar_cond2;
                                    ones(1,nGRP1_cond1) ...
                                    ones(1,nGRP1_cond2)*2];
                            elseif strcmp (grp, 'grp2')
                                cfg_stat.design = [%1:nHC_cond1+nHC_cond2;
                                    ones(1,nGRP2_cond1) ...
                                    ones(1,nGRP2_cond2)*2];
                            end

                            cfg_stat.design = [cfg_stat.design];
                            if strcmpi(covary,'yes')
                                age = eval(strcat(grp, '_ages'));
                                cfg_stat.design(2,:) = cell2mat ([age; age]);
                                cfg_stat.cvar = 2;
                            end

                            input10 = eval(strcat(grp, '_cond1'));
                            input20 = eval(strcat(grp, '_cond2'));
                            fprintf('Calculating stats: %u %u \n',length(input10), length(input20));
                            stat = ft_freqstatistics(cfg_stat, input10{:}, input20{:});
                            % use to look at uncorrected p-values (also uncomment cfg_stat.correctm =
                            % 'cluster' and comment out the sections that mask the clusters
                            ind = find(stat.prob<.05);
                            mask = stat.prob;
                            mask(:)=0;
                            mask(ind)=1;

                            inputAvg10 = eval(strcat(grp, 'a_cond1'));
                            inputAvg20 = eval(strcat(grp, 'a_cond2'));
                            pws_var = inputAvg10.powspctrm;
                            zmax                   = max(pws_var(:));
                            zmin                   = min(pws_var(:));
                            cfg_plot               = [];
                            cfg_plot.parameter     = 'powspctrm';
                            if isempty(zlimited)
                                cfg_plot.zlim          = [zmin zmax]
                            else
                                cfg_plot.zlim = [zlimited.min zlimited.max]
                            end

                            if zmax > 0
                                %Top part of graphs
                                subplot(3,1,1); ft_singleplotTFR(cfg_plot,inputAvg10);
                                title([grpName ' ' cond1 ' Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                subplot(3,1,2); ft_singleplotTFR(cfg_plot,inputAvg20);
                                title([grpName ' ' cond2 ' Grand Average: ' msr ]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                %cfg_stat               = [];
                                cfg_stat.parameter     = 'stat';

                                %Third section of graphs (the comparison
                                subplot(3,1,3); ft_singleplotTFR(cfg_stat,stat);
                                title(['T-statistic map: ' msr cluster]); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
                                hold on;
                                if strcmp(cluster,'mcc')
                                    if exist('stat.posclusters')
                                        pos_pvals=[stat.posclusters(:).prob];
                                        pind = find(pos_pvals<stat.cfg.alpha);
                                        pos = squeeze(ismember(stat.posclusterslabelmat, pind));
                                        %over plot a line with the cluster corrected p < .05 result
                                        contour(stat.time,stat.freq,squeeze(pos),[1 1],'w-','linewidth',3);
                                    end
                                    if exist('stat.negclusters')
                                        neg_pvals=[stat.negclusters(:).prob];
                                        nind = find(neg_pvals<stat.cfg.alpha);
                                        neg = squeeze(ismember(stat.negclusterslabelmat, nind));
                                        %over plot a line with the cluster corrected p < .05 result
                                        contour(stat.time,stat.freq,squeeze(neg),[1 1],'k-','linewidth',3);
                                    end
                                else
                                    contour(stat.time,stat.freq,squeeze(mask),[1 1],'w-','linewidth',2);
                                end
                                print(h,fid,'-depsc');
                                save(fid2,'stat');
                            end
                            close (h);
                        else
                            disp('Already plotted. Skipping.')

                        end
                    end
                end

            end

        end
    end
end
