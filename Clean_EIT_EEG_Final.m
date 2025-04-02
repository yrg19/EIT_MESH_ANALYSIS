%% Clean EEG and EIT spike data


clear 
clc
close all 

%% Setup packagaes 

setup = 1; 

if setup == 1
    addpath(genpath("/home/eit/Packages_and_tools/wave_clus-master")); 
    addpath("/home/eit/Packages_and_tools/Random_functions"); 
    run('/home/eit/Packages_and_tools/eidors-v3.10-ng/eidors/startup.m')
    cd('/home/yuval/MATLAB Add-Ons/Collections/EEGLAB')
    eeglab
    cd('/home/eit/Packages_and_tools/ScouseTom Software/Fast-neural-processing-master/biosign')
    biosig_installer
    savepath
end 

%% Load Data  

clear 
clc

s_num = 5; 
rec_name = 'EIT003'; %'eit'; 

cd('/home/eit/Documents/Spike_detection')

if ~exist(['E0',num2str(s_num),'/',rec_name], 'dir')

    mkdir(['E0',num2str(s_num),'/',rec_name])

end 

cd('/home/eit/Packages_and_tools/ScouseTom Software')

inpath = ['/home/yuval/E0',num2str(s_num),'/raw_data/']; 
addpath(inpath)
eeg_rec = pop_loadbv(inpath, [rec_name,'.vhdr']); 


outpath = ['/home/eit/Documents/Spike_detection/E0',num2str(s_num),'/', rec_name]; 

load(['/home/yuval/E0',num2str(s_num),'/raw_data/', rec_name, '_log.mat'])

stim = mk_stim_yuval(ExpSetup.Protocol,64,60e-6,1,32);

coord = readmatrix(['/home/yuval/E0',num2str(s_num),'/raw_data/electrodes_mm.txt']);

names = table2cell(readtable("electrode_names.xlsx")); 

%% Find Nearest Electrodes to ReReference 

k = 3; % Include the point itself and its nearest neighbor
[nearest, dist] = knnsearch(coord, coord, 'K', k);


% Number of clusters
numClusters = 8;

% Perform K-means clustering
[idx, centroids] = kmeans(coord, numClusters);

% `idx` contains the cluster assignment for each point.
% `centroids` contains the 3D coordinates of the cluster centers.

legendEntries = {};

% Example of how to group points

for n = 1:numClusters

    group{n} = [coord(idx == n, :), find(idx == n)];

    scatter3(coord(idx == n, 1), coord(idx == n, 2), coord(idx == n, 3), 50, 'filled')
    hold on 
    legendEntries{end+1} = ['Cluster ' num2str(n)]; % Add label for each line


end 

    legend(legendEntries)


%% Filter EEG for chosen channel 

cd(outpath)

oeeg = eeg_rec.data'; 

% Rereference 

reeg = oeeg - mean(oeeg, 2);
    
sr = eeg_rec.srate; 
reeg = reeg(1:length(reeg),:); % remove last second of data because some recordings have issues. 
chanOI = 22; %channel 10 as leading channel 

% downsample
df = 10;
Fs = sr/ df;
ceeg = reeg(1:df:length(reeg),:);

% filter EEG spikes

[b, a] = iirnotch(50 / (Fs/2) , 50 / (Fs/2) / 35); 
all_filt = filtfilt(b, a, double(ceeg)); 
filt_eeg = filtfilt(b, a, double(ceeg(:,chanOI)));

[b, a] = iirnotch(1700 / (Fs/2) , 1700 / (Fs/2) / 35); 
all_filt = filtfilt(b, a, all_filt);
filt_eeg = filtfilt(b, a, filt_eeg);

[b, a] = butter(1,[1 200]/ (Fs/2),'bandpass');
all_filt = filtfilt(b, a, all_filt); 
filt_eeg = filtfilt(b, a, filt_eeg);

plot_filt = all_filt(1: 10:length(all_filt),:); 

[~, sortedIdx] = sort(idx);
plot_filt = plot_filt(:,sortedIdx);

save('all_data_filtered','plot_filt');

% Time

Time= linspace(1,(length(filt_eeg)/Fs),length(filt_eeg));

plot(Time,filt_eeg)
hold on

x = all_filt(:,chanOI); 

plot(x)
legend({'Chan10 Filt', 'Ch10 from All Filt'})

%% Save to plot with SEEG 

cd(outpath)

% Create a basic EEGLAB structure
EEG = pop_importdata('data', plot_filt', 'setname', 'Clustered EEG', 'srate', 500);

% Reorder channels by clusters
[~, sortedIdx] = sort(idx); % Sort indices by cluster
EEG.data = EEG.data(sortedIdx, :); % Reorder the data

% Initialize chanlocs structure
numChannels = size(coord, 1);
chanlocs = struct('labels', [], 'X', [], 'Y', [], 'Z', []);

% Populate chanlocs
for i = 1:numChannels
    %chanlocs(i).labels = ['Ch' num2str(i) '_C' num2str(idx(sortedIdx(i)))]; % Assign channel label
    chanlocs(i).labels = names{sortedIdx(i)}; 
    chanlocs(i).X = coord(sortedIdx(i), 1); % X-coordinate
    chanlocs(i).Y = coord(sortedIdx(i), 2); % Y-coordinate
    chanlocs(i).Z = coord(sortedIdx(i), 3); % Z-coordinate
end

% Assign to EEG structure
EEG.chanlocs = chanlocs;


% Save the updated dataset for EEGLAB
EEG = eeg_checkset(EEG);
pop_saveset(EEG, 'filename', 'ClusteredEEG.set', 'filepath', outpath);

%% Define segments 

str_stp = []; stim_str = []; 

event_eeg = struct2table(eeg_rec.event); 

stims = contains(event_eeg.type, 'S  3');
stim_str = event_eeg.latency(stims); 
stim_str = stim_str/df; %downsample the triggers 
remove_segs = diff(stim_str); 
if any(remove_segs) < 3*Fs
    warning('Incorrectly labeled segment found')
    del = find(remove_segs < 3*Fs); 
    stim_str(del) = []; 
end 
seglength = median(diff(stim_str)); 

%sanity check the timings to make sure they overlap 
figure()
plot(ceeg(:,chanOI))
hold on 
xline(stim_str)


%create a 2 minute segment from before stimulation times 
str_stp(:,1) =stim_str; %starting times - 120000 because of artifact at recording start

% narrowband filter and downsample

[b, a] = butter(1,[1 40]/ (Fs/2),'bandpass');
all_filt = filtfilt(b, a, all_filt); 
filt_eeg = filtfilt(b, a, filt_eeg);

df2 = 20; 
filt_eeg = filt_eeg(1:df2:length(filt_eeg),:); 
all_filt = all_filt(1:df2:length(all_filt),:); 

Fs2 = Fs/ df2; 

str_stp = round(str_stp ./ df2); 
str_stp(:,2) = str_stp + 120*Fs2; 

% Double check still looks the same! 

figure()
plot(filt_eeg)
hold on 
plot(all_filt(:,chanOI))


% Sanity check all timings. 
figure()
plot(ceeg(1:df2:length(ceeg),chanOI))
hold on 
xline(str_stp(:,1), 'r', 'LineWidth',2)
xline(str_stp(:,2), 'k')


%% Amp detect spikes 

cd(outpath)

par = set_parameters_DEFAULT; 
par.segments_length = round(length(filt_eeg)/(Fs2*60)); % add 0.1 because its slightly bigger than 2 minutes and fixes a future bug 
par.sr = Fs2; 
par.stdmin = 5; 
par.detect_order = 0; 
par.detect_fmax = 100; 
par.detect_fmin = 1; 
par.sort_order = 0; 
par.sort_fmax = 100; 
par.sort_fmin = 1; 
par.ref_ms = 350; 
par.w_pre = round(0.15*Fs2); %100 ms before 
par.w_post = round(0.3*Fs2); % 150ms after 
par.detection = 'pos';
par.min_clus = 5;
par.interpolation = 'n'; 
par.cont_segment = 0; 


%% Get spikes 

for c = 1:numChannels
    
    chanOI = c; 
    
    data = squeeze(all_filt(:,chanOI));
    
    data(abs(data)> 10*std(data)) = 0; 
    
    % plot(data)
    
    save('whole_recording_leading_chan_EEG','data')
    
    inp = {'whole_recording_leading_chan_EEG.mat'}; 
    
    Get_spikes(inp,par);

    load('whole_recording_leading_chan_EEG_spikes.mat', 'spikes')

    % disp(['Channel ', names{c}, ' contains ', num2str(size(spikes,1)), ' spikes'])

    subplot(8,8,c)

    plot(spikes')
    ylim([-500 700])
    
    title([names{c}, ' - ', num2str(size(spikes,1))])

end 
%% Manual spikes removal

Time = linspace(-1*par.w_pre/Fs2*1000 ,par.w_post/Fs2*1000 , par.w_post + par.w_pre);

load('whole_recording_leading_chan_EEG_spikes.mat');

spikes_rod_chans = []; start_point= []; end_point = [];

chans_OI = nearest(chanOI, :); 

% Loop through each index and create a linspace for each start and end pair
for i = 1:length(idx)
    start_point = idx(i) - par.w_pre;
    end_point = idx(i) + par.w_post-1;
    spikes_rod_chans(i,:,:) = squeeze(all_filt(start_point:end_point,chans_OI)); %, 8, 7, 36, 35, 34, 33]));
end

removed = []; spk_keep = [];

for sp = 1:size(spikes,1)

    plot(Time, squeeze(spikes_rod_chans(sp,:,:)))
    hold on
    xline(0)
    hold off
    ylim([-400 800])
    title(['Spike Number ', num2str(sp)])

    answer = questdlg('Would you like keep this spike?', ...
        'Yes','No');


    % Handle response
    switch answer
        case 'Yes'
            spk_keep(sp) = 1;
        case 'No'
            spk_keep(sp) = 0;
            removed = [removed; sp];
        case 'Cancel'
            error('Manual spike removal terminated by user')
    end

end

spikes(~logical(spk_keep),:) = [];
idx(~logical(spk_keep)) = [];
index(~logical(spk_keep)) = [];

disp(['kept ', num2str(sp - length(removed)), ' out of ', num2str(sp), ' spikes'])

save(['whole_recording_leading_chan_EEG_spikes_clean.mat'], 'spikes', 'par','threshold','index','idx','sr_psegment','psegment', 'spk_keep');

%% Do Clustering 


inp = {['whole_recording_leading_chan_EEG_spikes_clean.mat']};

Do_clustering(inp, par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except eeg_rec rec_name outpath stim


%% Filter and demodulate EIT

oeeg = eeg_rec.data';
sr = eeg_rec.srate;

[b, a] = butter(1, [1400 2000]/(sr/2),'bandpass');

num_chans = size(oeeg,2); 

for c = 10%:num_chans 

    tic 

    sig = oeeg(:,c); 

    sig_filt = filtfilt(b, a, double(sig));

    sig_hilb = abs(hilbert(sig_filt)); 

    all_hilb(:,c) = single(sig_hilb); 

    disp(['Finished Channel ', num2str(c)])

    toc 
end 

%% Save demodulated signal 

cd(outpath) 
save(['demod_EIT_', rec_name], 'all_hilb', '-v7.3')


%% Downsample and filter again signal 

chanOI = 10; 

df = 50; 
Fs = sr / df; 

ds_hilb = double(all_hilb(1:df:size(all_hilb,1), :)); 

% filter again 
[b, a] = butter(1,[20 80]/(Fs/2),'bandpass');
filt_hilb = filtfilt(b,a,ds_hilb); 

subplot(3,1,1)
plot(all_hilb(:,chanOI))
title('Demodulated EIT')

subplot(3,1,2) 
plot(ds_hilb(:,chanOI))
title('DownSampled Demodulated EIT')

subplot(3,1,3)
plot(filt_hilb(:,chanOI))
title('Filtered, Downsampled Demodulated EIT')

clear oeeg sig_hilb sig_filt  


%% Split into segmenets 

event_eeg = struct2table(eeg_rec.event);

stims = contains(event_eeg.type, 'S  3');
stim_str = event_eeg.latency(stims);
remove_segs = diff(stim_str);
if any(remove_segs) < 3*sr
    warning('Incorrectly labeled segment found')
    del = find(remove_segs < 3*sr);
    stim_str(del) = [];
end 
seglength = median(diff(stim_str));


str_stp(:,1) = stim_str; %[sr; stim_str]; %starting times - 120000 because of artifact at recording start
str_stp = round(str_stp ./ (sr/Fs)); 
str_stp(:,2) = str_stp(:,1) + 120* Fs; 

%str_stp = str_stp(1:length(stim),:);


%sanity check start and stop times 
figure()
plot(filt_hilb(:,10),'LineWidth',2)
hold on 
xline(str_stp(:,1),'LineWidth', ...
    2,'Color','r')
legend({'start time'})
xline(str_stp(:,2),'LineWidth',0.3, 'Color','k')
legend({'stop time'})
xlabel('Samples')
ylabel('dZ')
title('Segement start and stop')



%% Find spikes per segment 

load('times_whole_recording_leading_chan_EEG_spikes_clean.mat')
load('whole_recording_leading_chan_EEG_spikes_clean.mat')

us = Fs / 250; 

timepoints = idx .* us; 

for i = 1:size(str_stp, 1)
    
    % Get start and stop time for the current epoch
    start_time = str_stp(i, 1);
    stop_time = str_stp(i, 2);
    
    % Find time points within the current epoch
    points_in_epochs{i} = timepoints(timepoints >= start_time & timepoints <= stop_time);

    clust_in_epochs{i} = cluster_class(find(timepoints >= start_time & timepoints <= stop_time),1);

end


%% Cut EIT timeseries into injection segments 

num_segs = length(str_stp)-1; 
num_chans = size(filt_hilb,2); 

EIT_all_seg_raw = nan(num_segs, 120*Fs, num_chans); 

tmpts = arrayfun(@colon, str_stp(1:41,1), str_stp(1:41,2), 'Uniform', false);
tmpts = cell2mat(tmpts);

hh = filt_hilb(tmpts,:);

hh = reshape(hh,num_segs,120*Fs+1,num_chans);

hh = shiftdim(hh,1); 

EIT_all_seg_raw = shiftdim(hh,2);

%% Cut unfiltered data into segements for standing voltage calc 

num_segs = length(str_stp)-1; 
num_chans = size(ds_hilb,2); 

unfiltered_segments = nan(num_segs, 120*Fs, num_chans); 

tmpts = arrayfun(@colon, str_stp(1:41,1), str_stp(1:41,2), 'Uniform', false);
tmpts = cell2mat(tmpts);

hh = ds_hilb(tmpts,:); 

hh = reshape(hh,num_segs,120*Fs+1,num_chans);

hh = shiftdim(hh,1); 

unfiltered_segments = shiftdim(hh,2);


%% Standing Voltage Values 

EIT_all_seg = nan(num_segs, 120*Fs+1, num_chans); 

for s= 1:num_segs 
    
    tmp = squeeze(EIT_all_seg_raw(s,:,:)); 
    
    un_seg = squeeze(unfiltered_segments(s,:,:)); 
    
    SV(s,:)= mean(un_seg(2*Fs:end-(2*Fs),:), 1); % removing injection artefact from calulation 
    
    tmp2 = (tmp ./ SV(s,:)) * 100; % ((tmp - SV(s,:)) ./ SV(s,:)) * 100; %
    
    EIT_all_seg_prct(s,:,:)= tmp2; 
end 

rm_sv =  prctile(SV(:), 6);
SV2 = SV; 
SV2(SV(:) < rm_sv) = NaN; 
EIT_all_seg_prct2 = EIT_all_seg_prct; 
 
[r, c] = ind2sub([size(SV,1),size(SV,2)], find(SV(:) < rm_sv));

EIT_all_seg_prct(r, :, c) = NaN; 

%% Plot segment with spike indices 

s = 1; 
ch = 10; 

plot(EIT_all_seg_prct(s,:,ch))

hold on 

points = points_in_epochs{s} - (str_stp(1,1) + (size(EIT_all_seg_prct,2) * (s-1))); 

plot(points, 2, 'o')

%% Split EIT segments into spikes 

all_EIT_spikes_prct= []; av_EIT_spks_prct =[];
all_EIT_spikes_raw= []; av_EIT_spks_raw =[];


for s = 1:length(points_in_epochs)/2

    tmps = []; points = []; r = []; tmpmat = []; 

    points = [points_in_epochs{s}, points_in_epochs{s+20}];

    points = points - (str_stp(1,1) + (size(EIT_all_seg,2) * (s-1))); 

    r = find(points - (2*Fs) < 50); 
    r = [r, find((119*Fs) - points < 0)]; 

    disp(length(r))

    points(r) =[]; 

    removed_spikes{s} = r; 
    
    tmps(:,1) = points - 0.2*Fs;
    tmps(:,2) = points + 0.5*Fs;

    %PERCENTAGE 
    tmpmat = squeeze(EIT_all_seg_prct(s,:,:)); 

    plot(tmpmat(:,10)); hold on; plot(points, 2,'o')

    for t = 1:length(tmps)

        b = tmpmat(tmps(t,1): tmps(t,2), :);

        all_EIT_spikes_prct{s}(t, :,:) = b; 

    end

    j = squeeze(mean(all_EIT_spikes_prct{s},1, 'omitnan')); 

    if s > 40
        s2 = s - 40;
    elseif s > 20
        s2 = s - 20;
    elseif s < 20
        s2 = s;
    end

    inj_chans = find(full(stim(s2).stim_pattern) ~=0); 

    inj_chans = [inj_chans; 32]; 

    j(:,inj_chans) = []; 

    av_EIT_spks_prct = [av_EIT_spks_prct, j];

    % RAW 
    tmpmat = squeeze(EIT_all_seg_raw(s,:,:)); 

    plot(tmpmat(:,10)); hold on; plot(points, 2,'o')

    for t = 1:length(tmps)

        b = tmpmat(tmps(t,1): tmps(t,2), :);

        all_EIT_spikes_raw{s}(t, :,:) = b; 

    end

    j = squeeze(mean(all_EIT_spikes_raw{s},1, 'omitnan')); 

    if s > 40
        s2 = s - 40;
    elseif s > 20
        s2 = s - 20;
    elseif s < 20
        s2 = s;
    end

    inj_chans = find(full(stim(s2).stim_pattern) ~=0); 

    inj_chans = [inj_chans; 32]; 

    j(:,inj_chans) = []; 

    av_EIT_spks_raw = [av_EIT_spks_raw, j];
   
end 

%% Divide spikes into clusters

clust_EIT_spikes_raw= []; clust_EIT_spikes_prct = []; 

for c = 1:length(unique(cluster_class(:,1)))

            
    % PRCT
    tmp =[];

    for s = 1:num_segs/2 

        clust_vec = [clust_in_epochs{s}; clust_in_epochs{s+20}];
        clust_vec(removed_spikes{s}) =[]; 

        id = find(clust_vec == c);
        
        j = all_EIT_spikes_prct{s}(id,:,:);     

        k = squeeze(mean(j,1, 'omitnan')); 

        if s > 40
            s2 = s - 40; 
        elseif s > 20
            s2 = s - 20; 
        elseif s < 20
            s2 = s; 
        end 

        inj_chans = find(full(stim(s2).stim_pattern) ~=0);

        inj_chans = [inj_chans; 32];

        k(:, inj_chans) = [];

        tmp = [tmp, k];

        

    end

    clust_EIT_spikes_prct(c,:,:) = tmp;
    
    % RAW 
    
     tmp =[];

    for s = 1:num_segs/2

        clust_vec = [clust_in_epochs{s}; clust_in_epochs{s+20}];
        clust_vec(removed_spikes{s}) =[]; 

        id = find(clust_vec == c);
        
        j = all_EIT_spikes_raw{s}(id,:,:);     

        k = squeeze(mean(j,1, 'omitnan')); 

        if s > 40
            s2 = s - 40; 
        elseif s > 20
            s2 = s - 20; 
        elseif s < 20
            s2 = s; 
        end 

        inj_chans = find(full(stim(s2).stim_pattern) ~=0);

        inj_chans = [inj_chans; 32];

        k(:, inj_chans) = [];

        tmp = [tmp, k];

        

    end

    clust_EIT_spikes_raw(c,:,:) = tmp;

end 




%% Save EIT data

cd(outpath)

save('EIT_spikes_new', 'av_EIT_spks_raw', 'all_EIT_spikes_raw','EIT_all_seg_raw', 'EIT_all_seg_raw', 'clust_EIT_spikes_raw',...
    'ds_hilb','all_EIT_spikes_prct', 'av_EIT_spks_prct','clust_EIT_spikes_prct','EIT_all_seg_prct', 'rm_sv','str_stp','stim','-v7.3'); 

