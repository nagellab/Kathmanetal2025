clear all
% close all

use_flies = [21:27, 30:39]; %cwoo
use_flies = [21:26, 30, 32:35, 37:39]; %cwoo
sdm = .5;

%% Set up
headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Desktop/';
cd(headpath)

allflies_filenames

for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    fly_struct{f} = syncstruct;
end


% fly_struct = struct_PFG;
scale = 1;

% Initialize matrix to hold data
mat_dff_off = NaN(1, 100*25); % 15 seconds of data
mat_dff_on = NaN(1, 100*25); % 15 seconds of data

count = 0;
odor_offsets = [];
odor_onsets = [];

for fly = 1:length(fly_struct)
    data = fly_struct{fly};

    % Find bump threshold from all trials
    all_dff = cat(1, data(arrayfun(@(t) isfield(t, 'bump_amp') && ~isempty(t.bump_amp) && any(~isnan(t.bump_amp)), data)).bump_amp);
    all_dff = all_dff(~isnan(all_dff));
    base_thr = nanmean(all_dff)+nanstd(all_dff)*sdm;


    use_trials = find(arrayfun(@(t) t.closed_loop && sum(diff(t.odor) == 1) == 1 && sum(diff(t.odor) == -1) == 1, data));% && all(isnan(t.plume_t)), data); %cwoo

    for trial = use_trials

        [trial_t, trial_dff, odoroni, odoroffi, trial_odor] = extract_struct_data(data(trial));



%         % Grab appropriate data from trial struct
%         trial_odor = data(trial).odor; % odor signal
%         trial_dff = data(trial).dff; % bump amplitude
%         trail_t = data(trial).tt; % time
    
        % Find odor start and end
        trial_odor_bin = find(trial_odor > 0);
        odor_start = trial_odor_bin(1);
        odor_end = trial_odor_bin(end);

        cutoff_offset = odor_end; % offset has to happen after end of odour
        cutoff_onset = odor_end; % onset has to happen before end of odour
    
        % Define filter parameters
        fs = 100; fc = 2; % Designs a band-pass filter (cutoff frequency = 2 Hz)
        [b, a] = butter(2,fc/fs); % Designs a band-pass filter
        trial_bump_filt = filtfilt(b,a,trial_dff);
    
        state_positive = sum(diff(trial_bump_filt(end-200:end)) > 0) / length(diff(trial_bump_filt(end-200:end))) == 1;
        state_negative = sum(diff(trial_bump_filt(end-200:end)) < 0) / length(diff(trial_bump_filt(end-200:end))) == 1;
        
        % Check if dff signal is ok    
        if ~(state_positive || state_negative)
            % Grab bump onset and offset times
%             [bump_onset_idx_filt, bump_offset_idx_filt] = bump_finder(data, trial, 0, scale); % bump onset/offset
            oni = []; offi = []; 
            [bumpi, bump_onset_idx_filt, bump_offset_idx_filt] = bumpfinder2(trial_dff, trial_t, 3, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration

            % Check delay from odor offset to bump offset
            bump_offset_time = 0;
            bump_onset_time = -90000;
            for bump = 1:length(bump_onset_idx_filt)
                bump_onset = bump_onset_idx_filt(bump);
                bump_offset = bump_offset_idx_filt(bump);

                if bump_offset > cutoff_offset && bump_onset < cutoff_onset && bump_onset > odor_start
                    bump_offset_time = trial_t(bump_offset) - trial_t(odor_end);
                    bump_onset_time = trial_t(bump_onset) - trial_t(odor_start);
                end
            end

            odor_offsets = [odor_offsets, bump_offset_time];
            odor_onsets = [odor_onsets, bump_onset_time];

            count = count + 1;
    
            trial_dff_off = trial_dff(odor_end-5*100+1:odor_end+20*100);
            trial_dff_on = trial_dff(odor_start-5*100+1:odor_start+20*100);
    
            mat_dff_off(count, :) = trial_dff_off;
            mat_dff_on(count, :) = trial_dff_on;
            
        end
    end
end

[a_sorted, a_order] = sort(odor_offsets, 'descend');
mat_dff_off = mat_dff_off(a_order,:);

mat_dff_off = mat_dff_off(a_sorted > 0, :);


t = -5:0.01:20;
t = t(1:end-1);

figure; hold on
clm = [0, max(max(mat_dff_off))];
imagesc(t, (1:size(mat_dff_off, 1)),mat_dff_off, clm)
colormap('gray')
% colormap(redblue(100))
% colormap redblue
xlim([-5,20])
xlabel('time to odor offset (s)')
ylim([0, size(mat_dff_off, 1)+1])
ylabel('Trial #')
cb = colorbar(); 
ylabel(cb,'Bump amplitude','FontSize',12,'Rotation',90)
hColorbar.Label.Position(1) = 10;
plot([0,0], [0.5, size(mat_dff_off, 1)-0.5], '--r', 'LineWidth',1.5)
title(['Fly: ', num2str(fly)])
caxis([0,3]);
caxis([0,1]);

%%
%%

function [tt, dff, oncrossi, offcrossi, odor, fps] = extract_struct_data(trl)

%     nonan = ~isnan(trl.calc_ts(:,1));
    nonan = ~isnan(trl.bump_amp(:,1));
    fps = cat(1,trl.fps);%fps = median(fps(~isnan(fps)));
    dff = trl.bump_amp(nonan, 1);
%     heading = smoothdata(trl.calc_heading(nonan), 'movmean', 10); %cwoo or 2pstim
    heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(trl.calc_windpos(nonan))), 'movmean', 100)));
    uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 1);

    speed = trl.calc_speed(nonan);
    x = -cumsum(speed.*sin(heading*pi/180))/fps;
    y = cumsum(speed.*cos(heading*pi/180))/fps;
    tt = trl.calc_ts(nonan);
    avel = smoothdata(trl.calc_deltaz(nonan), 'movmean', 150);
    abs_avel = smoothdata(abs(trl.calc_deltaz(nonan)), 'movmean', 50);
    uwv = cos(heading*pi/180).*speed;
    fvel = trl.calc_deltapitch(nonan)*9.52/2;
    fvel = smoothdata(fvel, 'movmean', 100); 
    %find stim times before trial sorting bc used to rule out plume trials
    wind = trl.wind(nonan);
    odor = trl.odor(nonan);


    oncrossi = find(diff(odor) == 1);
    offcrossi = find(diff(odor) == -1);
    if odor(1) == 1
        oncrossi = [1; oncrossi]; end
    if odor(end) == 1
        offcrossi = [offcrossi; length(odor)]; end
    
%     oncrossi = find(diff(odor > .5) == 1, 1, 'first');
%     offcrossi = find(diff(odor > .5) == -1, 1, 'last');

end