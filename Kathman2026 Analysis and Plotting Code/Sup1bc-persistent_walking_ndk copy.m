% all_flys = {struct_hdk, struct_PFG, struct_VT062617};
clear all
close all
% use_flies = [27, 32, 36, 38, 109:112]; %plumes
use_flies = [21:27, 30:39]; %cwoo


headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Desktop/';
cd(headpath)

allflies_filenames
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    use_trials = find(arrayfun(@(t) t.closed_loop && sum(diff(t.odor) == 1) == 1 && sum(diff(t.odor) == -1) == 1, syncstruct));% && all(isnan(t.plume_t)), data); %cwoo
    all_flys{1}{f} = syncstruct(use_trials);
end

all_scales = [1, 2, 1];

% fly_structs = struct_hdk;
% scale = 1;

% Filter
fs = 100; fc = 2; % Designs a band-pass filter (cutoff frequency = 2 Hz)
[b, a] = butter(2,fc/fs); % Designs a band-pass filter

array_goals = {};
array_goals_dir = {};
array_bumps = {};

array_goals_ctrl = {};
array_goals_dir_ctrl = {};

all_dur = [];
all_dur_ctrl = [];

count = 0;

for cell = 1:length(all_flys)
    fly_structs = all_flys{cell};
    scale = all_scales(cell);

    % Loop through flies
    for fly = 1:length(fly_structs)
    % for fly = 1:1
        
        goal_dur = [];
        bump_dur = [];
        goal_dir = [];

        goal_dur_ctrl = [];
        goal_dir_ctrl = [];

        % Grab fly data
        data = fly_structs{fly};

        %find good trials and bump threshold
        lengths = arrayfun(@(x) length(x.odor), data);
        use_t = [];
        allamp_thr = [];
        for t = 1:length(data)
            odoron = find(diff(data(t).odor) == 1);
            allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
            if length(data(t).odor) == mode(lengths) && length(odoron) == 1 && data(t).closed_loop == 1
                [data(t).tt, data(t).dff, data(t).odor, data(t).x, data(t).y, data(t).heading, data(t).fvel] = extract_struct_data(data(t));
                use_t = [use_t, t];
            end
        end
        base_thr = mean(allamp_thr)+std(allamp_thr)*0.5;

    
        % Loop through trials
        for trial = use_t
        % for trial = 3:3
            

            % Grab data
            trial_x = data(trial).x;
            trial_y = data(trial).y;
%             trial_upos = data(trial).pos_uw;
%             trial_head = data(trial).uwheading;
%             trial_fvel = data(trial).fvel;
            trial_dff = data(trial).dff;
            trial_odor = data(trial).odor;
            trial_tt = data(trial).tt;

%             % Grab data
%             [trial_tt, trial_dff, trial_odor, trial_x, trial_y, trial_head, trial_fvel] = extract_struct_data(data(trial));

   
            % Find odor onset and offset
            trial_odor_bin = find(trial_odor > 0);
            odor_start = trial_odor_bin(1);
            odor_end = trial_odor_bin(end);
    
            cutoff_offset = odor_end; % offset has to happen after end of odour
            cutoff_onset = odor_end; % onset has to happen before end of odour
    
            % Check if behavioural and dff data is ok
            % Determine if data is ok
            trial_bump_filt = filtfilt(b,a,trial_dff);
            state_positive = sum(diff(trial_bump_filt(end-200:end)) > 0) / length(diff(trial_bump_filt(end-200:end))) == 1;
            state_negative = sum(diff(trial_bump_filt(end-200:end)) < 0) / length(diff(trial_bump_filt(end-200:end))) == 1;
            state_behavior = max([max(diff(trial_x)), max(diff(trial_y))]) > 5;
            if ~(state_positive || state_negative || state_behavior)
                
                % Find persistent walking
                [stable_onset, stable_offset, goal_dir_i] = persist_walk_finder_v2(data, trial, 0);
                % Find bumps
%                 [bump_onset, bump_offset] = bump_finder(data, trial, 0, scale); % bump onset/offset
                [~, bump_onset, bump_offset] = bumpfinder2(trial_dff, trial_tt, 3, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration

                % Find persistent walking control (baseline)
                [stable_onset_ctrl, stable_offset_ctrl, goal_dir_ctrl_i] = persist_walk_finder_v2_ctrl(data, trial, 0);
    
                % Variables to track if goal / bump straddle odor offset
                state_goal = false;
                state_bump = false;
                
                
    
                %% Loop through bump periods
                for bump = 1:length(bump_onset)
                    % Grab data
                    bump_onset_i = bump_onset(bump);
                    bump_offset_i = bump_offset(bump);
    
                    if bump_offset_i > cutoff_offset && bump_onset_i < cutoff_onset
                        bump_dur = [bump_dur, trial_tt(bump_offset_i) - trial_tt(odor_end)];
                        goal_dur = [goal_dur, trial_tt(stable_offset) - trial_tt(odor_end)];
                        goal_dir = [goal_dir, goal_dir_i];
                        goal_dur_ctrl = [goal_dur_ctrl, trial_tt(stable_offset_ctrl) - trial_tt(1)];
                        goal_dir_ctrl = [goal_dir_ctrl, goal_dir_ctrl_i];
                        all_dur = [all_dur, trial_tt(stable_offset) - trial_tt(odor_end)];
                        all_dur_ctrl = [all_dur_ctrl, trial_tt(stable_offset_ctrl) - trial_tt(1)];
                        state_bump = true;
                    end                
                end
            end
        end
        count = count + 1;
        array_goals{1, count} = goal_dur;
        array_bumps{1, count} = bump_dur;
        array_goals_dir{1, count} = goal_dir;

        array_goals_ctrl{1, count} = goal_dur_ctrl;
        array_goals_dir_ctrl{1, count} = goal_dir_ctrl;
    end
end


%%
figure; hold on
dur_goal_all = [];
dur_bump_all = [];
dir_goal_all = [];
dur_goal_all_ctrl = [];
dir_goal_all_ctrl = [];


plot([0,30], [0, 30], '--k')
for fly = 1:length(array_goals)
    goal_dur = array_goals{fly};
    goal_dir = array_goals_dir{fly};
    bump_dur = array_bumps{fly};
    goal_dur_ctrl = array_goals_ctrl{fly};
    goal_dir_ctrl = array_goals_dir_ctrl{fly};

    dur_goal_all = [dur_goal_all, goal_dur];
    dur_bump_all = [dur_bump_all, bump_dur];
    dir_goal_all = [dir_goal_all, goal_dir];

    dur_goal_all_ctrl = [dur_goal_all_ctrl, goal_dur_ctrl];
    dir_goal_all_ctrl = [dir_goal_all_ctrl, goal_dir_ctrl];

    plot(goal_dur, bump_dur, 'o', 'Color', 'k')
    
end
xlim([0,30])
ylim([0,30])
xlabel('goal durations (s)')
ylabel('bump durations (s)')

dur_goal_all = dur_goal_all';
dur_bump_all = dur_bump_all';
dur_goal_all_ctrl = dur_goal_all_ctrl';

%%
dur_goal_all;
dir_goal_all;

% dur_goal_all(dur_goal_all > 30) = 30;

% figure; hold on
xplot = dur_goal_all' .* sin(dir_goal_all*pi/180);
yplot = dur_goal_all' .* cos(dir_goal_all*pi/180);

figure; hold on

% Plot the vectors
for i = 1:size(xplot,2)
    plot([0, xplot(i)], [0, yplot(i)], 'color', [90,90,90]/255)
end

% Mean direction vector
plot([0, mean(xplot)], [0, mean(yplot)], '-r', 'LineWidth',1.5)

% Polar plot setup
max_radius = 30; % or set to max(dur_goal_all)
radii = [10, 20, 30]; % example radii in seconds
angle_ticks = 0:30:330; % degrees for tick marks
label_angle = 165; % angle (in degrees) where to place circle labels

% Concentric circles as duration grid lines
theta = linspace(0, 2*pi, 500);
for r = radii
    % Draw circle
    plot(r*sin(theta), r*cos(theta), 'k:', 'LineWidth', 0.8)
    
    % Label each circle (at a fixed angle)
    angle_rad = deg2rad(label_angle);
    xlab = (r + 0.8) * sin(angle_rad);  % small offset outward
    ylab = (r + 0.8) * cos(angle_rad);
    text(xlab, ylab, sprintf('%d s', r), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 9)
end

% Degree ticks around circle
for ang = angle_ticks
    r_tick = max_radius * 1.05;
    x_outer = r_tick * sin(deg2rad(ang));
    y_outer = r_tick * cos(deg2rad(ang));
    x_inner = (r_tick - 1) * sin(deg2rad(ang));
    y_inner = (r_tick - 1) * cos(deg2rad(ang));
    plot([x_inner x_outer], [y_inner y_outer], 'k')
    text(x_outer * 1.05, y_outer * 1.05, sprintf('%d°', ang), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8)
end

% Formatting
xlim([-1.2 1.2] * max_radius)
ylim([-1.2 1.2] * max_radius)
axis equal
axis off
title('Direction and Duration Vectors (odor)')


% xlim([-30.5,30.5])
% ylim([-30.5,30.5])
% 
% x = -30:0.01:30;
% plot(x, sqrt(30^2 - x.^2), '-k')
% plot(x, -sqrt(30^2 - x.^2), '-k')
% 
% axis equal

%% Ctrls
dur_goal_all_ctrl;
dir_goal_all_ctrl;


dur_goal_all_ctrl(dur_goal_all_ctrl > 20) = 20;

% figure; hold on
xplot = dur_goal_all_ctrl' .* sin(dir_goal_all_ctrl*pi/180);
yplot = dur_goal_all_ctrl' .* cos(dir_goal_all_ctrl*pi/180);

figure; hold on

% Plot the vectors
for i = 1:size(xplot,2)
    plot([0, xplot(i)], [0, yplot(i)], 'color', [90,90,90]/255)
end

% Mean direction vector
plot([0, mean(xplot)], [0, mean(yplot)], '-r', 'LineWidth',1.5)


% Polar plot setup
max_radius = 30; % or set to max(dur_goal_all)
radii = [10, 20, 30]; % example radii in seconds
angle_ticks = 0:30:330; % degrees for tick marks
label_angle = 165; % angle (in degrees) where to place circle labels

% Concentric circles as duration grid lines
theta = linspace(0, 2*pi, 500);
for r = radii
    % Draw circle
    plot(r*sin(theta), r*cos(theta), 'k:', 'LineWidth', 0.8)
    
    % Label each circle (at a fixed angle)
    angle_rad = deg2rad(label_angle);
    xlab = (r + 0.8) * sin(angle_rad);  % small offset outward
    ylab = (r + 0.8) * cos(angle_rad);
    text(xlab, ylab, sprintf('%d s', r), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 9)
end

% Degree ticks around circle
for ang = angle_ticks
    r_tick = max_radius * 1.05;
    x_outer = r_tick * sin(deg2rad(ang));
    y_outer = r_tick * cos(deg2rad(ang));
    x_inner = (r_tick - 1) * sin(deg2rad(ang));
    y_inner = (r_tick - 1) * cos(deg2rad(ang));
    plot([x_inner x_outer], [y_inner y_outer], 'k')
    text(x_outer * 1.05, y_outer * 1.05, sprintf('%d°', ang), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8)
end

% Formatting
xlim([-1.2 1.2] * max_radius)
ylim([-1.2 1.2] * max_radius)
axis equal
axis off
title('Direction and Duration Vectors (baseline)')


%%
figure; hold on
prob_all = NaN(length(array_goals), 10);
tot_trials = 0;
for fly = 1:length(array_goals)
    data = array_goals{fly};
    [counts,centers] = hist(data, 3:3:30);
    prob_all(fly, :) = counts / sum(counts);
    tot_trials = tot_trials + length(array_goals{fly});
end

hist_mean = nanmean(prob_all);
hist_sem = nanstd(prob_all) / sqrt(size(prob_all, 1));

% plot(centers, hist_mean, 'color', [3, 165, 80]/255, 'LineWidth',2)
% plot(centers, hist_mean+hist_sem, 'color', [215,230,211]/255, 'LineWidth',2)
% plot(centers, hist_mean-hist_sem, 'color', [215,230,211]/255, 'LineWidth',2)
% 

% Define upper and lower bounds
upper = hist_mean + hist_sem;
lower = hist_mean - hist_sem;

% Create fill area (light green shade)
fill([centers, fliplr(centers)], ...
     [upper, fliplr(lower)], ...
     [215,230,211]/255, ...      % Light green shade
     'EdgeColor', 'none', ...
     'FaceAlpha', 1);            % 1 = opaque, <1 = transparent

% Plot mean line on top
plot(centers, hist_mean, 'color', [3, 165, 80]/255, 'LineWidth', 2)


xlabel('post-odor duration')
ylabel('probability')
% ylim([0,0.36])
xlim([0,30])

prob_all = NaN(length(array_goals_ctrl), 10);
for fly = 1:length(array_goals_ctrl)
    data = array_goals_ctrl{fly};
    [counts,centers] = hist(data, 3:3:30);
    prob_all(fly, :) = counts / sum(counts);
end

hist_mean = nanmean(prob_all);
hist_sem = std(prob_all) / sqrt(size(prob_all, 1));

plot(centers, hist_mean, 'color', 'k', 'LineWidth',2)
plot(centers, hist_mean+hist_sem, 'color', [111,111,111]/255, 'LineWidth',2)
plot(centers, hist_mean-hist_sem, 'color', [111,111,111]/255, 'LineWidth',2)

xlabel('post-odor duration')
ylabel('probability')
ylim([0,0.7])
xlim([0,30])
title(['nTrials=' num2str(tot_trials)])


%% aarons plot of same data

centers = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
counts_odor = [0.300751879699248, 0.270676691729323, 0.105263157894737, 0.0827067669172932, 0.105263157894737, 0.0676691729323308, 0.0225563909774436, 0.0451127819548872, 0, 0];
counts_base = [0.521212121212121, 0.151515151515152, 0.0727272727272727, 0.0181818181818182, 0.0242424242424242, 0, 0.00606060606060606, 0, 0, 0.0121212121212121];



xdata = centers';
ydata = counts_odor';


figure(5); hold on
tau1 = decayfit(xdata, ydata, 'm');
ydata = counts_base';
tau2 = decayfit(xdata, ydata, 'k');


xlabel('Straight walking durs (post-dor v base)')
ylabel('Probability')
legend('95% CI', 'Odor runs', 'Fit', '95% CI', 'Base runs', 'Fit')
title(['Exp Fit with Shaded 95% CI, post-odor tau=', num2str(tau1), ', base tau=', num2str(tau2), 'n=135'])


function tau = decayfit(xdata, ydata, cstring)
    % Create a fit type: A*exp(-x/tau) + C
    ft = fittype('A*exp(-x/tau) + C', ...
        'independent', 'x', 'coefficients', {'A', 'tau', 'C'});
    
    % Fit the model
    startPoints = [max(ydata)-min(ydata), 3, min(ydata)];
    [fitresult, gof] = fit(xdata, ydata, ft, 'StartPoint', startPoints);
    fitresult
    tau = fitresult.tau;

    % Evaluate the fitted curve and confidence intervals
    xfit = linspace(min(xdata), max(xdata), length(xdata));
    [yCI, yfit] = predint(fitresult, xfit, 0.95, 'functional');  % 'functional' gives CI on mean prediction
    
    % plot
    x_fill = [xfit, fliplr(xfit)];
    y_fill = [yCI(:,1)', fliplr(yCI(:,2)')];
    fill(x_fill, y_fill, cstring, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.15);   % lower = more transparent
    %     plot(cnt, N, 'Color', 'k', 'LineWidth', 2)
    plot(xdata, ydata, 'Color', cstring, 'LineWidth', 2)
    plot(xfit, yfit, ['--' cstring], 'LineW', 1)
end

function [tt, dff, odor, x, y, uwheading, fvel] = extract_struct_data(trl)

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