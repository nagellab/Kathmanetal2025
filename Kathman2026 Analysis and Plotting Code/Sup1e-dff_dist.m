%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2024 paper/sandbox/exploratory figs/bump postion/';
cd(headpath)

allflies_filenames

%%
% use_flies = [21:26,30, 33:35, 37:39]; sdm = 1;
use_flies = [21:27, 30, 32:39]; sdm = .5;
% use_flies = [41:51]; sdm = .5;

print_figs = 0;

fig_no = 1;
all_stim = cell(1,3);
all_t = 0;
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;

    %find use_trials and base_thresh for whole expt
    allamp_thr = []; use_trials = [];
    for t = 1:length(data)
        if ~isempty(data(t).bump_amp)
            nonan = ~isnan(data(t).bump_amp(:,1));
    
            %find stim times before trial sorting bc used to rule out plume trials
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)

            wind = data(t).wind(nonan);
            idyl = wind >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            windoni = idy(wind(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
    
            allamp(t) = quantile(data(t).bump_amp(~isnan(data(t).bump_amp(:,1))), .95); 
            allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
            if length(odoroni) == 1  && data(t).closed_loop %&& length(windoni) == 1 
                use_trials = [use_trials, t];
            end
        end
    end
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;

    if ~isempty(use_trials)
        for t = use_trials
            nonan = ~isnan(data(t).bump_amp(:,1));

            %% find stim times (again, dumb)
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            odoroni = odoroni + 10;
            odoroffi = idy(odor(idy(idy < length(odor))+1) < .5); %all pts where next pt is subthresh (aka downward crossing)
            odoroffii = odoroffi + 10; %compesate for odor delay (100ms)

            wind = data(t).wind(nonan);
            idyl = wind >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            windoni = idy(wind(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            windoni = windoni + 10;
            windoffi = idy(wind(idy(idy < length(wind))+1) < .5); %all pts where next pt is subthresh (aka downward crossing)
            windoffii = windoffi + 10; %compesate for odor delay (100ms)            

            %% find dff and behav metrics
            dff = data(t).bump_amp(nonan, 1);
            speed = data(t).calc_speed(nonan);
            fps = cat(1,data.fps);fps = median(fps(~isnan(fps))); %fictrac fs
            if data(t).closed_loop
                heading = data(t).calc_windpos(nonan);
            else
                heading = data(t).calc_heading(nonan);
            end
            heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
            uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 100); 
            tt = data(t).calc_ts(nonan);
            col = data(t).bump_col(nonan, 1);
            com = data(t).bump_com(nonan, 1);
            pos = smoothdata(data(t).bump_col(nonan), 'movmean', 100); %200);
            avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 200);
            abs_avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 50);
            uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
            fvel = data(t).calc_deltapitch(nonan)*9.52/2;
            fvel = smoothdata(fvel, 'movmean', 200); 

            %% find bump data
            im = data(t).dff_green(1:8,:);
            ts = data(t).dff_ts;
            fs = size(im, 2)/max(ts); %thor fs           
            [pos, uwpos, ~, thr] = bumpposfinder(im, ts, tt, base_thr, sdm); %thr is min of expt thresh (base_thr) and trial thresh

%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 3, thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 3, base_thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
            [bump, oni, offi, thrsh] = bumpfinder2(max(im), ts, 5, 3, thr); 
            bumpi = find(bump);
            
            %% meat
            oi = find(odor);
            do = dff(oi);
            all_stim{3} = [all_stim{3}; do];
%             all_stim{3} = [all_stim{3}, mean(do)];

            wi = find(wind & ~odor);
            dw = dff(wi);
            all_stim{2} = [all_stim{2}; dw];
%             all_stim{2} = [all_stim{2}, mean(dw)];

            bi = find(~wind & ~odor);
            db = dff(bi);
            all_stim{1} = [all_stim{1}; db];
%             all_stim{1} = [all_stim{1}, mean(db)];
            all_t = all_t + 1;
        end
    end
end
% cmap = lines(3);
cmap = {'k'; 'g'; 'm'};
figure(fig_no); clf; hold on
set(gcf, 'Position', [1 1 1020 550]); %[100 100 500 1000])
for s = 1:3

%     hist(all_stim{s}, 100)

    bnw = 0.005;
    edg = -.2:bnw:2.5;
    N = histcounts(all_stim{s}, edg, 'Normalization', 'probability'); 
    N = smoothdata(N, 'movmean', 10);
    plot(edg(2:end)-bnw/2, N, 'Color', cmap{s}, 'LineWidth', 2)
%     plot(edg(2:end)-bnw/2, N, 'Color', cmap(s, :), 'LineWidth', 2)
    set(gca,'YScale','log')
end

legend('No stim', 'Wind only', 'Wind+Odor')
xlabel('∆F/F')
ylabel('Probability')
title(['total_trls = ', num2str(all_t), '; nBlank=' num2str(size(all_stim{1},1)), ', nWind=' num2str(size(all_stim{2},1)), ', nOdor=' num2str(size(all_stim{3},1))])




