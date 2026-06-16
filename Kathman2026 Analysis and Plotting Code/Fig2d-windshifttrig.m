%% all fly
%%
clear all
% close all

headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/';
cd(headpath)

allflies_filenames

%%
use_flies = [116:120]; %[61:65];%, 109:112]; 
sdm = 0.6; 
print_figs = 0;
add2txt = '';

shift = 35; %s (can find in allo_pos, but these are all at i = 35, so hard coded it out of laziness)

figure(1); clf; hold on;
set(gcf, 'Position', [1 1 620 1550]);

figure(2); clf; hold on;
set(gcf, 'Position', [1 1 1420 550]);

figure(200); clf;
set(gcf, 'Position', [1 1 1420 950]);

k = 0;
j = 1;
ii = 0;
prevY1 = 1;
prevY2 = 1;
prevY3 = 1;
prevY4 = 1;

shiftdiff_pos = [];
shiftdiff_head = [];

%%
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;


    allamp_thr = []; use_trials = [];
    for t = 1:length(data)
        if ~isempty(data(t).bump_amp) && ~isempty(data(t).wind)
            %metrics to sort trials by
            nonan = ~isnan(data(t).bump_amp(:,1));

            %find stim times before trial sorting bc used to rule out plume trials
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
    
            allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
    %         if length(oncrossi) > 5
%             if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 && max(odor) > 0%cwoo
    %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
    %         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 && max(odor) > 0%owoo
                use_trials = [use_trials, t];
%             end
        end
    end
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;

    if ~isempty(use_trials)
        for t = use_trials
            nonan = ~isnan(data(t).bump_amp(:,1));
    
            %find stim times before trial sorting bc used to rule out plume trials
            odor = data(t).odor(nonan);
            wind = data(t).wind(nonan);
            if data(t).shift == 90
                allo_pos = data(t).allo_pos(nonan);
            end
            ampthresh = .5;%0.125;
            idyl = odor >= ampthresh; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(odor(idy-1)<ampthresh);
    

            if ~isempty(data(t).bump_amp) %closed loop, odor on, wind on, not plume, not bad/empty trial
    
                %% metrics to plot or compute by
                offcrossi = idy(odor(idy(idy < length(odor))+1)<ampthresh); %all pts where next pt is subthresh (aka downward crossing)
                offcrossi = offcrossi + 10; %compesate for odor delay (100ms)

                dff = data(t).bump_amp(nonan, 1);
                speed = data(t).calc_speed(nonan);
                fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%                 if data(t).closed_loop
%                     heading = data(t).calc_windpos(nonan);
%                 else
                    heading = data(t).calc_heading(nonan);
%                 end
                heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
                uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 1);   
                xx = -cumsum(speed.*sin(heading*pi/180))/fps;
                yy = cumsum(speed.*cos(heading*pi/180))/fps;                
                tt = data(t).calc_ts(nonan);
                col = data(t).bump_col(nonan, 1);
                com = data(t).bump_com(nonan, 1);
%                 uwpos = unwrap(data(t).bump_col(nonan);
                oldpos = smoothdata(data(t).bump_col(nonan), 'movmean', 200); %200);
                newpos = data(t).bump_uwspos(nonan);
%                 pos = data(t).bump_col(nonan);
                avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 200);
                abs_avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 50);
                uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
                fvel = data(t).calc_deltapitch(nonan)*9.52/2;
                fvel = smoothdata(fvel, 'movmean', 700); 

                im = data(t).dff_green(1:8,:);
                ts = data(t).dff_ts;
                fs = size(im, 2)/max(ts);

                    
                %% find pos
                %smooth data
                sim = [];
                for c = 1:8
                    sim(c,:) = smoothdata(im(c,:), 'movmean', 13);
                end
                for r = 1:size(im,2)
                    sim(:,r) = smoothdata(sim(:,r), 'movmean', 3);
                end

                %find bump amp (m), peak position (p), and set a threshold (thr) for this trial
                [m,p] = max(sim); 
                trl_thr = mean(m)+(sdm*std(m)); %std of single smoothed trial
                thr = min([trl_thr, base_thr]); %use smallest of single trial and expt thresholds
%                 thr = base_thr;

                %convert position to radians
                p = (p-1)*(2*pi)/7.5;

                %find bumps with thor times
                % bmp = find(m > thr); %basic threshold (without gap filler/min bumplength)
                [bump, oni, offi, thrsh] = bumpfinder2(m, ts, 5, 5, thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
%                 [bump, oni, offi, thrsh] = bumpfinder2(m, ts, 5, 3, thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
                bumpi = find(bump);
                [~, shifti]=min(abs(ts-shift)); %shift index in thor time
                obumps = find(offi > shifti & oni < shifti); %which on/off straddles shift
                [~, ficshifti]=min(abs(tt-shift)); %shift index in fictrac time


%                 %find bumps with fictrac times
%                 ficbumpi = [];
%                 for i = 1:length(oni)
%                     ficbumpi = [ficbumpi; find(tt >= ts(oni(i)) & tt <= ts(offi(i)))];
%                 end
                
                
                %make bump pos stay put when off
                pos = nan(1,length(m));
                if ~isempty(bumpi)
                    pos(bumpi) = p(bumpi);
                    for i = 1:length(m)
                        if i > bumpi(1) && isnan(pos(i))
                            pos(i) = pos(i-1);
                        end
                    end
                end

%                 %unwrap, smooth, rewrap, set start pos, convert to column dimensions
%                 spos = wrapTo2Pi(smoothdata(unwrap(pos), 'movmean', 15));
%                 spos(isnan(spos)) = pi;
%                 spos = (spos/(2*pi/7)) + 1;

                %unwrap, smooth, set start pos, convert to column dimensions
                spos_uw = smoothdata(unwrap(pos), 'movmean', 45);
                spos_uw(isnan(spos_uw)) = pi;
                spos_uw = (spos_uw/(2*pi/7)) + 1; 
%                 spos_uw = 



%                 k = bumpi >= shifti - 10 & bumpi <= shifti + 10; 
                if any(obumps) && abs(data(t).shift) > 0 && any(bumpi >= shifti - 10 & bumpi <= shifti - 3) %bumpi 10s before shift
%                     plot(ts, spos_uw - spos_uw(shifti))
                    
                    if data(t).shift > 0 %invert bump pos and heading for right shifts to normalize direction
                        spos_uw = -spos_uw;
                        uwheading = -uwheading;
                    end
                    spos_uw = spos_uw * 180/pi; %in deg
                    uwheading = uwheading * 180/pi; %in deg

                    bumpi = bumpi(bumpi >= oni(obumps) & bumpi <= offi(obumps));
                    ficbumpi = find(tt >= ts(oni(obumps)) & tt <= ts(offi(obumps)));
                    bumpi = bumpi(bumpi >= shifti - (fs*5) & bumpi <= shifti + (fs*2) + (fs*5));
                    ficbumpi = ficbumpi(ficbumpi >= ficshifti - (fps*5) & ficbumpi <= ficshifti + (fps*2) + (fps*5));
                    
                    if ismember(ii+1, [4, 6, 10, 12, 14, 15, 16, 17]) %replace w turnfinder
%                     if ismember(ii+1, [1, 2, 4, 6, 10, 12, 14, 15, 16, 17]) 
%                     if ismember(ii+1, [2, 4, 6, 10, 12, 14, 15, 16, 17, 21]) 
%                     if ismember(ii+1, [4, 6, 10, 12, 14, 15, 16, 17]) 


                    figure(1);
                    subplot(1,3,1); hold on
                    y = spos_uw(bumpi)-min(spos_uw(bumpi))+prevY1+10;
                    scatter(ts(bumpi), y, 10, 'filled')
                    prevY1 = max(y);

%                     subplot(1,4,2); hold on
%                     scatter(tt(ficbumpi), oldpos(ficbumpi)+prevY2+10, 10, 'k', 'filled')
%                     prevY2 = max(oldpos(ficbumpi)) + prevY2 + 10;

                    subplot(1,3,2); hold on
                    y = uwheading(ficbumpi)-min(uwheading(ficbumpi))+prevY3+(10 * 180/pi); %in deg;
                    scatter(tt(ficbumpi), y, 10, 'filled')
                    prevY3 = max(y);
%                     if ismember(ii+1, [2, 4, 6, 10, 12, 14, 15, 16, 17, 21]) 
                    if ismember(ii+1, [4, 6, 10, 12, 14, 15, 16, 17]) 
%                     if ismember(ii+1, [1, 2, 4, 6, 10, 12, 14, 15, 16, 17]) 
                        plot(43, y(end), '*r')
                    end

                    subplot(1,3,3); hold on
%                     y = abs_avel(ficbumpi)-min(fvel(ficbumpi))+prevY4+10;
                    y = fvel(ficbumpi)-min(fvel(ficbumpi))+prevY4+(10 * 180/pi);
                    scatter(tt(ficbumpi), y, 10, 'filled')
                    prevY4 = max(y);  

                    figure(100);
                    subplot(5, 5, j);
                    scatter(xx, yy, 10, 'k', 'filled')
                    xlim([-350 350])
                    ylim([-150 550])
                    title(['fly' num2str(f), ', t' num2str(t)])
                    j = j + 1;



                    figure(2); %clf
                    nrm = 1; %flag to normalize or not
                    normtime = 0.25; %sec before shift
                    subplot(1,3,1); hold on
                    cropi = bumpi(ts(bumpi) > 35-3); 
                    plot(ts(cropi)-35, spos_uw(cropi)-(spos_uw(shifti-round(normtime*fs))*nrm), 'k')
%                     plot(ts(bumpi)-35, spos_uw(bumpi)-(spos_uw(shifti-round(normtime*fs))*nrm), 'k')
                    plot([shift-35 shift-35], [-200 460], '--k')
                    plot([shift+2-35 shift+2-35], [-200 460], '--k')
                    ylim([-200 460])
                    ylabel(['Bump Pos(deg)'])
                    xlabel(['Time (s)'])
    
                    ficcropi = ficbumpi(tt(ficbumpi) > 35-3);
                    subplot(1,3,2); hold on
                    plot(tt(ficcropi)-35, uwheading(ficcropi)-(uwheading(ficshifti-round(normtime*fps))*nrm), 'k')
%                     plot(tt(ficbumpi)-35, uwheading(ficbumpi)-(uwheading(ficshifti-round(normtime*fps))*nrm), 'k')
                    plot([shift-35 shift-35], [-200 460], '--k')
                    plot([shift+2-35 shift+2-35], [-200 460], '--k')
                    ylim([-200 460])
                    ylabel(['Heading (deg)'])
                    xlabel(['Time (s)'])

                    ficropi = ficbumpi(tt(ficbumpi) > 35-3);
                    subplot(1,3,3); hold on
%                     plot(tt(ficcropi)-35, abs_avel(ficcropi)-(fvel(ficshifti-round(normtime*fps))*nrm), 'k')
                    plot(tt(ficbumpi)-35, fvel(ficbumpi)-(fvel(ficshifti-round(normtime*fps))*nrm), 'k')
                    plot([shift-35 shift-35], [-8 8], '--k')
                    plot([shift+2-35 shift+2-35], [-8 8], '--k')
                    ylim([-8 8])
                    ylabel(['Fwd Vel (mm/s)'])
                    xlabel(['Time (s)'])

                    %collect summary data
                    shiftdiff_pos = [shiftdiff_pos, nanmean(spos_uw(shifti+2:shifti+2+round(2*fs))) - nanmean(spos_uw(shifti-round(2*fs):shifti))];
                    shiftdiff_head = [shiftdiff_head, nanmean(uwheading(ficshifti+2:ficshifti+2+round(2*fps))) - nanmean(uwheading(ficshifti-round(2*fps):ficshifti))];

%                     shiftdiff_pos = [shiftdiff_pos, nanmean(spos_uw(shifti+2:shifti+2+round(2*fs))) - nanmean(spos_uw(shifti-round(2*fs):shifti))];
%                     shiftdiff_head = [shiftdiff_head, nanmean(uwheading(ficshifti+2:ficshifti+2+round(2*fps))) - nanmean(uwheading(ficshifti-round(2*fps):ficshifti))];

%                     shiftdiff_pos = [shiftdiff_pos, abs(nanmean(spos_uw(shifti+2:shifti+2+round(2*fs))) - nanmean(spos_uw(shifti-round(2*fs):shifti)))];
%                     shiftdiff_head = [shiftdiff_head, abs(nanmean(uwheading(ficshifti+2:ficshifti+2+round(2*fps))) - nanmean(uwheading(ficshifti-round(2*fps):ficshifti)))];


                    figure(200);
                    subplot(7,4,k+1); hold on
                    plot(tt(ficbumpi)-35, uwheading(ficbumpi)-(uwheading(ficshifti-round(normtime*fps))*nrm)./10, 'k')
                    plot(tt(ficbumpi(1:end-1))-35, diff(uwheading(ficbumpi)-(uwheading(ficshifti-round(normtime*fps))*nrm)), 'r')
                    if ismember(ii+1, [4, 6, 10, 12, 14, 15, 16, 17]) 
%                     if ismember(ii+1, [1, 2, 4, 6, 10, 12, 14, 15, 16, 17]) 
                        plot(9, 0.04, '*m')
                    end                    
                    xlim([-8 10])
%                     ylim([-.05 .05])
                    plot([0 0], ylim, '--k')
                    plot([2 2], ylim, '--k')                    

                    k = k + 1;
                    end
% keyboard

                    ii = ii + 1;
                end
            end
        end      
    end
end

figure(1);
subplot(1, 3, 1)
plot([35 35], [0 prevY1 + 10], '--k')
plot([37 37], [0 prevY1 + 10], '--k')
plot(tt,allo_pos*.1 + prevY1 + 15, 'k')
plot(tt,data(1).odor(nonan)*10 + prevY1 + 30, 'k')
ylim([0 prevY1+(42*180/pi)]) %304
xlim([27 45])
ylabel(['Bump Pos (columns)'])
xlabel(['Time (s)'])
title(['newpos, N = 5, n=' num2str(k)])

% subplot(1, 4, 2)
% plot([35 35], [0 prevY2], '--k')
% plot([37 37], [0 prevY2], '--k')
% plot(tt,data(1).odor(nonan)*10 + prevY2 + 10, 'k')
% % ylim([0 prevY + 30])
% title(['old pos, n=' num2str(ii)])

subplot(1, 3, 2)
plot([35 35], [0 prevY3 + 10], '--k')
plot([37 37], [0 prevY3 + 10], '--k')
plot(tt,allo_pos*.1 + prevY3 + 15, 'k')
plot(tt,data(1).odor(nonan)*10 + prevY3 + 30, 'k')
ylim([0 prevY3+42])
xlim([27 45])
ylabel(['Heading(deg)'])
xlabel(['Time (s)'])
title(['uwheading, N = 5, n=' num2str(k)])

subplot(1, 3, 3)
plot([35 35], [0 prevY4 + 10], '--k')
plot([37 37], [0 prevY4 + 10], '--k')
plot(tt,allo_pos*.1 + prevY4 + 15, 'k')
plot(tt,data(1).odor(nonan)*10 + prevY4 + 30, 'k')
ylim([0 prevY4+42])
xlim([27 45])
ylabel(['Fwd Vel (mm/s)'])
xlabel(['Time (s)'])
title(['fvel, N = 5, n=' num2str(k)]) 

if print_figs
            print_fig(1, gcf, [subpath 'allwindshifts_pos-head-fvel', add2txt], 'eps', 1); 
%     exportgraphics(gcf, [subpath 'fly' num2str(use_flies(f)) '_dff_trialsplits_plumetrials', add2txt, '.png'])
end

%
figure(3);clf;hold on
scatter(ones(1, length(shiftdiff_pos)) + (randn(1,length(shiftdiff_pos))/20), shiftdiff_pos, 15, 'k', 'filled')
scatter(ones(1, length(shiftdiff_head)) + 1 + (randn(1,length(shiftdiff_head))/20), shiftdiff_head, 15, 'k', 'filled')
[~, p] = ttest2(shiftdiff_pos, shiftdiff_head);
[~, p1] = ttest(shiftdiff_pos);
[~, p2] = ttest(shiftdiff_head);
y = [];
y = [shiftdiff_pos', shiftdiff_head'];
boxplot(y)
% sigstar([1 2], p)
xlim([0 3])
title(['m1=', num2str(mean(shiftdiff_pos), '%.3f') ', p1=', num2str(p1, '%.3f'),...
    ', m2=', num2str(mean(shiftdiff_head), '%.3f') , ', p2=', num2str(p2, '%.3f')])





