%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames_all

%%

% subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2023 paper/files/fig6/';
subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/';
subpath = '/Users/kathmn01/Desktop/';

use_flies{1} = [1,2,3,4,5];
use_flies{2} = [41:46, 49:51]; %47, 48 errors for some reason
use_flies{2} = [41:51]; %52g12
% use_flies{3} = [21:27, 30, 32:39];%[32, 38] onset/offset only flies; %[21:26, 33, 35:37, 39];%[21, 22, 24:26, 33, 35:37, 39]; %[21:26, 33, 37];%, 35:39]; %[21:26, 33, 36:38];%; %; %[6,8:10];%1:17; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging
% % use_flies{3} = [21:27 32:39];
% use_flies{3} = [21:26,30, 32:35, 37:39];
% use_flies{3} = 25;
% % use_flies{1} = [81, 83:86];
% % use_flies{1} = [141:146];
% % use_flies{1} = [81, 84:87];
% % % use_flies{1} = [66, 67, 68, 71:74];%ss00020
% % % use_flies{1} = 151:156;
% % use_flies{1} = [81, 84:87, 151:156]; 

add2txt{1} = '_21D07';
add2txt{2} = '_51G12_cwoo_all';
% add2txt{3} = '_62617_cwoo_nofly30bc10sstim(little resp)';
% add2txt{1} = '_hksplit_both';
ylm{1} = [.35 1.5];
ylm{2} = [.1 3.2];
ylm{3} = [.1 .8];

print_figs = 0;
plot_trials = 0;
pulses = 0;

fig_no = 100;
for l = 2; %3%1:length(use_flies)
%     close all
    allwint = [];
    allmet = [];
    %%%%% add filter to exclude flies with no bump at all or no (good) cwoo trials, currently doing manually %%%%%
    good_trl = [];
    h = figure(fig_no); clf; hold on
    % set(gcf, 'Position', [1 1 1420 850]); %[100 100 500 1000])
    moffs = []; mons = []; all_stim = [];
    
    
    
    
    
    
    
    for f = 1:length(use_flies{l})
%         if ismember(f, 6:11)
%             pulses = 1;
%         else
%             pulses = 0;
%         end
        pulses = 0;

        load(filename{use_flies{l}(f)})
        data = syncstruct;
        for s = 1%:2 %stim [odor, wind]
        jj = 1;
        ons = []; offs = [];
%         if use_flies{l}(f) == 22
%             use_trials = [1:26, 30:length(data)];
%         elseif use_flies{l}(f) == 37
%             use_trials = 8:length(data);
%         else
            use_trials = 1:length(data);
%         end
        for t = use_trials
            win1 = []; wint1 = [];
            if ~isempty(data(t).bump_amp)
                nonan = ~isnan(data(t).bump_amp(:,1));
                if s == 1
                    stimi = data(t).odor(nonan);
        
                    odor = data(t).odor(nonan);
                    if pulses && length(odor) > 4500
                        idyl = odor >= 0.5; %all superthresh pts
                        idy = find(idyl); %superthresh indices
                        idy = idy(idy>1); %all subthresh after 1
                        odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            
                        stimi = zeros(length(odor), 1);
%                         stimi(2485:4500) = 1;
                        stimi(odoroni(1):odoroni(1)+1500) = 1;
                    end
    
                else
                    stimi = data(t).wind(nonan);
                end
                ampthresh = .5;%0.125;
                idyl = stimi >= ampthresh; %all superthresh pts
                idy = find(idyl); %superthresh indices
                idy = idy(idy>1); %all subthresh after 1
                oncrossi = idy(stimi(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
        %             ontimes = tt(oncrossi);
%     if f == 8
%         keyboard; end
%                 if data(t).closed_loop && ~isempty(data(t).odor_on) && data(t).odor_on && data(t).wind_on ...
%                         && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %&& length(oncrossi) < 2 %cwoo
        %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
        %         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %owoo
%                 if ~isempty(data(t).bump_amp) && ~data(t).odor_on %odor off
                if data(t).closed_loop && ~isempty(data(t).odor_on) && data(t).odor_on && data(t).wind_on ...
                        && ~isempty(data(t).bump_amp)
                    if (~pulses && length(oncrossi) < 2) || pulses
                    good_trl = [good_trl, t+100*f];
    
                    dff = data(t).bump_amp(nonan, 1);
                    heading = data(t).calc_windpos(nonan);
                    speed = data(t).calc_speed(nonan);
                    fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
                    x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
                    y = cumsum(speed.*cos(heading*pi/180))/fps;
                    tt = data(t).calc_ts(nonan);
    %                 strait = data(t).straight(nonan, 1);
                    pos = smoothdata(data(t).bump_col(nonan, 1), 'movmean', 200);
    
    %                 oncrossi = 2470; %use this for blanks so has a time stamp to center window at
    
                    %trig data
                    if s == 1
                        win = [12, 27]; %viewing window, not whole trial
                    else
                        win = [5, 15]; 
                    end
                    met = dff;
                    for c = 1:length(oncrossi)
                        win1 = (oncrossi(c)-win(1)*fps):(oncrossi(c)+win(2)*fps);
                        win1 = win1(win1 < length(tt));
                        wint1 = tt(win1)-tt(oncrossi(c));
                        
                        if ~isempty(ons) && length(win1) ~= size(ons,1) %align to stim start
                            this_met = met(win1);
                            this_on = nan(size(ons,1),1);
                            this_on(end-length(this_met)+1:end) = this_met;
                            ons = [ons, this_on];
                        else
                            plot_wint = wint1;
                            ons = [ons, met(win1)];
                        end
                        if plot_trials
                            figure(f+2); hold on; 
                            subplot(2,1,1); hold on; plot(wint1, stimi(win1), 'k')
        %                     title(['fly#' num2str(use_flies{l}(f)) add2txt{l}], 'interpreter', 'none')
                            subplot(2,1,2); hold on; plot(wint1, met(win1), 'k')
    
                        end
                    end
                    jj = jj + 1;
                    end
                end
            end
        end
        title(['fly#' num2str(use_flies{l}(f)) ', n=' num2str(jj-1) ', ' add2txt{l}], 'interpreter', 'none')
        size(ons,2)
        if ~isempty(ons)
            moffs = [moffs, mean(offs,2)];
            mons = [mons, mean(ons,2)];
            
            if ~isempty(win1)
                figure(fig_no); hold on
                subplot(2,1,1); hold on; plot(wint1, stimi(win1), 'k')
            end
            figure(fig_no); hold on
            subplot(2,1,2); hold on; plot(plot_wint, smoothdata(mean(ons,2), 10), 'k'); 
            plot(wint1, stimi(win1), 'k')
            if plot_trials
                figure(f+2); hold on
                subplot(2,1,2); hold on; plot(plot_wint, smoothdata(mean(ons,2), 10), 'r', 'LineW', 2.5)
    %             exportgraphics(gcf, [subpath 'fly' num2str(use_flies{l}(f)) '_cwoo_mndffodorresp' add2txt{l}, '.' export_type])
%                 if f == 2
%                     ylim([0 1.75])
%                 end

                if f == 3 || f == 5
                    ylim([0 4.75])
                end
            end

            allwint = [allwint; plot_wint'];
            allmet = [allmet; smoothdata(mean(ons,2), 10)'];
        end
        end
    end
    figure(fig_no); hold on
    stimtype{1} = 'odor'; stimtype{2} = 'wind';
    subplot(2,1,1); ylabel(stimtype{s})
    subplot(2,1,2); hold on; plot(plot_wint, mean(mons,2), 'r', 'LineW', 2)
    ylim(ylm{l})
%     ylim([.1 1.5])
%     if l == 1 && f == 5
%         keyboard; end
    title(['N=' num2str(length(use_flies{l})) add2txt{l}], 'interpreter', 'none')


    if print_figs
        exportgraphics(gcf, [subpath 'odortrig' add2txt{l}, '.png'])
        print_fig(1, gcf, [subpath 'odortrig' add2txt{l}], 'eps', 1);
    end
    fig_no = fig_no + 1;
%     exportgraphics(gcf, [subpath 'odortrig' add2txt{l} '_N10_missing45bcerror', '.png'])
%     print_fig(1, gcf, [subpath 'odortrig' add2txt{l} '_N10_missing45bcerror'], 'eps', 1);

    %plot heatmap
    figure(50); clf;  hold on
    imagesc(allwint(1,:), 1:size(allmet,1), allmet, [0.1 1.5])
    plot([0 0], [0.5 size(allmet,1)+.5], '--k')
    plot([15 15], [0.5 size(allmet,1)+.5], '--k')
    colorbar
%     xlim([0 65])
    if print_figs
        exportgraphics(gcf, [subpath 'odortrig_heatmap' add2txt{l}, '.png'])
        print_fig(1, gcf, [subpath 'odortrig_heatmap' add2txt{l}], 'eps', 1);
    end

end

    



