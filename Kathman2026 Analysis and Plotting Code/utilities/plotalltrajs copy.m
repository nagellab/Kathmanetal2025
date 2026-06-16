%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames

%%
% subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2023 paper/files/fig6/';
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/exploratory figs/traj splits/';
use_flies = [21:27, 30:40]; sdm = 0.5;%[21:40, 91:93];%[21:26, 33, 37]; %[21:26,30, 32:35, 37:39]; %21:39;%[6,8:10];%1:17; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging
use_flies = [25, 26, 21, 23]; sdm = 0.5;
% use_flies = 41:51; sdm = 0; %sdm not used in colornorm
% use_flies = 103; sdm = 0.5;
% use_flies = [81, 84:87];
% use_flies = [81];
% % use_flies = 151:156;

pulses = 0;

overlay_trls = 1;
print_figs = 0;
add2txt = '_cwoo_caxis0d5';
%%

j = 1; %counts subplot (flies or trials depending on flag)
sp = [7,10];
sp = [8,3]; %[5, 4];
sp = [7,5];
sp = [5, 6];
% sp = [1,1];
if overlay_trls
    figure(1); clf;
    set(gcf, 'Position', [1 1 1250 850]); %[1 1 1450 850]); %[100 100 500 1000])
    sp = [4,5]; %[3,4 ]; %4,5
end


strait_cond = cell(1);
for f = 1:length(use_flies)
    if ~overlay_trls
        figure(f); clf
%         set(gcf, 'Position', [1 1 750 1950]); %[100 100 500 1000])
        set(gcf, 'Position', [1 1 1550 350]); %[100 100 500 1000])
        j = 1;
    else
%         subplot(1, 5, f)
        subplot(1, 6, f)
    end


    load(filename{use_flies(f)})
    data = syncstruct;

    jj = 1; %counts trials per fly
    ons = []; offs = [];


    allamp_thr = []; use_trials = [];
    for t = 1:length(data)
        nonan = ~isnan(data(t).bump_amp(:,1));

        %find stim times before trial sorting bc used to rule out plume trials
        odor = data(t).odor(nonan);
        idyl = odor >= 0.5; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        oncrossi = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)

        allamp(t) = quantile(data(t).bump_amp(~isnan(data(t).bump_amp(:,1))), .90); 
        allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 && max(odor) > 0%cwoo (closed loop, odor on, wind on, not plume, not bad/empty trial)
        if data(t).closed_loop% && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
            use_trials = [use_trials, t];
        end
    end
%     ylm = [-.05 nanmean(allamp)*1.65];
    ylm = [0 nanmax(allamp)];
    ylm = [0 2.52];
%     sdm = .35;
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;

%     if length(use_trials) > sp(1)*sp(2)
%         use_trials = use_trials(1:(sp(1)*sp(2)));
%     end
    for t = use_trials

        nonan = ~isnan(data(t).bump_amp(:,1));
        odor = data(t).odor(nonan);


        ampthresh = .5;%0.125;
        idyl = odor >= ampthresh; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        oncrossi = idy(odor(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)

        if pulses
                odor = zeros(length(odor), 1);
%             odor(2485:4500) = 1;
            odor(oncrossi(1):oncrossi(1)+1500) = 1;
        end

%         if data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %cwoo          %closed loop, odor on, wind on, not plume, not bad/empty trial
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
%         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp) %owoo
        if ~isempty(data(t).bump_amp)% && length(oncrossi) > 0 && length(oncrossi) < 3%cwoo

            dff = data(t).bump_amp(nonan, 1);
            if data(t).closed_loop
                heading = data(t).calc_windpos(nonan);
            else
                heading = data(t).calc_heading(nonan);
            end
            speed = data(t).calc_speed(nonan);
            fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
            x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
            y = cumsum(speed.*cos(heading*pi/180))/fps;
            tt = data(t).calc_ts(nonan);
%             strait = data(t).straight(nonan, 1);
%             odor = data(t).odor(nonan);
% 
%             if pulses
%                     odor = zeros(length(odor), 1);
%     %             odor(2485:4500) = 1;
%                 odor(odoroni(1):odoroni(1)+1500) = 1;
%             end


            wind = data(t).wind(nonan);
            fvel = smoothdata(abs(data(t).calc_deltapitch(nonan))*9.52/2, 'movmean', 100);
            uwv = cos(heading*pi/180).*speed;
            
            if max(odor) > 1
                goal = mean(heading(find(odor)));
            else
                goal = mean(heading);
            end
            goalv = cos((heading-goal)*pi/180).*speed;
            
            hold on
            stm = odor + wind; %make stim color code
            scatter(x, y, 6, [1; -stm(2:end)], 'filled') %plot traj colored by stim
            
% % %             subplot(sp(1), sp(2),j); 
% % %             hold on            
% % %             z = smoothdata(dff, 'movmean', 200);
% % %             scatter(x, y, 5, z, 'filled')
% % %             caxis(ylm)
% % % %             caxis([0 max(z)])
% % %             colormap jet
% % %             plot(x(1), y(1), 'rx', 'LineW', 1.1,'MarkerSize',4)
% % %             o = find(odor > 0.5);
% % %             w = find(wind > 0.5);
% % %             plot(x(o(1)), y(o(1)), 'kx', 'LineW', 1.2 ,'MarkerSize',2)
% % %             plot(x(o(end)), y(o(end)), 'kx', 'LineW', 1.2 ,'MarkerSize',2)
% % %             plot(x(w(1)), y(w(1)), 'kx', 'LineW', .1 ,'MarkerSize',2)
% % %             plot(x(w(end)), y(w(end)), 'kx', 'LineW', .1 ,'MarkerSize',2)
% % %             colorbar

            xlim([-425 425])
            ylim([-325 825])
%             scatter(x, y, 8, [1.5; odor(2:end)], 'filled')
            
            if ~overlay_trls
                if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %cwoo
                    fill([-300 -280 -280 -300], [-200 -200 400 400], 'k')
                elseif data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
                    fill([-300 -280 -280 -300], [-200 -200 400 400], 'g')
                elseif ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp) %owoo
                    fill([-300 -280 -280 -300], [-200 -200 400 400], 'r')
                end

                set(gca,'Yticklabel',[]) 
                set(gca,'Xticklabel',[])
                title(['t' num2str(t) ',uwv_odr=' num2str(mean(uwv(find(odor))), '%.1f'), ',clr=', num2str(ylm(2), '%.1f')], 'Interpreter', 'none')
%                 title(['t' num2str(t) ',uwv=' num2str(mean(uwv(find(odor))), '%.1f'), ',glv=' num2str(mean(goalv(find(odor))), '%.1f')])
                j = j + 1;
            end
            jj = jj + 1;

        end
    end
    sgtitle(['Fly ' num2str(use_flies(f))])
    if overlay_trls
        set(gca,'Yticklabel',[]) 
        set(gca,'Xticklabel',[])
        title(['fly #' num2str(use_flies(f)) ', n=' num2str(jj)])
        j = j + 1;
    else
        if print_figs
            print_fig(1, gcf, [subpath 'alltrial_traj_colornorm_fly' num2str(use_flies(f)) add2txt], 'eps', 1); 
%             exportgraphics(gcf, [subpath 'alltrial_traj_colornorm_fly' num2str(use_flies(f)), '.png'])
        end
    end

end
% if print_figs
%     print_fig(1, gcf, [subpath 'allflies_trajoverlay_cwoo_hACK'], 'eps', 1); 
%     exportgraphics(gcf, [subpath 'allflies_trajoverlay_cwoo_hACK', '.png'])
% end

% fly_splits = 1;
% figure(2); clf; hold on
% j = 1;
% sp = [4,5];
% for f = 1:length(use_flies)
%     if ~isempty(strait_cond{f})
%         if fly_splits
%             subplot(sp(1), sp(2), f); hold on
%         %     subplot(4,3,f); hold on
%             for i = 1:3
%                 if ~isempty(strait_cond{f})
%                 sz(j) = length(strait_cond{f}(:,i));
%                 scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, strait_cond{f}(:,i), 10, 'b', 'filled')
%                 end
%             end
%             [~,ps(1)] = ttest2(strait_cond{f}(:,1), strait_cond{f}(:,2));
%             [~,ps(2)] = ttest2(strait_cond{f}(:,2), strait_cond{f}(:,3));
%             [~,ps(3)] = ttest2(strait_cond{f}(:,1), strait_cond{f}(:,3));
% %             sigstar({[1,2], [2,3], [1,3]}, ps)
%             plot(1:3, nanmedian(strait_cond{f}), 'k', 'LineW', 2)
%             ylabel('Straightness index')
%             title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j))])
%             ylim([0.5 1.2]); xlim([0.5 3.5])
%         else
%             plot(1:3, nanmedian(strait_cond{f}), 'k', 'LineW', 2)
%             ylabel('Straightness index')
%             title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j))])
%             ylim([0.82 1]); xlim([0.5 3.5])
%         end
%         j = j + 1;
%     end
% end




