clear all
close all

%% Inputs
% use_flies = [21:27, 30:39]; sdm = 0.5; line = 'hAck';%for some reason cutting out 40?
% use_flies = 41:51; sdm = 0; line = 'FC1'; %sdm not used in colornorm
% use_flies = [21:27, 30:39, 41:51]; line = 'alllines';
use_flies = [21:27, 30:39]; line = '62617';
% use_flies = [41:51]; line = '52g12';
% use_flies = 25;

print_figs = 0;
add2txt = ['_' line '_cwoo'];

%% dir and filenames
headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2023 paper/sandbox/exploratory figs/';
cd(headpath)

allflies_filenames_all

%% Calc
all_uwv = [];
all_strt = [];
all_fvel = [];
all_speed = [];

N = 0;
figure(1); set(gcf, 'Position', [1 1 1300 1200]); %[100 100 500 1000])
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;



%% Pad straightness field to align by frame number if needed
% This fixes trials where data(t).straightness is stored only at sampled frames
% instead of being NaN-padded to full frame length.

mxfrms = 0;

% Find longest trial by frame number
for t = 1:length(data)
    if isfield(data, 'frame') && ~isempty(data(t).frame)
        mxfrms = max(mxfrms, max(data(t).frame));
    end
end

for t = 1:length(data)

    if isfield(data, 'straight') && isfield(data, 'frame') && ...
            ~isempty(data(t).straight) && ~isempty(data(t).frame)

        straightness_raw = data(t).straight;
        frame_raw = data(t).frame;

        % Force column vectors
        straightness_raw = straightness_raw(:);
        frame_raw = frame_raw(:);

        % Only pad if straightness matches the unpadded frame vector
        if length(straightness_raw) == length(frame_raw)

            padded_straightness = nan(mxfrms, 1);
            padded_straightness(frame_raw) = straightness_raw;

            data(t).straight = padded_straightness;

        end
    end
end




    muwv = [];
    mstrt = [];
    mfvel = []; 
    mspeed = [];

    %find use_trials and alltrial_dff data to set bump threshold
    use_trials = [];
    for t = 1:length(data)
        nonan = ~isnan(data(t).bump_amp(:,1));

        %find stim times before trial sorting bc used to rule out plume trials
        odori = data(t).odor(nonan);
        idyl = odori >= 0.5; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        odoroni = idy(odori(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
        if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(odoroni) == 1 && max(odori) > 0 && odori(end) == 0%cwoo
            odoroffi = idy(odori(idy+1) < 0.5); %all pts where next pt is subthresh (aka downward crossing)
            if length(odoroffi) == 1 
                if ~isempty(data(t).straight)
                use_trials = [use_trials, t];
                end
            end
        end
    end


    if ~isempty(use_trials) && length(use_trials) > 5
        for t = use_trials
            nonan = ~isnan(data(t).bump_amp(:,1));
            
            %define metrics
            dff = data(t).bump_amp(nonan, 1);
            speed = data(t).calc_speed(nonan);
            fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
            if data(t).closed_loop
                heading = data(t).calc_windpos(nonan);
            else
                heading = data(t).calc_heading(nonan);
            end
            heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 200)));
            tt = data(t).calc_ts(nonan);
            avel = data(t).calc_deltaz(nonan);
%             avel = smoothdata(avel, 'movmean', 500);
            fvel = data(t).calc_deltapitch(nonan)*9.52/2;
%             fvel = smoothdata(fvel, 'movmean', 30); 
            uwv = cos(heading*pi/180).*speed;
            x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
            y = cumsum(speed.*cos(heading*pi/180))/fps;
            strait = data(t).straight(nonan);

            %find stim times
            odor = data(t).odor(nonan);
            wind = data(t).wind(nonan);
            odori = find(odor)+10;
            windi = find(wind)+10;

%             if windi(1) > 100 && windi(end) < length(wind) %if wind isn't on whole time
%                 basei = [1:windi(1)];
%                 stimi = odori(1:length(basei));
%             end
            basei = 1:1000;
            stimi = odori(1):(odori(1)+1000);
            if nanmean(fvel) < 200 && nanmean(strait(basei)) > .5 && nanmean(strait(stimi)) > .5
                muwv = [muwv; [mean(uwv(basei)), mean(uwv(stimi))]];
                mstrt = [mstrt; [nanmean(strait(basei)), nanmean(strait(stimi))]];                
                mfvel = [mfvel; [nanmean(fvel(basei)), nanmean(fvel(stimi))]]; 
                mspeed = [mspeed; [nanmean(speed(basei)), nanmean(speed(stimi))]]; 
            end
        end
    end
    
    if ~isempty(muwv)
        all_uwv = [all_uwv; mean(muwv)];
        all_strt = [all_strt; mean(mstrt)];
        all_fvel = [all_fvel; mean(mfvel)];
        all_speed = [all_speed; mean(mspeed)];
        figure(1);
        if use_flies(f) < 40 %hAck
            subplot(1,4,1); hold on;
            plot([1 2], mean(muwv), 'k', 'LineWidth', 1)
            plot([1 2], mean(muwv), '.k', 'MarkerSize', 30)

            subplot(1,4,2); hold on;
            plot([1 2], mean(mstrt), 'k', 'LineWidth', 1)
            plot([1 2], mean(mstrt), '.k', 'MarkerSize', 30)

            subplot(1,4,3); hold on;
            plot([1 2], mean(mfvel), 'k', 'LineWidth', 1)
            plot([1 2], mean(mfvel), '.k', 'MarkerSize', 30)

            subplot(1,4,4); hold on;
            plot([1 2], mean(mspeed), 'k', 'LineWidth', 1)
            plot([1 2], mean(mspeed), '.k', 'MarkerSize', 30)
        else %FC1
            plot([1 2], mean(muwv), 'r', 'LineWidth', 1)
            plot([1 2], mean(muwv), '.r', 'MarkerSize', 30)
        end

        N = N + 1;
    end            


end

subplot(1,4,1); hold on;
plot([.6 1.4], [nanmean(all_uwv(:,1)) nanmean(all_uwv(:,1))], 'r', 'LineW', 4)
plot([1.6 2.4], [nanmean(all_uwv(:,2)) nanmean(all_uwv(:,2))], 'r', 'LineW', 4)
[~,p] = ttest2(all_uwv(:,1), all_uwv(:,2));
sigstar({[1,2]}, p)

title([line ' allflies, if red=FC1, N=', num2str(N), ', p=' num2str(p, '%.7f')])
ylabel('uwVel'); xlabel('base, odor')
xlim([0 3]); %ylim([])


subplot(1,4,2); hold on;
plot([.6 1.4], [nanmean(all_strt(:,1)) nanmean(all_strt(:,1))], 'r', 'LineW', 4)
plot([1.6 2.4], [nanmean(all_strt(:,2)) nanmean(all_strt(:,2))], 'r', 'LineW', 4)
[~,p] = ttest2(all_strt(:,1), all_strt(:,2));
sigstar({[1,2]}, p)

title([line ' allflies, if red=FC1, N=', num2str(N), ', p=' num2str(p, '%.7f')])
ylabel('straightness'); xlabel('base, odor')
xlim([0 3]); %ylim([])

subplot(1,4,3); hold on;
plot([.6 1.4], [nanmean(all_fvel(:,1)) nanmean(all_fvel(:,1))], 'r', 'LineW', 4)
plot([1.6 2.4], [nanmean(all_fvel(:,2)) nanmean(all_fvel(:,2))], 'r', 'LineW', 4)
[~,p] = ttest2(all_fvel(:,1), all_fvel(:,2));
sigstar({[1,2]}, p)

title([line ' allflies, if red=FC1, N=', num2str(N), ', p=' num2str(p, '%.7f')])
ylabel('fvel'); xlabel('base, odor')
xlim([0 3]); %ylim([])

subplot(1,4,4); hold on;
plot([.6 1.4], [nanmean(all_speed(:,1)) nanmean(all_speed(:,1))], 'r', 'LineW', 4)
plot([1.6 2.4], [nanmean(all_speed(:,2)) nanmean(all_speed(:,2))], 'r', 'LineW', 4)
[~,p] = ttest2(all_speed(:,1), all_speed(:,2));
sigstar({[1,2]}, p)

title([line ' allflies, if red=FC1, N=', num2str(N), ', p=' num2str(p, '%.7f')])
ylabel('speed'); xlabel('base, odor')
xlim([0 3]); %ylim([])



if print_figs
    print_fig(1, gcf, [subpath 'uwVel_baseVodor', add2txt], 'eps', 1); 
    exportgraphics(gcf, [subpath 'uwVel_baseVodor', add2txt, '.png'])
end


%% old (check why diff, gave slightly diff line graph, but doesn't run to check)
% %% all fly
% %%
% clear all
% close all
% 
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
% cd(headpath)
% 
% allflies_filenames
% 
% %%
% use_flies = [21:27, 30, 32:39]; sdm = 0.5; %N=15 (16 here, but 1 has no cwoo trials
% % use_flies = [41:51]; sdm = 1; %N = 10?
% 
% print_figs = 0;
% line_label = 'HACK';
% add2filename = '_new_test';
% 
% 
% subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/exploratory figs/';
% fig_no = 1;
% 
% uwv_cond = cell(1, length(use_flies));
% all_uwv = [];
% for f = 1:length(use_flies)
%     load(filename{use_flies(f)})
%     data = syncstruct;
% 
% 
%     for s = 1%:2 %stim [odor, wind]
%     jj = 1;
% 
%     muwv = []; mbase = [];
%     for t = 1:length(data)
% 
% 
%         nonan = ~isnan(data(t).bump_amp(:,1));
%         odor = data(t).odor(nonan);
%         ampthresh = .5;%0.125;
%         idyl = odor >= ampthresh; %all superthresh pts
%         idy = find(idyl); %superthresh indices
%         idy = idy(idy>1); %all subthresh after 1
%         oncrossi = idy(odor(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
% %             ontimes = tt(oncrossi);
% 
%         
% %         if data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %all non-plumes
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %cwoo
% %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
% %         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp)  && length(oncrossi) < 2 %owoo
% %         if ~isempty(data(t).bump_amp) %closed loop, odor on, wind on, not plume, not bad/empty trial
% 
% 
% %         if ~isempty(data(t).bump_amp)% && length(oncrossi) < 2
% 
%             dff = []; speed = []; x = []; y = []; tt = []; uwv = []; strait = [];
%             
%             dff = data(t).bump_amp(nonan, 1);
%             if data(t).closed_loop
%                 heading = smoothdata(data(t).calc_windpos(nonan), 'movmean', 100);
%             else
%                 heading = smoothdata(data(t).calc_heading(nonan), 'movmean', 100);
%             end
%             
% %             dff = data(t).bump_amp(nonan);
% %             heading = data(t).calc_windpos(nonan);
%             speed = data(t).calc_speed(nonan);
%             fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%             x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
%             y = cumsum(speed.*cos(heading*pi/180))/fps;
%             tt = data(t).calc_ts(nonan);
% %             strait = data(t).straight(nonan, 1).^3;
%             uwv = cos(heading*pi/180).*speed;
%             avel = abs(data(t).calc_deltaz(nonan));
%             fvel = data(t).calc_deltapitch(nonan)*9.52/2;
% 
% 
%             %straighness index
%             pth = []; strait = [];
%             for ii = 1:length(x)-1
%                 pth(ii) = sqrt((x(ii+1)-x(ii))^2 + (y(ii+1)-y(ii))^2);
%             end
%             pth = [pth, pth(end)];
% 
%             i = 1;
%             win = 150; %window size
%             for j = 1:length(x)-win
%                 strait(i) = sqrt((x(j+win)-x(j))^2 + (y(j+win)-y(j))^2)/sum(pth(j:j+win));
%                 strait(i) = strait(i)^2;
%                 i = i + 1;
%             end
%             straight = nan(length(tt),1);
%             straiti = round(win/2):(round(win/2)+length(strait)-1);
%             straight(straiti) = strait;
%             strait = []; strait = straight;
% 
% 
% 
%             wi = data(t).wind(nonan);
%             wii = find(wi);
%             if s == 1
%                 stimi = data(t).odor(nonan);
%             else
%                 stimi = data(t).wind(nonan);
%             end
%             stimii = find(stimi);
%             baseii = stimi(1:wii(1));
%             if ~isempty(baseii) && length(baseii) > 1
%                 muwv = [muwv mean(uwv(stimii))];
%                 mbase = [mbase mean(uwv(baseii))];
%             end
%         end
%     end
%     end
%     if ~isempty(muwv) && ~isempty(mbase)
%         all_uwv = [all_uwv; [mean(muwv) mean(mbase)]];
% %         plot()
%     end
% end
% 
% 
% % uwVel
% fly_splits = 0;
% figure(fig_no); clf; hold on
% % set(gcf, 'Position', [1 1 1420 850]);
% j = 1;
% sp = [4,4]; muwv = [];
% for f = 1:length(use_flies)
%     if ~isempty(uwv_cond{f})
%         if fly_splits
%             subplot(sp(1), sp(2), f); hold on
%         %     subplot(4,3,f); hold on
%             for i = 1:4
%                 if ~isempty(uwv_cond{f})
%                 sz(j) = length(uwv_cond{f}(:,i));
%                 scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, uwv_cond{f}(:,i), 10, 'k', 'filled')
%                 end
%             end
%             ps = []; 
%             [~,ps(1)] = ttest2(uwv_cond{f}(:,2), uwv_cond{f}(:,3));
%             [~,ps(2)] = ttest2(uwv_cond{f}(:,3), uwv_cond{f}(:,4));
%             sigstar({[2,3], [3,4]}, ps)
% 
%             plot(1:4, nanmedian(uwv_cond{f}), 'b', 'LineW', 2)
%             ylabel('uwVel')
%             title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j))])
%             ylim([-5 20]); xlim([0.5 4.5])
%         else
% 
%             muwv(1:2,j) = nanmedian(uwv_cond{f}(:,1:2));
%             muwv(3:4,j) = nanmedian(uwv_cond{f}(:,3:4));
%             muwv(5:6,j) = nanmedian(uwv_cond{f}(:,5:6));
% 
%             plot(1:2, nanmedian(uwv_cond{f}(:,1:2)), 'k', 'LineW', 2)
%             plot(3:4, nanmedian(uwv_cond{f}(:,3:4)), 'k', 'LineW', 2)
%             plot(5:6, nanmedian(uwv_cond{f}(:,5:6)), 'k', 'LineW', 2)
%             scatter([1:6], nanmedian(uwv_cond{f}(:,1:6)), 50, 'k', 'filled')
% %             scatter([1:6]+randn(1)/30, nanmedian(uwv_cond{f}(:,1:6)), 50, 'k', 'filled')
% 
%             ylabel('uwVel')
%             title(['N=' num2str(j)])
% %             ylim([0.82 1]); 
%             xlim([.5 6.5])
%         end
%         j = j + 1;
%     end
% end
% % ylim([0 10.5])
% ps = []; 
% [~,ps(1)] = ttest2(muwv(1,:), muwv(2,:));
% [~,ps(2)] = ttest2(muwv(3,:), muwv(4,:));
% [~,ps(3)] = ttest2(muwv(5,:), muwv(6,:));
% sigstar({[1,2], [3,4], [5,6]}, ps)
% 
% % exportgraphics(gcf, [subpath 'allflies_bumpcond_uwvel_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
% if print_figs
%     exportgraphics(gcf, [subpath 'bumpcond_uwvel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename, '.png'])
%     print_fig(1, gcf, [subpath 'bumpcond_uwvel_cwoo_', line_label, '_N' num2str(length(use_flies)) '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
% end
% 
% 
% 
