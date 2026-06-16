%% all fly
%%
clear all
% close all

% use_flies = [21:27, 30:39]; %vt062617
% use_flies = [21, 22, 24:26, 30, 31, 33:35, 37]; %vt062617 no nans
% use_flies = [41:51]; %52g12
% use_flies = [41, 42, 44:47, 51]; %52g12 no nans
use_flies = 25; %Figure 1d
use_flies = [22, 33]; %Figure 4d, vt062617
use_flies = [45, 51]; %Figure 4d, 52g12
% use_flies = 51;

cylm = [0 .6]; %1f
cylm = [0 .8; 0 0.8]; %vt062617 4d
cylm = [0 3; 0 4]; %52g12 4d



% use_flies = [21:27, 30:39]; %vt062617
% use_flies = [21, 22, 24:26, 30, 31, 33:35, 37]; %vt062617 no nans
% use_flies = [41:51]; %52g12
% use_flies = [41, 42, 44:47, 51]; %52g12 no nans


line_label = 'paper';
% track_met = 'pmove_odor'; metlim = [0 .9];
track_met = 'uwv_odor'; metlim = [0 10]; zscale = 1;
% track_met = 'fvel_odor'; metlim = [2 8];
% track_met = 'bmpmet'; metlim = [0 1];
% % % cylm = repmat([0.1, 2], length(use_flies), 1);

plot_heatmaps = 1;
plot_reg = 0;
plot_trialseries = 0;
plot_nullcomp = 0;
plot_corr = 0;

pulses = 0;

print_figs = 0;
add2txt = 'paper';

fig_no = 1;

%% load data and start figs
headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames
subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/exploratory figs/dffVbehave metric comparisons/';
subpath = '/Users/kathmn01/Desktop/';

if plot_heatmaps
    figure(fig_no); clf;
    set(gcf, 'Position', [1 1 1420 150]);
end
if plot_trialseries
    figure(fig_no + 1); clf;
    set(gcf, 'Position', [1 1 1420 850]);
end
if plot_reg
    figure(fig_no + 2); clf;
    set(gcf, 'Position', [1 1 900 850]);
end
if plot_nullcomp
    figure(fig_no + 3); clf;
    set(gcf, 'Position', [1 1 550 650]); xlim([0.5, 2.5]); ylim([-1 1])
end
if plot_corr
    figure(fig_no + 5); clf;
    set(gcf, 'Position', [1 1 1450 850]); xlim([0.5, 2.5]); ylim([-1 1])
end
%% collect data

i = 1;
all_wd = [];
all_bhvauto = [];
allbmp = []; alltrackm = [];
all_ccm = nan(18, 1);
all_nullccm = nan(18,1);
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;

    %find use_trials and alltrial_dff data to set bump threshold
    use_trials = [];
    for t = 1:length(data)
        if ~pulses
            nonan = ~isnan(data(t).bump_amp(:,1));
    
            %find stim times before trial sorting bc used to rule out plume trials
            odori = data(t).odor(nonan);
            idyl = odori >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odori(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ...
                    ~isempty(data(t).bump_amp) && length(odoroni) == 1 && max(odori) > 0 && odori(end) == 0%cwoo
    %         if length(odoroni) == 1
    %             odoroffi = idy(odori(idy+1) < 0.5); %all pts where next pt is subthresh (aka downward crossing)
                odoroffi = idy(odori(idy(idy < length(odori))+1)<0.5);
                if length(odoroffi) == 1 
                    use_trials = [use_trials, t];
                end
            end
        else
            use_trials = [use_trials, t];
        end
    end

    hm = [];
    trackm = [];
    bmp = [];
    for t = use_trials
        nonan = ~isnan(data(t).bump_amp(:,1));
        odor = data(t).odor(nonan);




        if pulses && length(odor) > 4500

            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)

            odor = zeros(length(odor), 1);
%             odor(2485:4500) = 1;
            odor(odoroni(1):odoroni(1)+1500) = 1;
        end

        %concat activity metric
        dff = [];
        dff= data(t).bump_amp(:,1);
%         dff = (dff-min(dff)); %baseline shifted
%             dff = (dff-min(dff))./max(dff); %normalized
        hm = [hm, dff];
        bmpmet = quantile(dff, .95);
        bmp = [bmp, quantile(dff, .95)];
%             bmp = [bmp, nanmean(dff)];
%             bmp = [bmp, nanmean(dff(find(odor)))];
%             bmp = [bmp, quantile(dff(find(odor)), .75)];


        %concat behavior metrics
        if data(t).closed_loop
            heading = data(t).calc_windpos(nonan);
        else
            heading = data(t).calc_heading(nonan);
        end
        heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 200)));
        speed = data(t).calc_speed(nonan);
        fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
        x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
        y = cumsum(speed.*cos(heading*pi/180))/fps;
        tt = data(t).calc_ts(nonan);
%         strait = data(t).straight(nonan);
        avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 500);
        fvel = abs(data(t).calc_deltapitch(nonan))*9.52/2;
%             fvel = smoothdata(abs(data(t).calc_deltapitch(nonan))*9.52/2, 'movmean', 20);
        uwv = cos(heading*pi/180).*speed;
        pmove = zeros(length(fvel), 1);
        pmove(fvel > 2) = 1;
        goal = mean(heading(find(odor)));
        goalv = cos((heading-goal)*pi/180).*speed;



        fvel_odor = fvel(find(odor));
        uwv_odor = uwv(find(odor));
        pmove_odor = pmove(find(odor));
        goalv_odor = goalv(find(odor));

% if f == 5
%     keyboard
% end

        met = eval(track_met);
%         if nanmean(met) < 50
%             trackm = [trackm, nanmean(met)]; %mean behavior for this trial
%         else
%             trackm = [trackm, nan]; %mean behavior for this trial
%         end

        trackm = [trackm, nanmean(met)]; %mean behavior for this trial
    end
%% Corr

    ccm = []; nullccm = [];

    x = bmp(~isnan(trackm)); %bump of trials with a behavior value
    y = trackm(~isnan(trackm)); %behavior of trials with behavior value
    ccm = corrcoef(x, y);

    randx = x(:,randperm(length(x)));
    randy = y(:,randperm(length(y)));
    nullccm = corrcoef(x, randy);

    if length(x) > 5
        all_ccm(i) = ccm(1,2);
        all_nullccm(i) = nullccm(1,2);
    end    

    normbmp = (bmp-min(bmp))/(max(bmp)-min(bmp)); %normalize within one fly
    normtrackm = (trackm-min(trackm))/(max(trackm)-min(trackm));

% if f == 5
%     keyboard
% end

%% heatmap

    if plot_heatmaps
        if ~isempty(hm)
%             figure(fig_no); subplot(3,4,i); hold on
            figure(fig_no); subplot(1,2,i); hold on
%             figure(fig_no); subplot(1,7,i); hold on
            title(['fly ' num2str(use_flies(f)) ', R=' num2str(ccm(1,2), '%.2f')])
    
            dff_ts = data(1).dff_ts;
            ft_ts = data(1).calc_ts + 3;
%             newnonan = ~isnan(ft_ts);
            
            % Build the trackm block (repeated 5 times) to prepend to heatmap
            if ~isempty(trackm)
                % Make sure trackm is a column vector
                trackm_col = trackm(:)*zscale*cylm(f,2)/metlim(2);  
                % Repeat each row across 5 columns
                trackm_block = repmat(trackm_col, 1, 500);  
                % Concatenate with hm (rows must match)
                combined_hm = [trackm_block, hm(320:end, :)'];
            
                % Create corresponding x values: 5 fake columns for trackm + ft_ts
                combined_x = [(ft_ts(321)-200*.01:.01:ft_ts(321)), ft_ts(320:end)'];
%                 combined_x = [(-8:-1), ft_ts(newnonan)'];
            else
                combined_hm = hm(newnonan,:);
                combined_x = ft_ts(newnonan);
            end
            
            % Plot the combined heatmap
            imagesc(combined_x, 1:size(combined_hm,1), combined_hm, cylm(f,:))
            colorbar
            
            % Axis limits
            xlim([-5 65])
            if size(hm,2) < 20
                ylim([0 20])
            end
% % % % %             dff_ts = data(1).dff_ts;
% % % % %             ft_ts = data(1).calc_ts + 3;
% % % % %             newnonan = ~isnan(ft_ts);
% % % % %             imagesc(ft_ts(newnonan), 1:size(hm,2), hm(newnonan,:)', cylm(f,:))
% % % % % %             imagesc(ft_ts(newnonan), 1:size(hm,2), hm(newnonan,:)', [0 max(max(hm(newnonan, :)))])
% % % % %             colorbar
% % % % %     
% % % % %     
% % % % %             xlim([0 65])
% % % % %             if size(hm,2) < 20
% % % % %                 ylim([0 20])
% % % % %             end
% % % % %     
% % % % %     %         trackm = []; trackm = bmp;
% % % % %             if ~isempty(trackm)
% % % % %                 h1 = gca;
% % % % %                 h2 = axes;
% % % % %                 subplot(4,5,i, h2); hold on
% % % % %         %         figure(1)
% % % % %                 xlim([0 65])
% % % % %                 if size(hm,2) < 20
% % % % %                     ylim([0 20])
% % % % %                 end
% % % % %                 
% % % % %                 imagesc(1:5, 1:length(trackm), repmat(trackm, 5, 1)', metlim)
% % % % %                 set(gca,'visible','off')
% % % % % 
% % % % %     %             colorbar('west')
% % % % %             end
% % % % % %                 plot([28 28], [0 30], '--k', 'LineW', 2)
% % % % % %                 plot([43 43], [0 30], '--k', 'LineW', 2)

                xlim([0 65])
                if size(hm,2) < 20
                    ylim([0 20])
                end
    
        end
        sgtitle([line_label ', ' track_met], 'interp', 'none')

    end


%% trial series
    if plot_trialseries
        figure(fig_no + 1); subplot(5,5,i); hold on
        title(['fly ' num2str(use_flies(f)) ', R=' num2str(ccm(1,2), '%.2f')])
    
        plot(1:length(bmp), normbmp, 'g')
        plot(1:length(bmp), normtrackm, 'k')
        scatter(1:length(bmp), normbmp, 5, 'g', 'filled')
        scatter(1:length(bmp), normtrackm, 5, 'k', 'filled')
        sgtitle([line_label ', ' track_met], 'interp', 'none')
    end


%% regression plots
    if plot_reg
        figure(fig_no + 2); subplot(5,5,i); hold on
    
        plot(normbmp, normtrackm , '.k', 'MarkerSize', 10)
    %     plot(x, y , '.k', 'MarkerSize', 10)
        xlabel('bump Q4'); ylabel(['mean ' track_met], 'interp', 'none')
        title(['fly ' num2str(use_flies(f)) ', R=' num2str(ccm(1,2), '%.2f')])
        xlim([-.2 1.2])
        ylim([-.2 1.2])
    %     xlim([0 .8]); 
    %     ylim([-6 12])  
        sgtitle([line_label ', ' track_met], 'interp', 'none')
    end

%% auto and x corr
    if plot_corr
        figure(fig_no + 5); subplot(4,5,i); hold on
        
        if ~isempty(normbmp) & ~isempty(normtrackm)
%             [bmpauto,l] = xcorr(normbmp-mean(normbmp));
            [behvauto, ll] = xcorr(normtrackm-mean(normtrackm));
            pad = zeros(1,81);
            mid = round(length(ll)/2);
            pad(41-(mid-1):41+(mid-1)) = behvauto;
            all_bhvauto = [all_bhvauto; pad];
            smbhv = smoothdata(behvauto, 'movmean', 5);
            ll = ll/100;

            [mx,mxi] = max(smbhv);
            [mn,mni] = min(smbhv);
            
            thr = mn+(mx-mn)/2;
            over = smbhv >= thr; %all superthresh pts
            overi = find(over); %superthresh indices
            overi = overi(overi>1); %all subthresh after 1
            oncrossi = overi(smbhv(overi-1)<thr);
            offcrossi = overi(smbhv(overi(overi < length(smbhv))+1)<thr); %all pts where next pt is subthresh (aka downward crossing)
%             offcrossi = overi(smbhv(overi+1)<thr); %all pts where next pt is subthresh (aka downward crossing)

            wdi = [];
            if ~isempty(oncrossi) && any(ll(oncrossi) < 0)
                ii = find(ll(oncrossi) < 0);
                wdi(1) = max(oncrossi(ii))-1; 
            end
            if ~isempty(offcrossi) && any(ll(offcrossi) > 0)
                ii = find(ll(offcrossi) > 0);
                wdi(2) = min(offcrossi(ii))+1; 
            end
            if ~isempty(wdi) && ~any(wdi == 0)
%                 plot(ll(wdi), smbhv(wdi), '*k'); 
                all_wd = [all_wd, diff(ll(wdi))];
            end

%             plot(l, bmpauto, 'g')
            plot(ll, behvauto, 'k')
            plot(ll, smbhv, 'r')
            xlim([-.40 .40])
            sgtitle(['Autocorr_', line_label ', ' track_met], 'interp', 'none')

%             plot(ll(oncrossi-1), smbhv(oncrossi-1), '*k')
%             plot(ll(offcrossi+1), smbhv(offcrossi+1), '*k')

%             plot(ll(mxi), smbhv(mxi), '*g')
%             plot(ll(mni), smbhv(mni), '*g')
            plot([ll(1) ll(end)], [thr thr], '--k')
%             [ht, cnt, wd] = findpeaks(smbhv,ll,'Annotate','extents');
%             figure(100); clf; findpeaks(smbhv, ll, 'Annotate', 'extents','MinPeakDistance',6, 'MinPeakProminence',.1);
%             ii = round(length(wd)/2)+1;
%             all_wd = [all_wd, wd(ii)];
            title(['fly ' num2str(use_flies(f)) ', wd=' num2str(diff(ll(wdi)))])
            print_fig(1, gcf, [subpath line_label '_splitflies_autocorrs_' track_met add2txt], 'eps', 1); 
        end

%         figure(fig_no + 6); subplot(5,5,i); hold on
%         title(['fly ' num2str(use_flies(f))])
%         if ~isempty(normbmp) & ~isempty(normtrackm)
%             bmpbehvxcorr = xcorr(normbmp, normtrackm);
%             plot(bmpbehvxcorr, 'k')
%         end
%         sgtitle(['XCorr_', line_label ', ' track_met], 'interp', 'none')
    end
    i = i + 1;
end

%% shuffle comp
    if plot_nullcomp
        figure(fig_no + 3); hold on
        for f = 1:size(all_ccm,1)    
            plot([1 2], [all_ccm(f) all_nullccm(f)], '.k', 'MarkerS', 20)
            plot([1 2], [all_ccm(f) all_nullccm(f)], 'k')
        end
    
        xlim([.5 2.5])
        ylabel(['CC of resp and ' track_met], 'interp', 'none')
        xlabel('1=CC, 2 = shuffleCC')
    
        [h,p] = ttest2(all_ccm, all_nullccm);
        title([line_label, ', BehavMet=' track_met ', p=', num2str(p,'%.3f')], 'interp', 'none')
    end

%%
plot_allfly_autocorr = 0;
if plot_allfly_autocorr
    figure(fig_no + 100); clf; hold on
    scatter(ones(1, length(all_wd))+randn(1, length(all_wd))/30, all_wd, 50, 'k')
    boxplot(all_wd)
    ylim([0 .20]);
    xlim([0 2])
    ylabel('peak half-height width (s)')
    title([line_label '_' track_met], 'Interpreter', 'none') 
    print_fig(1, gcf, [subpath line_label '_allflies_autocorr-peakwidth_' track_met add2txt], 'eps', 1); 
    
    figure(fig_no + 101); clf; hold on;
    lt =  -.4:.01:.4;
    m = mean(all_bhvauto);
    se = std(all_bhvauto)./sqrt(size(all_bhvauto,1));
    st = std(all_bhvauto);
    us = m + se;
    ls = m - se;
    fill([lt, lt(end:-1:1)], [us, ls(end:-1:1)], 'k', 'FaceAlpha',0.1)
    plot(lt, m, 'k', 'LineW', 4)
    ylim([-.5 2])
    xlabel('time lag (s)')
    title([line_label '_' track_met], 'Interpreter', 'none')
    print_fig(1, gcf, [subpath line_label '_allflies_meanautocorr_' track_met add2txt], 'eps', 1); 
end

if print_figs
    if plot_heatmaps
        figure(fig_no)
        print_fig(1, gcf, [subpath line_label '_allflies_heatmap_' track_met add2txt], 'eps', 1); 
    %     exportgraphics(gcf, [subpath 'hACK_allflies_dffheatmap_uwv_cwoo.png'])
        exportgraphics(gcf, [subpath line_label '_allflies_heatmap_' track_met add2txt '.png'])
    end
    if plot_reg
        figure(fig_no + 2)
        print_fig(1, gcf, [subpath line_label '_allflies_reg_' track_met add2txt], 'eps', 1); 
    %     exportgraphics(gcf, [subpath 'hACK_allflies_dffheatmap_uwv_cwoo.png'])
        exportgraphics(gcf, [subpath line_label '_allflies_reg_' track_met add2txt '.png'])
    end
    if plot_nullcomp
        figure(fig_no + 3)
        print_fig(1, gcf, [subpath line_label '_allflies_nullcomp_' track_met add2txt], 'eps', 1); 
    %     exportgraphics(gcf, [subpath 'hACK_allflies_dffheatmap_uwv_cwoo.png'])
        exportgraphics(gcf, [subpath line_label '_allflies_nullcomp_' track_met add2txt '.png'])
    end
end






