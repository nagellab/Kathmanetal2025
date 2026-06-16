%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames

%%
% use_flies = [21:27, 30, 32:39]; sdm = 0.5; %N=15 (16 here, but 1 has no cwoo trials
use_flies = [21:26, 30, 32:35, 37:39]; sdm = 0.5; %N=14 
% use_flies = 25;
% use_flies = [41:51]; sdm = 1; %N = 10?
% use_flies = [73]; sdm = .5; %[66, 67, 73, 74]
% use_flies = [73]; sdm = .5; %[66, 67, 73, 74]

% use_flies = [142:146]; sdm = .25; %splits

print_figs = 0;
line_label = 'fcsplit';
add2filename = '_test';

ttype = nan(100,20);

subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/exploratory figs/bump no bump/';
fig_no = 30;

for ss = sdm %[0 0.5 1 2]
    sdm = ss;

    strait_cond = cell(1, length(use_flies));
    uwv_cond = cell(1, length(use_flies));
    avel_cond = cell(1, length(use_flies));
    fvel_cond = cell(1, length(use_flies));
    muwv = [];
    allgdtrl = 0;
    for f = 1:length(use_flies)
        load(filename{use_flies(f)})
        data = syncstruct;

        vals = arrayfun(@(x) x.closed_loop, data)'; % extract field values into vector
        ttype(1:length(vals), f) = vals;
        numOnes = sum(vals == 1);
        numZeros = sum(vals == 0);
% 
        fprintf('%d, %d, %d, %.2f\n', f, numOnes, numZeros, numZeros/numOnes);
    
        allamp_thr = [];
        for i = 1:length(data)
            allamp(i) = quantile(data(i).bump_amp(~isnan(data(i).bump_amp(:,1))), .95); 
            allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
    %         if any(data(i).bump_amp(~isnan(data(i).bump_amp(:,1))) < -2)
    %             keyboard
    %         end
        end
    %     ylm = [-.05 nanmean(allamp)*1.5];
        base_thr = nanmean(allamp_thr)+nanstd(allamp_thr)*sdm;
        
    
    
        for s = 1%:2 %stim [odor, wind]
        jj = 1;
        ons = []; offs = [];
    %     figure(f); clf;
    %     set(gcf, 'Position', [1 1 1420 850]); %[100 100 500 1000])
        for t = 1:length(data)
    %         clearvars -except headpath filename data f t jj thresh use_flies strait_cond uwv_cond s muwv
    
    
            nonan = ~isnan(data(t).bump_amp(:,1));
            odor = data(t).odor(nonan);
            ampthresh = .5;%0.125;
            idyl = odor >= ampthresh; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(odor(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
    %             ontimes = tt(oncrossi);
    
            
    %         if data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %all non-plumes
            if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %cwoo
    %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
    %         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp)  && length(oncrossi) < 2 %owoo
    %         if ~isempty(data(t).bump_amp) %closed loop, odor on, wind on, not plume, not bad/empty trial
                if t < 400
                allgdtrl = allgdtrl + 1;
    %         if ~isempty(data(t).bump_amp)% && length(oncrossi) < 2
    
                dff = []; speed = []; x = []; y = []; tt = []; uwv = []; strait = [];
                
                dff = data(t).bump_amp(nonan, 1);
                if data(t).closed_loop
                    heading = smoothdata(data(t).calc_windpos(nonan), 'movmean', 100);
                else
                    heading = smoothdata(data(t).calc_heading(nonan), 'movmean', 100);
                end
                
    %             dff = data(t).bump_amp(nonan);
    %             heading = data(t).calc_windpos(nonan);
                speed = data(t).calc_speed(nonan);
                fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
                x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
                y = cumsum(speed.*cos(heading*pi/180))/fps;
                tt = data(t).calc_ts(nonan);
%                 strait = data(t).straight(nonan, 1).^3;
                uwv = cos(heading*pi/180).*speed;
                avel = abs(data(t).calc_deltaz(nonan));
                fvel = data(t).calc_deltapitch(nonan)*9.52/2;
    
    
                %straighness index
                pth = []; strait = [];
                for ii = 1:length(x)-1
                    pth(ii) = sqrt((x(ii+1)-x(ii))^2 + (y(ii+1)-y(ii))^2);
                end
                pth = [pth, pth(end)];
                
                % pth = sqrt(x.^2 + y.^2);
                i = 1;
                win = 150; %window size
                for j = 1:length(x)-win
                    strait(i) = sqrt((x(j+win)-x(j))^2 + (y(j+win)-y(j))^2)/sum(pth(j:j+win));
                    strait(i) = strait(i)^2;
                    i = i + 1;
                end
                straight = nan(length(tt),1);
                straiti = round(win/2):(round(win/2)+length(strait)-1);
                straight(straiti) = strait;
                strait = []; strait = straight;
    
    
    
    %             thresh = mean(dff)+var(dff)*5;
    %             if thresh < 0.2
    %                 thresh = 0.2; end
    %             if var(dff) < 0.0025
    %                 thresh = 0.3; end
    %             if thresh > 0.4
    %                 thresh = 0.4; end
    % 
    % 
    % %             straiti = data(t).straight(nonan);
    % %             if isnan(thresh(use_flies(f)))
    % %                 thresh(use_flies(f)) = 0.275; end
    %             bumpi = smoothdata(dff, 'movmean', 100) > thresh;
    
%                 bumpi = bumpfinder2(dff, tt, 3, 1, base_thr);
                bumpi = bumpfinder2(dff, tt, 3, 3, base_thr);
    
    %             figure(fig_no); hold on;
    %             plot(tt, dff, 'k')
    %             scatter(tt(find(bumpi)), dff(find(bumpi)), 'r', 'filled')
    %             plot([0, 60], [base_thr, base_thr], '--k')
    %             ylim([0, 2])
    %             fig_no = fig_no + 1;
    %             base_thr
    
                if s == 1
                    stimi = data(t).odor(nonan);
                else
                    stimi = data(t).wind(nonan);
                end
                stimii = find(stimi);
                muwv = [muwv mean(uwv(stimii))];
                if ~isempty(stimii)% & mean(uwv(stimii))>2.5 %filter trials by uwvel during odor
    
                    pre_bumpi = []; prenobumpi = []; odor_bumpi = []; odor_nobumpi = []; post_bumpi = []; post_nobumpi = [];
                    pre_bumpi = find(bumpi & ~stimi & tt < tt(stimii(1)));
                    pre_nobumpi = find(~bumpi & ~stimi & tt < tt(stimii(1)));
                    odor_bumpi = find(bumpi & stimi);
                    odor_nobumpi = find(~bumpi & stimi);
                    post_bumpi = find(bumpi & ~stimi & tt > tt(stimii(end)));
                    post_nobumpi = find(~bumpi & ~stimi & tt > tt(stimii(end)));
    
                    %straightness
%                     %%get rid of this conditional, but as is, only thing that changes by adding it is one fly (f24) and only at 1 std threshold (for some reason)
%                     if nanmean(strait(pre_bumpi)) > 0.5 %fix this, for some reason many of these are very low (unusually low) in F24, prob something about being close to beginning of trial?
%                         strait_cond{f}(jj,1) = nanmean(strait(pre_bumpi)); 
%                     else
%                         strait_cond{f}(jj,1) = nan;
%                     end
%     
%                     if nanmean(strait(pre_nobumpi)) > 0.5 %fix this, for some reason many of these are very low in F24, prob someething about being close to beginning of trial?
%                         strait_cond{f}(jj,2) = nanmean(strait(pre_nobumpi)); 
%                     else
%                         strait_cond{f}(jj,2) = nan;
%                     end

                    strait_cond{f}(jj,1) = nanmean(strait(pre_bumpi)); 
                    strait_cond{f}(jj,2) = nanmean(strait(pre_nobumpi));

                    strait_cond{f}(jj,3) = nanmean(strait(odor_bumpi));
                    strait_cond{f}(jj,4) = nanmean(strait(odor_nobumpi));
                    strait_cond{f}(jj,5) = nanmean(strait(post_bumpi));
                    strait_cond{f}(jj,6) = nanmean(strait(post_nobumpi));
    
                    %uwVel
                    uwv_cond{f}(jj,1) = nanmean(uwv(pre_bumpi));
                    uwv_cond{f}(jj,2) = nanmean(uwv(pre_nobumpi));
                    uwv_cond{f}(jj,3) = nanmean(uwv(odor_bumpi));
                    uwv_cond{f}(jj,4) = nanmean(uwv(odor_nobumpi)); 
                    uwv_cond{f}(jj,5) = nanmean(uwv(post_bumpi));
                    uwv_cond{f}(jj,6) = nanmean(uwv(post_nobumpi));  
    
                    %aVel
                    avel_cond{f}(jj,1) = nanmean(avel(pre_bumpi));
                    avel_cond{f}(jj,2) = nanmean(avel(pre_nobumpi));
                    avel_cond{f}(jj,3) = nanmean(avel(odor_bumpi));
                    avel_cond{f}(jj,4) = nanmean(avel(odor_nobumpi)); 
                    avel_cond{f}(jj,5) = nanmean(avel(post_bumpi));
                    avel_cond{f}(jj,6) = nanmean(avel(post_nobumpi));  
    
    
                    %fVel
                    fvel_cond{f}(jj,1) = nanmean(fvel(pre_bumpi));
                    fvel_cond{f}(jj,2) = nanmean(fvel(pre_nobumpi));
                    fvel_cond{f}(jj,3) = nanmean(fvel(odor_bumpi));
                    fvel_cond{f}(jj,4) = nanmean(fvel(odor_nobumpi)); 
                    fvel_cond{f}(jj,5) = nanmean(fvel(post_bumpi));
                    fvel_cond{f}(jj,6) = nanmean(fvel(post_nobumpi));  
    
        %             figure(1); hold on; j = 1;          
        %             for i = 1:4
        %                 if ~isempty(uwv_cond{f})
        %                 sz(j) = length(uwv_cond{f}(jj,i));
        %                 scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, uwv_cond{f}(:,i), 10, 'k', 'filled')
        %                 end
        %             end
        
                    jj = jj + 1;
                end
                end
            end
        end
        end
    
    %     subpath = '/Users/kathmn01/Desktop/all syncstructs/';
    % %     print_fig(1, gcf, [subpath 'all_trials_fly' num2str(use_flies(f))], 'png', 1); 
    %     exportgraphics(gcf, [subpath 'all_trials_fly' num2str(use_flies(f)), '.png'])
        
    end
    
    % % % %normalized distribution of uwv
    % % % figure(10); clf; hold on
    % % % hist(muwv_cwoo, 100)
    % % % figure(11); clf; hold on
    % % % hist(muwv_owoo, 100)
    % % % % 
    % % % % muwv_owoo = muwv;
    % % % N = [];
    % % % figure(100); clf; hold on
    % % % bnw = 0.5;
    % % % edg = -10:bnw:15;
    % % % N = histcounts(muwv_cwoo, edg, 'Normalization', 'probability'); 
    % % % N = smoothdata(N, 'movmean', 10);
    % % % plot(edg(2:end)-bnw/2, N, 'k', 'LineW', 2)
    % % % 
    % % % N = [];
    % % % N = histcounts(muwv_owoo, edg, 'Normalization', 'probability'); 
    % % % N = smoothdata(N, 'movmean', 10);
    % % % plot(edg(2:end)-bnw/2,N, 'r', 'LineW', 2)
    % % % 
    % % %     subpath = '/Users/kathmn01/Desktop/new figs/';
    % % % %     print_fig(1, gcf, [subpath 'all_trials_fly' num2str(use_flies(f))], 'png', 1); 
    % % %     exportgraphics(gcf, [subpath 'allfly_alltrl_uwvdist_cwoo-owoo_black-red' num2str(use_flies(f)), '.png'])
    
    % uwVel
    fly_splits = 0;
    figure(fig_no); clf; hold on
    % set(gcf, 'Position', [1 1 1420 850]);
    j = 1;
    sp = [4,4]; muwv = [];
    for f = 1:length(use_flies)
        if ~isempty(uwv_cond{f})
            if fly_splits
                subplot(sp(1), sp(2), f); hold on
            %     subplot(4,3,f); hold on
                for i = 1:4
                    if ~isempty(uwv_cond{f})
                    sz(j) = length(uwv_cond{f}(:,i));
                    scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, uwv_cond{f}(:,i), 10, 'k', 'filled')
                    end
                end
                ps = []; 
                [~,ps(1)] = ttest2(uwv_cond{f}(:,2), uwv_cond{f}(:,3));
                [~,ps(2)] = ttest2(uwv_cond{f}(:,3), uwv_cond{f}(:,4));
                sigstar({[2,3], [3,4]}, ps)
    
                plot(1:4, nanmedian(uwv_cond{f}), 'b', 'LineW', 2)
                ylabel('uwVel')
                title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j)), ', p='  num2str(ps)])
                ylim([-5 20]); xlim([0.5 4.5])
            else
    
                muwv(1:2,j) = nanmedian(uwv_cond{f}(:,1:2));
                muwv(3:4,j) = nanmedian(uwv_cond{f}(:,3:4));
                muwv(5:6,j) = nanmedian(uwv_cond{f}(:,5:6));
    
                plot(1:2, nanmedian(uwv_cond{f}(:,1:2)), 'k', 'LineW', 2)
                plot(3:4, nanmedian(uwv_cond{f}(:,3:4)), 'k', 'LineW', 2)
                plot(5:6, nanmedian(uwv_cond{f}(:,5:6)), 'k', 'LineW', 2)
                scatter([1:6], nanmedian(uwv_cond{f}(:,1:6)), 50, 'k', 'filled')
    %             scatter([1:6]+randn(1)/30, nanmedian(uwv_cond{f}(:,1:6)), 50, 'k', 'filled')
    
                ylabel('uwVel')
                title(['N=' num2str(j)])
    %             ylim([0.82 1]); 
                xlim([.5 6.5])
            end
            j = j + 1;
        end
    end
    % ylim([0 10.5])
    ps = []; 
    [~,ps(1)] = ttest2(muwv(1,:), muwv(2,:));
    [~,ps(2)] = ttest2(muwv(3,:), muwv(4,:));
    [~,ps(3)] = ttest2(muwv(5,:), muwv(6,:));
    sigstar({[1,2], [3,4], [5,6]}, ps)
    
    % exportgraphics(gcf, [subpath 'allflies_bumpcond_uwvel_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
    if print_figs
        exportgraphics(gcf, [subpath 'bumpcond_uwvel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename, '.png'])
        print_fig(1, gcf, [subpath 'bumpcond_uwvel_cwoo_', line_label, '_N' num2str(length(use_flies)) '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
    end
    
    %     exportgraphics(gcf, [subpath 'bumpcond_uwvel_cwoo_FC1_N' num2str(length(use_flies)), '_fly', num2str(use_flies(f)) '.png'])
    %     print_fig(1, gcf, [subpath 'bumpcond_uwvel_cwoo_FC1_N' num2str(length(use_flies)), '_fly', num2str(use_flies(f))], 'eps', 1); 
    
    %straightness
    fly_splits = 0;
    figure(fig_no+1); clf; hold on
    % set(gcf, 'Position', [1 1 1420 850]);
    j = 1;mstr = [];
    
    for f = 1:length(use_flies)
        if ~isempty(strait_cond{f})
            if fly_splits
                subplot(sp(1), sp(2), f); hold on
            %     subplot(4,3,f); hold on
                for i = 1:4
                    if ~isempty(strait_cond{f})
                    sz(j) = length(strait_cond{f}(:,i));
                    scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, strait_cond{f}(:,i), 10, 'k', 'filled')
                    end
                end
                ps = [];
                [~,ps(1)] = ttest2(strait_cond{f}(:,2), strait_cond{f}(:,3));
                [~,ps(2)] = ttest2(strait_cond{f}(:,3), strait_cond{f}(:,4));
                sigstar({[2,3], [3,4]}, ps)
    
                plot(1:4, nanmedian(strait_cond{f}), 'b', 'LineW', 2)
                ylabel('Straightness index')
                title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j)), ', p='  num2str(ps)])
                ylim([0.5 1.2]); xlim([0.5 4.5])
            else
                mstr(1:2,j) = nanmedian(strait_cond{f}(:,1:2));
                mstr(3:4,j) = nanmedian(strait_cond{f}(:,3:4));
                mstr(5:6,j) = nanmedian(strait_cond{f}(:,5:6));
    
                plot(1:2, nanmedian(strait_cond{f}(:,1:2)), 'k', 'LineW', 2)
                plot(3:4, nanmedian(strait_cond{f}(:,3:4)), 'k', 'LineW', 2)
                plot(5:6, nanmedian(strait_cond{f}(:,5:6)), 'k', 'LineW', 2)
                scatter([1:6], nanmedian(strait_cond{f}(:,1:6)), 50, 'k', 'filled')  
    
    
    
                ylabel('Straightness index')
                title(['N=' num2str(j), ', p='  num2str(ps)])
    %             ylim([0.82 1]); 
                xlim([.5 6.5])
            end
            j = j + 1;
        end
    end
    
    ps = []; 
    [~,ps(1)] = ttest2(mstr(1,:), mstr(2,:));
    [~,ps(2)] = ttest2(mstr(3,:), mstr(4,:));
    [~,ps(3)] = ttest2(mstr(5,:), mstr(6,:));
    sigstar({[1,2], [3,4], [5,6]}, ps)
    title(['N=' num2str(j), ', p='  num2str(ps)])
    
    % exportgraphics(gcf, [subpath 'allflies_bumpcond_strait_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
    if print_figs
        exportgraphics(gcf, [subpath 'bumpcond_strait_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename, '.png'])
        print_fig(1, gcf, [subpath 'bumpcond_strait_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
    end
    
    %     exportgraphics(gcf, [subpath 'bumpcond_strait_cwoo_FC1_N' num2str(length(use_flies)), '_fly', num2str(use_flies(f)) '.png'])
    %     print_fig(1, gcf, [subpath 'bumpcond_strait_cwoo_FC1_N' num2str(length(use_flies)) '_fly', num2str(use_flies(f))], 'eps', 1); 
    
    
    %avel
    fly_splits = 0;
    figure(fig_no+2); clf; hold on
    % set(gcf, 'Position', [1 1 1420 850]);
    j = 1;mstr = [];
    
    for f = 1:length(use_flies)
        if ~isempty(avel_cond{f})
            if fly_splits
                subplot(sp(1), sp(2), f); hold on
            %     subplot(4,3,f); hold on
                for i = 1:4
                    if ~isempty(avel_cond{f})
                    sz(j) = length(avel_cond{f}(:,i));
                    scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, avel_cond{f}(:,i), 10, 'k', 'filled')
                    end
                end
                ps = [];
                [~,ps(1)] = ttest2(avel_cond{f}(:,2), avel_cond{f}(:,3));
                [~,ps(2)] = ttest2(avel_cond{f}(:,3), avel_cond{f}(:,4));
                sigstar({[2,3], [3,4]}, ps)
    
                plot(1:4, nanmedian(avel_cond{f}), 'b', 'LineW', 2)
                ylabel('abs(aVel)')
                title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j)), ', p='  num2str(ps)])
                ylim([0.5 1.2]); xlim([0.5 4.5])
            else
                mstr(1:2,j) = nanmedian(avel_cond{f}(:,1:2));
                mstr(3:4,j) = nanmedian(avel_cond{f}(:,3:4));
                mstr(5:6,j) = nanmedian(avel_cond{f}(:,5:6));
    
                plot(1:2, nanmedian(avel_cond{f}(:,1:2)), 'k', 'LineW', 2)
                plot(3:4, nanmedian(avel_cond{f}(:,3:4)), 'k', 'LineW', 2)
                plot(5:6, nanmedian(avel_cond{f}(:,5:6)), 'k', 'LineW', 2)
                scatter([1:6], nanmedian(avel_cond{f}(:,1:6)), 50, 'k', 'filled')  
    
    
    
                ylabel('abs(aVel)')
                title(['N=' num2str(j)])
    %             ylim([0.82 1]); 
                xlim([.5 6.5])
            end
            j = j + 1;
        end
    end
    
    ps = []; 
    [~,ps(1)] = ttest2(mstr(1,:), mstr(2,:));
    [~,ps(2)] = ttest2(mstr(3,:), mstr(4,:));
    [~,ps(3)] = ttest2(mstr(5,:), mstr(6,:));
    sigstar({[1,2], [3,4], [5,6]}, ps)
    title(['N=' num2str(j), ', p='  num2str(ps)])
    
    % exportgraphics(gcf, [subpath 'allflies_bumpcond_strait_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
    if print_figs
        exportgraphics(gcf, [subpath 'bumpcond_avel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename '.png'])
        print_fig(1, gcf, [subpath 'bumpcond_avel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
    end
    
    %     exportgraphics(gcf, [subpath 'bumpcond_avel_cwoo_FC1_N' num2str(length(use_flies)), '_fly', num2str(use_flies(f)) '.png'])
    %     print_fig(1, gcf, [subpath 'bumpcond_avel_cwoo_FC1_N' num2str(length(use_flies)) '_fly', num2str(use_flies(f))], 'eps', 1); 
    
    
    %fvel
    fly_splits = 0;
    figure(fig_no+3); clf; hold on
    % set(gcf, 'Position', [1 1 1420 850]);
    j = 1;mstr = [];
    
    for f = 1:length(use_flies)
        if ~isempty(fvel_cond{f})
            if fly_splits
                subplot(sp(1), sp(2), f); hold on
            %     subplot(4,3,f); hold on
                for i = 1:4
                    if ~isempty(fvel_cond{f})
                    sz(j) = length(fvel_cond{f}(:,i));
                    scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, fvel_cond{f}(:,i), 10, 'k', 'filled')
                    end
                end
                ps = [];
                [~,ps(1)] = ttest2(fvel_cond{f}(:,2), fvel_cond{f}(:,3));
                [~,ps(2)] = ttest2(fvel_cond{f}(:,3), fvel_cond{f}(:,4));
                sigstar({[2,3], [3,4]}, ps)
    
                plot(1:4, nanmedian(fvel_cond{f}), 'b', 'LineW', 2)
                ylabel('fVel')
                title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j))])
                ylim([0.5 1.2]); xlim([0.5 4.5])
            else
                mstr(1:2,j) = nanmedian(fvel_cond{f}(:,1:2));
                mstr(3:4,j) = nanmedian(fvel_cond{f}(:,3:4));
                mstr(5:6,j) = nanmedian(fvel_cond{f}(:,5:6));
    
                plot(1:2, nanmedian(fvel_cond{f}(:,1:2)), 'k', 'LineW', 2)
                plot(3:4, nanmedian(fvel_cond{f}(:,3:4)), 'k', 'LineW', 2)
                plot(5:6, nanmedian(fvel_cond{f}(:,5:6)), 'k', 'LineW', 2)
                scatter([1:6], nanmedian(fvel_cond{f}(:,1:6)), 50, 'k', 'filled')  
    
    
    
                ylabel('fVel')
                title(['N=' num2str(j)])
    %             ylim([0.82 1]); 
                xlim([.5 6.5])
            end
            j = j + 1;
        end
    end
    
    ps = []; 
    [~,ps(1)] = ttest2(mstr(1,:), mstr(2,:));
    [~,ps(2)] = ttest2(mstr(3,:), mstr(4,:));
    [~,ps(3)] = ttest2(mstr(5,:), mstr(6,:));
    sigstar({[1,2], [3,4], [5,6]}, ps)
    title(['N=' num2str(j), ', p='  num2str(ps)])
    
    if print_figs
    % exportgraphics(gcf, [subpath 'allflies_bumpcond_strait_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
        exportgraphics(gcf, [subpath 'bumpcond_fvel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename '.png'])
        print_fig(1, gcf, [subpath 'bumpcond_fvel_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
    
    %     exportgraphics(gcf, [subpath 'bumpcond_avel_cwoo_FC1_N' num2str(length(use_flies)), '_fly', num2str(use_flies(f)) '.png'])
    %     print_fig(1, gcf, [subpath 'bumpcond_avel_cwoo_FC1_N' num2str(length(use_flies)) '_fly', num2str(use_flies(f))], 'eps', 1); 
    end
    fig_no = fig_no + 1;
end


%% tried to clean this code (below), but figures don't look the same. not sure why,  need to check!!!

% clear all
% close all
% 
% %% Inputs
% use_flies = [21:27, 30:39]; sdm = 0.5; line = 'hAck';%for some reason cutting out 40?
% % use_flies = 41:51; sdm = 0; line = 'FC1'; %sdm not used in colornorm
% % use_flies = [21:27, 30:39, 41:51]; line = 'alllines';
% 
% met = 'strait'; %'abs_avel';
% fly_splits = 0;
% 
% print_figs = 0;
% add2txt = ['_' line '_cwoo'];
% 
% %% dir and filenames
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
% subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2023 paper/sandbox/exploratory figs/';
% cd(headpath)
% 
% allflies_filenames
% 
% %% Calc
% fig_no = 10;
% for sdm = [0.5] %[0 0.5 1]
%     all_met = cell(1, length(use_flies));
%     for f = 1:length(use_flies)
%         load(filename{use_flies(f)})
%         data = syncstruct;
%     
%         jj = 1;
%     
%         %find use_trials and alltrial_dff data to set bump threshold
%         use_trials = [];
%         allamp_thr = [];
%         for t = 1:length(data)
%             nonan = ~isnan(data(t).bump_amp(:,1));
%     
%             %find stim times before trial sorting bc used to rule out plume trials
%             odori = data(t).odor(nonan);
%             idyl = odori >= 0.5; %all superthresh pts
%             idy = find(idyl); %superthresh indices
%             idy = idy(idy>1); %all subthresh after 1
%             odoroni = idy(odori(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
%             if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(odoroni) == 1 && max(odori) > 0 && odori(end) == 0%cwoo
%                 odoroffi = idy(odori(idy+1) < 0.5); %all pts where next pt is subthresh (aka downward crossing)
%                 if length(odoroffi) == 1 
%                     use_trials = [use_trials, t];
%                 end
%             end
%             allamp_thr = [allamp_thr; data(t).bump_amp(nonan, 1)];
%         end
%         base_thr = nanmean(allamp_thr)+nanstd(allamp_thr)*sdm;
%     
%     
%         if ~isempty(use_trials) && length(use_trials) > 5
%             for t = use_trials
%                 nonan = ~isnan(data(t).bump_amp(:,1));
%                 
%                 %define metrics
%                 dff = data(t).bump_amp(nonan, 1);
%                 speed = data(t).calc_speed(nonan);
%                 fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%                 if data(t).closed_loop
%                     heading = data(t).calc_windpos(nonan);
%                 else
%                     heading = data(t).calc_heading(nonan);
%                 end
%                 heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 200)));
%                 tt = data(t).calc_ts(nonan);
%                 avel = data(t).calc_deltaz(nonan);
%     %             avel = smoothdata(avel, 'movmean', 500);
%                 abs_avel = abs(avel);
%                 fvel = data(t).calc_deltapitch(nonan)*9.52/2;
%     %             fvel = smoothdata(fvel, 'movmean', 30); 
%                 uwv = cos(heading*pi/180).*speed;
%                 x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
%                 y = cumsum(speed.*cos(heading*pi/180))/fps;
%                 strait = data(t).straight(nonan);
%                 old_strait = strait;
% 
%                 %need to check why this is diff than bump_stats straigt met!!!!!!
%                 [strait] = calc_straightness(x, y, 150); %seems to be lower magnitude than bump_stat version, but maybe less errors
%                 
% 
% %                 % check straightness metric
% %                 figure;
% %                 set(gcf, 'Position', [1 1 1420 850]); %[100 100 500 1000])
% %                 subplot(1,2,1); scatter(x, y, 30, old_strait, 'filled'); title('old calc')
% %                 xlim([-200 200]); ylim([-200 250]); caxis([.75 1])
% %                 subplot(1,2,2); scatter(x, y, 30, strait, 'filled'); title('new calc')
% %                 xlim([-200 200]); ylim([-200 250]); caxis([.75 1])
%     
%                 %find stim times
%                 odor = data(t).odor(nonan);
%                 wind = data(t).wind(nonan);
%                 odori = find(odor)+10;
%                 windi = find(wind)+10;
%     
% 
%                 bump = bumpfinder2(dff, tt, 3, 2, base_thr);
%     
%                 pre_bumpi = []; prenobumpi = []; odor_bumpi = []; odor_nobumpi = []; post_bumpi = []; post_nobumpi = [];
%     
%                 %conditions (bump state, odor state, before/after odor)
%                 pre_bumpi = find(bump & ~odor & tt < tt(odori(1)));
%                 pre_nobumpi = find(~bump & ~odor & tt < tt(odori(1)));
%                 odor_bumpi = find(bump & odor);
%                 odor_nobumpi = find(~bump & odor);
%                 post_bumpi = find(bump & ~odor & tt > tt(odori(end)));
%                 post_nobumpi = find(~bump & ~odor & tt > tt(odori(end)));
%     
%     
%                 this_met = eval(met);
%                 all_met{f}(jj,1) = nanmean(this_met(pre_bumpi));
%                 all_met{f}(jj,2) = nanmean(this_met(pre_nobumpi));
%                 all_met{f}(jj,3) = nanmean(this_met(odor_bumpi));
%                 all_met{f}(jj,4) = nanmean(this_met(odor_nobumpi)); 
%                 all_met{f}(jj,5) = nanmean(this_met(post_bumpi));
%                 all_met{f}(jj,6) = nanmean(this_met(post_nobumpi)); 
%     
%                 jj = jj + 1;
%             end
%         end
%     end
%     
%     
%     figure(fig_no); clf; hold on
%     % set(gcf, 'Position', [1 1 1420 850]);
%     j = 1;
%     sp = [4,4]; mnmet = [];
%     for f = 1:length(use_flies)
%         if ~isempty(all_met{f})
%             if fly_splits
%                 subplot(sp(1), sp(2), f); hold on
%             %     subplot(4,3,f); hold on
%                 for i = 1:4
%                     if ~isempty(all_met{f})
%                     sz(j) = length(all_met{f}(:,i));
%                     scatter(ones(1,sz(j)).*i+randn(1,sz(j))/40, all_met{f}(:,i), 10, 'k', 'filled')
%                     end
%                 end
%                 ps = []; 
%                 [~,ps(1)] = ttest2(all_met{f}(:,2), all_met{f}(:,3));
%                 [~,ps(2)] = ttest2(all_met{f}(:,3), all_met{f}(:,4));
%                 sigstar({[2,3], [3,4]}, ps)
%     
%                 plot(1:4, nanmedian(all_met{f}), 'b', 'LineW', 2)
%                 ylabel(met, 'interp', 'none')
%                 title(['fly #' num2str(use_flies(f)) ', n=' num2str(sz(j))])
%                 ylim([-5 20]); xlim([0.5 4.5])
%             else
%     
%                 mnmet(1:2,j) = nanmedian(all_met{f}(:,1:2));
%                 mnmet(3:4,j) = nanmedian(all_met{f}(:,3:4));
%                 mnmet(5:6,j) = nanmedian(all_met{f}(:,5:6));
%     
%                 plot(1:2, nanmedian(all_met{f}(:,1:2)), 'k', 'LineW', 2)
%                 plot(3:4, nanmedian(all_met{f}(:,3:4)), 'k', 'LineW', 2)
%                 plot(5:6, nanmedian(all_met{f}(:,5:6)), 'k', 'LineW', 2)
%                 scatter([1:6], nanmedian(all_met{f}(:,1:6)), 50, 'k', 'filled')
%     %             scatter([1:6]+randn(1)/30, nanmedian(all_met{f}(:,1:6)), 50, 'k', 'filled')
%     
%                 ylabel(met, 'interp', 'none')
%                 title(['Line=' line ', sdm= ' num2str(sdm) ', N=' num2str(j)], 'interp', 'none')
%     %             ylim([0.82 1]); 
%                 xlim([.5 6.5])
%             end
%             j = j + 1;
%         end
%     end
%     % ylim([0 10.5])
%     ps = []; 
%     [~,ps(1)] = ttest2(mnmet(1,:), mnmet(2,:));
%     [~,ps(2)] = ttest2(mnmet(3,:), mnmet(4,:));
%     [~,ps(3)] = ttest2(mnmet(5,:), mnmet(6,:));
%     sigstar({[1,2], [3,4], [5,6]}, ps)
%     
%     % exportgraphics(gcf, [subpath 'allflies_bumpcond_uwvel_cwoo_62617_N' num2str(length(use_flies)), '.eps'])
%     if print_figs
%         exportgraphics(gcf, [subpath 'bumpcond_' met '_cwoo_', line_label, '_N' num2str(length(use_flies)), '_threshmeanand' num2str(sdm) 'std' add2filename, '.png'])
%         print_fig(1, gcf, [subpath 'bumpcond_' met '_cwoo_', line_label, '_N' num2str(length(use_flies)) '_threshmeanand' num2str(sdm) 'std' add2filename], 'eps', 1); 
%     end
%     fig_no = fig_no + 1;
% end

