
clear all
close all
plot_pooled_curves = 1;
print_figs = 0;
add2txt = 'hksplit_both_matfit';
afit = 0; %manual fitting method

use_flies = [27, 32, 36, 38, 109:112]; %plumes
use_flies = [21:27, 30:39]; %cwoo
use_flies = [21:26, 30, 32:35, 37:39]; %cwoo
use_flies = [21:26, 33:35, 37:39]; %cwoo, no 30
fig_flyid = 1:14;
% use_flies = 21;

% use_flies = [81, 84:87];
% use_flies = 151:156;
% use_flies = [81, 84:87, 151:156]; %hk splits

sdm = .5;
fig_no = 10;
%% Set up
headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2024 paper/sandbox/exploratory figs/';
subpath = '/Users/kathmn01/Desktop/';
cd(headpath)

allflies_filenames

all_durs = [];
all_durs_odor = [];
all_durs_offset = [];
all_durs_onset = [];

perFlyOffsets = cell(length(use_flies), 1);   % holds offset durs per fly
total_bumps_per_fly = zeros(length(use_flies),1);
%%
all_t = 0;
for f = 1:length(use_flies)

    load(filename{use_flies(f)})
    data = syncstruct;

    % Find trials 
%     use_trials = [];
%     for t = 1:length(data) % or use "use_trials"
% %         if ~isempty(data(t).bump_amp) && data(t).closed_loop
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && sum(diff(data(t).odor) == 1) == 1%cwoo
%             use_trials = [use_trials, t];
%         end
%     end


    use_trials = find(arrayfun(@(t) t.closed_loop && sum(diff(t.odor) == 1) == 1 && sum(diff(t.odor) == -1) == 1, data));% && all(isnan(t.plume_t)), data); %cwoo
%     use_trials = arrayfun(@(t) t.closed_loop && sum(diff(t.odor) == 1) > 2, data);% && all(isnan(t.plume_t)), data); %plume
    all_t = all_t + length(use_trials);

    % Find bump threshold from all trials
    all_dff = cat(1, data(arrayfun(@(t) isfield(t, 'bump_amp') && ~isempty(t.bump_amp) && any(~isnan(t.bump_amp)), data)).bump_amp);
    all_dff = all_dff(~isnan(all_dff));
    base_thr = nanmean(all_dff)+nanstd(all_dff)*sdm;


    for t = use_trials

        %extract trial values
        [tt, dff, odoroni, odoroffi, odor, fps] = extract_struct_data(data(t));
        oni = []; offi = []; 
        [bump, oni, offi] = bumpfinder2(dff, tt, 3, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration

% Count ALL detected bumps in this trial (before inclusion filtering)
        total_bumps_per_fly(f) = total_bumps_per_fly(f) + length(offi);

        for i = 1:length(oni)
            all_durs = [all_durs; (offi(i)-oni(i))/fps];
%             if oni(i) > odoroni && oni(i) < odoroffi
%                 all_durs_odor = [all_durs_odor; (offi(i)-oni(i))/100];
%                 all_durs_onset = [all_durs_onset; (offi(i)-odoroni)/100];
                if oni(i) < odoroffi && offi(i) > odoroffi && oni(i) > odoroni

                    this_dur = (offi(i)-oni(i))/fps;
                    this_onset = (offi(i)-odoroni)/fps;
                    this_offset = (offi(i)-odoroffi)/fps;

                    all_durs_odor = [all_durs_odor; this_dur];
                    all_durs_onset = [all_durs_onset; this_onset];
                    all_durs_offset = [all_durs_offset; this_offset];
                    
                    perFlyOffsets{f} = [perFlyOffsets{f}, this_offset];
                end
%             end
        end

    end
end


%% histograms
if plot_pooled_curves
% bump dist (pooled)
bnw = 2;
edg = 2:bnw:80;
cnt = edg(2:end)-bnw/2;
smth = 1;
ylm = [0 .3];

cmap = lines;
figure(fig_no); clf; hold on
set(gcf, 'Position', [1 1 1520 550]); %[100 100 500 1000])


%% all durs
subplot(1,3,1); hold on
N = histcounts(all_durs, edg, 'Normalization', 'probability'); 
sN = smoothdata(N, 'movmean', smth);

% % plot aaron's data set to find his tau
% figure; hold on % aarons data
% cnt = [3,6,9,12,15,18,21,24,27,30];
% % N = [0.3055, 0.2102, 0.1348, 0.0958, 0.0382, 0.0462, 0.0462, 0.0680, 0.0107, 0.0445];
% N = [0.5947, 0.2028, 0.0777, 0.0492, 0.0229, 0.0137, 0.0115, 0.0022, 0.0063, 0.019];


if afit
    %aaron's fit (confounded by a max value)
    [tau, yfit, ci] = fit_exponential_decay(cnt, N, 0, max(N)*1.1);
    plot(cnt, sN, 'Color', 'k', 'LineWidth', 2)
    plot(cnt, yfit, 'Color', 'r', 'LineWidth', 2)
    plot(cnt, ci(1,:), '--k')
    plot(cnt, ci(2,:), '--k')
    title(['tau=', num2str(tau) ', n = ' num2str(length(all_durs)) ', n_odor = ' num2str(length(all_durs_odor))], 'Interpreter', 'none')

else

    xdata = cnt';
    ydata = N';
    
    % Create a fit type: A*exp(-x/tau) + C
    ft = fittype('A*exp(-x/tau) + C', ...
        'independent', 'x', 'coefficients', {'A', 'tau', 'C'});
    
    % Fit the model
    startPoints = [max(ydata)-min(ydata), 3, min(ydata)];
    [fitresult, gof] = fit(xdata, ydata, ft, 'StartPoint', startPoints);
    fitresult

    % Evaluate the fitted curve and confidence intervals
    xfit = linspace(min(xdata), max(xdata), length(xdata));
    [yCI, yfit] = predint(fitresult, xfit, 0.95, 'functional');  % 'functional' gives CI on mean prediction
    
    % plot
    x_fill = [xfit, fliplr(xfit)];
    y_fill = [yCI(:,1)', fliplr(yCI(:,2)')];
    fill(x_fill, y_fill, [0.8 0.8 0.8], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5);   % lower = more transparent
%     plot(cnt, N, 'Color', 'k', 'LineWidth', 2)
    plot(cnt, sN, 'Color', 'k', 'LineWidth', 2)
    plot(xfit, yfit, 'r', 'LineW', 1)
    
    
    N = histcounts(all_durs_odor, edg, 'Normalization', 'probability'); 
    sN = smoothdata(N, 'movmean', smth);
    plot(cnt, sN, 'Color', 'b', 'LineWidth', 2)
    
    ylim(ylm)
    xlabel('Bump durs (all v odor)')
    ylabel('Probability')
    title(['Exp Fit with Shaded 95% CI, tau=', num2str(fitresult.tau), ', n=', num2str(length(all_durs))])
    legend('95% CI', 'AllBumps', 'Fit', 'OdorBumps')
end

%% onset durs
subplot(1,3,3); hold on
N = histcounts(all_durs_onset, edg, 'Normalization', 'probability'); 
% plot(cnt, N, 'Color', 'k', 'LineWidth', 2)
sN = smoothdata(N, 'movmean', smth);
% N = sN * max(N)/max(sN);
plot(cnt, sN, 'Color', 'k', 'LineWidth', 2)
%     set(gca,'YScale','log')
ylim(ylm)
xlabel('Offset from odor on')


%% offset durs

edg = 0:bnw:80;
cnt = edg(2:end)-bnw/2;

subplot(1,3,2); hold on
N = histcounts(all_durs_offset, edg, 'Normalization', 'probability'); 
% plot(cnt, N, 'Color', 'k', 'LineWidth', 2)
sN = smoothdata(N, 'movmean', smth);
% N = sN * max(N)/max(sN);

if afit
    %aaron's fit exp
    [tau, yfit, ci] = fit_exponential_decay(cnt, N, 0, max(N)*1.1);
    
    plot(cnt, sN, 'Color', 'k', 'LineWidth', 2)
    plot(cnt, yfit, 'Color', 'r', 'LineWidth', 2)
    plot(cnt, ci(1,:), '--k')
    plot(cnt, ci(2,:), '--k')
    title(['tau=', num2str(tau) ', n = ' num2str(length(all_durs_offset))], 'Interpreter', 'none')
else
    xdata = cnt';
    ydata = N';
    
    % Create a fit type: A*exp(-x/tau) + C
    ft = fittype('A*exp(-x/tau) + C', ...
        'independent', 'x', 'coefficients', {'A', 'tau', 'C'});
    
    % Fit the model
    startPoints = [max(ydata)-min(ydata), 3, min(ydata)];
    [fitresult, gof] = fit(xdata, ydata, ft, 'StartPoint', startPoints);
    fitresult
    % Evaluate the fitted curve and confidence intervals
    xfit = linspace(min(xdata), max(xdata), length(xdata));
    [yCI, yfit] = predint(fitresult, xfit, 0.95, 'functional');  % 'functional' gives CI on mean prediction
    
    % plot
    x_fill = [xfit, fliplr(xfit)];
    y_fill = [yCI(:,1)', fliplr(yCI(:,2)')];
    fill(x_fill, y_fill, [0.8 0.8 0.8], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5);   % lower = more transparent
    plot(cnt, sN, 'Color', 'k', 'LineWidth', 2)
    plot(xfit, yfit, 'r', 'LineW', 1)
    
    ylim(ylm)
    xlabel('Offset from odoroff')
    ylabel('Probability')
    legend('95% CI', 'OdorBumps', 'Fit')
    title(['Exp Fit with Shaded 95% CI, tau=', num2str(fitresult.tau), ', n=', num2str(length(all_durs_offset))])
end
end
%% Per-fly scatter of bump offset durations (rows = flies, ordered by mean)

% Remove empty flies and compute means
flyMeans = nan(length(perFlyOffsets),1);
for k = 1:length(perFlyOffsets)
    v = perFlyOffsets{k};
    v = v(~isnan(v));
    perFlyOffsets{k} = v;                   % clean back
    if ~isempty(v)
        flyMeans(k) = mean(v, 'omitnan');
    end
end
valid = ~cellfun(@isempty, perFlyOffsets);

if any(valid)
    % Order flies by mean offset duration 
    [~, ordWithinValid] = sort(flyMeans(valid), 'ascend');
    idxValid = find(valid);
    ord = idxValid(ordWithinValid);

    % Build y positions (row per fly)
    nF = numel(ord);
    yPad = 0.35;                % row half-height for visual padding
    
    fig_no = fig_no + 1;
    figure(fig_no); clf; hold on
    set(gcf, 'Position', [100 100 430 300])

    % Scatter each fly’s offsets on its row
    for r = 1:nF
        k = ord(r);
        x = perFlyOffsets{k};
        y = r * ones(size(x));
        scatter(x, y, 20, 'filled')  % dots
    end

    % Overplot the mean for each fly as a short vertical tick
    for r = 1:nF
        k = ord(r);
        mu = mean(perFlyOffsets{k}, 'omitnan');
        plot([mu mu], [r - yPad, r + yPad], 'k-', 'LineWidth', 2)
    end


    % Y tick labels = fly IDs from use_flies
    yticks(1:nF);
%     ylabels = arrayfun(@(k) sprintf('Fly %d', fig_flyid(k)), ord, 'UniformOutput', false);
    ylabels = arrayfun(@(k) sprintf('Fly %d', use_flies(k)), ord, 'UniformOutput', false);
    yticklabels(ylabels)

%     % Gridlines between rows
%     for r = 0.5:1:(nF+0.5)
%         yline(r, ':', 'Color', [0.8 0.8 0.8])
%     end

    % Nice x-limits
    xmax = max(0.5, max(cellfun(@(v) max(v), perFlyOffsets(valid))));
    xlim([0, xmax*1.05])
    ylim([0.5, nF + 0.5])


    % One-way ANOVA across flies for offset durations
    
    % Flatten all data and make a group vector
    allVals = [];
    allGroups = [];
    for k = 1:length(perFlyOffsets)
        v = perFlyOffsets{k};
        if ~isempty(v)
            allVals = [allVals, v];                         % append durations
            allGroups = [allGroups, repmat(k, 1, length(v))]; % append fly IDs
        end
    end
    
    % Convert to columns (ANOVA likes column vectors)
    allVals = allVals(:);
    allGroups = allGroups(:);
    
    % Run one-way ANOVA
    p = anova1(allVals, allGroups);
%     [p, tbl, stats] = anova1(allVals, allGroups);
%     title('One-way ANOVA of bump offset durations across flies');
    fprintf('One-way ANOVA p-value across flies: %.6g\n', p);

    figure(fig_no);
        % Cosmetics
    xlabel('Bump offset duration (s)')
    ylabel('Fly (ordered by mean offset duration)')
    title(sprintf('Per-fly bump offset durations (n=%d flies, %d bumps, p=%.4f)', ...
                  nF, numel(cell2mat(perFlyOffsets(ord)')), p))


%     box on
else
    warning('No valid per-fly offset durations to plot.');
end
 %%

% % 
% % a0 = max(N);               % starting amplitude
% % tau_guess = (cnt(end)-cnt(1)) / 2;  % crude estimate of decay timescale
% % r0 = 1 / tau_guess;
% % 
% % % Define exponential model
% % ft = fittype('1/a*exp(-x/a)', 'independent', 'x', 'dependent', 'y');
% % % ft = fittype('a*exp(-r*x)', 'independent', 'x', 'dependent', 'y');
% % 
% % [fitresult, gof, output] = fit(cnt', N', ft, 'StartPoint', [a0, r0]);
% % % plot(fitresult, cnt, N, 'r', 'LineWidth', 2)
% % 
% % % Confidence intervals
% % ci = confint(fitresult, 0.95);  % 95% CI for parameters a and b
% % 
% % xfit = linspace(min(cnt), max(cnt), 100);
% % % [yfit, delta] = predint(fitresult, cnt, 0.95, 'functional');  % 'functional' = confidence on mean
% % [yfit, delta] = predint(fitresult, xfit, 0.95, 'functional');  % 'functional' = confidence on mean
% % 
% % 
% % % set(gca,'YScale','log')
% % % ylim(ylm)
% % xlabel('Offset from odoroff')
% % legend('data', 'exp fit')
% % title(['tau=' num2str(tau)])%, ', CI=' num2str(ci)])


%%
if print_figs
    print_fig(1, gcf, [subpath 'bumpdurdist_cwoo_hAck_' add2txt], 'eps', 1); 
    exportgraphics(gcf, [subpath 'bumpdurdist_cwoo_hAck_' add2txt, '.png'])
end
%%

function [tt, dff, oncrossi, offcrossi, odor, fps] = extract_struct_data(trl)

    nonan = ~isnan(trl.calc_ts(:,1));
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




%% alt lines

%     % Find bump threshold from selected trials
%     all_dff = [];
%     for t = 1:size(data,1) % or use "use_trials"
%         if ~isempty(data(t).bump_amp) && data(t).closed_loop
%             dff = data(t).bump_amp;
%             all_dff = [all_dff; dff(~isnan(dff))];  % concatenate
%         end
%     end







% % 
% % a0 = max(N);               % starting amplitude
% % tau_guess = (cnt(end)-cnt(1)) / 2;  % crude estimate of decay timescale
% % r0 = 1 / tau_guess;
% % 
% % % Define exponential model
% % % ft = fittype('a*exp(b*x)', 'independent', 'x', 'dependent', 'y');
% % % ft = fittype('a*exp(-r*x)', 'independent', 'x', 'dependent', 'y');
% % ft = fittype('a*exp(b*x)', 'independent', 'x', 'dependent', 'y');
% % 
% % [fitresult, gof, output] = fit(cnt', N', ft);
% % % [fitresult, gof, output] = fit(cnt', N', ft, 'StartPoint', [a0, r0]);
% % plot(fitresult, cnt, N, 'r')
% % % 
% % % % Confidence intervals
% % % ci = confint(fitresult, 0.95);  % 95% CI for parameters a and b
% % % 
% % % xfit = linspace(min(cnt), max(cnt), 100);
% % % % [yfit, delta] = predint(fitresult, cnt, 0.95, 'functional');  % 'functional' = confidence on mean
% % % [yfit, delta] = predint(fitresult, xfit, 0.95, 'functional');  % 'functional' = confidence on mean
% % % 
% % % % Plot
% % % % figure; hold on
% % % % plot(xdata, ydata, 'bo')
% % % % plot(xfit, yfit, 'g-', 'LineWidth', 2)
% % % plot(xfit, yfit(:,1) + delta, 'k--')
% % % plot(xfit, yfit(:,1) - delta, 'k--')
% % % % legend('Data', 'Fit', '95% CI')
