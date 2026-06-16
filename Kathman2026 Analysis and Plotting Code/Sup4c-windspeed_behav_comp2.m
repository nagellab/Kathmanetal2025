%% windspeed_behavior_metrics_odor_vs_wholetrial_AVEL
% Same as your original script, except:
%   - replaces Forward velocity metric with mean |angular velocity| (deg/s)
%
% Metrics plotted vs wind speed [8 16 33]:
%   1) Heading SD (deg)
%   2) Upwind velocity uwv (mm/s): cosd(heading).*speed
%   3) |Angular velocity| (deg/s): abs(smooth(calc_deltaz*(180/pi)))
%
% Makes TWO sets of figures:
%   Set A: metrics computed DURING ODOR (onIdx:offIdx)
%   Set B: metrics computed over WHOLE TRIAL (all samples after nonan mask)
%
% Plot style:
%   - gray dots = trials (jittered)
%   - black dot+line = fly mean ± SEM at each wind
% Text above each wind:
%   - nT = number of trials
%   - nF = number of flies contributing
% Stats in title:
%   - Fly-averaged ANOVA p (primary)
%   - Trial-level Spearman rho/p (exploratory)

clear; close all;

%% ---------------- USER PARAMS ----------------
use_flies   = [21:27, 31:39];   % <-- edit
keep_winds  = [8 16 33];        % only include these wind speeds
heading_smooth_win = 100;       % samples for smoothing heading
avel_smooth_win    = 150;       % samples for smoothing angular velocity
min_odor_samples   = 20;        % reject trials with too-short odor windows (odor-only set)

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath);
allflies_filenames               % defines filename{flyId}

%% -------- Figure save directory --------
saveDir = '/Users/kathmn01/Desktop/windspeed_behavior_figs/';
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

%% ---------------- COLLECT TRIAL-LEVEL METRICS ----------------
% We'll collect in two "modes": odor-only and whole-trial
M_odor  = init_metric_store();
M_whole = init_metric_store();

nChecked = 0;
nKept = 0;

for ff = 1:numel(use_flies)
    flyId = use_flies(ff);

    Sload = load(filename{flyId}, 'syncstruct');
    data = Sload.syncstruct;

    % single-odor-pulse trials
    use_trials = find(arrayfun(@(t) t.closed_loop && ...
        sum(diff(t.odor)==1)==1 && sum(diff(t.odor)==-1)==1, data));
    if isempty(use_trials), continue; end

    for t = use_trials
        nChecked = nChecked + 1;

        % wind speed filter
        if ~isfield(data(t),'wind_speed') || isempty(data(t).wind_speed) || isnan(data(t).wind_speed)
            continue;
        end
        w = data(t).wind_speed;
        if ~ismember(w, keep_winds), continue; end

        % extract signals
        [tt, odor, heading_deg, speed, uwv, abs_avel] = extract_trial_signals(data(t), heading_smooth_win, avel_smooth_win);
        if isempty(tt) || isempty(odor) || isempty(heading_deg) || isempty(speed) || isempty(abs_avel)
            continue;
        end

        % whole-trial indices (all samples)
        idxWhole = 1:numel(tt);

        % odor indices
        [onIdx, offIdx] = find_single_pulse_onoff(odor);
        hasOdorWindow = ~isempty(onIdx) && ~isempty(offIdx) && offIdx > onIdx;
        if hasOdorWindow
            idxOdor = onIdx:offIdx;
        else
            idxOdor = [];
        end

        % ---- Compute WHOLE-TRIAL metrics (always)
        yH_whole   = heading_sd_deg(heading_deg(idxWhole));
        yUWV_whole = mean(uwv(idxWhole), 'omitnan');
        yA_whole   = mean(abs_avel(idxWhole), 'omitnan');

        if all(isfinite([yH_whole yUWV_whole yA_whole]))
            M_whole = append_metrics(M_whole, flyId, w, yH_whole, yUWV_whole, yA_whole);
        end

        % ---- Compute ODOR metrics (only if odor exists and long enough)
        if hasOdorWindow && numel(idxOdor) >= min_odor_samples
            yH_odor   = heading_sd_deg(heading_deg(idxOdor));
            yUWV_odor = mean(uwv(idxOdor), 'omitnan');
            yA_odor   = mean(abs_avel(idxOdor), 'omitnan');

            if all(isfinite([yH_odor yUWV_odor yA_odor]))
                M_odor = append_metrics(M_odor, flyId, w, yH_odor, yUWV_odor, yA_odor);
            end
        end

        nKept = nKept + 1;
    end
end

fprintf('Checked %d single-pulse trials; processed %d trials at winds [%s].\n', ...
    nChecked, nKept, num2str(keep_winds));

% Ensure we only plot winds actually present in the collected data
windVals_odor  = sort(unique(M_odor.trialWind));
windVals_whole = sort(unique(M_whole.trialWind));

if isempty(windVals_odor)
    warning('No ODOR-window trials retained (maybe odor too short / missing).');
end
if isempty(windVals_whole)
    error('No WHOLE-trial data retained. Check filters/use_flies/winds.');
end

%% ---------------- SET A: ODOR ONLY (3 FIGURES) ----------------
if ~isempty(windVals_odor)
    make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialHeadingSD, windVals_odor, ...
        'Heading SD during odor (deg)', sprintf('Heading SD vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);

    make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialUWV, windVals_odor, ...
        'Mean upwind velocity during odor (mm/s)', sprintf('Upwind velocity vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);

    make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialAbsAvel, windVals_odor, ...
        'Mean |angular velocity| during odor (deg/s)', sprintf('|ang vel| vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);
end

%% ---------------- SET B: WHOLE TRIAL (3 FIGURES) ----------------
make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialHeadingSD, windVals_whole, ...
    'Heading SD whole trial (deg)', sprintf('Heading SD vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);

make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialUWV, windVals_whole, ...
    'Mean upwind velocity whole trial (mm/s)', sprintf('Upwind velocity vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);

make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialAbsAvel, windVals_whole, ...
    'Mean |angular velocity| whole trial (deg/s)', sprintf('|ang vel| vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);

%% ===================== HELPERS =====================

function M = init_metric_store()
    M.trialFly = [];
    M.trialWind = [];
    M.trialHeadingSD = [];
    M.trialUWV = [];
    M.trialAbsAvel = [];
end

function M = append_metrics(M, flyId, w, hSD, uwv, absAvel)
    M.trialFly(end+1,1) = flyId;
    M.trialWind(end+1,1) = w;
    M.trialHeadingSD(end+1,1) = hSD;
    M.trialUWV(end+1,1) = uwv;
    M.trialAbsAvel(end+1,1) = absAvel;
end

function y = heading_sd_deg(heading_deg)
    % unwrap heading and compute SD in degrees
    uw = unwrap(deg2rad(heading_deg(:)));
    y = std(rad2deg(uw), 0, 'omitnan');
end

function make_summary_plot(trialFly, trialWind, y, windVals, ylab, mainTitle, saveDir)
    % Trial-level Spearman (exploratory)
    [rho, pSpearman] = corr(trialWind, y, 'Type','Spearman', 'Rows','complete');

    % Fly-averaged ANOVA (primary)
    flyMeans = [];
    windGroups = [];
    for wi = 1:numel(windVals)
        w = windVals(wi);
        selW = (trialWind == w);
        fliesW = unique(trialFly(selW));
        for f = fliesW(:)'
            flyMeans(end+1,1) = mean(y(selW & trialFly==f), 'omitnan'); %#ok<AGROW>
            windGroups(end+1,1) = w; %#ok<AGROW>
        end
    end
    if numel(unique(windGroups)) >= 2 && numel(flyMeans) >= 3
        pAnova = anova1(flyMeans, windGroups, 'off');
    else
        pAnova = NaN;
    end

    % Plot
    figure('Color','w'); clf;
    set(gcf, 'Position', [200 200 950 450]);
    ax = axes; hold(ax,'on'); box(ax,'on');

    rng(0);
    jit = 0.25;

    nTrialByWind = zeros(numel(windVals),1);
    nFlyByWind   = zeros(numel(windVals),1);
    flyMuByWind  = nan(numel(windVals),1);
    flySEMByWind = nan(numel(windVals),1);

    for wi = 1:numel(windVals)
        w = windVals(wi);
        sel = (trialWind == w);

        yy = y(sel);
        xx = w + jit*(rand(sum(sel),1)-0.5);
        plot(xx, yy, '.', 'Color', [0.65 0.65 0.65], 'MarkerSize', 12);

        nTrialByWind(wi) = sum(sel);

        % Fly means at this wind
        fliesW = unique(trialFly(sel));
        muFly = nan(numel(fliesW),1);
        for k = 1:numel(fliesW)
            f = fliesW(k);
            muFly(k) = mean(y(sel & trialFly==f), 'omitnan');
        end

        nFlyByWind(wi) = numel(muFly);
        flyMuByWind(wi) = mean(muFly, 'omitnan');
        if numel(muFly) >= 2
            flySEMByWind(wi) = std(muFly, 'omitnan') / sqrt(numel(muFly));
        else
            flySEMByWind(wi) = NaN;
        end
    end

    errorbar(windVals, flyMuByWind, flySEMByWind, 'ko-', ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    % annotate counts (same as your original)
    yl = ylim(ax);
    yText = yl(2) - 0.05*(yl(2)-yl(1));
    for wi = 1:numel(windVals)
        txt = sprintf('nT=%d, nF=%d', nTrialByWind(wi), nFlyByWind(wi));
        text(windVals(wi), yText, txt, 'HorizontalAlignment','center', 'FontSize', 9);
    end

    xlabel('Wind speed');
    ylabel(ylab);
    xticks(windVals);
    grid(ax,'off');
    box(ax,'off');
    set(ax,'TickDir','out');

    title(ax, { ...
        mainTitle, ...
        sprintf('Fly-avg ANOVA p=%s,   Trial Spearman \\rho=%s (p=%s)', ...
            format_p(pAnova), fmt_num(rho,2), format_p(pSpearman)) ...
        });

    %% ---- Save figure ----
    safeName = regexprep(mainTitle, '[^a-zA-Z0-9]', '_');
    figFile = fullfile(saveDir, [safeName '.fig']);
    pngFile = fullfile(saveDir, [safeName '.png']);
    savefig(gcf, figFile);
    exportgraphics(gcf, pngFile, 'Resolution', 300);

    hold(ax,'off');
end

function [tt, odor, heading_deg, speed, uwv, abs_avel] = extract_trial_signals(trl, heading_smooth_win, avel_smooth_win)
    tt = []; odor = []; heading_deg = []; speed = []; uwv = []; abs_avel = [];

    if ~isfield(trl,'calc_ts') || isempty(trl.calc_ts), return; end
    nonan = ~isnan(trl.calc_ts(:,1));

    tt   = trl.calc_ts(nonan);
    odor = trl.odor(nonan);

    % speed
    if isfield(trl,'calc_speed') && ~isempty(trl.calc_speed)
        speed = trl.calc_speed(nonan);
    else
        speed = nan(size(tt));
    end

    % heading
    if isfield(trl,'calc_windpos') && ~isempty(trl.calc_windpos)
        h = trl.calc_windpos(nonan);      % degrees
        h = deg2rad(h);
        h = unwrap(h);
        h = smoothdata(h,'movmean',heading_smooth_win);
        heading_deg = wrapTo180(rad2deg(h));
    else
        heading_deg = nan(size(tt));
    end

    % upwind velocity (as in your original)
    uwv = cosd(heading_deg(:)) .* speed(:);

    % |angular velocity| (deg/s)
    if isfield(trl,'calc_deltaz') && ~isempty(trl.calc_deltaz)
        avel = trl.calc_deltaz(nonan) * (180/pi);
        avel = smoothdata(avel, 'movmean', avel_smooth_win);
        abs_avel = abs(avel);
    else
        abs_avel = nan(size(tt));
    end
end

function [onIdx, offIdx] = find_single_pulse_onoff(odor)
    onIdx = []; offIdx = [];
    odor = odor(:);

    oncrossi  = find(diff(odor)==1);
    offcrossi = find(diff(odor)==-1);

    if ~isempty(odor) && odor(1)==1
        oncrossi = [1; oncrossi(:)];
    else
        oncrossi = oncrossi(:);
    end
    if ~isempty(odor) && odor(end)==1
        offcrossi = [offcrossi(:); numel(odor)];
    else
        offcrossi = offcrossi(:);
    end

    if isempty(oncrossi) || isempty(offcrossi), return; end
    onIdx  = oncrossi(1);
    offIdx = offcrossi(end);
end

function s = format_p(p)
    if isnan(p)
        s = 'NaN';
    elseif p < 1e-3
        s = sprintf('%.1e', p);
    else
        s = sprintf('%.3g', p);
    end
end

function s = fmt_num(x, nd)
    if isnan(x), s = 'NaN'; return; end
    s = sprintf(['%.' num2str(nd) 'f'], x);
end



% %% windspeed_behav_comp_avel
% % Behavior vs wind speed (8,16,33 cm/s)
% % Metrics:
% %   - Upwind velocity
% %   - |Angular velocity|
% %   - Heading SD
% % Computed over:
% %   - Whole trial
% %   - Odor period only
% 
% clear; close all;
% 
% %% ---------------- USER PARAMS ----------------
% use_flies = [21:27, 31:39];
% 
% windKeep = [8 16 33];
% 
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
% cd(headpath);
% allflies_filenames
% 
% 
% %% ---------------- COLLECT DATA ----------------
% 
% ALL = struct();
% 
% for ff = 1:numel(use_flies)
% 
%     flyId = use_flies(ff);
% 
%     Sload = load(filename{flyId},'syncstruct');
%     data = Sload.syncstruct;
% 
%     % single pulse trials
%     use_trials = find(arrayfun(@(t) t.closed_loop && ...
%         sum(diff(t.odor)==1)==1 && sum(diff(t.odor)==-1)==1, data));
% 
%     if isempty(use_trials), continue; end
% 
% 
%     for t = use_trials
% 
%         if ~isfield(data(t),'wind_speed') || isnan(data(t).wind_speed)
%             continue;
%         end
% 
%         wind = data(t).wind_speed;
% 
%         if ~ismember(wind, windKeep)
%             continue;
%         end
% 
%         S = extract_struct_data(data(t));
%         if isempty(S.tt), continue; end
% 
%         % odor indices
%         if isempty(S.odorOn) || isempty(S.odorOff)
%             continue;
%         end
% 
% 
%         %% ----------- Metrics -----------
% 
%         % Upwind velocity
%         uwv = cosd(S.heading) .* S.speed;
% 
%         % |angular velocity|
%         aavel = abs(S.avel);
% 
%         % Heading SD
%         hsd = movstd(unwrap(deg2rad(S.heading)), 150, 'omitnan');
%         hsd = rad2deg(hsd);
% 
% 
%         %% ----------- Whole trial -----------
% 
%         WT.uwv  = mean(uwv,'omitnan');
%         WT.avel = mean(aavel,'omitnan');
%         WT.hsd  = mean(hsd,'omitnan');
% 
% 
%         %% ----------- Odor only -----------
% 
%         odIdx = S.odorOn:S.odorOff;
% 
%         OD.uwv  = mean(uwv(odIdx),'omitnan');
%         OD.avel = mean(aavel(odIdx),'omitnan');
%         OD.hsd  = mean(hsd(odIdx),'omitnan');
% 
% 
%         %% ----------- Store -----------
% 
%         key = wind_key(wind);
% 
%         if ~isfield(ALL,key)
% 
%             ALL.(key) = struct( ...
%                 'fly',[], ...
%                 'WT',struct('uwv',[],'avel',[],'hsd',[]), ...
%                 'OD',struct('uwv',[],'avel',[],'hsd',[]) );
%         end
% 
% 
%         ALL.(key).fly(end+1,1) = flyId;
% 
%         ALL.(key).WT.uwv(end+1,1)  = WT.uwv;
%         ALL.(key).WT.avel(end+1,1) = WT.avel;
%         ALL.(key).WT.hsd(end+1,1)  = WT.hsd;
% 
%         ALL.(key).OD.uwv(end+1,1)  = OD.uwv;
%         ALL.(key).OD.avel(end+1,1) = OD.avel;
%         ALL.(key).OD.hsd(end+1,1)  = OD.hsd;
% 
%     end
% end
% 
% 
% %% ---------------- PLOT ----------------
% 
% metrics = {'uwv','avel','hsd'};
% names   = {'Upwind velocity','|Angular velocity|','Heading SD'};
% units   = {'mm/s','deg/s','deg'};
% 
% periods = {'WT','OD'};
% pNames  = {'Whole trial','Odor only'};
% 
% 
% for pi = 1:2
% 
%     per = periods{pi};
% 
%     for mi = 1:3
% 
%         met = metrics{mi};
% 
%         figure('Color','w'); clf;
%         set(gcf,'Position',[200 200 700 450]);
% 
%         ax = axes; hold(ax,'on'); box(ax,'on');
% 
%         windVals = windKeep;
% 
%         trialWind = [];
%         trialVal  = [];
%         trialFly  = [];
% 
% 
%         %% Collect trial-level
% 
%         for w = windVals
% 
%             key = wind_key(w);
% 
%             if ~isfield(ALL,key), continue; end
% 
%             v = ALL.(key).(per).(met);
%             f = ALL.(key).fly;
% 
%             trialVal  = [trialVal; v];
%             trialWind = [trialWind; repmat(w,numel(v),1)];
%             trialFly  = [trialFly; f];
% 
%         end
% 
% 
%         %% Plot trials (gray)
% 
%         rng(0);
%         jit = 0.12;
% 
%         for wi = 1:numel(windVals)
% 
%             w = windVals(wi);
% 
%             sel = trialWind==w;
% 
%             x = w + jit*(rand(sum(sel),1)-0.5);
%             y = trialVal(sel);
% 
%             plot(x,y,'.','Color',[.6 .6 .6],'MarkerSize',12);
% 
%         end
% 
% 
%         %% Fly means
% 
%         muFly = nan(numel(windVals),1);
%         semFly = nan(numel(windVals),1);
% 
%         for wi = 1:numel(windVals)
% 
%             w = windVals(wi);
% 
%             sel = trialWind==w;
% 
%             flies = unique(trialFly(sel));
% 
%             m = nan(numel(flies),1);
% 
%             for k = 1:numel(flies)
% 
%                 f = flies(k);
% 
%                 m(k) = mean(trialVal(sel & trialFly==f),'omitnan');
% 
%             end
% 
%             muFly(wi) = mean(m,'omitnan');
% 
%             if numel(m)>1
%                 semFly(wi) = std(m,'omitnan')/sqrt(numel(m));
%             end
% 
%         end
% 
% 
%         errorbar(windVals,muFly,semFly,'ko-','LineWidth',1.5,...
%             'MarkerFaceColor','k','MarkerSize',6);
% 
% 
%         %% Stats
% 
%         [rho,pS] = corr(trialWind,trialVal,'Type','Spearman','Rows','complete');
% 
%         flyMeans = [];
%         flyGroup = [];
% 
%         for wi = 1:numel(windVals)
% 
%             w = windVals(wi);
% 
%             sel = trialWind==w;
% 
%             flies = unique(trialFly(sel));
% 
%             for k = 1:numel(flies)
% 
%                 f = flies(k);
% 
%                 flyMeans(end+1,1) = mean(trialVal(sel & trialFly==f),'omitnan'); %#ok<AGROW>
%                 flyGroup(end+1,1) = w; %#ok<AGROW>
% 
%             end
%         end
% 
%         if numel(unique(flyGroup))>=2
%             pA = anova1(flyMeans,flyGroup,'off');
%         else
%             pA = NaN;
%         end
% 
% 
%         %% Title
% 
%         title(ax,{ ...
%             sprintf('%s vs wind speed (%s)',names{mi},pNames{pi}), ...
%             sprintf('ANOVA p=%s, Spearman \\rho=%.2f (p=%s)', ...
%             format_p(pA),rho,format_p(pS))});
% 
% 
%         xlabel('Wind speed (cm/s)');
%         ylabel(sprintf('%s (%s)',names{mi},units{mi}));
% 
%         xticks(windVals);
%         set(ax,'TickDir','out');
% 
%     end
% end
% 
% 
% %% ================= HELPERS =================
% 
% function S = extract_struct_data(trl)
% 
%     S = struct('tt',[],'speed',[],'avel',[],'heading',[], ...
%                'odor',[],'odorOn',[],'odorOff',[]);
% 
%     if ~isfield(trl,'calc_ts'), return; end
% 
%     nonan = ~isnan(trl.calc_ts(:,1));
% 
%     S.tt = trl.calc_ts(nonan);
% 
%     S.speed = sqrt( ...
%         trl.calc_deltapitch(nonan).^2 + ...
%         trl.calc_deltaroll(nonan).^2 );
% 
%     S.avel = trl.calc_deltaz(nonan)*(180/pi);
%     S.avel = smoothdata(S.avel,'movmean',150);
% 
%     S.heading = wrapTo180(rad2deg( ...
%         smoothdata(unwrap(deg2rad(trl.calc_windpos(nonan))), ...
%         'movmean',100)));
% 
%     S.odor = trl.odor(nonan);
% 
%     on = find(diff(S.odor)==1);
%     off = find(diff(S.odor)==-1);
% 
%     if ~isempty(S.odor) && S.odor(1)==1
%         on = [1;on];
%     end
% 
%     if ~isempty(S.odor) && S.odor(end)==1
%         off = [off;numel(S.odor)];
%     end
% 
%     if ~isempty(on),  S.odorOn = on(1); end
%     if ~isempty(off), S.odorOff = off(end); end
% 
% end
% 
% 
% function key = wind_key(w)
%     key = sprintf('w_%d',round(w));
% end
% 
% 
% function s = format_p(p)
% 
%     if isnan(p)
%         s = 'NaN';
%     elseif p<1e-3
%         s = sprintf('%.1e',p);
%     else
%         s = sprintf('%.3g',p);
%     end
% end




% % % %% windspeed_behavior_metrics_odor_vs_wholetrial
% % % % Plots behavior metrics vs wind speed for winds [8 16 33], using single-odor-pulse trials:
% % % %   1) Heading SD (deg)
% % % %   2) Upwind velocity uwv (mm/s): cosd(heading).*speed
% % % %   3) Forward velocity fvel (mm/s): calc_deltapitch*9.52/2
% % % %
% % % % Makes TWO sets of figures:
% % % %   Set A: metrics computed DURING ODOR (onIdx:offIdx)
% % % %   Set B: metrics computed over WHOLE TRIAL (all samples after nonan mask)
% % % %
% % % % Plot style:
% % % %   - gray dots = trials (jittered)
% % % %   - black dot+line = fly mean ± SEM at each wind
% % % % Stats in title:
% % % %   - Fly-averaged ANOVA p (primary)
% % % %   - Trial-level Spearman rho/p (exploratory)
% % % 
% % % clear; close all;
% % % 
% % % %% ---------------- USER PARAMS ----------------
% % % use_flies   = [21:27, 31:39];   % <-- edit
% % % keep_winds  = [8 16 33];        % only include these wind speeds
% % % heading_smooth_win = 100;       % samples for smoothing heading
% % % min_odor_samples   = 20;        % reject trials with too-short odor windows (odor-only set)
% % % 
% % % headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
% % % cd(headpath);
% % % allflies_filenames               % defines filename{flyId}
% % % 
% % % %% -------- Figure save directory --------
% % % saveDir = '/Users/kathmn01/Desktop/windspeed_behavior_figs/';
% % % 
% % % if ~exist(saveDir, 'dir')
% % %     mkdir(saveDir);
% % % end
% % % 
% % % %% ---------------- COLLECT TRIAL-LEVEL METRICS ----------------
% % % % We'll collect in two "modes": odor-only and whole-trial
% % % M_odor  = init_metric_store();
% % % M_whole = init_metric_store();
% % % 
% % % nChecked = 0;
% % % nKept = 0;
% % % 
% % % for ff = 1:numel(use_flies)
% % %     flyId = use_flies(ff);
% % % 
% % %     Sload = load(filename{flyId}, 'syncstruct');
% % %     data = Sload.syncstruct;
% % % 
% % %     % single-odor-pulse trials
% % %     use_trials = find(arrayfun(@(t) t.closed_loop && ...
% % %         sum(diff(t.odor)==1)==1 && sum(diff(t.odor)==-1)==1, data));
% % %     if isempty(use_trials), continue; end
% % % 
% % %     for t = use_trials
% % %         nChecked = nChecked + 1;
% % % 
% % %         % wind speed filter
% % %         if ~isfield(data(t),'wind_speed') || isempty(data(t).wind_speed) || isnan(data(t).wind_speed)
% % %             continue;
% % %         end
% % %         w = data(t).wind_speed;
% % %         if ~ismember(w, keep_winds), continue; end
% % % 
% % %         % extract signals
% % %         [tt, odor, heading_deg, speed, fvel, uwv] = extract_trial_signals(data(t), heading_smooth_win);
% % %         if isempty(tt) || isempty(odor) || isempty(heading_deg) || isempty(speed) || isempty(fvel)
% % %             continue;
% % %         end
% % % 
% % %         % whole-trial indices (all samples)
% % %         idxWhole = 1:numel(tt);
% % % 
% % %         % odor indices
% % %         [onIdx, offIdx] = find_single_pulse_onoff(odor);
% % %         hasOdorWindow = ~isempty(onIdx) && ~isempty(offIdx) && offIdx > onIdx;
% % %         if hasOdorWindow
% % %             idxOdor = onIdx:offIdx;
% % %         else
% % %             idxOdor = [];
% % %         end
% % % 
% % %         % ---- Compute WHOLE-TRIAL metrics (always)
% % %         yH_whole = heading_sd_deg(heading_deg(idxWhole));
% % % %         yUWV_whole = mean_upwind_velocity(heading_deg(idxWhole), speed(idxWhole));
% % %         yUWV_whole = mean(uwv(idxWhole), 'omitnan');
% % % 
% % %         yF_whole = mean(fvel(idxWhole), 'omitnan');
% % % 
% % %         if all(isfinite([yH_whole yUWV_whole yF_whole]))
% % %             M_whole = append_metrics(M_whole, flyId, w, yH_whole, yUWV_whole, yF_whole);
% % %         end
% % % 
% % %         % ---- Compute ODOR metrics (only if odor exists and long enough)
% % %         if hasOdorWindow && numel(idxOdor) >= min_odor_samples
% % %             yH_odor = heading_sd_deg(heading_deg(idxOdor));
% % % %             yUWV_odor = mean_upwind_velocity(heading_deg(idxOdor), speed(idxOdor));
% % %             yUWV_odor = mean(uwv(idxOdor), 'omitnan');
% % % 
% % %             yF_odor = mean(fvel(idxOdor), 'omitnan');
% % % 
% % %             if all(isfinite([yH_odor yUWV_odor yF_odor]))
% % %                 M_odor = append_metrics(M_odor, flyId, w, yH_odor, yUWV_odor, yF_odor);
% % %             end
% % %         end
% % % 
% % %         nKept = nKept + 1;
% % %     end
% % % end
% % % 
% % % fprintf('Checked %d single-pulse trials; processed %d trials at winds [%s].\n', ...
% % %     nChecked, nKept, num2str(keep_winds));
% % % 
% % % % Ensure we only plot winds actually present in the collected data
% % % windVals_odor  = sort(unique(M_odor.trialWind));
% % % windVals_whole = sort(unique(M_whole.trialWind));
% % % 
% % % if isempty(windVals_odor)
% % %     warning('No ODOR-window trials retained (maybe odor too short / missing).');
% % % end
% % % if isempty(windVals_whole)
% % %     error('No WHOLE-trial data retained. Check filters/use_flies/winds.');
% % % end
% % % 
% % % %% ---------------- SET A: ODOR ONLY (3 FIGURES) ----------------
% % % if ~isempty(windVals_odor)
% % %     make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialHeadingSD, windVals_odor, ...
% % %         'Heading SD during odor (deg)', sprintf('Heading SD vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);
% % % 
% % %     make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialUWV, windVals_odor, ...
% % %         'Mean upwind velocity during odor (mm/s)', sprintf('Upwind velocity vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);
% % % 
% % %     make_summary_plot(M_odor.trialFly, M_odor.trialWind, M_odor.trialFwdV, windVals_odor, ...
% % %         'Mean forward velocity during odor (mm/s)', sprintf('Forward velocity vs wind (ODOR only),  winds=%s', num2str(keep_winds)), saveDir);
% % % end
% % % 
% % % %% ---------------- SET B: WHOLE TRIAL (3 FIGURES) ----------------
% % % make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialHeadingSD, windVals_whole, ...
% % %     'Heading SD whole trial (deg)', sprintf('Heading SD vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);
% % % 
% % % make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialUWV, windVals_whole, ...
% % %     'Mean upwind velocity whole trial (mm/s)', sprintf('Upwind velocity vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);
% % % 
% % % make_summary_plot(M_whole.trialFly, M_whole.trialWind, M_whole.trialFwdV, windVals_whole, ...
% % %     'Mean forward velocity whole trial (mm/s)', sprintf('Forward velocity vs wind (WHOLE trial), winds=%s', num2str(keep_winds)), saveDir);
% % % 
% % % %% ===================== HELPERS =====================
% % % 
% % % function M = init_metric_store()
% % %     M.trialFly = [];
% % %     M.trialWind = [];
% % %     M.trialHeadingSD = [];
% % %     M.trialUWV = [];
% % %     M.trialFwdV = [];
% % % end
% % % 
% % % function M = append_metrics(M, flyId, w, hSD, uwv, fwd)
% % %     M.trialFly(end+1,1) = flyId;
% % %     M.trialWind(end+1,1) = w;
% % %     M.trialHeadingSD(end+1,1) = hSD;
% % %     M.trialUWV(end+1,1) = uwv;
% % %     M.trialFwdV(end+1,1) = fwd;
% % % end
% % % 
% % % function y = heading_sd_deg(heading_deg)
% % %     % unwrap heading and compute SD in degrees
% % %     uw = unwrap(deg2rad(heading_deg(:)));
% % %     y = std(rad2deg(uw), 0, 'omitnan');
% % % end
% % % 
% % % function mu = mean_upwind_velocity(heading_deg, speed)
% % %     mu = mean(cosd(heading_deg(:)) .* speed(:), 'omitnan');
% % % end
% % % 
% % % function make_summary_plot(trialFly, trialWind, y, windVals, ylab, mainTitle, saveDir)
% % %     % Trial-level Spearman (exploratory)
% % %     [rho, pSpearman] = corr(trialWind, y, 'Type','Spearman', 'Rows','complete');
% % % 
% % %     % Fly-averaged ANOVA (primary)
% % %     flyMeans = [];
% % %     windGroups = [];
% % %     for wi = 1:numel(windVals)
% % %         w = windVals(wi);
% % %         selW = (trialWind == w);
% % %         fliesW = unique(trialFly(selW));
% % %         for f = fliesW(:)'
% % %             flyMeans(end+1,1) = mean(y(selW & trialFly==f), 'omitnan'); %#ok<AGROW>
% % %             windGroups(end+1,1) = w; %#ok<AGROW>
% % %         end
% % %     end
% % %     if numel(unique(windGroups)) >= 2 && numel(flyMeans) >= 3
% % %         pAnova = anova1(flyMeans, windGroups, 'off');
% % %     else
% % %         pAnova = NaN;
% % %     end
% % % 
% % %     % Plot
% % %     figure('Color','w'); clf;
% % %     set(gcf, 'Position', [200 200 950 450]);
% % %     ax = axes; hold(ax,'on'); box(ax,'on');
% % % 
% % %     rng(0);
% % %     jit = 0.25;
% % % 
% % %     nTrialByWind = zeros(numel(windVals),1);
% % %     nFlyByWind   = zeros(numel(windVals),1);
% % %     flyMuByWind  = nan(numel(windVals),1);
% % %     flySEMByWind = nan(numel(windVals),1);
% % % 
% % %     for wi = 1:numel(windVals)
% % %         w = windVals(wi);
% % %         sel = (trialWind == w);
% % % 
% % %         yy = y(sel);
% % %         xx = w + jit*(rand(sum(sel),1)-0.5);
% % %         plot(xx, yy, '.', 'Color', [0.65 0.65 0.65], 'MarkerSize', 12);
% % % 
% % %         nTrialByWind(wi) = sum(sel);
% % % 
% % %         % Fly means at this wind
% % %         fliesW = unique(trialFly(sel));
% % %         muFly = nan(numel(fliesW),1);
% % %         for k = 1:numel(fliesW)
% % %             f = fliesW(k);
% % %             muFly(k) = mean(y(sel & trialFly==f), 'omitnan');
% % %         end
% % % 
% % %         nFlyByWind(wi) = numel(muFly);
% % %         flyMuByWind(wi) = mean(muFly, 'omitnan');
% % %         if numel(muFly) >= 2
% % %             flySEMByWind(wi) = std(muFly, 'omitnan') / sqrt(numel(muFly));
% % %         else
% % %             flySEMByWind(wi) = NaN;
% % %         end
% % %     end
% % % 
% % %     errorbar(windVals, flyMuByWind, flySEMByWind, 'ko-', ...
% % %         'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 6);
% % % 
% % %     % annotate counts
% % %     yl = ylim(ax);
% % %     yText = yl(2) - 0.05*(yl(2)-yl(1));
% % %     for wi = 1:numel(windVals)
% % %         txt = sprintf('nT=%d, nF=%d', nTrialByWind(wi), nFlyByWind(wi));
% % %         text(windVals(wi), yText, txt, 'HorizontalAlignment','center', 'FontSize', 9);
% % %     end
% % % 
% % %     xlabel('Wind speed');
% % %     ylabel(ylab);
% % %     xticks(windVals);
% % %     grid(ax,'off');
% % %     box(ax, 'off')
% % %     set(ax,'TickDir','out');
% % % 
% % %     title(ax, { ...
% % %         mainTitle, ...
% % %         sprintf('Fly-avg ANOVA p=%s,   Trial Spearman \\rho=%s (p=%s)', ...
% % %             format_p(pAnova), fmt_num(rho,2), format_p(pSpearman)) ...
% % %         });
% % % 
% % %     %% ---- Save figure ----
% % %     safeName = regexprep(mainTitle, '[^a-zA-Z0-9]', '_');
% % %     
% % %     figFile = fullfile(saveDir, [safeName '.fig']);
% % %     pngFile = fullfile(saveDir, [safeName '.png']);
% % %     
% % %     savefig(gcf, figFile);
% % %     exportgraphics(gcf, pngFile, 'Resolution', 300);
% % % 
% % % 
% % % 
% % %     hold(ax,'off');
% % % end
% % % 
% % % function [tt, odor, heading_deg, speed, fvel, uwv] = extract_trial_signals(trl, heading_smooth_win)
% % %     tt = []; odor = []; heading_deg = []; speed = []; fvel = [];
% % % 
% % %     if ~isfield(trl,'calc_ts') || isempty(trl.calc_ts), return; end
% % %     nonan = ~isnan(trl.calc_ts(:,1));
% % % 
% % %     tt   = trl.calc_ts(nonan);
% % %     odor = trl.odor(nonan);
% % % 
% % %     % speed
% % %     if isfield(trl,'calc_speed') && ~isempty(trl.calc_speed)
% % %         speed = trl.calc_speed(nonan);
% % %     else
% % %         speed = nan(size(tt));
% % %     end
% % % 
% % %     % forward vel
% % %     if isfield(trl,'calc_deltapitch') && ~isempty(trl.calc_deltapitch)
% % %         fvel = trl.calc_deltapitch(nonan)*9.52/2;
% % %         fvel = smoothdata(fvel, 'movmean', 150);
% % %     else
% % %         fvel = nan(size(tt));
% % %     end
% % % 
% % %     % heading
% % %     if isfield(trl,'calc_windpos') && ~isempty(trl.calc_windpos)
% % %         % calc_windpos is in degrees
% % %         h = trl.calc_windpos(nonan);
% % %     
% % %         % unwrap → smooth → wrap
% % %         h = deg2rad(h);
% % %         h = unwrap(h);
% % %         h = smoothdata(h,'movmean',100);
% % %         heading_deg = wrapTo180(rad2deg(h));
% % %     else
% % %         heading_deg = nan(size(tt));
% % %     end
% % % 
% % %     uwv = cosd(heading_deg(:)) .* speed(:);
% % % 
% % % end
% % % 
% % % function [onIdx, offIdx] = find_single_pulse_onoff(odor)
% % %     onIdx = []; offIdx = [];
% % %     odor = odor(:);
% % % 
% % %     oncrossi  = find(diff(odor)==1);
% % %     offcrossi = find(diff(odor)==-1);
% % % 
% % %     if ~isempty(odor) && odor(1)==1
% % %         oncrossi = [1; oncrossi(:)];
% % %     else
% % %         oncrossi = oncrossi(:);
% % %     end
% % %     if ~isempty(odor) && odor(end)==1
% % %         offcrossi = [offcrossi(:); numel(odor)];
% % %     else
% % %         offcrossi = offcrossi(:);
% % %     end
% % % 
% % %     if isempty(oncrossi) || isempty(offcrossi), return; end
% % %     onIdx  = oncrossi(1);
% % %     offIdx = offcrossi(end);
% % % end
% % % 
% % % function s = format_p(p)
% % %     if isnan(p)
% % %         s = 'NaN';
% % %     elseif p < 1e-3
% % %         s = sprintf('%.1e', p);
% % %     else
% % %         s = sprintf('%.3g', p);
% % %     end
% % % end
% % % 
% % % function s = fmt_num(x, nd)
% % %     if isnan(x), s = 'NaN'; return; end
% % %     s = sprintf(['%.' num2str(nd) 'f'], x);
% % % end
