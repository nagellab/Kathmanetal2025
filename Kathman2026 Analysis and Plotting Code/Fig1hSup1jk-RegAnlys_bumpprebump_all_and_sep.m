%% RegAnlys_all_and_sep
% Detect EXACTLY ONE set of bumps/pre-bumps, then:
%   1) Plot ALL bumps pooled (like RegAnlys_indiscrim_ReviewerReq)
%   2) Split the SAME bumps into epochs by bump ONSET: pre / during / post odor,
%      and plot each epoch (paired bump vs pre-bump)
%   3) Summary-only plots: epoch-colored scatters of BUMPS only (one figure per metric)
%
% Key principle: detection + inclusion (valid pre-window, finite x/y) happens ONCE,
% then the resulting bump list is reused for all plotting modes.

clear; close all;

%% ---------------- USER PARAMS ----------------
use_flies = [21:26, 30, 32:35, 37:39];  % your set
sdm      = 0.5;                         % SD multiplier for bumpfinder threshold
odorThr  = 0.5;                         % bumptrig odor threshold crossing
minOn    = 3;                           % bumpfinder2 minOn (sec)
maxOff   = 3;                           % bumpfinder2 maxOff (sec)
preN     = 500;                         % samples immediately before bump onset for "pre" window

% Where your fly files live (must contain allflies_filenames.m that defines filename{flyId})
headpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/';
cd(headpath);
allflies_filenames

% Output folder (optional)
saveFigs = false;
outdir   = '/Users/kathmn01/Desktop/revision figs/';
if saveFigs && ~exist(outdir,'dir'); mkdir(outdir); end

%% ---------------- METRICS ----------------
% Each metric defines:
%  - getY:     how to get the raw sample-by-sample metric vector
%  - statfun:  how to summarize that metric within a bump window (and within pre-window)
metrics = {
    struct( ...
        'name','Abs ang vel', ...
        'yLabel','|ang vel| (deg/s)', ...
        'getY', @(S) abs(S.avel), ...
        'statfun', @(v) mean(v,'omitnan'), ...
        'forceYLim', [0 110])

    struct( ...
        'name','Heading SD', ...
        'yLabel','heading SD (deg)', ...
        'getY', @(S) S.heading, ...
        'statfun', @(h) rad2deg(std(unwrap(deg2rad(h)), 'omitnan')), ...
        'forceYLim', [])
};

%% ---------------- MAIN ----------------
for mi = 1:numel(metrics)
    M = metrics{mi};

    % Collect ONCE for this metric
    R = collect_bumps_and_prebumps(use_flies, filename, sdm, odorThr, minOn, maxOff, preN, M.getY, M.statfun);

    % ---- 1) ALL bumps pooled (bump + pre) ----
    fig1 = figure('Name', ['ALL pooled - ' M.name], 'Color','w');

    muYb_all = mean(R.all.yb,'omitnan');
    muYp_all = mean(R.all.yp,'omitnan');

    % pooled paired test: bump vs pre (same as annotate_pair)
    good_all = isfinite(R.all.yb) & isfinite(R.all.yp);
    n_all = nnz(good_all);
    if n_all >= 2
        [~, p_all] = ttest(R.all.yb(good_all), R.all.yp(good_all));
    else
        p_all = NaN;
    end
    
    ttl = sprintf('%s, ALL bumps pooled, flies=%d, bumps=%d, paired t-test: N=%d, p=%s', ...
        M.name, numel(unique(R.all.flyUsed)), numel(R.all.xb), n_all, format_p(p_all));

%     ttl = sprintf('%s, ALL bumps pooled, flies=%d, bumps=%d, meanY: bump=%.3g, pre=%.3g', ...
%         M.name, numel(unique(R.all.flyUsed)), numel(R.all.xb), muYb_all, muYp_all);

%     ttl = sprintf('%s, ALL bumps pooled, flies=%d, bumps=%d', ...
%         M.name, numel(unique(R.all.flyUsed)), numel(R.all.xb));
    pH = plot_scatter_with_left_marginals(R.all.xb, R.all.yb, R.all.xp, R.all.yp, ...
        'TitleStr', ttl, 'XLabel','mean bump amp (dF/F)', 'YLabel', M.yLabel, 'NBins', 30);
    if ~isempty(M.forceYLim); ylim(pH.axS, M.forceYLim); end

    if saveFigs
        saveas(fig1, fullfile(outdir, sprintf('RegAnlys_ALL_%s_%d.fig', sanitize_name(M.name), minOn)));
    end

    % ---- 2) SAME bumps split by epoch; ONE FIGURE PER EPOCH (each with marginals) ----
    epochOrder = {'pre','during','post'};
    epochTitle = {'Pre-odor','During odor','Post-odor'};

    for ei = 1:3
        key = epochOrder{ei};
        D = R.phase.(key);
        figE = figure('Name', sprintf('%s - %s', epochTitle{ei}, M.name), 'Color','w');

        muYb = mean(D.yb,'omitnan');
        muYp = mean(D.yp,'omitnan');
        
        ttlE = sprintf('%s, %s, flies=%d, bumps=%d, meanY: bump=%.3g, pre=%.3g', ...
            M.name, epochTitle{ei}, numel(unique(D.flyUsed)), numel(D.xb), muYb, muYp);


%         ttlE = sprintf('%s, %s, flies=%d, bumps=%d', ...
%             M.name, epochTitle{ei}, numel(unique(D.flyUsed)), numel(D.xb));

        pE = plot_scatter_with_left_marginals(D.xb, D.yb, D.xp, D.yp, ...
            'TitleStr', ttlE, 'XLabel','mean bump amp (dF/F)', 'YLabel', M.yLabel, 'NBins', 30);
        if ~isempty(M.forceYLim); ylim(pE.axS, M.forceYLim); end

        if saveFigs
            saveas(figE, fullfile(outdir, sprintf('RegAnlys_EPOCH_%s_%s_%d.fig', key, sanitize_name(M.name), minOn)));
        end
    end

    % ---- 3) Summary categorical scatters: 6 columns (pre/dur/post) x (pre-window vs bump-window) ----
    fig3 = figure('Name', ['Summary 6-column - ' M.name], 'Color','w');
    plot_categorical_phase_pairs( ...
        R.phase.pre.yp,    R.phase.pre.yb, ...
        R.phase.during.yp, R.phase.during.yb, ...
        R.phase.post.yp,   R.phase.post.yb, ...
        'MetricName', M.name, 'YLabel', M.yLabel);
    if ~isempty(M.forceYLim); ylim(M.forceYLim); end

    if saveFigs
        saveas(fig3, fullfile(outdir, sprintf('RegAnlys_SUMMARY_%s_%d.fig', sanitize_name(M.name), minOn)));
    end

    % Per-fly counts printout
    disp('--- Per-fly counts (detected vs plotted) ---');
    disp(R.perFly);
end

%% ====================== COLLECTION ======================

function R = collect_bumps_and_prebumps(use_flies, filename, sdm, odorThr, minOn, maxOff, preN, getYfun, statfun)
% Returns:
%  R.all   : pooled xb/yb/xp/yp across all epochs
%  R.phase : struct with fields pre/during/post, each with xb/yb/xp/yp/flyUsed
%  R.perFly: table with selected trial counts and bump counts

    % outputs
    keys = {'pre','during','post'};
    for kk = 1:numel(keys)
        k = keys{kk};
        R.phase.(k) = struct('xb',[],'yb',[],'xp',[],'yp',[],'flyUsed',[]);
    end
    R.all = struct('xb',[],'yb',[],'xp',[],'yp',[],'flyUsed',[]);

    nF = numel(use_flies);
    nTrialsSelected        = zeros(nF,1);
    nTrialsWithDetected    = zeros(nF,1);
    nDetected_pre          = zeros(nF,1);
    nDetected_during       = zeros(nF,1);
    nDetected_post         = zeros(nF,1);
    nDetected_total        = zeros(nF,1);
    nPlotted_pre           = zeros(nF,1);
    nPlotted_during        = zeros(nF,1);
    nPlotted_post          = zeros(nF,1);
    nPlotted_total         = zeros(nF,1);

    for fi = 1:nF
        flyId = use_flies(fi);

        Sload = load(filename{flyId}, 'syncstruct');
        data  = Sload.syncstruct;

        % --- trial selection + threshold pool (MATCH bumptrig crossing logic) ---
        use_trials = [];
        allamp_thr = [];

        for t = 1:numel(data)
            if ~bumptrig_trial_ok(data(t), odorThr)
                continue;
            end

            nonan = ~isnan(data(t).bump_amp(:,1));
            allamp_thr = [allamp_thr; data(t).bump_amp(nonan,1)]; %#ok<AGROW>
            use_trials(end+1) = t; %#ok<AGROW>
        end

        nTrialsSelected(fi) = numel(use_trials);
        if isempty(use_trials) || isempty(allamp_thr) || all(~isfinite(allamp_thr))
            continue;
        end

        base_thr = mean(allamp_thr,'omitnan') + std(allamp_thr,'omitnan')*sdm;

        trialHadDetected = false(size(use_trials));

        % --- loop trials ---
        for ui = 1:numel(use_trials)
            t = use_trials(ui);

            S = extract_struct_for_reganlys(data(t), odorThr);
            if isempty(S.tt) || isempty(S.dff) || isempty(S.odorOn) || isempty(S.odorOff)
                continue;
            end

            Y = getYfun(S);
            if numel(Y) ~= numel(S.dff)
                continue;
            end

            [~, oni, offi] = bumpfinder2(S.dff, S.tt, maxOff, minOn, base_thr);
            if isempty(oni)
                continue;
            end
            trialHadDetected(ui) = true;

            for bi = 1:numel(oni)
                on  = oni(bi);
                off = offi(bi);
                if off <= on, continue; end

                % epoch by ONSET only
                phaseKey = phase_from_onset(on, S.odorOn, S.odorOff);

                % count detected (even if not plot-valid)
                nDetected_total(fi) = nDetected_total(fi) + 1;
                nDetected_pre(fi)    = nDetected_pre(fi)    + strcmp(phaseKey,'pre');
                nDetected_during(fi) = nDetected_during(fi) + strcmp(phaseKey,'during');
                nDetected_post(fi)   = nDetected_post(fi)   + strcmp(phaseKey,'post');

%                 % pre window validity
%                 preStart = max(1, on - preN);
%                 preEnd   = on - 1;
%                 if preEnd <= preStart
%                     continue;
%                 end

                % pre window validity (MATCHED to bump length)
                bumpLen  = off - on + 1;
                preEnd   = on - 1;
                preStart = preEnd - bumpLen + 1;   % == on - bumpLen
                
                % require full-length pre window
                if preStart < 1 || preEnd <= preStart
                    continue;
                end

                x_b = mean(S.dff(on:off), 'omitnan');
                y_b = statfun(Y(on:off));

                x_p = mean(S.dff(preStart:preEnd), 'omitnan');
                y_p = statfun(Y(preStart:preEnd));

                if ~(isfinite(x_b) && isfinite(y_b) && isfinite(x_p) && isfinite(y_p))
                    continue;
                end

                % append (PHASE)
                R.phase.(phaseKey).xb(end+1,1) = x_b; %#ok<AGROW>
                R.phase.(phaseKey).yb(end+1,1) = y_b; %#ok<AGROW>
                R.phase.(phaseKey).xp(end+1,1) = x_p; %#ok<AGROW>
                R.phase.(phaseKey).yp(end+1,1) = y_p; %#ok<AGROW>
                R.phase.(phaseKey).flyUsed(end+1,1) = flyId; %#ok<AGROW>

                % append (ALL)
                R.all.xb(end+1,1) = x_b; %#ok<AGROW>
                R.all.yb(end+1,1) = y_b; %#ok<AGROW>
                R.all.xp(end+1,1) = x_p; %#ok<AGROW>
                R.all.yp(end+1,1) = y_p; %#ok<AGROW>
                R.all.flyUsed(end+1,1) = flyId; %#ok<AGROW>

                % plotted counts
                nPlotted_total(fi) = nPlotted_total(fi) + 1;
                nPlotted_pre(fi)    = nPlotted_pre(fi)    + strcmp(phaseKey,'pre');
                nPlotted_during(fi) = nPlotted_during(fi) + strcmp(phaseKey,'during');
                nPlotted_post(fi)   = nPlotted_post(fi)   + strcmp(phaseKey,'post');
            end
        end

        nTrialsWithDetected(fi) = nnz(trialHadDetected);
    end

    R.perFly = table( ...
        use_flies(:), nTrialsSelected, nTrialsWithDetected, ...
        nDetected_pre, nDetected_during, nDetected_post, nDetected_total, ...
        nPlotted_pre, nPlotted_during, nPlotted_post, nPlotted_total, ...
        'VariableNames', {'fly','nTrialsSelected','nTrialsWithDetectedBumps', ...
                          'nDetected_pre','nDetected_during','nDetected_post','nDetected_total', ...
                          'nPlotted_pre','nPlotted_during','nPlotted_post','nPlotted_total'});
end

%% ====================== TRIAL / EXTRACTION HELPERS ======================
function ok = bumptrig_trial_ok(trl, odorThr)
% Mirrors bumptrig_offonly_pulse_cntrlcomp trial acceptance.
    ok = false;

    if ~isfield(trl,'closed_loop') || ~trl.closed_loop, return; end
    if ~isfield(trl,'odor_on')     || ~trl.odor_on,     return; end
    if ~isfield(trl,'wind_on')     || ~trl.wind_on,     return; end
    if isfield(trl,'plume_t') && any(~isnan(trl.plume_t)), return; end
    if ~isfield(trl,'bump_amp') || isempty(trl.bump_amp), return; end
    if ~isfield(trl,'odor')     || isempty(trl.odor),     return; end

    nonan = ~isnan(trl.bump_amp(:,1));
    od = trl.odor(nonan);
    if isempty(od) || all(~isfinite(od)), return; end

    [odoroni, odoroffi] = odor_crossings_like_bumptrig(od, odorThr);
    if numel(odoroni) ~= 1, return; end
    if numel(odoroffi) ~= 1, return; end

    if ~(max(od) > 0), return; end
    if ~(od(end) < odorThr), return; end  % bumptrig often expects it to end at 0; this is robust

    ok = true;
end

function [odoroni, odoroffi] = odor_crossings_like_bumptrig(od, thr)
% Replicates bumptrig crossing method:
%   idy = find(od>=thr); idy=idy(idy>1);
%   odoroni  = idy(od(idy-1)<thr);
%   odoroffi = idy(od(idy+1)<thr);
    idy = find(od >= thr);
    idy = idy(idy > 1);

    odoroni = idy(od(idy-1) < thr);

    idy2 = idy(idy < numel(od));
    odoroffi = idy2(od(idy2+1) < thr);
end

function S = extract_struct_for_reganlys(trl, odorThr)
% Extract dff/tt and vectors, and compute odorOn/odorOff indices using bumptrig crossings.
    S = struct('dff',[],'tt',[],'avel',[],'heading',[],'odor',[],'odorOn',[],'odorOff',[]);

    if ~isfield(trl,'calc_ts') || isempty(trl.calc_ts), return; end
    if ~isfield(trl,'bump_amp') || isempty(trl.bump_amp), return; end

    nonan = ~isnan(trl.bump_amp(:,1));
    if ~any(nonan), return; end

    S.dff = trl.bump_amp(nonan,1);
    S.tt  = trl.calc_ts(nonan);

    % angular velocity
    if isfield(trl,'calc_deltaz') && ~isempty(trl.calc_deltaz)
        S.avel = smoothdata(trl.calc_deltaz(nonan) * (180/pi), 'movmean', 150);
    else
        S.avel = nan(size(S.dff));
    end

    % heading (windpos-based)
    if isfield(trl,'calc_windpos') && ~isempty(trl.calc_windpos)
        S.heading = wrapTo180(rad2deg(smoothdata( ...
            unwrap(deg2rad(trl.calc_windpos(nonan))), 'movmean', 100)));
    else
        S.heading = nan(size(S.dff));
    end

    if ~isfield(trl,'odor') || isempty(trl.odor), return; end
    S.odor = trl.odor(nonan);

    [odoroni, odoroffi] = odor_crossings_like_bumptrig(S.odor, odorThr);
    if numel(odoroni) ~= 1 || numel(odoroffi) ~= 1
        return;
    end
    S.odorOn  = odoroni(1);
    S.odorOff = odoroffi(1);
end

function key = phase_from_onset(onIdx, odorOnIdx, odorOffIdx)
    if onIdx < odorOnIdx
        key = 'pre';
    elseif onIdx <= odorOffIdx
        key = 'during';
    else
        key = 'post';
    end
end

function s = sanitize_name(str)
    s = regexprep(str, '[^\w]+', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
end

%% ====================== PLOTTING ======================
function p = plot_scatter_with_left_marginals(xb, yb, xp, yp, varargin)
    ip = inputParser;
    ip.addParameter('TitleStr','', @ischar);
    ip.addParameter('XLabel','',   @ischar);
    ip.addParameter('YLabel','',   @ischar);
    ip.addParameter('ColorBump',[0 0 0], @(c) isnumeric(c)&&numel(c)==3);
    ip.addParameter('ColorPre', [0 0 1], @(c) isnumeric(c)&&numel(c)==3);
    ip.addParameter('NBins', 30, @(n) isnumeric(n)&&isscalar(n)&&n>1);
    ip.addParameter('Norm', 'pdf', @(s) ischar(s) || isstring(s));
    ip.parse(varargin{:});
    P = ip.Results;

    leftMargin = 0.08; gap = 0.02;
    wMarg = 0.10; wScat = 0.62;
    bottom = 0.12; height = 0.78;

    posB = [leftMargin, bottom, wMarg, height];
    posP = [leftMargin + wMarg + gap, bottom, wMarg, height];
    posS = [leftMargin + 2*wMarg + 2*gap, bottom, wScat, height];

    axB = axes('Position', posB); hold(axB,'on'); box(axB,'off'); set(axB,'TickDir','out');
    axP = axes('Position', posP); hold(axP,'on'); box(axP,'off'); set(axP,'TickDir','out');
    axS = axes('Position', posS); hold(axS,'on'); box(axS,'off'); set(axS,'TickDir','out');

    scatter(axS, xb, yb, 12, 'filled', 'MarkerFaceColor', P.ColorBump, 'MarkerFaceAlpha', 0.8);
    scatter(axS, xp, yp, 12, 'filled', 'MarkerFaceColor', P.ColorPre,  'MarkerFaceAlpha', 0.6);
    xlabel(axS, P.XLabel); ylabel(axS, P.YLabel);
    title(axS, P.TitleStr);

    % hist edges
    yAll = [yb(:); yp(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        edges = linspace(0, 1, P.NBins+1);
    else
        edges = linspace(min(yAll), max(yAll), P.NBins+1);
    end

    histogram(axB, yb, edges, 'Orientation','horizontal', 'Normalization',char(P.Norm), ...
        'EdgeColor','none', 'FaceColor',P.ColorBump, 'FaceAlpha',0.45);
    histogram(axP, yp, edges, 'Orientation','horizontal', 'Normalization',char(P.Norm), ...
        'EdgeColor','none', 'FaceColor',P.ColorPre,  'FaceAlpha',0.45);

    linkaxes([axB axP axS], 'y');
    axP.YTickLabel = [];
    axB.XTick = [];
    axP.XTick = [];
    title(axB,'bump','FontWeight','normal');
    title(axP,'pre','FontWeight','normal');

    p = struct('axB',axB,'axP',axP,'axS',axS);
end

function plot_epoch_panel(PH, ttl, ylab)
% Single-axis epoch panel: bump vs pre in the same scatter axis (no marginals).
    if isempty(PH.xb)
        text(0.5,0.5,'(no plotted points)','HorizontalAlignment','center');
        axis off;
        title(ttl);
        return;
    end

    hold on; box off; set(gca,'TickDir','out');
    scatter(PH.xb, PH.yb, 14, 'filled', 'MarkerFaceAlpha', 0.85); % bumps
    scatter(PH.xp, PH.yp, 18, '^', 'filled', 'MarkerFaceAlpha', 0.65); % pre

    xlabel('mean bump amp (dF/F)');
    ylabel(ylab);
    title(sprintf('%s (n=%d)', ttl, numel(PH.xb)));
    legend({'bump','pre'}, 'Location','best');
end

function plot_summary_epoch_scatter(phase, ylab, metricName)
% One axes: bumps only, colored by epoch.
    hold on; box off; set(gca,'TickDir','out');

    epochOrder = {'pre','during','post'};
    epochLabel = {'pre','during','post'};

    for ei = 1:3
        k = epochOrder{ei};
        if isempty(phase.(k).xb), continue; end
        scatter(phase.(k).xb, phase.(k).yb, 14, 'filled', 'MarkerFaceAlpha', 0.85);
    end

    xlabel('mean bump amp (dF/F)');
    ylabel(ylab);
    title(sprintf('%s, bumps only, colored by epoch', metricName));
    legend(epochLabel, 'Location','best');
end

function plot_categorical_phase_pairs(yPre_pre, yPre_bump, yDur_pre, yDur_bump, yPost_pre, yPost_bump, varargin)

    p = inputParser;
    p.addParameter('MetricName','', @ischar);
    p.addParameter('YLabel','', @ischar);
    p.addParameter('ColorPre',[0 0 1], @(c) isnumeric(c)&&numel(c)==3);
    p.addParameter('ColorBump',[0 0 0], @(c) isnumeric(c)&&numel(c)==3);
    p.parse(varargin{:});
    P = p.Results;

    hold on; box on;

    x_pre_pre   = 1; x_pre_bump  = 2;
    x_dur_pre   = 4; x_dur_bump  = 5;
    x_post_pre  = 7; x_post_bump = 8;

    jitter = @(n) (rand(n,1)-0.5)*0.18;

    scatter(x_pre_pre  + jitter(numel(yPre_pre)),  yPre_pre,  14, 'filled', 'MarkerFaceColor', P.ColorPre,  'MarkerFaceAlpha',0.65);
    scatter(x_pre_bump + jitter(numel(yPre_bump)), yPre_bump, 14, 'filled', 'MarkerFaceColor', P.ColorBump, 'MarkerFaceAlpha',0.75);

    scatter(x_dur_pre  + jitter(numel(yDur_pre)),  yDur_pre,  14, 'filled', 'MarkerFaceColor', P.ColorPre,  'MarkerFaceAlpha',0.65);
    scatter(x_dur_bump + jitter(numel(yDur_bump)), yDur_bump, 14, 'filled', 'MarkerFaceColor', P.ColorBump, 'MarkerFaceAlpha',0.75);

    scatter(x_post_pre  + jitter(numel(yPost_pre)),  yPost_pre,  14, 'filled', 'MarkerFaceColor', P.ColorPre,  'MarkerFaceAlpha',0.65);
    scatter(x_post_bump + jitter(numel(yPost_bump)), yPost_bump, 14, 'filled', 'MarkerFaceColor', P.ColorBump, 'MarkerFaceAlpha',0.75);

    draw_mean_bar(x_pre_pre,   yPre_pre);
    draw_mean_bar(x_pre_bump,  yPre_bump);
    draw_mean_bar(x_dur_pre,   yDur_pre);
    draw_mean_bar(x_dur_bump,  yDur_bump);
    draw_mean_bar(x_post_pre,  yPost_pre);
    draw_mean_bar(x_post_bump, yPost_bump);

    set(gca,'XLim',[0 9]);
    set(gca,'XTick',[1.5 4.5 7.5]);
    set(gca,'XTickLabel',{'pre-odor','during odor','post-odor'});
    ylabel(P.YLabel);
    title(sprintf('%s: behavior metric (pre-bump vs bump) by odor phase', P.MetricName), 'Interpreter','none');

    yAll = [yPre_pre(:); yPre_bump(:); yDur_pre(:); yDur_bump(:); yPost_pre(:); yPost_bump(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll), return; end
    yMin = min(yAll); yMax = max(yAll);
    pad  = 0.08 * (yMax - yMin + eps);
    ylim([yMin - 0.05*(yMax-yMin+eps), yMax + 3*pad]);

        % ---- mean labels above each column ----
    yMeanText = yMax + 2.2*pad;   % position for mean labels
    
    label_column_mean(x_pre_pre,   yPre_pre,   yMeanText, 'mean=%.3g', P.ColorPre);
    label_column_mean(x_pre_bump,  yPre_bump,  yMeanText, 'mean=%.3g', P.ColorBump);
    
    label_column_mean(x_dur_pre,   yDur_pre,   yMeanText, 'mean=%.3g', P.ColorPre);
    label_column_mean(x_dur_bump,  yDur_bump,  yMeanText, 'mean=%.3g', P.ColorBump);
    
    label_column_mean(x_post_pre,  yPost_pre,  yMeanText, 'mean=%.3g', P.ColorPre);
    label_column_mean(x_post_bump, yPost_bump, yMeanText, 'mean=%.3g', P.ColorBump);

    annotate_pair(x_pre_pre,  x_pre_bump,  yPre_bump,  yPre_pre,  yMax + 1*pad);
    annotate_pair(x_dur_pre,  x_dur_bump,  yDur_bump,  yDur_pre,  yMax + 1*pad);
    annotate_pair(x_post_pre, x_post_bump, yPost_bump, yPost_pre, yMax + 1*pad);

    legend({'pre-bump','bump'}, 'Location','best');
end

function draw_mean_bar(x, y)
    y = y(isfinite(y));
    if isempty(y), return; end
    m = mean(y,'omitnan');
    hw = 0.22;
    line([x-hw x+hw], [m m], 'Color','r', 'LineWidth', 2);
end

function annotate_pair(x1, x2, yBump, yPre, yText)
    good = isfinite(yBump) & isfinite(yPre);
    n = nnz(good);
    if n >= 2
        [~, p] = ttest(yBump(good), yPre(good));
    else
        p = NaN;
    end
    xc = mean([x1 x2]);
    text(xc, yText, sprintf('N=%d, p=%s', n, format_p(p)), ...
        'HorizontalAlignment','center', 'FontSize',10);

    line([x1 x1 x2 x2], [yText, yText-0.15, yText-0.15, yText], ...
        'Color',[0 0 0], 'LineWidth',1);
end

function s = format_p(p)
    if isnan(p), s = 'NaN';
    elseif p < 1e-3, s = sprintf('%.1e', p);
    else, s = sprintf('%.3g', p);
    end
end

function label_column_mean(x, y, yText, fmt, colorRGB)
    y = y(isfinite(y));
    if isempty(y), return; end
    m = mean(y,'omitnan');
    text(x, yText, sprintf(fmt, m), ...
        'HorizontalAlignment','center', ...
        'FontSize', 10, ...
        'Color', colorRGB);
end
