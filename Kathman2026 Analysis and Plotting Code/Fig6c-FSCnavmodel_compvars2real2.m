% 

%% FSCnavmodel_compvars2real_MEAN_MEDIAN_unitCI_REALGATEONLY_v2.m
% REAL ONLY:
%   1) Use pre-odor frames pre_start_idx:odor_on (odoroni(1))
%   2) Smooth fvel/avel with movmean(150) (inside extract_struct_data)
%   3) Gate trials FIRST (trial-level gate on that pre-odor segment)
%   4) THEN keep samples where fvel >= 0, applied to BOTH fvel and abs(avel)
%
% MODEL:
%   - NO gating
%   - NO sample filtering (uses ALL samples; v can be negative)
%
% Outputs (NO bootstrap):
%   - Mean + Median (pooled raw; visual only)
%   - Mean unit-CI (t-based) for model and real separately + overlap flag
%   - Median unit-CI (order-stat) for model and real separately + overlap flag

clear all
close all

%% ================= USER SETTINGS =================
% Real data path and fly list
headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath)
allflies_filenames             % expects `filename{}` cell array

use_flies  = [21:27, 30:39];   % cwoo

% Segment selection: pre_start_idx : odor_on_index
pre_start_idx = 100;

% Trial-level walking gate for REAL (applied BEFORE any sample filtering)
use_walk_gate   = false;
walk_gate_speed = .1;           % mm/s
walk_gate_frac  = 0.4;        % keep trial if fraction(fvel>walk_gate_speed) > walk_gate_frac

% Angular metric
use_abs_avel = true;           % compare abs(avel)

% Plot bins (visualization only)
% IMPORTANT: include negatives for v so MODEL plotting is not implicitly filtered.
bnw_v = 0.35;  vedg = -10:bnw_v:20;   vcnt = vedg(2:end) - bnw_v/2;
bnw_a = 6;     aedg = 0:bnw_a:325;    acnt = aedg(2:end) - bnw_a/2;

% Model params
ntrials = 50;
samp    = 1050;
T = 3;
avvars = [.7 0.3 0.05 0.3];

% CI level
alpha = 0.05;  % 95%

%% ================= REAL: per-fly pooled samples (REAL gating only) =================
real_by_fly = build_real_by_fly_REALgateOnly(filename, use_flies, pre_start_idx, ...
                                             use_walk_gate, walk_gate_speed, walk_gate_frac, ...
                                             use_abs_avel);

% pooled raw samples (for plotting + pooled markers)
real_v_all = vertcat(real_by_fly.v{:});
real_a_all = vertcat(real_by_fly.a{:});

% visualization histograms (NOT used for stats)
realPv = histcounts(real_v_all, vedg, 'Normalization','probability');
realPa = histcounts(real_a_all, aedg, 'Normalization','probability');

% independent units (REAL = flies): one mean/median per fly (computed on pooled per-fly samples)
real_v_mean_units   = cellfun(@(x) mean(x,'omitnan'),   real_by_fly.v);
real_v_median_units = cellfun(@(x) median(x,'omitnan'), real_by_fly.v);
real_a_mean_units   = cellfun(@(x) mean(x,'omitnan'),   real_by_fly.a);
real_a_median_units = cellfun(@(x) median(x,'omitnan'), real_by_fly.a);

real_v_mean_units   = real_v_mean_units(isfinite(real_v_mean_units));
real_v_median_units = real_v_median_units(isfinite(real_v_median_units));
real_a_mean_units   = real_a_mean_units(isfinite(real_a_mean_units));
real_a_median_units = real_a_median_units(isfinite(real_a_median_units));

%% ================= MODEL: per-trial samples (NO gating, NO filtering) =================
model_v_trials = cell(ntrials,1);
model_a_trials = cell(ntrials,1);

modelPv_trials = zeros(ntrials, numel(vcnt));
modelPa_trials = zeros(ntrials, numel(acnt));

for i = 1:ntrials
    res = FSCnavmodel_15Hz(samp, 'pulse', '', T, avvars, 0);

    vraw = res.v(:);                  % signed mm/s (NO filtering)
    araw = rad2deg(res.a(:));         % signed deg/s

    % MODEL: use all finite samples (no gate, no v>=0 filter)
    v_use = vraw(isfinite(vraw));

    if use_abs_avel
        a_use = abs(araw);
    else
        a_use = araw;
    end
    a_use = a_use(isfinite(a_use));

    model_v_trials{i} = v_use;
    model_a_trials{i} = a_use;

    % visualization histograms per trial (bins include negative v)
    modelPv_trials(i,:) = histcounts(v_use, vedg, 'Normalization','probability');
    modelPa_trials(i,:) = histcounts(a_use, aedg, 'Normalization','probability');
end

% visualization distributions (mean across trials)
modelPv = mean(modelPv_trials,1); modelPv = modelPv / max(eps, sum(modelPv));
modelPa = mean(modelPa_trials,1); modelPa = modelPa / max(eps, sum(modelPa));

% pooled raw samples (for visual markers)
model_v_all = vertcat(model_v_trials{:});
model_a_all = vertcat(model_a_trials{:});

% independent units (MODEL = trials): one mean/median per trial
model_v_mean_units   = cellfun(@(x) mean(x,'omitnan'),   model_v_trials);
model_v_median_units = cellfun(@(x) median(x,'omitnan'), model_v_trials);
model_a_mean_units   = cellfun(@(x) mean(x,'omitnan'),   model_a_trials);
model_a_median_units = cellfun(@(x) median(x,'omitnan'), model_a_trials);

model_v_mean_units   = model_v_mean_units(isfinite(model_v_mean_units));
model_v_median_units = model_v_median_units(isfinite(model_v_median_units));
model_a_mean_units   = model_a_mean_units(isfinite(model_a_mean_units));
model_a_median_units = model_a_median_units(isfinite(model_a_median_units));

%% ================= VISUAL MARKERS (pooled raw): mean + median =================
real_v_mean   = mean(real_v_all,'omitnan');
real_v_median = median(real_v_all,'omitnan');
real_a_mean   = mean(real_a_all,'omitnan');
real_a_median = median(real_a_all,'omitnan');

model_v_mean   = mean(model_v_all,'omitnan');
model_v_median = median(model_v_all,'omitnan');
model_a_mean   = mean(model_a_all,'omitnan');
model_a_median = median(model_a_all,'omitnan');

%% ================= Unit-wise CIs (NO bootstrap; NO diff CI) =================
% Mean unit-CIs
[real_v_mean_est,  real_v_mean_ci]  = mean_ci_t(real_v_mean_units,  alpha);
[model_v_mean_est, model_v_mean_ci] = mean_ci_t(model_v_mean_units, alpha);

[real_a_mean_est,  real_a_mean_ci]  = mean_ci_t(real_a_mean_units,  alpha);
[model_a_mean_est, model_a_mean_ci] = mean_ci_t(model_a_mean_units, alpha);

meanCI_overlap_v = intervals_overlap(real_v_mean_ci, model_v_mean_ci);
meanCI_overlap_a = intervals_overlap(real_a_mean_ci, model_a_mean_ci);

% Median unit-CIs
[real_v_med_est,  real_v_med_ci]  = median_ci_orderstat(real_v_median_units,  alpha);
[model_v_med_est, model_v_med_ci] = median_ci_orderstat(model_v_median_units, alpha);

[real_a_med_est,  real_a_med_ci]  = median_ci_orderstat(real_a_median_units,  alpha);
[model_a_med_est, model_a_med_ci] = median_ci_orderstat(model_a_median_units, alpha);

medianCI_overlap_v = intervals_overlap(real_v_med_ci, model_v_med_ci);
medianCI_overlap_a = intervals_overlap(real_a_med_ci, model_a_med_ci);

%% ================= PLOTS =================
figure(1); clf
set(gcf,'Position',[1 1 1520 550])

% ---- Forward velocity ----
subplot(1,2,1); hold on
plot(vcnt, modelPv, 'k', 'LineWidth', 2)
plot(vcnt, realPv,  'r', 'LineWidth', 2)
set(gca,'YScale','log')
xlabel('forward velocity (mm/s)')
ylabel('log Probability')
legend('model','real','Location','best')

% pooled mean/median markers (visual only)
xline(model_v_median, '-',  'LineWidth', 2, 'Color', 'k');
xline(model_v_mean,   '--', 'LineWidth', 2, 'Color', 'k');
xline(real_v_median,  '-',  'LineWidth', 2, 'Color', 'r');
xline(real_v_mean,    '--', 'LineWidth', 2, 'Color', 'r');

title(sprintf(['Forward v\n' ...
    'MEAN unit-CI model %.3f [%.3f,%.3f] | real %.3f [%.3f,%.3f] | overlap=%d\n' ...
    'MEDIAN unit-CI model %.3f [%.3f,%.3f] | real %.3f [%.3f,%.3f] | overlap=%d\n' ...
    'Pooled mean/med: model %.3f/%.3f | real %.3f/%.3f'], ...
    model_v_mean_est, model_v_mean_ci(1), model_v_mean_ci(2), ...
    real_v_mean_est,  real_v_mean_ci(1),  real_v_mean_ci(2),  meanCI_overlap_v, ...
    model_v_med_est,  model_v_med_ci(1),  model_v_med_ci(2), ...
    real_v_med_est,   real_v_med_ci(1),   real_v_med_ci(2),   medianCI_overlap_v, ...
    model_v_mean, model_v_median, real_v_mean, real_v_median))

% ---- Angular |ω| ----
subplot(1,2,2); hold on
plot(acnt, modelPa, 'k', 'LineWidth', 2)
plot(acnt, realPa,  'r', 'LineWidth', 2)
set(gca,'YScale','log')
xlabel('angular velocity magnitude |deg/s|')
ylabel('log Probability')
legend('model','real','Location','best')

% pooled mean/median markers (visual only)
xline(model_a_median, '-',  'LineWidth', 2, 'Color', 'k');
xline(model_a_mean,   '--', 'LineWidth', 2, 'Color', 'k');
xline(real_a_median,  '-',  'LineWidth', 2, 'Color', 'r');
xline(real_a_mean,    '--', 'LineWidth', 2, 'Color', 'r');

title(sprintf(['Angular |ω|\n' ...
    'MEAN unit-CI model %.3f [%.3f,%.3f] | real %.3f [%.3f,%.3f] | overlap=%d\n' ...
    'MEDIAN unit-CI model %.3f [%.3f,%.3f] | real %.3f [%.3f,%.3f] | overlap=%d\n' ...
    'Pooled mean/med: model %.3f/%.3f | real %.3f/%.3f'], ...
    model_a_mean_est, model_a_mean_ci(1), model_a_mean_ci(2), ...
    real_a_mean_est,  real_a_mean_ci(1),  real_a_mean_ci(2),  meanCI_overlap_a, ...
    model_a_med_est,  model_a_med_ci(1),  model_a_med_ci(2), ...
    real_a_med_est,   real_a_med_ci(1),   real_a_med_ci(2),   medianCI_overlap_a, ...
    model_a_mean, model_a_median, real_a_mean, real_a_median))

%% ================= Console summary =================
disp('--- SUMMARY ---')
fprintf('Units: REAL=%d flies | MODEL=%d trials\n', numel(real_v_mean_units), numel(model_v_mean_units));
fprintf('REAL: frames pre_start_idx:odoroni(1), with movmean(150) smoothing for fvel/avel\n');
fprintf('REAL: trial gate use=%d (frac(fvel>%.2f) > %.2f) computed on UNFILTERED pre-odor segment\n', ...
    use_walk_gate, walk_gate_speed, walk_gate_frac);
fprintf('REAL: after gate, sample mask fvel>=0 applied to BOTH fvel and avel (then abs(avel))\n');
fprintf('MODEL: no gate, no sample mask (uses all finite v and a samples)\n');

%% ===================== LOCAL FUNCTIONS =====================

function tf = intervals_overlap(ciA, ciB)
    a1 = ciA(1); a2 = ciA(2);
    b1 = ciB(1); b2 = ciB(2);
    tf = ~(a2 < b1 || b2 < a1);
end

function [mu, ci] = mean_ci_t(x, alpha)
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);
    mu = mean(x);
    if n < 2
        ci = [NaN NaN];
        return
    end
    se = std(x, 0) / sqrt(n);
    tcrit = tinv(1 - alpha/2, n-1);
    ci = [mu - tcrit*se, mu + tcrit*se];
end

function [med, ci] = median_ci_orderstat(x, alpha)
    x = x(:);
    x = x(isfinite(x));
    x = sort(x);
    n = numel(x);
    med = median(x, 'omitnan');
    if n < 3
        ci = [NaN NaN];
        return
    end
    k = binoinv(alpha/2, n, 0.5) + 1;
    lo = max(1, k);
    hi = min(n, n - k + 1);
    if lo > hi
        ci = [x(1) x(end)];
    else
        ci = [x(lo) x(hi)];
    end
end

function real_by_fly = build_real_by_fly_REALgateOnly(filename, use_flies, pre_start_idx, ...
                                                      use_walk_gate, walk_gate_speed, walk_gate_frac, ...
                                                      use_abs_avel)
% REAL ONLY:
% - Use frames pre_start_idx : odoroni(1)
% - Gate trials FIRST on that segment (no sample filtering)
% - Then keep samples where fvel >= 0 (mask applied to BOTH fvel and avel)
% - Then abs(avel) if requested

    real_by_fly.v = cell(numel(use_flies),1);
    real_by_fly.a = cell(numel(use_flies),1);

    total_trials = 0;
    kept_trials  = 0;

    for fi = 1:numel(use_flies)
        fidx = use_flies(fi);
        load(filename{fidx}) %#ok<LOAD>
        data = syncstruct;

        use_trials = find(arrayfun(@(t) t.closed_loop && ...
            sum(diff(t.odor) == 1) == 1 && sum(diff(t.odor) == -1) == 1, data));

        vcat = [];
        acat = [];

        for t = use_trials
            total_trials = total_trials + 1;

            [~, fvel, avel, ~, odoroni, ~, ~, ~] = extract_struct_data(data(t));

            if isempty(odoroni)
                continue
            end

            od = odoroni(1);
            if od <= pre_start_idx
                continue
            end

            seg = pre_start_idx : od;

            % --- TRIAL-LEVEL GATE FIRST (no sample filtering here) ---
            if use_walk_gate
                frac_fast = nanmean(fvel(seg) > walk_gate_speed);
                if frac_fast <= walk_gate_frac
                    continue
                end
            end

            % --- Sample mask AFTER gate: keep fvel>=0 for BOTH metrics ---
            fseg = fvel(seg);
            aseg = avel(seg);

            mask = isfinite(fseg) & isfinite(aseg) & (fseg >= 0);
            fseg = fseg(mask);
            aseg = aseg(mask);

            if isempty(fseg) || isempty(aseg)
                continue
            end

            if use_abs_avel
                aseg = abs(aseg);
            end

            kept_trials = kept_trials + 1;

            vcat = [vcat; fseg]; %#ok<AGROW>
            acat = [acat; aseg]; %#ok<AGROW>
        end

        real_by_fly.v{fi} = vcat;
        real_by_fly.a{fi} = acat;
    end

    fprintf('REAL gate+filter summary: kept_trials=%d / total_candidate_trials=%d\n', kept_trials, total_trials);
end

function [tt, fvel, avel, speed, oncrossi, offcrossi, odor, fps] = extract_struct_data(trl)
% Smooth exactly like your behavior script (movmean 150 for both)
    nonan = ~isnan(trl.calc_ts(:,1));
    fps = cat(1,trl.fps);

    tt = trl.calc_ts(nonan);

    avel = trl.calc_deltaz(nonan) * (180/pi);
    avel = smoothdata(avel, 'movmean', 150);

    fvel = trl.calc_deltapitch(nonan)*9.52/2;
    fvel = smoothdata(fvel, 'movmean', 150);

    odor = trl.odor(nonan);

    oncrossi  = find(diff(odor) == 1);
    offcrossi = find(diff(odor) == -1);

    if odor(1) == 1
        oncrossi = [1; oncrossi];
    end
    if odor(end) == 1
        offcrossi = [offcrossi; length(odor)];
    end

    speed = trl.calc_speed(nonan); %#ok<NASGU>
end

%