%% heatmap_pbump_alltrials_cwoo_stimaligned_withMeanPID
% Heatmap of pbump aligned to odor onset (t=0), plus mean pbump trace below.
% Adds mean PID trace (t,mnpid) loaded from meanpid.mat ABOVE the heatmap.

clear; close all;

%% ---------------- USER SETTINGS ----------------
use_flies = [21:26, 30, 32:35, 37:39];   % <-- edit
sdm  = 0.5;                   % bump threshold SD multiplier (per fly)
minOn = 3;                    % bumpfinder2 minon (sec)
maxOff = 3;                   % bumpfinder2 maxoff (sec)

tStart = -20;                 % seconds relative to odor onset
tEnd   =  35;
dtBin  = 0.05;

keepOnlyTrialsWithBumps = true;

meanpid_path = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/meanpid.mat';

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath);
allflies_filenames            % defines filename{flyId}

%% ---------------- LOAD MEAN PID ----------------
S = load(meanpid_path, 't', 'mnpid');
t_pid = S.t(:);
mnpid = S.mnpid(:);

%% ---------------- COMMON TIME AXIS (-20..35) ----------------
tEdges  = (tStart:dtBin:tEnd)';        % edges
tCommon = tEdges(1:end-1) + dtBin/2;   % centers for plotting
nBins   = numel(tCommon);

pbumpMat  = [];     % trials x nBins

nTotalTrialsChecked = 0;
nTrialsKept = 0;

%% ---------------- LOOP FLIES ----------------
for ff = 1:numel(use_flies)
    flyId = use_flies(ff);

    Sload = load(filename{flyId}, 'syncstruct');
    data = Sload.syncstruct;

    % Single odor pulse trials (exactly one on and one off)
    use_trials = find(arrayfun(@(t) t.closed_loop && ...
        sum(diff(t.odor) == 1) == 1 && sum(diff(t.odor) == -1) == 1, data));
    if isempty(use_trials), continue; end

    % Fly-specific bump threshold from all bump_amp
    hasBA = arrayfun(@(t) isfield(t,'bump_amp') && ~isempty(t.bump_amp) && any(~isnan(t.bump_amp)), data);
    if ~any(hasBA), continue; end
    all_dff = cat(1, data(hasBA).bump_amp);
    all_dff = all_dff(~isnan(all_dff));
    if isempty(all_dff), continue; end
    base_thr = nanmean(all_dff) + nanstd(all_dff)*sdm;

    % Loop trials
    for t = use_trials
        nTotalTrialsChecked = nTotalTrialsChecked + 1;

        [dff, tt, odor] = extract_trial_basic(data(t));
        if isempty(dff) || isempty(tt) || isempty(odor), continue; end

        odorOnIdx = find_odor_on(odor);
        if isempty(odorOnIdx), continue; end

        % pbump from bumpfinder2 (first output preferred)
        [pbump, oni, offi] = bumpfinder2(dff, tt, maxOff, minOn, base_thr); %#ok<ASGLU>
        if isempty(pbump) || numel(pbump) ~= numel(tt)
            pbump = zeros(size(tt));
            for k = 1:numel(oni)
                a = max(1, oni(k));
                b = min(numel(pbump), offi(k));
                if b >= a
                    pbump(a:b) = 1;
                end
            end
        end
        pbump = double(pbump(:) > 0);

        if keepOnlyTrialsWithBumps && ~any(pbump)
            continue;
        end

        % Align time to odor onset
        t0 = tt(odorOnIdx);
        trel = tt(:) - t0;

        % Initialize row as zeros everywhere (fill missing with 0)
        row = zeros(1, nBins);

        % Assign samples into bins; bin becomes 1 if ANY sample has pbump=1
        idx = floor((trel - tEdges(1)) / dtBin) + 1;  % 1-based bin index
        valid = idx >= 1 & idx <= nBins & isfinite(idx);

        idx = idx(valid);
        vals = pbump(valid);

        for s = 1:numel(idx)
            b = idx(s);
            if vals(s) == 1
                row(b) = 1;
            end
        end

        pbumpMat(end+1, :) = row; %#ok<AGROW>
        nTrialsKept = nTrialsKept + 1;
    end
end

fprintf('Checked %d single-pulse trials across %d flies.\n', nTotalTrialsChecked, numel(use_flies));
fprintf('Kept %d trials (keepOnlyTrialsWithBumps=%d).\n', nTrialsKept, keepOnlyTrialsWithBumps);

if isempty(pbumpMat)
    error('No trials kept. Check filters, thresholds, or use_flies.');
end

%% ---------------- MEAN TRACE ----------------
meanPb = mean(pbumpMat, 1);

%% ---------------- PLOT (PID top, heatmap middle, mean bottom) ----------------
figure('Color','w'); clf;
set(gcf, 'Position', [100 100 1200 900]);

ax0 = subplot(3,1,1);
plot(t_pid, mnpid, 'r', 'LineWidth', 2);
hold on; xline(0, 'k-', 'LineWidth', 1.2); hold off;
ylabel('Mean PID');
title(sprintf('Mean PID + pbump heatmap aligned to odor onset,  flies=%d, trials=%d', ...
    numel(use_flies), size(pbumpMat,1)));
set(ax0,'TickDir','out'); box off;
xlim([tStart tEnd]);

ax1 = subplot(3,1,2);
imagesc(tCommon, 1:size(pbumpMat,1), pbumpMat);
colormap(ax1, gray);
caxis([0 1]);
ylabel('Trial');
hold on; xline(0, 'r-', 'LineWidth', 1.5); hold off;
set(ax1, 'TickDir','out');
xlim([tCommon(1) tCommon(end)]);
box off

ax2 = subplot(3,1,3);
plot(tCommon, meanPb, 'k', 'LineWidth', 2);
hold on; xline(0, 'r-', 'LineWidth', 1.5); hold off;
xlabel('Time from odor onset (s)');
ylabel('Mean pbump');
set(ax2, 'TickDir','out');
xlim([tCommon(1) tCommon(end)]);
ylim([0 1]);
box off

linkaxes([ax0 ax1 ax2], 'x');

%% ===================== HELPERS =====================

function [dff, tt, odor] = extract_trial_basic(trl)
    dff = []; tt = []; odor = [];
    if ~isfield(trl,'calc_ts') || isempty(trl.calc_ts), return; end
    nonan = ~isnan(trl.calc_ts(:,1));
    tt   = trl.calc_ts(nonan);
    odor = trl.odor(nonan);
    if isfield(trl,'bump_amp') && ~isempty(trl.bump_amp)
        dff = trl.bump_amp(nonan,1);
    end
end

function odorOnIdx = find_odor_on(odor)
    odor = odor(:);
    oncrossi = find(diff(odor)==1);
    if ~isempty(odor) && odor(1)==1
        oncrossi = [1; oncrossi(:)];
    else
        oncrossi = oncrossi(:);
    end
    if isempty(oncrossi)
        odorOnIdx = [];
    else
        odorOnIdx = oncrossi(1);
    end
end
