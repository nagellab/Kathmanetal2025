function fig_no = plottrials_flex(data, savepath, f, which_trials, ...
                                  dff_yaxis, dff_coloraxis, base_thr, ...
                                  add2txt, print_flag, file_type, ...
                                  trajcmap, fig_no)
% plottrials_flex
%
% fig_no = plottrials_flex(data, savepath, f, which_trials, ...
%                          dff_yaxis, dff_coloraxis, base_thr, ...
%                          add2txt, print_flag, file_type, ...
%                          trajcmap, fig_no)
%
% Inputs:
%   data           : struct array with fields:
%                    .bump_amp, .dff_ts, .calc_ts, .calc_heading, .calc_windpos,
%                    .calc_deltaz, .calc_deltapitch, .calc_speed, .odor, .wind,
%                    .plume_x, .plume_y, .dff_green, .filename, .thorname, .fps,
%                    .closed_loop, .odor_on, .wind_on
%   savepath       : folder path for output (only used if print_flag ~= 0)
%   f              : fly index (used in titles and filenames)
%   which_trials   : vector of trial indices (if empty → all trials)
%   dff_yaxis      : [ymin ymax] for dF/F trace axis (if [] → default [-0.1 0.8])
%   dff_coloraxis  : [cmin cmax] for dF/F colormap (if [] → default [0 0.4])
%   base_thr       : baseline threshold override passed to bumpfinder2
%   add2txt        : extra text appended to figure title and filename
%   print_flag     : 0 = do not save; 1 = save using print_fig
%   file_type      : 'png', 'eps', 'fig' or char matrix of these (rows)
%   trajcmap       : colormap for trajectory (e.g. jet/parula/custom)
%   fig_no         : starting figure number (returned updated)
%
% Output:
%   fig_no         : updated figure number (if print_flag == 0)

    % --- configuration flags ---
    no_dffcolor = 0;   % if 1: trajectory not colored by dF/F
    surf        = 0;   % if 1: use surface plot for trajectory coloring
    bumppos     = 0;   % if 1: overlay bump positions on dF/F image

    loop_labels = {'Wind Open Loop', 'Wind Closed Loop'};

    % --- handle defaults / optional inputs ---
    if nargin < 4 || isempty(which_trials)
        which_trials = 1:numel(data);
    end

    if nargin < 5 || isempty(dff_yaxis)
        yl_ca = [-0.1 0.8];
    else
        yl_ca = dff_yaxis;
    end

    if nargin < 6 || isempty(dff_coloraxis)
        cl_ca = [0 0.4];
    else
        cl_ca = dff_coloraxis;
    end

    if nargin < 7 || isempty(base_thr)
        base_thr = 0.3;    % sensible default if not provided
    end

    if nargin < 8 || isempty(add2txt)
        add2txt = '';
    end

    if nargin < 9 || isempty(print_flag)
        print_flag = 0;
    end

    if nargin < 10 || isempty(file_type)
        file_type = 'png';
    end

    if nargin < 11 || isempty(trajcmap)
        trajcmap = parula; % default trajectory colormap
    end

    if nargin < 12 || isempty(fig_no)
        fig_no = 1;
    end

    smth_fun = 'movmean';
    columns  = 1:8;     % columns for dff_green

    % --- determine experiment duration: longest trial dff_ts ---
    exp_dur = 0;
    for t = 1:numel(data)
        if ~isempty(data(t).dff_ts)
            exp_dur = max(exp_dur, data(t).dff_ts(end));
        end
    end
    exp_dur = exp_dur + 1;

    % subplot layout
    m      = 5;  % rows
    n      = 4;  % cols
    layers = 0;  % single-layer by default

    % rows/columns for timeseries subplots
    s  = [1,2;  5,6;  9,10; 13,14; 17,18];
    st = [3:4, 7:8, 11:12, 15:16, 19:20];  % trajectory panel indices

    % If FPS is stored per trial, use median over data as global estimate
    fps_all = cat(1, data.fps);
    fps     = median(fps_all(~isnan(fps_all)));

    % --- main loop over trials ---
    for t = which_trials

        % skip trials with empty bump_amp
        if isempty(data(t).bump_amp)
            continue
        end

        % find odor onset crossings (used by some commented filters)
        nonan  = ~isnan(data(t).bump_amp(:,1));
        odor   = data(t).odor(nonan, 1);
        ampthresh = 0.5;
        idyl   = odor >= ampthresh;
        idy    = find(idyl);
        idy    = idy(idy > 1);
        oncrossi = idy(odor(idy-1) < ampthresh); %#ok<NASGU>

        % --- title fragments ---
        tt = data(t).calc_ts;
        sp = strfind(data(t).filename, '-'); %#ok<NASGU>
        us = strfind(data(t).filename, '_'); %#ok<NASGU>

        trli   = strfind(data(t).thorname, '_corr');
        trlstr = data(t).thorname(trli-3:trli-1);

        figure(fig_no); clf
        set(gcf, 'Position', [100 100 1200 500]);

        sgtitle(sprintf('F%d, t%d, Fly%d_%s, Trial %s - %s', ...
            f, t, f, data(t).filename(sp(1)-6:sp(1)-1), trlstr, ...
            loop_labels{data(t).closed_loop + 1}), ...
            'Interpreter', 'none');

        % =========================
        % Row 1 – heading / windpos
        % =========================
        if ~data(t).closed_loop
            heading = data(t).calc_heading;
        else
            heading = data(t).calc_windpos;
        end

        heading = wrapTo180(rad2deg( ...
            smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));

        % also define an unwrapped, smoothed heading if needed later
        uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 100);

        subplot(m, n, s(1,:)); hold on;
        ylabel('Heading (deg)');

        plot(data(t).calc_ts, heading, 'k');
        yl = [-180 180];
        yline(0, 'k:', 'LineWidth', 0.1);

        plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, 0.06, 'k');
        plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, 0.5, 'm');

        ylim(yl);
        xlim([3 exp_dur]);

        % =========================
        % Row 2 – angular velocity
        % =========================
        subplot(m, n, s(2,:)); hold on;
        ylabel('aVel (rad/s)');

        avel = smoothdata(data(t).calc_deltaz, smth_fun, 80);
        plot(data(t).calc_ts, avel, 'k');

        yl = [-7 7];
        yline(0, 'k:', 'LineWidth', 0.1);

        plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, 0.06, 'k');
        plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, 0.5, 'm');

        ylim(yl);
        xlim([3 exp_dur]);

        % =========================
        % Row 3 – forward velocity
        % =========================
        subplot(m, n, s(3,:)); hold on;
        ylabel('fVel (mm/s)');

        fvel = data(t).calc_deltapitch * 9.52 / 2;
        fvel = smoothdata(fvel, smth_fun, 80);

        plot(data(t).calc_ts, fvel, 'k');

        yl = [-8 17];
        yline(0, 'k:', 'LineWidth', 0.1);

        plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, 0.06, 'k');
        plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, 0.5, 'm');

        ylim(yl);
        xlim([3 exp_dur]);

        % ==================================
        % Row 4 – bump amplitude (dF/F trace)
        % ==================================
        if ~layers
            subplot(m, n, s(4,:)); hold on;
            ylabel('dF/F Green');

            dff = data(t).bump_amp(:,1);
            dffi = ~isnan(dff);

            plot(data(t).calc_ts, dff, 'k');

            plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl_ca, 0.06, 'k');
            plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl_ca, 0.5, 'm');

            ylim(yl_ca);
            xlim([3 exp_dur]);

            % ==================================
            % Row 5 – dF/F image (columns x time)
            % ==================================
            subplot(m, n, s(5,:)); hold on;
            ylabel('dF/F Green');

            sigma_a = 0.5;  % spatial smoothing across columns
            sigma_b = 0.5;  % temporal smoothing across time

            im   = data(t).dff_green(columns,:);
            smim = imgaussfilt(im, [sigma_a, sigma_b]);

            imagesc(data(t).dff_ts, 1:numel(columns), smim, cl_ca);
            colormap('gray');

            ylim([0.5 numel(columns)+0.5]);
            xlim([3 exp_dur]);
            colorbar('east');

            % Optional bump position markers on the image
            if bumppos
                nonan   = ~isnan(data(t).bump_amp(:,1));
                bumpdff = data(t).bump_amp(nonan, 1);
                tt_b    = data(t).calc_ts(nonan);

                bump = bumpfinder2(bumpdff, tt_b, 5, 3, base_thr);
                bumpi = find(bump);

                col = data(t).bump_col(nonan); %#ok<NASGU>
                com = data(t).bump_com(nonan); %#ok<NASGU>
                pos = smoothdata(data(t).bump_col(nonan), 'movmean', 200);

                scatter(tt_b(bumpi), pos(bumpi), 4, 'g', 'filled');
                scatter(tt_b(bumpi), heading(bumpi)*8/360 + 4.5, 2, 'r', 'filled');
            end

        else
            % multi-layer case (looking at nonstand ROIs, but not actively used)
            for l = 1:2
                % dF/F trace per layer
                subplot(m, n, s(4 + 2*(l-1), :)); hold on;
                ylabel('dF/F Green');

                dff  = data(t).bump_amp(:,l);
                dffi = ~isnan(dff);

                plot(data(t).calc_ts, dff, 'k');
                plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl_ca, 0.06, 'k');
                plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl_ca, 0.5, 'm');

                ylim(yl_ca);
                xlim([3 exp_dur]);

                % dF/F image per layer
                subplot(m, n, s(5 + 2*(l-1), :)); hold on;
                ylabel('dF/F Green');

                im = data(t).dff_diff(columns(l,:), :);
                imagesc(data(t).dff_ts, 1:numel(columns(1,:)), im, cl_ca);

                colormap('gray');
                ylim([0.5 numel(columns(1,:))+0.5]);
                xlim([3 exp_dur]);
                colorbar('east');
            end
        end

        % =========================
        % Trajectory panel (big plot)
        % =========================
        dff  = data(t).bump_amp(:,1);
        dffi = ~isnan(dff);

        dffts   = data(t).calc_ts(dffi);
        dffodor = data(t).odor(dffi);
        dffwind = data(t).wind(dffi);

        o = dffts(dffodor > 0.5); % odor times
        w = dffts(dffwind > 0.5); % wind times

        if isempty(o), o = [-1 -1]; end
        if isempty(w), w = [-1 -1]; end

        odor = data(t).odor(dffi);
        wind = data(t).wind(dffi); %#ok<NASGU>

        yd = [-75 410];
        xd = [-range(yd)/2, range(yd)/2];

        % smoothed dF/F for trajectory coloring
        yd_dff = smoothdata(dff(dffi), smth_fun, 200);

        % calc XY, plume or non-plume 
        if ~all(isnan(data(t).plume_x))
            % plume trials: use plume coordinates
            x = data(t).plume_x(dffi);
            y = data(t).plume_y(dffi);

            % re-center and flip to match plotting frame
            y = (y-406) * -1;
            x = (x-108) * -1;

            % optional: plume background (turned off to avoid repeated h5 read)
            % plume = h5read('/Users/.../10302017_10cms_bounded_2.h5', ...
            %                '/dataset2', [1 1 1], [216 406 1]);
        else
            % non-plume trials: integrate velocity to position by heading
            heading_i = heading(dffi);
            speed_i   = data(t).calc_speed(dffi);

            x = -cumsum(speed_i .* sin(heading_i*pi/180)) / fps;
            y =  cumsum(speed_i .* cos(heading_i*pi/180)) / fps;

%             % non-plume trials: integrate velocity to position by 2dim vel
%             dt   = median(diff(data(t).calc_ts));
%             
%             dRotX = data(t).calc_deltapitch(dffi);        % rad per frame, e.g. forward/back
%             dRotY = data(t).calc_deltaroll(dffi);        % rad per frame, e.g. lateral
%             
%             v_fwd = (dRotX * 9.52) / dt;           % mm/s
%             v_lat = (dRotY * 9.52) / dt;           % mm/s
% 
%             theta  = deg2rad(heading(dffi));              % world heading
%             v_fwd  = data(t).fwd_vel(dffi);               % mm/s, body forward
%             v_lat  = data(t).lat_vel(dffi);               % mm/s, body lateral
%             
%             vx = v_fwd .* cos(theta) - v_lat .* sin(theta);  % world X velocity
%             vy = v_fwd .* sin(theta) + v_lat .* cos(theta);  % world Y velocity
%             
%             x = cumsum(vx * dt);  % integrate to position
%             y = cumsum(vy * dt);



        end

        % --- create the trajectory subplot with two overlaid axes if plume ---
        subplot(m, n, st); hold on;
        set(gca, 'YTickLabel', [], 'XTickLabel', []);

        if ~all(isnan(data(t).plume_x))
            % two axes trick: ax1 (background), ax2 (trajectory)
            ax1 = axes; hold on;
            subplot(m, n, st, ax1);

            colormap(ax1, 'gray');
            xlim(ax1, xd);
            ylim(ax1, yd);

            % rectangular plume region outline
            plot([-107, 108, 108, -107, -107], [0, 0, 405, 405, 0], ...
                 'k-', 'LineWidth', 0.5);
            plot([-107, 0, 108], [0, 405, 0], 'k-', 'LineWidth', 0.5);

            ax2 = axes; hold on;
            subplot(m, n, st, ax2);

        else
            ax1 = gca;
            ax2 = ax1;
        end

        % --- trajectory drawing on ax2 ---
        if no_dffcolor
            plot(ax2, x, y, 'r', 'LineWidth', 2);
            if ~all(isnan(data(t).plume_x))
                scatter(ax2, x(odor > 0.5), y(odor > 0.5), 5, 'w', 'filled');
            end
        else
            if surf
                z  = zeros(size(x));
                yq = yd_dff;
                surface(ax2, [x x], [y y], [z z], [yq yq], ...
                        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
                colorbar(ax2);
                colormap(ax2, 'jet');
                caxis(ax2, cl_ca);
            else                                            
                %%% typical traj plotting method %%%
                scatter(ax2, x, y, 150, yd_dff, 'filled');



                if ~all(isnan(data(t).plume_x))
                    odor_idx = find(odor);
                    scatter(ax2, x(odor_idx)-12, y(odor_idx), 8, 'm', 'filled');
                end
                colormap(ax2, trajcmap);
                caxis(ax2, cl_ca);
            end
        end

        % entry position
        plot(ax2, 0, 406, 'xr', 'MarkerSize', 5);

        % global trajectory limits
        xlim(ax2, [-375 375]);
        ylim(ax2, [-200 550]);
        colorbar(ax2);

        if ~all(isnan(data(t).plume_x))
            xlim(ax2, [-240 240]);
            ylim(ax2, [-175 410]);

            % link background and trajectory axes
            linkaxes([ax1, ax2]);

            % hide top axes (trajectory overlay)
            ax2.Visible = 'off';
            ax2.XTick   = [];
            ax2.YTick   = [];
            set([ax1, ax2], 'Position', [0.543 0.11 0.355 0.815]);
        else
            % mark start
            plot(ax2, x(1), y(1), '*m', 'LineWidth', 3);

            % try to mark odor/wind on/off on trajectory
            [~, ~, idx] = unique(round(abs(dffts - o(1))));
            oo = find(idx == 1);
            if ~isempty(oo) && oo(1) <= numel(x)
                plot(ax2, x(oo(1)), y(oo(1)), 'xk', 'LineWidth', 2);
            end

            [~, ~, idx] = unique(round(abs(dffts - o(end))));
            oo = find(idx == 1);
            if ~isempty(oo) && oo(1) <= numel(x)
                plot(ax2, x(oo(1)), y(oo(1)), 'xk', 'LineWidth', 2);
            end

            [~, ~, idx] = unique(round(abs(dffts - w(1))));
            ww = find(idx == 1);
            if ~isempty(ww) && ww(1) <= numel(x)
                plot(ax2, x(ww(1)), y(ww(1)), 'xk', 'LineWidth', 1);
            end

            [~, ~, idx] = unique(round(abs(dffts - w(end))));
            ww = find(idx == 1);
            if ~isempty(ww) && ww(1) <= numel(x)
                plot(ax2, x(ww(1)), y(ww(1)), 'xk', 'LineWidth', 1);
            end
        end

        % --- save / advance figure number ---
        fig_no = print_fig(print_flag, gcf, ...
            [savepath 'fly' num2str(f) '_t' trlstr add2txt], ...
            file_type, fig_no);

    end % for t

end


% =========================
% Helper: stimfills
% =========================
function plot_stimfills(ax_handle, ts, binarized, ylim, facealpha, colorstr)
%plot_stimfills(gca, stim_times(time series of stim signal, not just ons/offs), stim_values(zeros/ones of when stim is off/on), y range of boxes, eg. [0 3], transparency, eg. .06)

    ts = ts(~isnan(ts));
    binarized = binarized(~isnan(binarized));
    binarized = make_row(binarized);
    on = find(diff(binarized)>0);
    off = find(diff(binarized)<0);
    if binarized(1) 
        on = [1 on]; end
    if binarized(end)
        off = [off length(binarized)]; end
    for s = 1:length(on)
        g = fill(ax_handle, [ts(on(s)), ts(off(s)), ts(off(s)), ts(on(s))], [ylim(2) ylim(2) ylim(1) ylim(1)], colorstr);
        set( g, 'edgecolor', 'none', 'facealpha', facealpha );
    end
end

% =========================
% Helper: print_fig
% =========================
function fig_no = print_fig(print_flag, fig_handle, fig_filename, file_type, fig_no)
% print_fig(print_flag, fig_handle, fig_filename, file_type, fig_no)
% print_flag: 1 → save, 0 → do not save (increment fig_no)
% file_type : 'eps', 'png', 'fig' or multiple as rows (e.g. ['eps'; 'png'])

    if print_flag
        for i = 1:size(file_type, 1)
            this_type = strtrim(file_type(i,:));

            if strcmp(this_type, 'eps')
                print_type = 'epsc';
                print_mode = '-painters';
            elseif strcmp(this_type, 'png')
                print_type = 'png';
                print_mode = '';
            elseif strcmp(this_type, 'fig')
                savefig(fig_handle, fig_filename);
                continue
            else
                error('Undefined filetype: %s', this_type);
            end

            print(fig_handle, ['-d' print_type], print_mode, ...
                  [fig_filename '.' this_type]);
        end
    else
        fig_no = fig_no + 1;
    end
end

% =========================
% Helper: bumpfinder2 (old/trial version; needs updating)
% =========================
function [bump, oni, offi, thresh] = bumpfinder2(dff, tt, maxoff, minon, base_thr)
% bumpfinder2
%
% [bump, oni, offi, thresh] = bumpfinder2(dff, tt, maxoff, minon, base_thr)
%
% dff      : non-NaN Ca trace
% tt       : time vector (same length as dff)
% maxoff   : max gap (sec) between bump segments to fill
% minon    : minimum bump duration (sec)
% base_thr : baseline threshold used for some edge cases
%
% bump     : logical vector (same length as dff) marking bump epochs
% oni      : sample indices of bump onsets
% offi     : sample indices of bump offsets
% thresh   : amplitude threshold used on smoothed dff

    fs = 1 / median(diff(tt));                % sampling frequency
    sdff = smoothdata(dff, 'movmean', round(fs));  % 1 s smoothing

    % --- segment trace into Ca "states" using change points ---
    ipt = findchangepts(sdff);                % positions of change points
    ipt = [1; ipt; numel(sdff)];

    bases = zeros(numel(ipt)-1, 1);
    for i = 2:numel(ipt)
        bases(i-1) = mean(sdff(ipt(i-1):ipt(i)));
    end

    % --- threshold between Ca states ---
    if numel(bases) > 1
        [~, c] = kmeans(bases', 2);
        thresh = min(c) + abs(diff(c)) / 2;
    else
        thresh = mean(sdff) + 2.5 * std(sdff);
    end

    % floor threshold
    if thresh < 0.2
        thresh = 0.2;
    end

    if numel(bases) == 2 && abs(diff(bases)) < 2.5 * std(sdff)
        % low separation between states → use provided baseline threshold
        thresh = base_thr;
    end

    % initial suprathreshold samples
    bumpi = find(sdff > thresh);

    % fill gaps smaller than maxoff (in seconds)
    bump = zeros(numel(sdff), 1);
    bump(bumpi) = 1;

    sample_gap_max = round(maxoff * fs);
    gaps = diff(bumpi);
    idx  = find(gaps > 1 & gaps < sample_gap_max);

    for g = 1:numel(idx)
        gap_idx = bumpi(idx(g))+1 : bumpi(idx(g)+1)-1;
        bump(gap_idx) = 1;
    end

    bumpi = find(bump);

    % --- find on/off crossings in logical bump vector ---
    if ~isempty(bumpi)
        ampthresh = 0.5;

        % ON crossings
        idyl    = bump >= ampthresh;
        idy     = find(idyl);
        idy     = idy(idy > 1);
        oncrossi  = idy(bump(idy-1) < ampthresh);

        % OFF crossings
        idyl    = bump <= ampthresh;
        idy     = find(idyl);
        idy     = idy(idy > 1);
        offcrossi = idy(bump(idy-1) > ampthresh);

        % handle edge cases (no full on/off pairs)
        if isempty(oncrossi) && numel(offcrossi) == 1
            oncrossi = 1;
        end
        if isempty(offcrossi) && numel(oncrossi) == 1
            offcrossi = numel(bump);
        end
        if isempty(oncrossi) && isempty(offcrossi)
            oncrossi  = 1;
            offcrossi = numel(bump);
        end

        if offcrossi(1) < oncrossi(1)
            oncrossi = [1; oncrossi];
        end
        if oncrossi(end) > offcrossi(end)
            offcrossi = [offcrossi; numel(bump)];
        end

        % only keep bumps longer than minon seconds
        minon_samples = round(minon * fs);
        keep = find((offcrossi - oncrossi) > minon_samples);

        bump  = zeros(numel(sdff), 1);
        for k = 1:numel(keep)
            bump(oncrossi(keep(k)) : offcrossi(keep(k))) = 1;
        end

        bumpi = find(bump);
        oni   = oncrossi(keep);
        offi  = offcrossi(keep);
    else
        oni  = [];
        offi = [];
    end
end





%%
% function fig_no = plottrials_flex(data, savepath, f, which_trials, dff_yaxis, dff_coloraxis, base_thr, add2txt, print_flag, file_type, trajcmap, fig_no )
% % if which_trials is blank, plot all trials
% % dff axes are 2 element vectors of lower and upper boundries (if empty: dff_yaxis = [-.1 .8], dff_coloraxis = [0 0.4];
% % add2txt goes into title and filename (with fly date, experiment# and trial#)
% % print_flag = 0 won't print figure and doesn't need a savepath
% % file_type is 'png' or 'eps'
%     
% 
% %     savepath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/Janelia odor nav 2023/';
%     no_dffcolor = 0;
%     surf = 0;
%     bumppos = 0;
% 
% 
%     loop = {'Wind Open Loop' , 'Wind Closed Loop'};
% %     fig_no = 10;
%     if isempty(dff_yaxis)
%         yl_ca = [-.1 .8]; %[-0.1 1.3]; %[0 1.25]; %[0 0.7][-0.15 .75]; %[-0.4 0.6]; %[-.5 .5]; %[-.5 1]; %[-1 1.5];    
%     else
%         yl_ca = dff_yaxis;
%     end
% 
%     if isempty(dff_coloraxis)
%         cl_ca = [0 .4]; %[0.05 .4]; %[0.1 .4]; %[0 .8]; %[0 0.5]
%     else        
%         cl_ca = dff_coloraxis;
%     end
% 
% 
%     smth_fun = 'movmean'; %'lowess'; %'gaussian'; %'movmean'
%     columns = 1:8;
%     
%     %set expt dur as longest trial in experiment (so if big diffs in trial lengths, will notice)
%     exp_dur = 0;
%     for t = 1:length(data)
%         if ~isempty(data(t).dff_ts)
%             exp_dur = max([data(t).dff_ts(end), exp_dur]);
%         end
%     end
%     exp_dur = exp_dur + 1;
% 
%     m = 5; n = 4;
%     layers = 0;
%     s = [1,2;5,6;9,10;13,14;17,18];
%     st = [3:4, 7:8, 11:12, 15:16, 19:20];
% 
%     if isempty(which_trials)
%         which_trials = 1:length(data);
%     end
%     for t = which_trials %[9, 18], [25 27 29 31 37 41 45 47]
% 
%         nonan = ~isnan(data(t).bump_amp(:,1));
%         odor = data(t).odor(nonan, 1);
%         wind = find(data(t).wind(nonan, 1));
%         ampthresh = .5;%0.125;
%         idyl = odor >= ampthresh; %all superthresh pts
%         idy = find(idyl); %superthresh indices
%         idy = idy(idy>1); %all subthresh after 1
%         oncrossi = idy(odor(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
% %             ontimes = tt(oncrossi);
%       
% 
%         
% %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %cwoo          %closed loop, odor on, wind on, not plume, not bad/empty trial
% %         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
% %         if ~data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 %owoo
% %         if ~isempty(data(t).bump_amp) && ~data(t).odor_on %odor off
%         if ~isempty(data(t).bump_amp)% && length(oncrossi) < 2 %all trials
% 
% 
%            
%                     if data(t).dff_ts(end) >= 65
%                         exp_dur = data(t).dff_ts(end) + 1;
%                     end
%                     exp_dur = 68;
% 
%                     %split dff traces into individual layers
% %                     if size(data(t).bump_amp,2) == 2
% %                         m = 7; layers = 1; columns = [1:8; 9:16];
% %                         s = [1,2;5,6;9,10;13,14;17,18;21,22;25,26];
% %                         st = [3:4, 7:8, 11:12, 15:16, 19:20, 23:24, 27:28];
% %                     end
%                     tt = data(t).calc_ts;
%                     sp = strfind(data(t).filename, '-');
%                     us = strfind(data(t).filename, '_');
% 
%                     trli = strfind(data(t).thorname, '_corr');
%                     trlstr = data(t).thorname(trli-3:trli-1);
% 
%                     figure(fig_no); clf 
%                     set(gcf, 'Position', [100 100 1200 500]); %[100 100 500 1000])
%                     sgtitle(['F', num2str(f), ',t', num2str(t), ', Fly' num2str(f) '_' data(t).filename(sp(1)-6:sp(1)-1) ', Trial ' trlstr ' - ' loop{data(t).closed_loop + 1}], 'Interpreter', 'none')
%             
%                     %Row 1 - heading
%                     if ~data(t).closed_loop
%                         heading = data(t).calc_heading;
%                     else
%                         heading = data(t).calc_windpos;
%                     end
% %                     heading = data(t).calc_heading; %use for wind shift trials
% 
%                 heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 50)));
% %                 heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
% %                 uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 50);
%                 uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 100);               
%                 
% 
%                     subplot(m, n, s(1,:)); hold on; ylabel('Heading(deg)'); 
% %                     plot(data(t).calc_ts, uwheading, 'color', [0 0 0])
% %                     yl = [-pi, 3*pi];
%                     plot(data(t).calc_ts, heading, 'color', [0 0 0])
% %                     plot(data(t).calc_ts, smoothdata(heading, smth_fun, 25), 'color', [0 0 0])
%                     yl = [-180 180];
%                     plot([0 exp_dur], [0 0], 'k:', 'Linewidth', 0.1)
%                     plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, .06, 'k')
%                     plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, .5, 'm')
%                     plot([35 35], ylim, 'k'); plot([37 37], ylim, 'k')
%                     ylim(yl); 
%                     xlim([0 exp_dur])  
% %                     if length(oncrossi > 2)
% %                         xlim([tt(wind(1)) tt(wind(end))])  
% %                     end
%                     xlim([3 exp_dur])  
%             
%             
%                     %Row 2 - avel or uwvel or straightness
%                     subplot(m, n, s(2,:)); hold on; 
% 
% %                     %uwv
% %                     speed = data(t).calc_speed;
% %                     uwv = cos(heading*pi/180).*speed;
% %                     plot(data(t).calc_ts, smoothdata(uwv, smth_fun, 80), 'color', [0 0 0]);
% %                     ylabel('uwVel (mm/s)'); yl = [-14 20];
% 
%                     %avel
%                     plot(data(t).calc_ts, smoothdata(data(t).calc_deltaz, smth_fun, 80), 'color', [0 0 0]);
%                     ylabel('aVel (rad/s)'); yl = [-7 7];
%                     plot([0 exp_dur], [0 0], 'k:', 'Linewidth', 0.1)
%                     plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, .06, 'k')
%                     plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, .5, 'm')
%                     plot([35 35], ylim, 'k'); plot([37 37], ylim, 'k')
%                     ylim(yl); xlim([0 exp_dur])
% %                     if length(oncrossi > 2)
% %                         xlim([tt(wind(1)) tt(wind(end))])  
% %                     end
%                     xlim([3 exp_dur])  
%             
%             %         %Straightness
%             %         subplot(m, n, s(2,:)); hold on; ylabel('straighness'); yl = [0.5 1];
%             %         plot(data(t).calc_ts, smoothdata(data(t).straight, smth_fun, 10), 'color', [0 0 0])
%             %         plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, .06)
%             %         plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, .5)
%             % %         ylim(yl); xlim([0 exp_dur])
%             
%                     %Row 3 - fvel or xpos
%                     subplot(m, n, s(3,:)); hold on; ylabel('fVel (mm/s)'); yl = [-8 17]; %[-8 17];
%                     fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%                     fvel = data(t).calc_deltapitch*9.52/2;%*fps; %%%added for 31a11 only bc forgot to scale to fps in processing%%%
%                     plot(data(t).calc_ts, smoothdata(fvel, smth_fun, 80), 'color', [0 0 0])
%                     plot([0 exp_dur], [0 0], 'k:', 'Linewidth', 0.1)
%                     plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, .06, 'k')
%                     plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, .5, 'm')
%                     plot([35 35], ylim, 'k'); plot([37 37], ylim, 'k')
%                     ylim(yl); xlim([0 exp_dur])   %[-4 4]; %[-6 6];
% %                     if length(oncrossi > 2)
% %                         xlim([tt(wind(1)) tt(wind(end))])  
% %                     end 
%                     xlim([3 exp_dur])  
% 
%             %         %Xpos
%             %         subplot(m, n, s(3,:)); hold on; ylabel('xpos (mm)'); yl = [-125 125];
%             %         plot(data(t).calc_ts, smoothdata(data(t).plume_x-107, smth_fun, 1), 'color', [0 0 0])
%             %         plot([0 exp_dur], [0 0], 'k:', 'Linewidth', 0.1)
%             %         plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl, .06)
%             %         plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl, .5)
%             %         ylim(yl); xlim([0 exp_dur])               
%             
%             
%             
%                     if ~layers
%                         %Row 4 - bump amp
%                         subplot(m, n, s(4,:)); hold on; ylabel('dF/F Green');
%                         dff = data(t).bump_amp(:,1);
%                         dffi = ~isnan(data(t).bump_amp(:,1));
%                         plot(data(t).calc_ts, dff, 'k')
%                         plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl_ca, .06, 'k')
%                         plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl_ca, .5, 'm')
%                         plot([35 35], ylim, 'k'); plot([37 37], ylim, 'k')
%                         ylim(yl_ca); xlim([0 exp_dur])
% %                         if length(oncrossi > 2)
% %                             xlim([tt(wind(1)) tt(wind(end))])  
% %                         end
%                         xlim([3 exp_dur])  
%                 
%                 
%                         %Row 5 - cmap
%                         subplot(m, n, s(5,:)); hold on; ylabel('dF/F Green');
%                 
%                         sigma_a=.5; %smoothing over columns
%                         sigma_b=.5; %smoothing over time
%                 %         im = data(t).dff_red(columns,:);
%                         im = data(t).dff_green(columns,:);
% %                         im = data(t).dff_diff(columns,:);
%                         smim=imgaussfilt(im,[sigma_a, sigma_b]);
%                 
%                 
% %                         imagesc(data(t).dff_ts, 1:length(columns), smim, [.1 .5])
%                         imagesc(data(t).dff_ts, 1:length(columns), smim, cl_ca)
%                         colormap('gray')
%                         ylim([0.5 length(columns)+0.5]); xlim([0 exp_dur])
%                         colorbar('east')
%                         xlim([3 exp_dur])  
% %                         if length(oncrossi > 2)
% %                             xlim([tt(wind(1)) tt(wind(end))])  
% %                         end
% %                         bumppos = 1;
%                         if bumppos
%                             nonan = ~isnan(data(t).bump_amp(:,1));
%                             bumpdff = data(t).bump_amp(nonan, 1);
%                             tt = data(t).calc_ts(nonan);
%                             bump = bumpfinder2(bumpdff, tt, 5, 3, base_thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
%                             bumpi = find(bump);
% 
%                 
%                 %             tt = data(t).calc_ts(~isnan(data(t).calc_ts));
%                             col = data(t).bump_col(nonan);
%                             com = data(t).bump_com(nonan);
%                             pos = smoothdata(data(t).bump_col(nonan), 'movmean', 200);
%                             scatter(tt(bumpi), pos(bumpi), 4, 'g', 'filled')
%                             scatter(tt(bumpi), heading(bumpi)*8/360+4.5, 2, 'r', 'filled')
% %                             scatter(tt(bumpi), wrapTo2Pi(uwheading(bumpi))*1.1+1, 2, 'r', 'filled')
%                         end
%                     else
%                         for l = 1:2
%                             subplot(m, n, s(4+l*2-2,:)); hold on; ylabel('dF/F Green');
%                             dff = data(t).bump_amp(:,l);
%                             dffi = ~isnan(data(t).bump_amp(:,l));
%                             plot(data(t).calc_ts, dff, 'k')
%                             plot_stimfills(gca, data(t).calc_ts, data(t).wind, yl_ca, .06, 'k')
%                             plot_stimfills(gca, data(t).calc_ts, data(t).odor, yl_ca, .5, 'm')
%                             ylim(yl_ca); xlim([0 exp_dur])
%                     
%                     
%                             %Row 5 - cmap
%                             subplot(m, n, s(5+l*2-2,:)); hold on; ylabel('dF/F Green');
%                     
%                             sigma_a=.1; %smoothing over columns
%                             sigma_b=.1; %smoothing over time
%                     %         im = data(t).dff_red(columns,:);
%                             im = data(t).dff_diff(columns(l,:),:);
%             %                 im = data(t).dff_green(columns(l,:),:);
%                     %         smim=imgaussfilt(im,[sigma_a, sigma_b]);
%                     
%                             imagesc(data(t).dff_ts, 1:length(columns(1,:)), im, cl_ca)
%                             colormap('gray')
%                             ylim([0.5 length(columns(1,:))+0.5]); xlim([0 exp_dur]) 
%                             colorbar('east')
% %                             if length(oncrossi > 2)
% %                                 xlim([tt(wind(1)) tt(wind(end))])  
% %                             end
%                             xlim([3 exp_dur])  
%                         end
% 
%                     end
%             
%                     %Trajectory
%                     dffts = data(t).calc_ts(dffi);
%                     dffodor = data(t).odor(dffi);
%                     dffwind = data(t).wind(dffi);
%                     o = dffts(dffodor > 0.5); %find stim onset/offset times
%                     w = dffts(dffwind > 0.5); %find stim onset/offset times
%                     odor = data(t).odor(dffi); %data(t).odor(~isnan(data(t).odor));
%                     wind = data(t).wind(dffi);
%                     if isempty(o)
%                         o = [-1 -1]; end
%                     if isempty(w)
%                         w = [-1 -1]; end
%             
%             
%                     yd = [-75 410]; xd = [-range(yd)/2 range(yd)/2];
%                     subplot(m,n,st); hold on
%                     set(gca,'YTickLabel',[], 'XTickLabel', []);
%                     if ~all(isnan(data(t).plume_x))
%                         ax1 = axes; hold on
%                         subplot(m,n,st,ax1)
% %                         g = fill(gca, [xd(1) xd(2) xd(2) xd(1)], [yd(1) yd(1) yd(2) yd(2)], 'k');
% % %                         g = fill(gca, [-300 500 500 -300], [-300 -300 500 500], 'k');
%                         
%                         plume = h5read('/Users/kathmn01/Dropbox (NYU Langone Health)/Data/10302017_10cms_bounded_2.h5', '/dataset2', [1 1 1], [216, 406, 1]);
%                         x = data(t).plume_x; x = x(dffi); %x = x(~isnan(x));
%                         y = data(t).plume_y; y = y(dffi); %y = y(~isnan(y));
% 
%                         y = (y-406).*-1;
%                         x = (x-108).*-1;
% 
%                         pix2mm = 0.74; mm2pix = 1/pix2mm;
% %                         imagesc(ax1, 1:216, 1:406, plume')
%                         colormap(ax1,'gray')
%             
%                         xlim(xd); ylim(yd);
% %                         xlim([108-275 108+275]); ylim([-100 450]);
%                         clearvars plume
% %                         add2txt = '_plume';
%                     else
%                         fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%                         heading = heading(dffi); %heading(~isnan(heading));
%                         speed = data(t).calc_speed(dffi); %data(t).calc_speed(~isnan(data(t).calc_speed));
%                         x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
%                         y = cumsum(speed.*cos(heading*pi/180))/fps;
% %                         add2txt = '';
%                     end
%                     
%                     
% %                     figure; plot(x, y, 'r', 'LineW', 3)
% %                     xlim([-250 250]); ylim([-200 400])
% % %                     keyboard
% 
% 
%                     yd = smoothdata(dff(dffi), smth_fun, 200);
%             
%                     ax2 = axes; hold on
%                     subplot(m,n,st,ax2)
% 
% %                     no_dffcolor = 1;
% %                     surf = 0;
%                     if no_dffcolor
%                         plot(x, y, 'r', 'Linewidth', 2)
%                         if ~all(isnan(data(t).plume_x))
%                             scatter(ax2, x(find(odor)), y(find(odor)), 5, 'white', 'filled')%'.', 'Linew', 2)
%                         end
%                     else
%                         if surf
%                             z = zeros(size(x));
%                             yq = yd;
%                             surface([x,x],[y,y],[z,z],[yq,yq], 'facecol','no','edgecol','interp','linew',3);
%                             colorbar;colormap(ax2, 'jet');
%                             caxis([cl_ca])
%                         else
% %                             scatter(ax2, x(find(wind)), y(find(wind)), 30, yd(find(wind)), 'filled')
%                             scatter(ax2, x, y, 150, yd, 'filled')
% %                             plot(ax2, x, y, 'w') %plot a line behind
% %                             scatter(ax2, x, y, 15, yd, 'filled')
%                             if ~all(isnan(data(t).plume_x))
%                                 scatter(ax2, x(find(odor))-12, y(find(odor)), 8, 'm', 'filled')%'.', 'Linew', 2)
% % %                                 scatter(ax2, x(find(odor))-8, y(find(odor)), 8, 'white', 'filled')%'.', 'Linew', 2)
% %                                 scatter(ax2, x(find(odor))-12, y(find(odor)), 8, 'white', 'filled')%'.', 'Linew', 2)
% % %                                 scatter(ax2, x(find(odor))-15, y(find(odor)), 8, 'white', 'filled')%'.', 'Linew', 2)
%                             end
% %                             colormap(ax2,'jet')
% %                             colormap(ax2,'parula')
%                             colormap(ax2,trajcmap)
%                             caxis(cl_ca)
%                         end
%                     end
%                     if ~all(isnan(data(t).plume_x))
%                         plot([-107, 108, 108, -107, -107], [0, 0, 405, 405, 0], 'k-', 'LineWidth', 0.5);
%                         plot([-107, 0, 108], [0, 405, 0], 'k-', 'LineWidth', 0.5)
% %                         plot([-107, 108, 108, -107, -107], [0, 0, 405, 405, 0], 'w-', 'LineWidth', 0.5);
% %                         plot([-107, 0, 108], [0, 405, 0], 'w-', 'LineWidth', 0.5)
%                     end
%                     plot(0, 406, 'xr', 'MarkerS', 5)
% 
%                     xlim(ax2,[-375 375]); ylim([-200 550]); %%%%%%%%%%%%%%%%%%%nonplume trials %%%%%%%%%%%%%%%
% %                     xlim(ax2,[-325 325]); ylim([-200 450]); %%%%%%%%%%%%%%%%%%%nonplume trials %%%%%%%%%%%%%%%
% %                     xlim(ax2,[-300 200]); ylim([-250 250]);
%                     colorbar(ax2);
%             
%                     if ~all(isnan(data(t).plume_x))
% %                         xlim(ax2, [108-275 108+275]); ylim([-100 450]);
% %                         xlim(ax2, [xd(1) xd(2)]); ylim([yd(1) yd(2)]);
%                         xlim(ax2, [-240 240]); ylim([-175 410]);
%                         %%Link them together
%                         linkaxes([ax1,ax2])
%             
%                         %%Hide the top axes
%                         ax2.Visible = 'off';
%                         ax2.XTick = [];
%                         ax2.YTick = [];
%                         set([ax1,ax2],'Position',[.543 .11 .355 .815]);
%             
%                     else
%                         hold on
%                         plot(ax2, x(1), y(1), '*m', 'Linewidth', 3)
%             
%             
%                         [~, ~, idx] = unique(round(abs(dffts - o(1))));
%                         oo = find(idx==1);
%                         if oo(1) <= length(x)
%                             plot(x(oo(1)), y(oo(1)), 'xk', 'Linewidth', 2); end
%                 
%                         [~, ~, idx] = unique(round(abs(dffts - o(end))));
%                         oo = find(idx==1);
%                         if oo(1) <= length(x)
%                             plot(x(oo(1)), y(oo(1)), 'xk', 'Linewidth', 2); end
%                 
%                         [~, ~, idx] = unique(round(abs(dffts - w(1))));
%                         ww = find(idx==1);   
%                         if ww(1) <= length(x)
%                             plot(x(ww(1)), y(ww(1)), 'xk', 'Linewidth', 1); end
%                 
%                         [~, ~, idx] = unique(round(abs(dffts - w(end))));
%                         ww = find(idx==1); 
%                         if ww(1) <= length(x)
%                             plot(x(ww(1)), y(ww(1)), 'xk', 'Linewidth', 1); end
%             
%                     end
%             
%                     sp = strfind(data(t).filename, '-');
%                     us = strfind(data(t).filename, '_');
%                     fig_no = print_fig(print_flag, gcf, [savepath 'fly' num2str(f) '_t' trlstr add2txt], file_type, fig_no); 
% %                     fig_no = print_fig(print_flag, gcf, [savepath data(t).filename(sp(1)-6:sp(1)-1) '_' data(t).filename(1:us(1)-1) ...
% %                         '_trial' num2str(t) add2txt], file_type, 1);  
% 
% 
%         end
%     end
%     
% end
% 
% 
% 
% function fig_no = print_fig(print_flag, fig_handle, fig_filename, file_type, fig_no)
% % print_fig(print_flag, fig_handle, fig_filename, file_type, fig_no)
% % print_flag binary
% % fig_filename is full file path and name
% % file_type: 'eps' or 'png' or multiple types as rows (e.g. ['eps'; 'png'])
% % fig_no = fig_no + 1 if print_flag = 0
% if print_flag
%     for i = 1:size(file_type, 1) 
%         if strcmp(file_type(i,:), 'eps')
%             print_type = 'epsc'; print_mode = '-painters';
%         elseif strcmp(file_type(i,:), 'png')
%             print_type = 'png'; print_mode = '';
%         elseif strcmp(file_type(i,:), 'fig')
%             savefig(fig_handle, fig_filename)
%         else
%             error('Undefined filetype')
%         end
%         if ~strcmp(file_type(i,:), 'fig')
%             print(fig_handle, ['-d' print_type], print_mode, [fig_filename '.' file_type(i,:)])
%         end
%     end
% else
%     fig_no = fig_no + 1;
% end
%     
% end
% 
% 
% % function [bump, oni, offi] = bumpfinder2(dff, tt, maxoff, minon)
% % % dff needs to be nonaned
% % 
% % %             ipt = []; bases = [];bumpi = [];thresh = [];
% %         sdff = smoothdata(dff, 'movmean', 100);
% % 
% %         %find different Ca states for this trial
% %         ipt = findchangepts(sdff,'Statistic','std','MinThreshold',2000); 
% %         ipt = [1; ipt; length(sdff)];
% %         for i = 2:length(ipt)
% %             bases(i-1) = mean(sdff(ipt(i-1):ipt(i)));
% %         end
% % 
% %         % find separation between Ca states (midpoint of kmeans)
% %         if length(bases) > 1   
% %             [idx, c] = kmeans(bases', 2);
% %             thresh = min(c)+abs(diff(c))/2;
% % %                 thresh = mean(bases);
% %         else
% %             thresh = mean(sdff)+std(sdff)*2.5; %set some boundries for trials with no bump
% %         end
% % 
% %         %set some boundries for trials with no bump
% %         if thresh < 0.2  
% %             thresh = 0.2; end
% %         if length(bases) == 2 && abs(diff(bases)) < std(sdff)*2.5% && thresh <= 0.2
% %             thresh = 0.3;
% %         end
% % %             thresh = 0.27; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manual theshhold if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         bumpi = find(sdff > thresh);
% % 
% %         %fill in gaps <300 samples wide
% %         gaps = diff(bumpi);idx = find(gaps > 1 & gaps < maxoff); 
% %         bump = zeros(length(sdff), 1);
% %         bump(bumpi) = 1;
% %         for g = 1:length(idx)
% %             this_gap = [];
% %             this_gap = bumpi(idx(g))+1:bumpi(idx(g)+1)-1;
% %             bump(this_gap) = 1;
% %         end
% %         bumpi = find(bump);
% % 
% %     if ~isempty(bumpi)
% %         ampthresh = .5;%0.125;
% %         idyl = bump >= ampthresh; %all superthresh pts
% %         idy = find(idyl); %superthresh indices
% %         idy = idy(idy>1); %all subthresh after 1
% %         oncrossi = idy(bump(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
% % 
% % 
% %         ampthresh = 0.5; %set thresh
% %         idyl = []; idy = []; 
% %         idyl = bump <= ampthresh; %all subthresh pts
% %         idy = find(idyl); %superthresh indices
% %         idy = idy(idy>1); %all subthresh after 1
% %         offcrossi = idy(bump(idy-1)>ampthresh); %all pts where prev pt is superthresh (aka crossing)
% % 
% %         if isempty(oncrossi) & length(offcrossi) == 1
% %             oncrossi = 1; end
% %         if isempty(offcrossi) & length(oncrossi) == 1
% %             offcrossi = length(bump); end
% %         if isempty(oncrossi) & isempty(offcrossi)
% %             oncrossi = 1; offcrossi = length(bump); end
% % 
% %         if offcrossi(1) < oncrossi(1)
% %             oncrossi = [1;oncrossi]; end
% %         if oncrossi(end) > offcrossi(end)
% %             offcrossi = [offcrossi; length(bump)]; end
% % 
% %         keep = find((offcrossi - oncrossi) > minon);
% %         bump = zeros(length(sdff), 1);
% % %             bump(bumpi) = 1;
% %         for k = 1:length(keep)
% %             bump(oncrossi(keep(k)):offcrossi(keep(k))) = 1;
% %         end
% %         bumpi = find(bump);
% %         oni = oncrossi(keep);
% %         offi = offcrossi(keep);
% %         ont = tt(oni);
% %         offt = tt(oni);
% %     else
% %         oni = [];
% %         offi = [];
% %     end
% % end
% 
% function [bump, oni, offi, thresh] = bumpfinder2(dff, tt, maxoff, minon, base_thr)
% % dff needs to be nonaned
%             fs = 1/median(diff(tt));
% %             ipt = []; bases = [];bumpi = [];thresh = [];
%             sdff = smoothdata(dff, 'movmean', round(fs)); %1 sec
% %             sdff = smoothdata(dff, 'movmean', 100); %hack
% 
%             %find different Ca states for this trial
% %             ipt = findchangepts(sdff,'Statistic','std','MinThreshold',round(fs*20)); %hack
% %             ipt = findchangepts(sdff,'Statistic','std','MinThreshold',2000);%hack!!!!!!!!!!!!!!!!!!!
%             ipt = findchangepts(sdff); %fc1!!!!!!!!!  %%%% i think this is now obsolete bc it doesn't show change pts anymore? check!!!
% %             ipt = [1, ipt, length(sdff)];
%             ipt = [1; ipt; length(sdff)];
%             for i = 2:length(ipt)
%                 bases(i-1) = mean(sdff(ipt(i-1):ipt(i)));
%             end
% 
%             % find separation between Ca states (midpoint of kmeans)
%             if length(bases) > 1   
%                 [idx, c] = kmeans(bases', 2);
%                 thresh = min(c)+abs(diff(c))/2;
% %                 thresh = mean(bases);
%             else
%                 thresh = mean(sdff)+std(sdff)*2.5; %set some boundries for trials with no bump
%             end
% 
%             %set some boundries for trials with no bump
%             if thresh < 0.2  
%                 thresh = 0.2; end
%             if length(bases) == 2 && abs(diff(bases)) < std(sdff)*2.5% && thresh <= 0.2
% %                 thresh = 0.3;
% %                 thresh = mean(sdff)+ std(sdff)*2.5;
% %                 thresh = mean(sdff)+ std(sdff)*1.5;
%                 thresh = base_thr;
%             end
% %             thresh = 0.23; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manual theshhold if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             bumpi = find(sdff > thresh);
% %             bumpi = sdff > thresh(use_flies(f));
% %             figure(2); clf; hold on; plot(tt, sdff); plot([0 65], [thresh, thresh], 'm')
% 
%             %fill in gaps <300 samples wide
%             gaps = diff(bumpi);idx = find(gaps > 1 & gaps < maxoff*round(fs)); 
%             bump = zeros(length(sdff), 1);
%             bump(bumpi) = 1;
%             for g = 1:length(idx)
%                 this_gap = [];
%                 this_gap = bumpi(idx(g))+1:bumpi(idx(g)+1)-1;
%                 bump(this_gap) = 1;
%             end
%             bumpi = find(bump);
% %             plot(tt(bumpi), sdff(bumpi), '*g')
% 
%         if ~isempty(bumpi)
%             ampthresh = .5;%0.125;
%             idyl = bump >= ampthresh; %all superthresh pts
%             idy = find(idyl); %superthresh indices
%             idy = idy(idy>1); %all subthresh after 1
%             oncrossi = idy(bump(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
% %             ontimes = tt(oncrossi);
% 
% %             for c = 1:length(oncrossi)
% %                 win1 = (oncrossi(c)-win(1)*fps):(oncrossi(c)+win(2)*fps);
% %                 wint1 = tt(win1)-tt(oncrossi(c));
% %                 ons = [ons, met(win1)];
% %                 subplot(2,2,1); hold on; plot(wint1, bumpi(win1), 'k')
% %                 subplot(2,2,3); hold on; plot(wint1, met(win1), 'k')
% %             end
% 
%             ampthresh = 0.5; %set thresh
%             idyl = []; idy = []; 
%             idyl = bump <= ampthresh; %all subthresh pts
%             idy = find(idyl); %superthresh indices
%             idy = idy(idy>1); %all subthresh after 1
%             offcrossi = idy(bump(idy-1)>ampthresh); %all pts where prev pt is superthresh (aka crossing)
% %             offtimes = tt(offcrossi);
% 
% %             for c = 1:length(offcrossi)
% %                 win2 = (offcrossi(c)-win(1)*fps):(offcrossi(c)+win(2)*fps);
% %                 wint2 = tt(win2)-tt(offcrossi(c));
% %                 offs = [offs, met(win2)];
% %                 subplot(2,2,2); hold on; plot(wint2, stimi(win2), 'k')
% %                 subplot(2,2,4); hold on; plot(wint2, met(win2), 'k')
% %             end
%             if isempty(oncrossi) & length(offcrossi) == 1
%                 oncrossi = 1; end
%             if isempty(offcrossi) & length(oncrossi) == 1
%                 offcrossi = length(bump); end
%             if isempty(oncrossi) & isempty(offcrossi)
%                 oncrossi = 1; offcrossi = length(bump); end
% 
%             if offcrossi(1) < oncrossi(1)
%                 oncrossi = [1;oncrossi]; end
%             if oncrossi(end) > offcrossi(end)
%                 offcrossi = [offcrossi; length(bump)]; end
% 
%             keep = find((offcrossi - oncrossi) > minon*round(fs)); %gets rid of bumps less than minon duration
%             bump = zeros(length(sdff), 1);
% %             bump(bumpi) = 1;
%             for k = 1:length(keep)
%                 bump(oncrossi(keep(k)):offcrossi(keep(k))) = 1;
%             end
%             bumpi = find(bump);
%             oni = oncrossi(keep);
%             offi = offcrossi(keep);
%             ont = tt(oni);
%             offt = tt(oni);
%         else
%             oni = [];
%             offi = [];
%         end
% end