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