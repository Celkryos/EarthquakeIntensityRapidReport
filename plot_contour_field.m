function [hC, hScat] = plot_contour_field( ...
    LonGrid, LatGrid, Z_field, meta, lon_sta, lat_sta, varargin)

    lon_sta = lon_sta(:);
    lat_sta = lat_sta(:);

    p = inputParser;
    addParameter(p, 'UseLogForColorbar', true, ...
        @(x)islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    useLogColorbar = logical(p.Results.UseLogForColorbar);

    hold_state = ishold;
    hold on;

    % ---- 1. 准备用于绘图的 Z_plot ----
    if meta.useLog
        Z_plot = log10(Z_field);   % 插值结果转回 log10 空间
    else
        Z_plot = Z_field;
    end

    % ---- 2. 用 Z_plot 的真实范围重新生成 levels ----
    mask = isfinite(Z_plot);
    zmin = min(Z_plot(mask));
    zmax = max(Z_plot(mask));

    % 保留“等值线数量”不变（沿用 meta.levels_interp 的长度）
    nLevels = max(2, numel(meta.levels_interp));
    levels  = linspace(zmin, zmax, nLevels);

    % ---- 3. 画 filled 等值线 ----
    [~, hC] = contourf(LonGrid, LatGrid, Z_plot, levels, ...
                       'LineStyle', 'none');
    colormap(jet);

    hScat = plot(lon_sta, lat_sta, 'k^', 'MarkerSize', 5);

    axis equal tight;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');

    % ---- 4. colorbar 刻度的处理 ----
    cb = colorbar;

    if useLogColorbar
        if meta.useLog
            % colorbar 显示物理值（10^log）
            cb.Ticks = levels;
            cb.TickLabels = arrayfun(@(x)sprintf('%.3g', 10.^x), ...
                                     levels, 'UniformOutput', false);
            ylabel(cb, 'Value (physical, log10-uniform levels)');
        else
            cb.Ticks = levels;
            ylabel(cb, 'Value');
        end
    else
        % 直接用插值空间
        cb.Ticks = levels;
        if meta.useLog
            ylabel(cb, 'log_{10}(Value)');
        else
            ylabel(cb, 'Value');
        end
    end

    if ~hold_state
        hold off;
    end
end
