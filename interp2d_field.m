function [LonGrid, LatGrid, Z_field, meta] = interp2d_field( ...
    lon_sta, lat_sta, val_sta, varargin)
% INTERP2D_FIELD - 在经纬度平面上对台站标量场做 2D 插值，插值采用自然邻域法
%
% [LonGrid, LatGrid, Z_field, meta] = interp2d_field(lon_sta, lat_sta, val_sta, ...)
%
% 输入：
%   lon_sta  - n×1 或 1×n，经度（度）
%   lat_sta  - n×1 或 1×n，纬度（度）
%   val_sta  - n×1 或 1×n，台站标量值（PGA/PGV/烈度等）
%
% 可选参数（Name-Value）：
%   'UseLog'     - 是否在 log10 空间插值（默认 false）。
%                  true  时：对 val_sta 取 log10 再插值，输出 Z_field 为原尺度。
%                  false 时：直接在原尺度插值。
%   'GridStep'   - 网格步长（度，默认 0.02）。
%   'LonRange'   - [lon_min lon_max]，若不给则自动按观测值范围。
%   'LatRange'   - [lat_min lat_max]，若不给则自动按观测值范围。
%   'Margin'     - 在自动范围基础上各方向额外扩展的边界（度，默认 0）。
%   'NumLevels'  - 自动生成等高线层数（默认 12）。
%   'Levels'     - 手动指定插值空间中的等值线层（向量）。
%                  对 UseLog=true 时是 log10 后的层。
%
% 输出：
%   LonGrid, LatGrid - 网格化经纬度（meshgrid 格式）
%   Z_field          - 插值后的场（物理尺度；UseLog=true 时已经 10.^log）
%   meta             - 结构体，包含：
%       .useLog          - 是否使用对数插值
%       .val_raw_min/max - 去掉 NaN 后观测值的 min/max（原尺度）
%       .interp_min/max  - 插值空间中的 min/max（log10 或线性）
%       .lon_range       - [lon_min lon_max] 实际使用范围
%       .lat_range       - [lat_min lat_max] 实际使用范围
%       .levels_interp   - 插值空间中的等高线层（log10 或线性）
%       .levels_plot     - 对应物理尺度中的层（=levels_interp 或 10.^levels_interp）
%   - 使用 scatteredInterpolant(...,'natural','none')：
%       * 空间插值：Natural Neighbor
%       * 凸包外返回 NaN，台站为基准
%

    % ---------- 输入整理 ----------
    lon_sta = lon_sta(:);
    lat_sta = lat_sta(:);
    val_sta = val_sta(:);

    if ~(numel(lon_sta) == numel(lat_sta) && numel(lat_sta) == numel(val_sta))
        error('interp2d_field:SizeMismatch', ...
              'lon_sta, lat_sta, val_sta 长度必须相同。');
    end

    % ---------- 参数解析 ----------
    p = inputParser;
    p.FunctionName = 'interp2d_field';

    addParameter(p, 'UseLog', false, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'GridStep', 0.02, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'LonRange', [], @(x)isnumeric(x) && numel(x)==2);
    addParameter(p, 'LatRange', [], @(x)isnumeric(x) && numel(x)==2);
    addParameter(p, 'Margin', 0.0, @(x)isnumeric(x) && isscalar(x) && x>=0);
    addParameter(p, 'NumLevels', 12, @(x)isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'Levels', [], @(x)isnumeric(x) && isvector(x));

    parse(p, varargin{:});

    useLog    = logical(p.Results.UseLog);
    dGrid     = p.Results.GridStep;
    lonRange  = p.Results.LonRange;
    latRange  = p.Results.LatRange;
    marginDeg = p.Results.Margin;
    nLevels   = p.Results.NumLevels;
    levels_in = p.Results.Levels;

    % ---------- 先处理观测值中的 NaN/Inf ----------
    mask_valid = isfinite(lon_sta) & isfinite(lat_sta) & isfinite(val_sta);
    lon_sta = lon_sta(mask_valid);
    lat_sta = lat_sta(mask_valid);
    val_sta = val_sta(mask_valid);

    if isempty(lon_sta)
        error('interp2d_field:NoValidData', '有效观测点为空。');
    end

    meta = struct();
    meta.useLog        = useLog;
    meta.val_raw_min   = min(val_sta);
    meta.val_raw_max   = max(val_sta);

    % ---------- 插值空间中的数据 ----------
    if useLog
        mask_pos = val_sta > 0;
        if ~all(mask_pos)
            warning('interp2d_field:NonPositiveValues', ...
                'UseLog=true 但存在非正值，将被丢弃。');
            lon_sta = lon_sta(mask_pos);
            lat_sta = lat_sta(mask_pos);
            val_sta = val_sta(mask_pos);
        end
        val_interp = log10(val_sta);
    else
        val_interp = val_sta;
    end

    if isempty(val_interp)
        error('interp2d_field:NoValidInterpData', ...
              '用于插值的有效数据为空（可能全部为非正或 NaN）。');
    end

    meta.interp_min = min(val_interp);
    meta.interp_max = max(val_interp);

    % ---------- 网格范围 ----------
    if isempty(lonRange)
        lon_min = min(lon_sta) - marginDeg;
        lon_max = max(lon_sta) + marginDeg;
    else
        lon_min = lonRange(1);
        lon_max = lonRange(2);
    end

    if isempty(latRange)
        lat_min = min(lat_sta) - marginDeg;
        lat_max = max(lat_sta) + marginDeg;
    else
        lat_min = latRange(1);
        lat_max = latRange(2);
    end

    if ~(lon_max > lon_min && lat_max > lat_min)
        error('interp2d_field:InvalidRange', '经纬度范围设置有误。');
    end

    meta.lon_range = [lon_min lon_max];
    meta.lat_range = [lat_min lat_max];

    % ---------- 生成网格 ----------
    lon_vec = lon_min : dGrid : lon_max;
    lat_vec = lat_min : dGrid : lat_max;

    [LonGrid, LatGrid] = meshgrid(lon_vec, lat_vec);

    % ---------- 构造插值器（Natural Neighbor + 凸包外 NaN） ----------
    F = scatteredInterpolant(lon_sta, lat_sta, val_interp, ...
                             'natural', 'none');

    Z_interp = F(LonGrid, LatGrid);  % 插值空间中的值

    % ---------- 变回物理尺度 ----------
    if useLog
        Z_field = 10.^Z_interp;
    else
        Z_field = Z_interp;
    end

    % ---------- 自动或手动等高线层 ----------
    if ~isempty(levels_in)
        % 用户手动指定的层：认为是插值空间中的值
        levels_interp = sort(levels_in(:).');
    else
        vmin = meta.interp_min;
        vmax = meta.interp_max;
        if vmin == vmax
            % 防止全常数时 levels 退化
            dv = max(abs(vmin), 1);
            vmin = vmin - 0.1*dv;
            vmax = vmax + 0.1*dv;
        end
        levels_interp = linspace(vmin, vmax, nLevels);
    end

    meta.levels_interp = levels_interp;

    if useLog
        meta.levels_plot = 10.^levels_interp;  % 物理空间中的等值线（PGV/PGA）
    else
        meta.levels_plot = levels_interp;      % 烈度等线性量
    end
    fprintf("好了");
end
