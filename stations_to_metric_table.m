function metric_table = stations_to_metric_table(stations, event_info)
% STATIONS_TO_METRIC_TABLE
% 将单个事件的 stations 结构体转换为 table。
%
% 每一行 = 一个事件-台站样本。
% 会自动收集所有 pga_T..._h / pgv_T..._h 字段。
% 对某些台站缺失的指标字段，用 NaN 填充。

    % ---------- 空输入处理 ----------
    if isempty(stations) || isempty(fieldnames(stations))
        metric_table = table();
        return;
    end

    names = fieldnames(stations);

    % ---------- 先筛选有效台站，并收集所有指标字段 ----------
    valid_names = {};
    all_metric_fields = {};

    for i = 1:numel(names)
        sta_name = names{i};
        s = stations.(sta_name);

        if ~isfield(s, 'is_valid') || ~s.is_valid
            continue;
        end

        valid_names{end+1, 1} = sta_name; %#ok<AGROW>

        fns = fieldnames(s);
        is_metric = ~cellfun(@isempty, ...
            regexp(fns, '^(pga|pgv)_T\d+_\d+_h$', 'once'));

        all_metric_fields = [all_metric_fields; fns(is_metric)]; %#ok<AGROW>
    end

    if isempty(valid_names)
        metric_table = table();
        return;
    end

    all_metric_fields = unique(all_metric_fields, 'stable');

    % ---------- 构造统一字段模板 ----------
    template = struct();

    % 事件信息
    template.group_m = NaN;
    template.event_n = NaN;
    template.event_code = "";

    % 台站信息
    template.station_code = "";
    template.station_lat = NaN;
    template.station_long = NaN;
    template.dist_km = NaN;
    template.noise_ratio = NaN;
    template.magnitude = NaN;

    % 指标字段，全部先预置为 NaN
    for k = 1:numel(all_metric_fields)
        template.(all_metric_fields{k}) = NaN;
    end

    rows = repmat(template, numel(valid_names), 1);

    % ---------- 逐台站填表 ----------
    for i = 1:numel(valid_names)
        sta_name = valid_names{i};
        s = stations.(sta_name);

        row = template;

        % 事件信息
        if isfield(event_info, 'm')
            row.group_m = event_info.m;
        end

        if isfield(event_info, 'n')
            row.event_n = event_info.n;
        end

        if isfield(event_info, 'event_code')
            row.event_code = string(event_info.event_code);
        end

        % 台站信息
        row.station_code = string(sta_name);

        if isfield(s, 'station_lat')
            row.station_lat = s.station_lat;
        end

        if isfield(s, 'station_long')
            row.station_long = s.station_long;
        end

        if isfield(s, 'dist_km')
            row.dist_km = s.dist_km;
        end

        if isfield(s, 'noise_ratio')
            row.noise_ratio = s.noise_ratio;
        end

        if isfield(s, 'magnitude')
            row.magnitude = s.magnitude;
        end

        % 指标字段
        for k = 1:numel(all_metric_fields)
            fn = all_metric_fields{k};

            if isfield(s, fn)
                val = s.(fn);

                if isnumeric(val) && isscalar(val) && isfinite(val)
                    row.(fn) = val;
                elseif isnumeric(val) && isscalar(val)
                    row.(fn) = val;
                else
                    row.(fn) = NaN;
                end
            end
        end

        rows(i) = row;
    end

    % ---------- 转 table ----------
    metric_table = struct2table(rows);
end