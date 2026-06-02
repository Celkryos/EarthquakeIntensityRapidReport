function result = process_one_event(event_info, T_list, fc_high, varargin)
% PROCESS_ONE_EVENT
% 处理一个地震事件：
%   1. 依次读取 kik / knet ascii 目录
%   2. 初步基线矫正
%   3. EW/NS/UD 配对
%   4. SNR 筛除
%   5. 多频段滤波
%   6. 积分
%   7. 水平合成 PGA / PGV
%   8. 输出 stations 和 metric_table
%
% 输入：
%   event_info.m
%   event_info.n
%   event_info.event_code
%   event_info.data_dirs      cell，例如 {kik_ascii_dir, knet_ascii_dir}
%   event_info.networks       cell，例如 {'kik','knet'}
%
% 输出 result:
%   result.adata
%   result.stations
%   result.metric_table
%   result.unpaired
%   result.rejected_info
%   result.event_info

    p = inputParser;
    addRequired(p, 'event_info');
    addRequired(p, 'T_list', @isnumeric);
    addRequired(p, 'fc_high', @isnumeric);
    addParameter(p, 'CheckWindow', 2.0, @isnumeric);
    addParameter(p, 'NoiseThreshold', 0.3, @isnumeric);
    parse(p, event_info, T_list, fc_high, varargin{:});

    if numel(T_list) ~= numel(fc_high)
        error('T_list 和 fc_high 长度必须一致。');
    end

    fprintf('\n============================================================\n');
    fprintf('开始处理事件：G%d - E%d - %s\n', ...
        event_info.m, event_info.n, event_info.event_code);
    fprintf('============================================================\n');

    %% 1. 读取 kik / knet
    adata = {};

    for k = 1:numel(event_info.data_dirs)
        folder_path = event_info.data_dirs{k};

        if isfield(event_info, 'networks') && numel(event_info.networks) >= k
            net = event_info.networks{k};
        else
            net = sprintf('net%d', k);
        end

        if ~isfolder(folder_path)
            fprintf('跳过不存在的 %s 目录：%s\n', net, folder_path);
            continue;
        end

        fprintf('\n--- 读取 %s: %s ---\n', net, folder_path);
        tmp = batch_read_earthquake_data(folder_path);

        for i = 1:numel(tmp)
            if isempty(tmp{i}), continue; end

            % 添加事件来源信息，后面查错和统计会很方便
            tmp{i}.source_network = net;
            tmp{i}.group_m = event_info.m;
            tmp{i}.event_n = event_info.n;
            tmp{i}.event_code = event_info.event_code;
        end

        adata = [adata; tmp(:)]; %#ok<AGROW>
    end

    if isempty(adata)
        warning('事件 G%d-E%d 没有读取到任何记录。', event_info.m, event_info.n);

        result = struct();
        result.adata = {};
        result.stations = struct();
        result.metric_table = table();
        result.unpaired = [];
        result.rejected_info = [];
        result.event_info = event_info;
        return;
    end

    %% 2. 初步矫正
    adata = chubu(adata);

    %% 3. 记录匹配
    [stations, unpaired] = pair_earthquake_stations(adata);

    %% 4. 剔除低信噪比记录
    [adata, stations, rejected_info] = check_truncated_records( ...
        adata, stations, ...
        'CheckWindow', p.Results.CheckWindow, ...
        'NoiseThreshold', p.Results.NoiseThreshold);

    %% 5. 多频段滤波
    adata = multiband_filter_controller(adata, stations, T_list, fc_high);

    %% 6. 积分
    acc_fields = cell(size(T_list));
    for kT = 1:numel(T_list)
        T_str = sprintf('T%.3fs', T_list(kT));
        acc_fields{kT} = matlab.lang.makeValidName(['acc_' T_str]);
    end

    for i = 1:numel(adata)
        if isempty(adata{i}), continue; end

        if ~isfield(adata{i}, 'is_valid') || ~adata{i}.is_valid
            continue;
        end

        % 仅对水平分量积分，跳过 U-D
        if isfield(adata{i}, 'direction')
            d = upper(char(adata{i}.direction));
            if contains(d, 'U')
                continue;
            end
        end

        adata{i} = acc2vel(adata{i}, acc_fields);
    end

    %% 7. 合成水平时程并生成 pga/pgv 字段
    stations = hcsc(adata, stations, T_list);

    %% 8. 转成 table
    metric_table = stations_to_metric_table(stations, event_info);

    %% 9. 打包输出
    result = struct();
    result.adata = adata;
    result.stations = stations;
    result.metric_table = metric_table;
    result.unpaired = unpaired;
    result.rejected_info = rejected_info;
    result.event_info = event_info;

    fprintf('\n事件处理完成：G%d-E%d-%s，有效样本行数 = %d\n', ...
        event_info.m, event_info.n, event_info.event_code, height(metric_table));
end