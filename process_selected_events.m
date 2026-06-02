function batch_result = process_selected_events(event_dirs, selected_groups, T_list, fc_high, varargin)
% PROCESS_SELECTED_EVENTS
% 批量处理一个或多个组内的所有事件，并合并 metric_table。
%
% 输入：
%   event_dirs       - build_earthquake_read_dirs 生成的目录索引
%   selected_groups  - 例如 5 或 [3 5 7]
%   T_list
%   fc_high
%
% 输出 batch_result:
%   .metric_all
%   .event_results
%   .failed_events
%   .selected_groups

    p = inputParser;
    addRequired(p, 'event_dirs');
    addRequired(p, 'selected_groups', @isnumeric);
    addRequired(p, 'T_list', @isnumeric);
    addRequired(p, 'fc_high', @isnumeric);

    addParameter(p, 'SaveEachEvent', false, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'SaveFolder', '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'KeepFullResult', false, @(x)islogical(x) || isnumeric(x));

    addParameter(p, 'CheckWindow', 2.0, @isnumeric);
    addParameter(p, 'NoiseThreshold', 0.3, @isnumeric);

    parse(p, event_dirs, selected_groups, T_list, fc_high, varargin{:});

    save_each_event = logical(p.Results.SaveEachEvent);
    keep_full_result = logical(p.Results.KeepFullResult);
    save_folder = char(p.Results.SaveFolder);

    if save_each_event && isempty(save_folder)
        save_folder = fullfile(pwd, 'processed_events');
    end

    if save_each_event && ~isfolder(save_folder)
        mkdir(save_folder);
    end

    metric_all = table();

    event_results = struct([]);
    failed_events = struct('group_m', {}, 'event_n', {}, ...
                           'event_code', {}, 'message', {});

    n_done = 0;
    n_failed = 0;

    group_list = [event_dirs.m];

    fprintf('\n============================================================\n');
    fprintf('开始批量处理 selected_groups = %s\n', mat2str(selected_groups));
    fprintf('============================================================\n');

    for ig = 1:numel(selected_groups)
        m = selected_groups(ig);

        iG = find(group_list == m, 1);
        if isempty(iG)
            warning('未找到组 G%d，跳过。', m);
            continue;
        end

        events = event_dirs(iG).events;

        fprintf('\n#################### 处理 G%d，共 %d 个事件 ####################\n', ...
            m, numel(events));

        for n = 1:numel(events)
            event_info = events(n);

            try
                result = process_one_event( ...
                    event_info, T_list, fc_high, ...
                    'CheckWindow', p.Results.CheckWindow, ...
                    'NoiseThreshold', p.Results.NoiseThreshold);

                T_event = result.metric_table;

                if ~isempty(T_event) && height(T_event) > 0
                    metric_all = append_metric_table(metric_all, T_event);
                else
                    warning('G%d-E%d-%s 的 metric_table 为空。', ...
                        event_info.m, event_info.n, event_info.event_code);
                end

                n_done = n_done + 1;

                if keep_full_result
                    event_results(n_done).event_info = result.event_info;
                    event_results(n_done).stations = result.stations;
                    event_results(n_done).rejected_info = result.rejected_info;
                    event_results(n_done).unpaired = result.unpaired;
                    event_results(n_done).metric_table = result.metric_table;
                else
                    event_results(n_done).event_info = result.event_info;
                    event_results(n_done).metric_table = result.metric_table;
                end

                if save_each_event
                    out_name = sprintf('G%d_E%d_%s_metric.mat', ...
                        event_info.m, event_info.n, event_info.event_code);
                    out_path = fullfile(save_folder, out_name);

                    metric_table = result.metric_table; 
                    event_summary = result.event_info; 
                    rejected_info = result.rejected_info; 
                    unpaired = result.unpaired; 

                    save(out_path, 'metric_table', 'event_summary', ...
                         'rejected_info', 'unpaired');
                end

            catch ME
                n_failed = n_failed + 1;

                failed_events(n_failed).group_m = event_info.m;
                failed_events(n_failed).event_n = event_info.n;
                failed_events(n_failed).event_code = event_info.event_code;
                failed_events(n_failed).message = ME.message;

                fprintf(2, '\n处理失败：G%d-E%d-%s\n原因：%s\n', ...
                    event_info.m, event_info.n, event_info.event_code, ME.message);

                continue;
            end
        end
    end

    batch_result = struct();
    batch_result.metric_all = metric_all;
    batch_result.event_results = event_results;
    batch_result.failed_events = failed_events;
    batch_result.selected_groups = selected_groups;

    fprintf('\n============================================================\n');
    fprintf('批量处理完成。\n');
    fprintf('成功事件数：%d\n', n_done);
    fprintf('失败事件数：%d\n', n_failed);
    fprintf('合并样本数：%d\n', height(metric_all));
    fprintf('============================================================\n');
end