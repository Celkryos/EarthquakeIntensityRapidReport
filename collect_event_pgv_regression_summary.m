function event_regression_summary = collect_event_pgv_regression_summary(metric_folder, varargin)
% COLLECT_EVENT_PGV_REGRESSION_SUMMARY
% 对 metric_all_G1.mat ~ metric_all_G9.mat 中的每个地震事件分别做 PGV log-log 回归，
% 并提取 T=1s, 2s, 5s 的斜率和截距。
%
% 输出表每行 = 一个地震事件。
%
% 依赖：
%   analysis_metric_table.m

    p = inputParser;
    addRequired(p, 'metric_folder', @(x)ischar(x) || isstring(x));
    addParameter(p, 'FilePattern', 'metric_all_G*.mat', @(x)ischar(x) || isstring(x));
    addParameter(p, 'TargetPeriods', [1 2 5], @isnumeric);
    addParameter(p, 'BasisSuffix', 'T0_020_h', @(x)ischar(x) || isstring(x));
    addParameter(p, 'SaveOutput', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'OutputName', 'event_pgv_regression_summary.mat', @(x)ischar(x) || isstring(x));
    parse(p, metric_folder, varargin{:});

    metric_folder = char(metric_folder);
    file_pattern = char(p.Results.FilePattern);
    target_periods = p.Results.TargetPeriods(:)';
    basis_suffix = char(p.Results.BasisSuffix);

    files = dir(fullfile(metric_folder, file_pattern));

    if isempty(files)
        error('没有找到匹配文件：%s', fullfile(metric_folder, file_pattern));
    end

    % 按文件名排序，保证 G1, G2, ... 尽量顺序处理
    [~, idx_sort] = sort({files.name});
    files = files(idx_sort);

    rows = table();

    fprintf('\n============================================================\n');
    fprintf('开始逐事件 PGV 回归汇总，共找到 %d 个 mat 文件。\n', numel(files));
    fprintf('目标周期：%s s\n', mat2str(target_periods));
    fprintf('============================================================\n');

    for iFile = 1:numel(files)
        mat_path = fullfile(files(iFile).folder, files(iFile).name);

        fprintf('\n--- 读取文件 %d/%d: %s ---\n', ...
            iFile, numel(files), files(iFile).name);

        S = load(mat_path);

        if isfield(S, 'metric_all')
            metric_all = S.metric_all;
        elseif isfield(S, 'batch_result') && isfield(S.batch_result, 'metric_all')
            metric_all = S.batch_result.metric_all;
        else
            warning('文件 %s 中没有 metric_all，跳过。', files(iFile).name);
            continue;
        end

        if isempty(metric_all) || height(metric_all) == 0
            warning('文件 %s 的 metric_all 为空，跳过。', files(iFile).name);
            continue;
        end

        event_codes = unique(string(metric_all.event_code), 'stable');

        for iEvent = 1:numel(event_codes)
            event_code = event_codes(iEvent);

            mask_event = string(metric_all.event_code) == event_code;
            metric_event = metric_all(mask_event, :);

            if height(metric_event) < 5
                warning('事件 %s 样本数太少，跳过。N=%d', event_code, height(metric_event));
                continue;
            end

            group_m = mode(metric_event.group_m);
            event_n = mode(metric_event.event_n);

            fprintf('  分析 G%d-E%d-%s, N=%d\n', ...
                group_m, event_n, event_code, height(metric_event));

            try
                result = analysis_metric_table(metric_event, ...
                    'Groups', group_m, ...
                    'AnalysisType', 'PGV', ...
                    'BasisSuffix', basis_suffix, ...
                    'MakePlots', false, ...
                    'ColorByEvent', false);

                summary = result.summary;

                row = make_empty_event_row(event_code, group_m, event_n, target_periods);

                for kT = 1:numel(target_periods)
                    T0 = target_periods(kT);

                    idx = find(abs(summary.target_T_sec - T0) < 1e-9, 1);

                    if isempty(idx)
                        warning('事件 %s 缺少 T=%.3fs 的回归结果。', event_code, T0);
                        continue;
                    end

                    suffix = period_suffix(T0);

                    row.(['slope_' suffix])     = summary.slope_a(idx);
                    row.(['intercept_' suffix]) = summary.intercept_b(idx);
                    row.(['R2_' suffix])        = summary.R2(idx);
                    row.(['N_' suffix])         = summary.N(idx);
                end

                rows = [rows; row]; %#ok<AGROW>

            catch ME
                warning('事件 %s 回归失败：%s', event_code, ME.message);

                row = make_empty_event_row(event_code, group_m, event_n, target_periods);
                row.failed = true;
                row.fail_message = string(ME.message);

                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end

    event_regression_summary = rows;

    fprintf('\n============================================================\n');
    fprintf('逐事件 PGV 回归汇总完成。\n');
    fprintf('事件数：%d\n', height(event_regression_summary));
    fprintf('============================================================\n');

    if logical(p.Results.SaveOutput)
        out_path = fullfile(metric_folder, char(p.Results.OutputName));
        save(out_path, 'event_regression_summary');
        fprintf('已保存：%s\n', out_path);
    end
end


function row = make_empty_event_row(event_code, group_m, event_n, target_periods)
% 创建一行空结果，震级 magnitude 先填 NaN，后续手动补。

    row = table();

    row.event_code = string(event_code);
    row.group_m = group_m;
    row.event_n = event_n;
    row.magnitude = NaN;

    for k = 1:numel(target_periods)
        suffix = period_suffix(target_periods(k));

        row.(['slope_' suffix]) = NaN;
        row.(['intercept_' suffix]) = NaN;
        row.(['R2_' suffix]) = NaN;
        row.(['N_' suffix]) = NaN;
    end

    row.failed = false;
    row.fail_message = "";
end


function s = period_suffix(T0)
% 生成字段后缀：
%   1   -> T1s
%   2   -> T2s
%   5   -> T5s
%   0.5 -> T0_5s

    if abs(T0 - round(T0)) < 1e-9
        s = sprintf('T%ds', round(T0));
    else
        s = sprintf('T%.3fs', T0);
        s = strrep(s, '.', '_');
    end
end