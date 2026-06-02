%% select_events_by_group_eventn.m
% 从总表 metric_all 中，按照预设的 (group_m, event_n) 列表提取事件，
% 并拼成一个新的 table：metric_selected

%% ===================== 配置区 =====================

% 每一行是一个事件：[组号m, 组内事件号n]
target_pairs1 = [
    1, 2
    %5, 3
    6,3
    6,8
    7, 6
    7,11
];  % >= 7.0
target_pairs2 = [
    1, 1
    1,3
    5, 2
    5,6
    6,4
    6,5
    7, 1
    7,5
    7,7
    7,9
    7,10
]; % 6.5-6.9
target_pairs3 = [
    2,1
    3,1
    4,1
    4,2
    5,4
    5,5   
    6,1
    6,2
    6,6
    6,7
    6,9
    6,10
    7,2
    7,3
    7,4
    7,8
    8,1
    8,2
    8,3
];  % >= 6.0
target_pairs0=[
    5,3
];

target_pairs=target_pairs0;
% 是否打印每个事件的信息
verbose = true;

%% ===================== 提取区 =====================

metric_selected = table();

selected_info = table( ...
    [], [], strings(0,1), [], ...
    'VariableNames', {'group_m', 'event_n', 'event_code', 'N_rows'});

for i = 1:size(target_pairs, 1)

    m = target_pairs(i, 1);
    n = target_pairs(i, 2);

    mask = metric_all.group_m == m & metric_all.event_n == n;

    T_event = metric_all(mask, :);

    if isempty(T_event) || height(T_event) == 0
        warning('没有找到 G%d-E%d 对应的数据。', m, n);
        continue;
    end

    % 检查这个 m,n 是否只对应一个 event_code
    ev_codes = unique(string(T_event.event_code));

    if numel(ev_codes) > 1
        warning('G%d-E%d 对应了多个 event_code，请检查：%s', ...
            m, n, strjoin(ev_codes, ', '));
    end

    % 追加到总表
    metric_selected = [metric_selected; T_event]; %#ok<AGROW>

    % 记录摘要
    selected_info = [selected_info; table( ...
        m, n, strjoin(ev_codes, ', '), height(T_event), ...
        'VariableNames', {'group_m', 'event_n', 'event_code', 'N_rows'})]; %#ok<AGROW>

    if verbose
        fprintf('已提取 G%d-E%d：%s，N = %d\n', ...
            m, n, strjoin(ev_codes, ', '), height(T_event));
    end
end

fprintf('\n==================== 提取完成 ====================\n');
fprintf('目标事件数：%d\n', size(target_pairs, 1));
fprintf('成功提取事件数：%d\n', height(selected_info));
fprintf('合并后总行数：%d\n', height(metric_selected));

disp(selected_info);




%%
result_selected_pgv = analysis_metric_table(metric_selected, ...
    'Groups', [], ...
    'AnalysisType', 'PGV', ...
    'BasisSuffix', 'T0_020_h', ...
    'MakePlots', true, ...
    'ColorByEvent', true);