% 想挑出的事件编号
target_event_n = [3 4 5 8];   % [n1 n2 ...]

% 筛选这些 event_n 对应的所有行
mask = ismember(metric_all.event_n, target_event_n);

% 塞进新的 table
metric_small = metric_all(mask, :);

%%

metric_all=metric_all_clean;
%% 单个地震斜率与震级以及T的关系
metric_folder = 'D:\FAST';  % 8 个 mat 文件的目录

event_regression_summary = collect_event_pgv_regression_summary(metric_folder);

%% ↑的图
%event_regression_summary=e11;
figure;
hold on; box on; grid on;
egg = event_regression_summary([23:25, 28:31], :);

%egg=event_regression_summary;
x = egg.magnitude;

y1 = egg.slope_T1s;
y2 = egg.slope_T2s;
y5 = egg.slope_T5s;

%y1 = egg.R2_T1s;
%y2 = egg.R2_T2s;
%y5 = egg.R2_T5s;
Y = {y1, y2, y5};
period_labels = {'T = 1 s', 'T = 2 s', 'T = 5 s'};
markers = {'o', 's', '^'};

% 用 lines 生成 3 组颜色，保证散点和拟合线颜色一致
cmap = lines(3);

stat_text = strings(3, 1);

for k = 1:3
    y = Y{k};
    mask = isfinite(x) & isfinite(y);

    % 散点
    scatter(x(mask), y(mask), 42, ...
        markers{k}, ...
        'MarkerEdgeColor', cmap(k, :), ...
        'DisplayName', period_labels{k});

    % 线性拟合：slope = alpha * M + beta
    mdl = fitlm(x(mask), y(mask));

    alpha = mdl.Coefficients.Estimate(2);
    beta  = mdl.Coefficients.Estimate(1);
    pval  = mdl.Coefficients.pValue(2);
    R2    = mdl.Rsquared.Ordinary;

    % Pearson r
    [rval, ~] = corr(x(mask), y(mask), 'Type', 'Pearson');

    xfit = linspace(min(x(mask)), max(x(mask)), 100)';
    yfit = predict(mdl, xfit);

    plot(xfit, yfit, '-', ...
        'Color', cmap(k, :), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s fit', period_labels{k}));

    stat_text(k) = sprintf('%s: r=%.3f, R^2=%.3f, p=%.3g', ...
        period_labels{k}, rval, R2, pval);
end

xlabel('Magnitude');
ylabel('Slope');
title('Slope vs Magnitude in G7 (part)');

legend('Location', 'best');

% 在图左上角写统计量
xlim_now = xlim;
ylim_now = ylim;

text_x = xlim_now(1) + 0.03 * range(xlim_now);
text_y = ylim_now(2) - 0.06 * range(ylim_now);

text(text_x, text_y, strjoin(stat_text, newline), ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'Margin', 6, ...
    'FontSize', 10);


%% 增加震级字段

%% add_magnitude_to_metric_all.m
% 给合并后的 metric_all 添加 magnitude 字段
% 顺序：
% G1: E1 E2 E3
% G2: E1
% G3: E1
% G4: E1 E2
% G5: E2 E3 E4 E5 E6   % 注意 G5 从 E2 开始
% G6: E1-E10
% G7: E1-E11
% G8: E1-E3

%% ---------- 配置区 ----------

group_counts = [3, 1, 1, 2, 5, 10, 11, 3];

start_event_n = ones(1, 8);
start_event_n(5) = 2;   % 第五组从 event_n = 2 开始

mag_values = [
    6.6
    7.1
    6.6
    6.4
    6.4
    6.1
    6.2
    6.5
    7.6
    6.1
    6.0
    6.6
    6.2
    6.2
    7.3
    6.9
    6.8
    6.3
    6.1
    7.4
    6.3
    6.0
    6.5
    6.2
    6.2
    6.0
    6.9
    7.5
    6.6
    6.0
    6.9
    6.7
    7.5
    6.0
    6.2
    6.1
];

%% ---------- 构造 group_m / event_n / magnitude 对照表 ----------

n_events = sum(group_counts);

if numel(mag_values) ~= n_events
    error('震级数量不匹配：给了 %d 个，但按 group_counts 应有 %d 个。', ...
          numel(mag_values), n_events);
end

group_list = zeros(n_events, 1);
event_list = zeros(n_events, 1);

idx = 0;
for g = 1:numel(group_counts)
    n0 = start_event_n(g);
    n1 = n0 + group_counts(g) - 1;

    for n = n0:n1
        idx = idx + 1;
        group_list(idx) = g;
        event_list(idx) = n;
    end
end

magnitude_map = table( ...
    group_list, event_list, mag_values, ...
    'VariableNames', {'group_m', 'event_n', 'magnitude'});

disp(magnitude_map);

%% ---------- 写入 metric_all ----------

metric_all.magnitude = NaN(height(metric_all), 1);

for i = 1:height(magnitude_map)
    m = magnitude_map.group_m(i);
    n = magnitude_map.event_n(i);
    mag = magnitude_map.magnitude(i);

    mask = metric_all.group_m == m & metric_all.event_n == n;

    if ~any(mask)
        warning('metric_all 中没有找到 G%d-E%d。', m, n);
        continue;
    end

    metric_all.magnitude(mask) = mag;
end

%% ---------- 检查 ----------

missing_mask = isnan(metric_all.magnitude);

if any(missing_mask)
    warning('仍有 %d 行没有匹配到 magnitude。', sum(missing_mask));

    missing_events = unique(metric_all(missing_mask, ...
        {'group_m', 'event_n', 'event_code'}), 'rows');

    disp('未匹配到震级的事件：');
    disp(missing_events);
else
    fprintf('magnitude 添加完成，所有 %d 行均已匹配震级。\n', height(metric_all));
end

%% ---------- 可选：按事件检查 ----------

event_check = groupsummary(metric_all, ...
    {'group_m', 'event_n', 'event_code', 'magnitude'}, ...
    'numel', 'station_code');

disp(event_check);