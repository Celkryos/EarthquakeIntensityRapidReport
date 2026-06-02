function result = analysis_metric_table(metric_all, varargin)
% ANALYSIS_METRIC_TABLE
% 对一个或多个组的 metric_all 做 pooled log-log 回归分析。
%
% 回归形式固定为：
%   log10(Y) = a * log10(X) + b
%
% 转回物理空间：
%   Y = C * X^a,  C = 10^b
%
% 其中：
%   Y = 最短周期基准字段，例如 pgv_T0_020_h
%   X = 各长周期字段，例如 pgv_T0_100_h, pgv_T0_500_h, ...
%
% 用法示例：
%   result = analysis_metric_table(metric_all, ...
%       'Groups', 1, ...
%       'AnalysisType', 'PGV', ...
%       'BasisSuffix', 'T0_020_h');
%
%   result.summary

%% ---------- 参数 ----------
p = inputParser;
addRequired(p, 'metric_all', @istable);

addParameter(p, 'Groups', [], @isnumeric);  % 组编号
addParameter(p, 'AnalysisType', 'PGV', @(x)ischar(x) || isstring(x));
addParameter(p, 'BasisSuffix', 'T0_020_h', @(x)ischar(x) || isstring(x));
addParameter(p, 'MinN', 5, @isnumeric);
addParameter(p, 'MakePlots', true, @(x)islogical(x) || isnumeric(x));
addParameter(p, 'ColorByEvent', true, @(x)islogical(x) || isnumeric(x));

parse(p, metric_all, varargin{:});

groups = p.Results.Groups;
analysis_type = upper(char(string(p.Results.AnalysisType)));
basis_suffix = char(string(p.Results.BasisSuffix));
min_n = p.Results.MinN;
make_plots = logical(p.Results.MakePlots);
color_by_event = logical(p.Results.ColorByEvent);

%% ---------- 基本检查 ----------
if isempty(metric_all) || height(metric_all) == 0
    error('metric_all 为空。');
end

vars = metric_all.Properties.VariableNames;

if ~ismember('group_m', vars)
    error('metric_all 中缺少 group_m 字段。');
end

if ~ismember('event_code', vars)
    error('metric_all 中缺少 event_code 字段。');
end

%% ---------- 按组筛选 ----------
mask_group = true(height(metric_all), 1);

if ~isempty(groups)
    mask_group = ismember(metric_all.group_m, groups);
end

Tdata = metric_all(mask_group, :);

if height(Tdata) == 0
    error('筛选后没有数据。Groups = %s', mat2str(groups));
end

%% ---------- 字段模式 ----------
switch analysis_type
    case 'PGA'
        basis_field = ['pga_' basis_suffix];
        target_pattern = '^pga_T\d+_\d+_h$';
    case 'PGV'
        basis_field = ['pgv_' basis_suffix];
        target_pattern = '^pgv_T\d+_\d+_h$';
    otherwise
        error('AnalysisType 必须是 PGA 或 PGV。');
end

vars = Tdata.Properties.VariableNames;

if ~ismember(basis_field, vars)
    error('metric_all 中不存在基准字段：%s', basis_field);
end

match_idx = ~cellfun(@isempty, regexp(vars, target_pattern, 'once'));
target_fields = vars(match_idx);

% 排除基准字段自身
target_fields = setdiff(target_fields, {basis_field}, 'stable');

if isempty(target_fields)
    error('没有找到可用于回归的目标字段。');
end

%% ---------- 按周期排序 ----------
T_vals = nan(numel(target_fields), 1);

for k = 1:numel(target_fields)
    T_vals(k) = parse_period_from_field(target_fields{k});
end

[T_vals, order] = sort(T_vals);
target_fields = target_fields(order);

%% ---------- 基准 Y ----------
Y0 = Tdata.(basis_field);

valid_basis = isfinite(Y0) & Y0 > 0;

if sum(valid_basis) < min_n
    error('基准字段有效样本数太少：%s, N=%d', basis_field, sum(valid_basis));
end

%% ---------- 打印概要 ----------
if isempty(groups)
    group_label = 'ALL';
else
    group_label = ['G' strjoin(string(groups), '+G')];
end

% ---------- 全局事件颜色映射 ----------
event_list = unique(string(Tdata.event_code), 'stable');

% 如果希望按事件编号排序，而不是按表中出现顺序，用下面这一行替代上面一行：
% event_list = sort(unique(string(Tdata.event_code)));

event_cmap = lines(numel(event_list));

event_color_table = table( ...
    event_list(:), ...
    event_cmap(:,1), ...
    event_cmap(:,2), ...
    event_cmap(:,3), ...
    'VariableNames', {'event_code', 'R', 'G', 'B'});

fprintf('\n事件颜色映射：\n');
disp(event_color_table);
fprintf('\n============================================================\n');
fprintf('开始 pooled log-log 回归分析\n');
fprintf('Group      : %s\n', group_label);
fprintf('Type       : %s\n', analysis_type);
fprintf('Y basis    : %s\n', basis_field);
fprintf('Events     : %d\n', numel(event_list));
fprintf('Rows       : %d\n', height(Tdata));
fprintf('Targets    : %d\n', numel(target_fields));
fprintf('============================================================\n');

%% ---------- 结果容器 ----------
n_targets = numel(target_fields);

summary = table();
models = cell(n_targets, 1);

R2_history = nan(n_targets, 1);
N_history = nan(n_targets, 1);

%% ---------- Figure 1 ----------
fig_scatter = [];

if make_plots
    fig_scatter = figure( ...
        'Name', ['Pooled Regression - ' analysis_type], ...
        'Color', 'w', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.08 0.62 0.82]);

    n_cols = ceil(sqrt(n_targets));
    n_rows = ceil(n_targets / n_cols);
end

%% ---------- 循环回归 ----------
for k = 1:n_targets
    x_field = target_fields{k};
    X0 = Tdata.(x_field);

    mask = valid_basis & isfinite(X0) & X0 > 0;

    X = X0(mask);
    Y = Y0(mask);

    event_code_this = string(Tdata.event_code(mask));

    if numel(X) < min_n
        warning('%s 有效样本数不足，跳过。N=%d', x_field, numel(X));
        continue;
    end

    XX = log10(X);
    YY = log10(Y);

    mdl = fitlm(XX, YY);

    b = mdl.Coefficients.Estimate(1);    % 截距
    a = mdl.Coefficients.Estimate(2);    % 斜率
    r2 = mdl.Rsquared.Ordinary;
    rmse = mdl.RMSE;
    C = 10.^b;

    models{k} = mdl;
    R2_history(k) = r2;
    N_history(k) = numel(X);

    n_events = numel(unique(event_code_this));

    log_formula = sprintf('log10(Y) = %.4f log10(X) %+ .4f', a, b);
    physical_formula = sprintf('Y = %.4g X^{%.4f}', C, a);

    new_row = table( ...
        string(group_label), ...
        string(analysis_type), ...
        string(basis_field), ...
        string(x_field), ...
        T_vals(k), ...
        numel(X), ...
        n_events, ...
        a, ...
        b, ...
        C, ...
        r2, ...
        rmse, ...
        string(log_formula), ...
        string(physical_formula), ...
        'VariableNames', { ...
        'group_label', ...
        'analysis_type', ...
        'basis_field', ...
        'target_field', ...
        'target_T_sec', ...
        'N', ...
        'N_events', ...
        'slope_a', ...
        'intercept_b', ...
        'coef_C', ...
        'R2', ...
        'RMSE_log10', ...
        'formula_log10', ...
        'formula_physical'});

    summary = [summary; new_row]; %#ok<AGROW>

    %% ---------- 散点图 ----------
    if make_plots
        figure(fig_scatter);
        subplot(n_rows, n_cols, k);
        hold on; box on; grid on;


        if color_by_event
            % 使用全局 event_list / event_cmap，保证所有子图颜色一致
            h_legend = gobjects(numel(event_list), 1);

            for ie = 1:numel(event_list)
                ev = event_list(ie);
                mm = event_code_this == ev;

                % 先放一个 NaN 虚拟点，保证 legend 永远完整
                h_legend(ie) = scatter(NaN, NaN, 18, event_cmap(ie, :), ...
                    'filled', ...
                    'DisplayName', char(ev));

                if any(mm)
                    scatter(XX(mm), YY(mm), 10, event_cmap(ie, :), ...
                        'filled', ...
                        'MarkerFaceAlpha', 0.45, ...
                        'MarkerEdgeAlpha', 0.45, ...
                        'HandleVisibility', 'off');
                end
            end

            % 只在第一个子图显示图例，避免每个子图都挤满
            if k == 1
                legend(h_legend, 'Location', 'best', 'FontSize', 7);
            end

        else
            scatter(XX, YY, 10, 'filled', ...
                'MarkerFaceAlpha', 0.45, ...
                'MarkerEdgeAlpha', 0.45);
        end


        x_min = min(XX);
        x_max = max(XX);
        x_range = x_max - x_min;

        if x_range <= 0
            x_range = max(abs(x_min), 1) * 0.1;
        end

        x_fit = linspace(x_min, x_max, 100)';
        y_fit = predict(mdl, x_fit);

        plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5);

        title(sprintf('T = %.3f s', T_vals(k)), ...
            'FontSize', 10, 'FontWeight', 'bold');

        xlabel(['log_{10}(' strrep(x_field, '_', '\_') ')'], 'FontSize', 8);
        ylabel(['log_{10}(' strrep(basis_field, '_', '\_') ')'], 'FontSize', 8);

        y_min = min(YY);
        y_max = max(YY);
        y_range = y_max - y_min;

        if y_range <= 0
            y_range = max(abs(y_min), 1) * 0.1;
        end

        text_x = x_min + 0.04 * x_range;
        text_y = y_max - 0.12 * y_range;

        eq_str = sprintf('a=%.3f, b=%+.3f\nR^2=%.4f, N=%d', ...
            a, b, r2, numel(X));

        text(text_x, text_y, eq_str, ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none', ...
            'FontSize', 8, ...
            'Margin', 1);

        xlim([x_min - 0.08*x_range, x_max + 0.08*x_range]);
        ylim([y_min - 0.08*y_range, y_max + 0.08*y_range]);

        hold off;
    end
end

if make_plots && ~isempty(fig_scatter)
    figure(fig_scatter);
    sgtitle(sprintf('%s pooled log-log regression: %s, Y=%s', ...
        group_label, analysis_type, strrep(basis_field, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% ---------- Figure 2: R2 vs Period ----------
fig_r2 = [];

if make_plots
    valid_r2 = isfinite(T_vals) & isfinite(R2_history);

    if any(valid_r2)
        fig_r2 = figure( ...
            'Name', ['R2 vs Period - ' analysis_type], ...
            'Color', 'w', ...
            'Units', 'normalized', ...
            'Position', [0.70 0.15 0.28 0.42]);

        plot(T_vals(valid_r2), R2_history(valid_r2), ...
            '-o', 'LineWidth', 2, 'MarkerSize', 6);

        grid on; box on;

        xlabel('Period T (s) [X]');
        ylabel('R^2');
        title(sprintf('%s pooled %s, Y=%s', ...
            group_label, analysis_type, strrep(basis_field, '_', '\_')));

        for i = find(valid_r2(:))'
            text(T_vals(i), R2_history(i), ...
                sprintf(' %.4f', R2_history(i)), ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', 9);
        end

        ylim([max(0, min(R2_history(valid_r2)) - 0.05), 1.05]);
        xlim([0, max(T_vals(valid_r2)) * 1.1]);
    end
end

%% ---------- 输出 ----------
result = struct();
result.summary = summary;
result.models = models;
result.group_label = group_label;
result.groups = groups;
result.analysis_type = analysis_type;
result.basis_field = basis_field;
result.target_fields = target_fields;
result.target_T_sec = T_vals;
result.R2_history = R2_history;
result.N_history = N_history;
result.event_list = event_list;
result.event_cmap = event_cmap;
result.event_color_table = event_color_table;
result.fig_scatter = fig_scatter;
result.fig_r2 = fig_r2;

fprintf('\n--- 回归完成 ---\n');
disp(summary(:, {'target_field', 'target_T_sec', 'N', 'N_events', ...
    'slope_a', 'intercept_b', 'coef_C', 'R2'}));
end


function T_sec = parse_period_from_field(field_name)
% 从字段名中解析周期：
%   pgv_T0_020_h -> 0.020
%   pgv_T5_000_h -> 5.000

tok = regexp(field_name, 'T(\d+)_(\d+)', 'tokens', 'once');

if isempty(tok)
    T_sec = NaN;
    return;
end

int_part = str2double(tok{1});
dec_str = tok{2};
dec_part = str2double(dec_str) / 10^numel(dec_str);

T_sec = int_part + dec_part;
end