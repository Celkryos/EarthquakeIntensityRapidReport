function event_slope_table = extract_event_random_slopes(lme_rs_result, varargin)
% EXTRACT_EVENT_RANDOM_SLOPES
% 从事件随机斜率 LME 模型中提取每个事件的条件斜率，并画 slope vs magnitude。
%
% 输入：
%   lme_rs_result:
%       fit_lme_pgv_event_random_slope 输出的结果结构体
%
% 输出：
%   event_slope_table:
%       每一行 = 一个事件 × 一个周期
%
% 事件条件斜率：
%   slope_event = beta_logX ...
%               + beta_logX_mag * magnitude_c_event ...
%               + random_slope_logX

p = inputParser;
addRequired(p, 'lme_rs_result', @isstruct);
addParameter(p, 'MakePlot', true, @(x)islogical(x) || isnumeric(x));
parse(p, lme_rs_result, varargin{:});

make_plot = logical(p.Results.MakePlot);

models = lme_rs_result.models;
fit_tables = lme_rs_result.fit_tables;
summary = lme_rs_result.summary;

event_slope_table = table();

for k = 1:numel(models)

    lme = models{k};
    Tfit = fit_tables{k};

    if isempty(lme) || isempty(Tfit)
        continue;
    end

    T_sec = summary.T_sec(k);

    beta_logX = summary.beta_logX(k);
    beta_logX_mag = summary.beta_logX_mag(k);

    % ---------- 每个事件的基本信息 ----------
    event_vars = {'event_id', 'group_m', 'event_n', ...
        'event_code', 'magnitude', 'magnitude_c'};

    event_info = unique(Tfit(:, event_vars), 'rows');

    event_info.event_id_str = string(event_info.event_id);
    event_info.event_code = string(event_info.event_code);

    % group_m 之前是 categorical，这里转成数值便于后续看表
    event_info.group_m_num = str2double(string(event_info.group_m));

    % ---------- 随机效应 ----------
    % 你的 MATLAB 版本中 randomEffects 单输出是数值向量，
    % 所以这里必须用三输出形式。
    [B, Bnames, stats] = randomEffects(lme); %#ok<ASGLU>

    % Bnames 有些版本是 table，有些旧版本可能是 dataset
    if ~istable(Bnames)
        try
            Bnames = dataset2table(Bnames);
        catch
            error('Bnames 不是 table，且无法转换为 table。请检查 randomEffects(lme) 的输出格式。');
        end
    end

    fprintf('\nT = %.3f s randomEffects Bnames 字段：\n', T_sec);
    disp(Bnames.Properties.VariableNames');

    % 找到 Name / Level 字段
    name_var  = find_var_name(Bnames, {'Name'});
    level_var = find_var_name(Bnames, {'Level'});

    if isempty(name_var) || isempty(level_var)
        error('Bnames 表中没有找到 Name / Level 字段。请检查表头。');
    end

    re_name  = string(Bnames.(name_var));
    re_level = string(Bnames.(level_var));
    re_est   = B;

    % 提取 logX 随机斜率项
    idx_slope = re_name == "logX";

    re_slope = table();
    re_slope.event_id_str = re_level(idx_slope);
    re_slope.random_slope_logX = re_est(idx_slope);

    % 有些版本的 randomEffects 可能把截距叫做 '(Intercept)' 或 'Intercept'
    idx_intercept = re_name == "(Intercept)" | re_name == "Intercept";

    re_intercept = table();
    re_intercept.event_id_str = re_level(idx_intercept);
    re_intercept.random_intercept = re_est(idx_intercept);

    % ---------- 合并事件信息和随机斜率 ----------
    T_event = innerjoin(event_info, re_slope, 'Keys', 'event_id_str');

    if ~isempty(re_intercept)
        T_event = outerjoin(T_event, re_intercept, ...
            'Keys', 'event_id_str', ...
            'MergeKeys', true);
    else
        T_event.random_intercept = NaN(height(T_event), 1);
    end

    % ---------- 计算条件事件斜率 ----------
    T_event.T_sec = repmat(T_sec, height(T_event), 1);
    T_event.beta_logX = repmat(beta_logX, height(T_event), 1);
    T_event.beta_logX_mag = repmat(beta_logX_mag, height(T_event), 1);

    T_event.fixed_slope_at_M = beta_logX + beta_logX_mag .* T_event.magnitude_c;

    T_event.event_slope = T_event.fixed_slope_at_M + T_event.random_slope_logX;

    % 整理变量顺序
    keep_vars = {'T_sec', 'group_m_num', 'event_n', 'event_code', ...
        'magnitude', 'magnitude_c', ...
        'beta_logX', 'beta_logX_mag', ...
        'fixed_slope_at_M', ...
        'random_slope_logX', ...
        'event_slope', ...
        'random_intercept'};

    T_event = T_event(:, keep_vars);

    T_event.Properties.VariableNames{'group_m_num'} = 'group_m';

    event_slope_table = [event_slope_table; T_event]; %#ok<AGROW>
end

fprintf('\n================ 事件随机斜率提取完成 ================\n');
disp(event_slope_table);

if make_plot
    plot_event_random_slopes(event_slope_table);
end
end


function plot_event_random_slopes(event_slope_table)
% 画 event_slope vs magnitude

periods = unique(event_slope_table.T_sec, 'stable');

figure('Color', 'w');
hold on; box on; grid on;

cmap = lines(numel(periods));
markers = {'o', 's', '^', 'd', 'v', '>'};

stat_text = strings(numel(periods), 1);

for k = 1:numel(periods)
    T0 = periods(k);

    mask_T = event_slope_table.T_sec == T0;

    x = event_slope_table.magnitude(mask_T);
    y = event_slope_table.event_slope(mask_T);

    mask = isfinite(x) & isfinite(y);

    if sum(mask) < 3
        warning('T=%.3fs 有效事件数不足，跳过拟合。', T0);
        continue;
    end

    mk = markers{min(k, numel(markers))};

    scatter(x(mask), y(mask), 48, ...
        mk, ...
        'MarkerEdgeColor', cmap(k, :), ...
        'LineWidth', 1.1, ...
        'DisplayName', sprintf('T = %.0f s event slope', T0));

    mdl = fitlm(x(mask), y(mask));

    alpha = mdl.Coefficients.Estimate(2);
    beta  = mdl.Coefficients.Estimate(1);
    pval  = mdl.Coefficients.pValue(2);
    R2    = mdl.Rsquared.Ordinary;

    [rval, ~] = corr(x(mask), y(mask), 'Type', 'Pearson');

    xfit = linspace(min(x(mask)), max(x(mask)), 100)';
    yfit = predict(mdl, xfit);

    plot(xfit, yfit, '-', ...
        'Color', cmap(k, :), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('T = %.0f s fit', T0));

    stat_text(k) = sprintf('T = %.0f s: r=%.3f, R^2=%.3f, p=%.3g', ...
        T0, rval, R2, pval);
end

xlabel('Magnitude');
ylabel('Event conditional slope');
title('Event random-slope estimates vs Magnitude');

legend('Location', 'best');

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
end


function var_name = find_var_name(T, candidates)
% 在 table T 中寻找候选变量名

var_name = '';

vars = string(T.Properties.VariableNames);

for i = 1:numel(candidates)
    idx = find(vars == string(candidates{i}), 1);
    if ~isempty(idx)
        var_name = char(vars(idx));
        return;
    end
end
end