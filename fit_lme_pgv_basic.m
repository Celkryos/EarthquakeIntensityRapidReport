function lme_result = fit_lme_pgv_basic(metric_all, varargin)
% FIT_LME_PGV_BASIC
% PGV 混合效应模型（基础版，暂时不用了）：
%
%   模型：
%   logY ~ logX * magnitude_c + group_m + (1 | event_id)
%
% 其中：
%   logY = log10(pgv_T0_020_h)
%   logX = log10(pgv_T1_000_h / pgv_T2_000_h / pgv_T5_000_h)
%   magnitude_c = magnitude - mean(event-level magnitude)

    %% ---------- 参数 ----------
    p = inputParser;
    addRequired(p, 'metric_all', @istable);

    addParameter(p, 'TargetPeriods', [1 2 5], @isnumeric);
    addParameter(p, 'BasisField', 'pgv_T0_020_h', @(x)ischar(x) || isstring(x));
    addParameter(p, 'UseEventLevelMagnitudeMean', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'MakePlots', true, @(x)islogical(x) || isnumeric(x));

    parse(p, metric_all, varargin{:});

    target_periods = p.Results.TargetPeriods;
    basis_field = char(p.Results.BasisField);
    use_event_level_mag_mean = logical(p.Results.UseEventLevelMagnitudeMean);
    make_plots = logical(p.Results.MakePlots);

    %% ---------- 基本检查 ----------
    required_vars = {'group_m', 'event_n', 'event_code', 'station_code', ...
                     'magnitude', basis_field};

    for i = 1:numel(required_vars)
        if ~ismember(required_vars{i}, metric_all.Properties.VariableNames)
            error('metric_all 缺少字段：%s', required_vars{i});
        end
    end

    if any(isnan(metric_all.magnitude))
        warning('metric_all.magnitude 中存在 NaN，拟合时会自动剔除。');
    end

    %% ---------- 预处理 ----------
    Tbase = metric_all;

    % 构造唯一事件 ID，避免不同组内 event_n 重复
    Tbase.event_id = categorical( ...
        "G" + string(Tbase.group_m) + ...
        "_E" + string(Tbase.event_n) + ...
        "_" + string(Tbase.event_code));

    Tbase.group_m = categorical(Tbase.group_m);
    Tbase.station_code = categorical(string(Tbase.station_code));

    % 震级中心化
    if use_event_level_mag_mean % 平均震级按事件计算，而非记录
        event_mag_table = unique(Tbase(:, {'event_id', 'magnitude'}), 'rows');
        mag_mean = mean(event_mag_table.magnitude, 'omitnan');
    else
        mag_mean = mean(Tbase.magnitude, 'omitnan');
    end

    Tbase.magnitude_c = Tbase.magnitude - mag_mean;

    fprintf('\n============================================================\n');
    fprintf('开始混合效应模型：logY ~ logX * magnitude_c + group_m + (1 | event_id)\n');
    fprintf('Basis Y: %s\n', basis_field);
    fprintf('Magnitude center = %.4f\n', mag_mean);
    fprintf('Target periods: %s s\n', mat2str(target_periods));
    fprintf('============================================================\n');

    %% ---------- 容器 ----------
    nT = numel(target_periods);

    models = cell(nT, 1);
    fit_tables = cell(nT, 1);
    summary = table();

    %% ---------- 逐周期拟合 ----------
    for k = 1:nT
        T0 = target_periods(k);
        x_field = pgv_field_name(T0);

        if ~ismember(x_field, Tbase.Properties.VariableNames)
            warning('缺少字段 %s，跳过 T=%.3fs。', x_field, T0);
            continue;
        end

        Tfit = Tbase;

        y_raw = Tfit.(basis_field);
        x_raw = Tfit.(x_field);

        mask = isfinite(y_raw) & isfinite(x_raw) & ...
               y_raw > 0 & x_raw > 0 & ...
               isfinite(Tfit.magnitude_c);

        Tfit = Tfit(mask, :);

        Tfit.logY = log10(Tfit.(basis_field));
        Tfit.logX = log10(Tfit.(x_field));

        if height(Tfit) < 20
            warning('T=%.3fs 有效样本太少，跳过。N=%d', T0, height(Tfit));
            continue;
        end

        n_events = numel(unique(Tfit.event_id));
        n_groups = numel(unique(Tfit.group_m));
        n_stations = numel(unique(Tfit.station_code));

        fprintf('\n--- 拟合 T = %.3f s ---\n', T0);
        fprintf('N rows = %d, events = %d, groups = %d, stations = %d\n', ...
            height(Tfit), n_events, n_groups, n_stations);

        %% ---------- 拟合 ----------
        formula = 'logY ~ logX * magnitude_c + group_m + (1 | event_id)';

        lme = fitlme(Tfit, formula);

        models{k} = lme;
        fit_tables{k} = Tfit;

        coef = lme.Coefficients;

        % 可选：检查 MATLAB 实际给各项起了什么名字
        fprintf('模型系数项名称：\n');
        disp(coef.Name);

        beta_logX = get_coef_flexible(coef, {'logX'});
        beta_mag  = get_coef_flexible(coef, {'magnitude_c'});

        % 不同 MATLAB 版本可能是 logX:magnitude_c 或 magnitude_c:logX
        beta_inter = get_coef_flexible(coef, ...
            {'logX:magnitude_c', 'magnitude_c:logX'}, ...
            {'logX', 'magnitude_c'});

        [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme);
        [r2_cond, r2_fixed] = pseudo_r2_lme(lme, Tfit);

        new_row = table( ...
            T0, ...
            string(x_field), ...
            height(Tfit), ...
            n_events, ...
            n_groups, ...
            n_stations, ...
            mag_mean, ...
            beta_logX.Estimate, beta_logX.SE, beta_logX.pValue, ...
            beta_mag.Estimate, beta_mag.SE, beta_mag.pValue, ...
            beta_inter.Estimate, beta_inter.SE, beta_inter.pValue, ...
            r2_cond, ...
            r2_fixed, ...
            aic_val, ...
            bic_val, ...
            loglik_val, ...
            'VariableNames', { ...
                'T_sec', ...
                'x_field', ...
                'N_rows', ...
                'N_events', ...
                'N_groups', ...
                'N_stations', ...
                'magnitude_center', ...
                'beta_logX', 'SE_logX', 'p_logX', ...
                'beta_magnitude_c', 'SE_magnitude_c', 'p_magnitude_c', ...
                'beta_logX_mag', 'SE_logX_mag', 'p_logX_mag', ...
                'pseudo_R2_conditional', ...
                'pseudo_R2_fixed', ...
                'AIC', ...
                'BIC', ...
                'LogLikelihood'});

        summary = [summary; new_row]; %#ok<AGROW>

        disp(lme);

        fprintf('关键项：logX:magnitude_c = %.6f, p = %.4g\n', ...
            beta_inter.Estimate, beta_inter.pValue);

        %% ---------- 可选画图 ----------
        if make_plots
            figure('Color', 'w', 'Name', sprintf('LME residuals T=%.3fs', T0));
            plotResiduals(lme, 'fitted');
            title(sprintf('LME residuals, T = %.3f s', T0));
            grid on;
        end
    end

    %% ---------- 输出 ----------
    lme_result = struct();
    lme_result.models = models;
    lme_result.fit_tables = fit_tables;
    lme_result.summary = summary;
    lme_result.formula = 'logY ~ logX * magnitude_c + group_m + (1 | event_id)';
    lme_result.basis_field = basis_field;
    lme_result.target_periods = target_periods;
    lme_result.magnitude_center = mag_mean;

    fprintf('\n============================================================\n');
    fprintf('混合效应模型拟合完成。\n');
    fprintf('重点看 beta_logX_mag 与 p_logX_mag。\n');
    fprintf('如果 beta_logX_mag > 0 且显著，表示震级越大，log-log 斜率越大。\n');
    fprintf('============================================================\n');

    disp(summary);
end


function fn = pgv_field_name(T0)
% 周期 -> metric_all 字段名
%   1 -> pgv_T1_000_h
%   2 -> pgv_T2_000_h
%   5 -> pgv_T5_000_h

    T_str = sprintf('T%.3f', T0);
    T_str = strrep(T_str, '.', '_');
    fn = ['pgv_' T_str '_h'];
end


function out = get_coef_flexible(coef_table, exact_names, contains_terms)
% 安全提取系数。
%
% exact_names:
%   优先按完整名称匹配，例如 {'logX:magnitude_c','magnitude_c:logX'}
%
% contains_terms:
%   如果 exact_names 找不到，则找同时包含这些词的项，例如 {'logX','magnitude_c'}

    if nargin < 3
        contains_terms = {};
    end

    names = string(coef_table.Name);
    idx = [];

    % 1. 精确匹配
    for i = 1:numel(exact_names)
        idx = find(names == string(exact_names{i}), 1);
        if ~isempty(idx)
            break;
        end
    end

    % 2. 包含匹配
    if isempty(idx) && ~isempty(contains_terms)
        mask = true(size(names));

        for i = 1:numel(contains_terms)
            mask = mask & contains(names, string(contains_terms{i}));
        end

        idx = find(mask, 1);
    end

    if isempty(idx)
        warning('没有找到系数项：%s', strjoin(string(exact_names), ' / '));

        out = struct();
        out.Estimate = NaN;
        out.SE = NaN;
        out.pValue = NaN;
    else
        out = struct();
        out.Estimate = coef_table.Estimate(idx);
        out.SE = coef_table.SE(idx);
        out.pValue = coef_table.pValue(idx);
    end
end


% function [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme)
% % 安全读取 AIC / BIC / LogLikelihood。
% % 不同 MATLAB 版本不一定支持 lme.AIC 这种直接写法。
% 
%     aic_val = NaN;
%     bic_val = NaN;
%     loglik_val = NaN;
% 
%     try
%         crit = lme.ModelCriterion;
% 
%         if isprop_or_field(crit, 'AIC')
%             tmp = crit.AIC;
%             aic_val = tmp(1);
%         end
% 
%         if isprop_or_field(crit, 'BIC')
%             tmp = crit.BIC;
%             bic_val = tmp(1);
%         end
% 
%         if isprop_or_field(crit, 'LogLikelihood')
%             tmp = crit.LogLikelihood;
%             loglik_val = tmp(1);
%         end
% 
%     catch
%         % 某些旧版本可能没有 ModelCriterion
%         try
%             loglik_val = lme.LogLikelihood;
%         catch
%         end
%     end
% end
function [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme)
% 安全读取 AIC / BIC / LogLikelihood

    aic_val = NaN;
    bic_val = NaN;
    loglik_val = NaN;

    try
        crit = lme.ModelCriterion;

        % 如果是 table
        if istable(crit)
            aic_val = crit{1, 'AIC'};
            bic_val = crit{1, 'BIC'};
            loglik_val = crit{1, 'LogLikelihood'};
            return;
        end

        % 如果是对象 / dataset / 其他支持点索引的类型
        try
            tmp = crit.AIC;
            aic_val = tmp(1);
        catch
        end

        try
            tmp = crit.BIC;
            bic_val = tmp(1);
        catch
        end

        try
            tmp = crit.LogLikelihood;
            loglik_val = tmp(1);
        catch
        end

    catch
        % 某些旧版本可能没有 ModelCriterion
        try
            loglik_val = lme.LogLikelihood;
        catch
        end
    end
end

function tf = isprop_or_field(x, name)
    tf = false;

    try
        tf = isprop(x, name) || isfield(x, name);
    catch
        try
            x.(name); %#ok<VUNUS>
            tf = true;
        catch
            tf = false;
        end
    end
end


function [r2_cond, r2_fixed] = pseudo_r2_lme(lme, Tfit)
% 给 mixed model 做一个简单 pseudo-R2。
%
% pseudo_R2_conditional:
%   用包含随机效应的 fitted 值计算。
%
% pseudo_R2_fixed:
%   尝试用仅固定效应预测计算。
%   如果当前 MATLAB 版本不支持 Conditional=false，则返回 NaN。

    y = Tfit.logY;
    sst = sum((y - mean(y, 'omitnan')).^2);

    r2_cond = NaN;
    r2_fixed = NaN;

    if sst <= 0
        return;
    end

    % 条件预测：含随机效应
    try
        yhat_cond = fitted(lme);
        sse_cond = sum((y - yhat_cond).^2);
        r2_cond = 1 - sse_cond / sst;
    catch
    end

    % 固定效应预测：不含随机效应
    try
        yhat_fixed = predict(lme, Tfit, 'Conditional', false);
        sse_fixed = sum((y - yhat_fixed).^2);
        r2_fixed = 1 - sse_fixed / sst;
    catch
    end
end