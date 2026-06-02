function ratio_result = fit_lme_pgv_ratio_basic(metric_all, varargin)
% FIT_LME_PGV_RATIO_BASIC
% 对 宽频/低频 PGV 比值 做混合效应模型分析。
%
% 默认模型：
%   log_ratio ~ logX_c * magnitude_c + group_m + (1 + logX_c | event_id)
%
% 其中：
%   log_ratio = log10(PGV_broad) - log10(PGV_low)
%             = log10(PGV_broad ./ PGV_low)
%
%   PGV_broad 默认是 pgv_T0_020_h，对应 0.1–50 Hz 宽频带
%   PGV_low   是 pgv_T1_000_h / pgv_T2_000_h / pgv_T5_000_h
%
% 重点看：
%   beta_logX_c:
%       在平均震级下，宽频/低频比值是否随低频 PGV 幅值变化
%
%   beta_magnitude_c:
%       在平均 logX 下，震级是否影响宽频/低频比值
%
%   beta_logX_mag:
%       震级是否改变 logX_c 对 log_ratio 的影响
%
%   geom_ratio:
%       宽频 PGV / 低频 PGV 的几何平均倍数
%
%   scatter_factor_1sigma:
%       ratio 的 1 sigma 乘法离散因子

    %% ---------- 参数 ----------
    p = inputParser;
    addRequired(p, 'metric_all', @istable);

    addParameter(p, 'TargetPeriods', [1 2 5], @isnumeric);

    % 宽频字段。不要叫 HighFreqField，因为 pgv_T0_020_h 是 0.1–50 Hz 宽频带
    addParameter(p, 'BroadBandField', 'pgv_T0_020_h', ...
        @(x)ischar(x) || isstring(x));

    addParameter(p, 'UseEventLevelMagnitudeMean', true, ...
        @(x)islogical(x) || isnumeric(x));

    addParameter(p, 'FitMethod', 'ML', ...
        @(x)ischar(x) || isstring(x));

    addParameter(p, 'MakePlots', true, ...
        @(x)islogical(x) || isnumeric(x));

    parse(p, metric_all, varargin{:});

    target_periods = p.Results.TargetPeriods;
    broad_field = char(p.Results.BroadBandField);
    use_event_level_mag_mean = logical(p.Results.UseEventLevelMagnitudeMean);
    fit_method = char(p.Results.FitMethod);
    make_plots = logical(p.Results.MakePlots);

    %% ---------- 基本检查 ----------
    required_vars = {'group_m', 'event_n', 'event_code', ...
                     'station_code', 'magnitude', broad_field};

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

    % 构造唯一事件 ID，避免不同 group 内 event_n 重复
    Tbase.event_id = categorical( ...
        "G" + string(Tbase.group_m) + ...
        "_E" + string(Tbase.event_n) + ...
        "_" + string(Tbase.event_code));

    Tbase.group_m = categorical(Tbase.group_m);
    Tbase.station_code = categorical(string(Tbase.station_code));

    % 震级中心化：默认每个事件只贡献一个震级，避免台站多的事件权重过大
    if use_event_level_mag_mean
        [G_event, ~] = findgroups(Tbase.event_id);
        event_mag = splitapply(@(x) mean(x, 'omitnan'), ...
                               Tbase.magnitude, G_event);
        mag_mean = mean(event_mag, 'omitnan');
    else
        mag_mean = mean(Tbase.magnitude, 'omitnan');
    end

    Tbase.magnitude_c = Tbase.magnitude - mag_mean;

    %% ---------- 模型公式 ----------
    formula = 'log_ratio ~ logX_c * magnitude_c + group_m + (1 + logX_c | event_id)';

    fprintf('\n============================================================\n');
    fprintf('开始 宽频/低频 ratio 混合效应模型\n');
    fprintf('Formula: %s\n', formula);
    fprintf('Broad-band field: %s\n', broad_field);
    fprintf('Magnitude center = %.4f\n', mag_mean);
    fprintf('Target periods: %s s\n', mat2str(target_periods));
    fprintf('FitMethod: %s\n', fit_method);
    fprintf('============================================================\n');

    %% ---------- 容器 ----------
    nT = numel(target_periods);

    models = cell(nT, 1);
    fit_tables = cell(nT, 1);
    coef_tables = cell(nT, 1);
    covariance_tables = cell(nT, 1);
    summary = table();

    %% ---------- 逐周期拟合 ----------
    for k = 1:nT
        T0 = target_periods(k);
        low_field = pgv_field_name(T0);

        if ~ismember(low_field, Tbase.Properties.VariableNames)
            warning('缺少字段 %s，跳过 T=%.3fs。', low_field, T0);
            continue;
        end

        Tfit = Tbase;

        broad_raw = Tfit.(broad_field);
        low_raw   = Tfit.(low_field);

        mask = isfinite(broad_raw) & isfinite(low_raw) & ...
               broad_raw > 0 & low_raw > 0 & ...
               isfinite(Tfit.magnitude_c);

        Tfit = Tfit(mask, :);

        Tfit.log_broad = log10(Tfit.(broad_field));
        Tfit.logX      = log10(Tfit.(low_field));

        % 核心响应变量：宽频/低频 PGV 比值
        Tfit.log_ratio = Tfit.log_broad - Tfit.logX;

        % 中心化 logX，使 magnitude_c 主效应更好解释
        logX_center = mean(Tfit.logX, 'omitnan');
        Tfit.logX_c = Tfit.logX - logX_center;

        if height(Tfit) < 20
            warning('T=%.3fs 有效样本太少，跳过。N=%d', T0, height(Tfit));
            continue;
        end

        n_events = numel(unique(Tfit.event_id));
        n_groups = numel(unique(Tfit.group_m));
        n_stations = numel(unique(Tfit.station_code));

        fprintf('\n--- 拟合 ratio, T = %.3f s ---\n', T0);
        fprintf('Low-frequency field: %s\n', low_field);
        fprintf('N rows = %d, events = %d, groups = %d, stations = %d\n', ...
            height(Tfit), n_events, n_groups, n_stations);
        fprintf('logX center = %.6f\n', logX_center);

        %% ---------- 拟合 LME ----------
        try
            lme = fitlme(Tfit, formula, 'FitMethod', fit_method);
        catch ME
            warning('T=%.3fs 模型拟合失败：%s', T0, ME.message);
            continue;
        end

        models{k} = lme;
        fit_tables{k} = Tfit;
        coef_tables{k} = lme.Coefficients;

        try
            covariance_tables{k} = lme.CovarianceParameters;
        catch
            covariance_tables{k} = [];
        end

        coef = lme.Coefficients;

        %% ---------- 提取关键系数 ----------
        beta_logX = get_coef_flexible(coef, {'logX_c'});
        beta_mag  = get_coef_flexible(coef, {'magnitude_c'});

        beta_inter = get_coef_flexible(coef, ...
            {'logX_c:magnitude_c', 'magnitude_c:logX_c'}, ...
            {'logX_c', 'magnitude_c'});

        [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme);
        [r2_cond, r2_fixed] = pseudo_r2_lme(lme, Tfit, 'log_ratio');

        %% ---------- ratio 描述统计 ----------
        ratio_mean = mean(Tfit.log_ratio, 'omitnan');
        ratio_std  = std(Tfit.log_ratio, 'omitnan');

        % 原始尺度解释
        geom_ratio = 10.^ratio_mean;
        scatter_factor_1sigma = 10.^ratio_std;

        %% ---------- 输出表 ----------
        new_row = table( ...
            T0, ...
            string(low_field), ...
            height(Tfit), ...
            n_events, ...
            n_groups, ...
            n_stations, ...
            mag_mean, ...
            logX_center, ...
            ratio_mean, ...
            ratio_std, ...
            geom_ratio, ...
            scatter_factor_1sigma, ...
            beta_logX.Estimate, ...
            beta_logX.SE, ...
            beta_logX.pValue, ...
            beta_mag.Estimate, ...
            beta_mag.SE, ...
            beta_mag.pValue, ...
            beta_inter.Estimate, ...
            beta_inter.SE, ...
            beta_inter.pValue, ...
            r2_cond, ...
            r2_fixed, ...
            aic_val, ...
            bic_val, ...
            loglik_val, ...
            'VariableNames', { ...
                'T_sec', ...
                'low_field', ...
                'N_rows', ...
                'N_events', ...
                'N_groups', ...
                'N_stations', ...
                'magnitude_center', ...
                'logX_center', ...
                'mean_log_ratio', ...
                'std_log_ratio', ...
                'geom_ratio', ...
                'scatter_factor_1sigma', ...
                'beta_logX_c', ...
                'SE_logX_c', ...
                'p_logX_c', ...
                'beta_magnitude_c', ...
                'SE_magnitude_c', ...
                'p_magnitude_c', ...
                'beta_logX_mag', ...
                'SE_logX_mag', ...
                'p_logX_mag', ...
                'pseudo_R2_conditional', ...
                'pseudo_R2_fixed', ...
                'AIC', ...
                'BIC', ...
                'LogLikelihood'});

        summary = [summary; new_row]; %#ok<AGROW>

        disp(lme);

        fprintf('关键项：logX_c = %.6f, p = %.4g\n', ...
            beta_logX.Estimate, beta_logX.pValue);

        fprintf('关键项：magnitude_c = %.6f, p = %.4g\n', ...
            beta_mag.Estimate, beta_mag.pValue);

        fprintf('关键项：logX_c:magnitude_c = %.6f, p = %.4g\n', ...
            beta_inter.Estimate, beta_inter.pValue);

        fprintf('几何平均 ratio: broad/low = %.4f\n', geom_ratio);
        fprintf('1 sigma 乘法离散因子: x/÷ %.4f\n', scatter_factor_1sigma);

        %% ---------- 可选图 ----------
        if make_plots
            % 1. 残差图
            figure('Color', 'w', ...
                   'Name', sprintf('Ratio LME residuals T=%.3fs', T0));
            plotResiduals(lme, 'fitted');
            title(sprintf('Ratio LME residuals, T = %.3f s', T0));
            grid on;

            % 2. log_ratio vs magnitude，按 group 着色
            figure('Color', 'w', ...
                   'Name', sprintf('log ratio vs magnitude T=%.3fs', T0));
            hold on; box on; grid on;

            g = categories(Tfit.group_m);
            cmap = lines(numel(g));

            for ig = 1:numel(g)
                mask_g = Tfit.group_m == g{ig};

                scatter(Tfit.magnitude(mask_g), Tfit.log_ratio(mask_g), ...
                    12, cmap(ig, :), 'filled', ...
                    'MarkerFaceAlpha', 0.35, ...
                    'MarkerEdgeAlpha', 0.35, ...
                    'DisplayName', ['G' char(g{ig})]);
            end

            xlabel('Magnitude');
            ylabel(sprintf('log_{10}(%s / %s)', ...
                strrep(broad_field, '_', '\_'), ...
                strrep(low_field, '_', '\_')));

            title(sprintf('log ratio vs Magnitude, T = %.3f s', T0));
            legend('Location', 'bestoutside');
            yline(0, 'k--', 'LineWidth', 1);
            hold off;

            % 3. log_ratio vs logX，按 group 着色
            figure('Color', 'w', ...
                   'Name', sprintf('log ratio vs logX T=%.3fs', T0));
            hold on; box on; grid on;

            for ig = 1:numel(g)
                mask_g = Tfit.group_m == g{ig};

                scatter(Tfit.logX(mask_g), Tfit.log_ratio(mask_g), ...
                    12, cmap(ig, :), 'filled', ...
                    'MarkerFaceAlpha', 0.35, ...
                    'MarkerEdgeAlpha', 0.35, ...
                    'DisplayName', ['G' char(g{ig})]);
            end

            xlabel(sprintf('log_{10}(%s)', strrep(low_field, '_', '\_')));
            ylabel(sprintf('log_{10}(%s / %s)', ...
                strrep(broad_field, '_', '\_'), ...
                strrep(low_field, '_', '\_')));

            title(sprintf('log ratio vs logX, T = %.3f s', T0));
            legend('Location', 'bestoutside');
            yline(0, 'k--', 'LineWidth', 1);
            hold off;

            % 4. log_ratio vs logX_c，按 magnitude 着色
            figure('Color', 'w', ...
                   'Name', sprintf('log ratio vs logX_c colored by M T=%.3fs', T0));
            scatter(Tfit.logX_c, Tfit.log_ratio, ...
                14, Tfit.magnitude, 'filled', ...
                'MarkerFaceAlpha', 0.45, ...
                'MarkerEdgeAlpha', 0.45);
            box on; grid on;
            colorbar;
            xlabel(sprintf('centered log_{10}(%s)', strrep(low_field, '_', '\_')));
            ylabel(sprintf('log_{10}(%s / %s)', ...
                strrep(broad_field, '_', '\_'), ...
                strrep(low_field, '_', '\_')));
            title(sprintf('log ratio vs logX_c colored by Magnitude, T = %.3f s', T0));
            yline(0, 'k--', 'LineWidth', 1);
        end
    end

    %% ---------- 输出 ----------
    ratio_result = struct();
    ratio_result.models = models;
    ratio_result.fit_tables = fit_tables;
    ratio_result.coef_tables = coef_tables;
    ratio_result.covariance_tables = covariance_tables;
    ratio_result.summary = summary;
    ratio_result.formula = formula;
    ratio_result.broad_field = broad_field;
    ratio_result.target_periods = target_periods;
    ratio_result.magnitude_center = mag_mean;
    ratio_result.fit_method = fit_method;

    fprintf('\n============================================================\n');
    fprintf('ratio 混合效应模型拟合完成。\n');
    fprintf('重点看 beta_logX_c、beta_magnitude_c、beta_logX_mag。\n');
    fprintf('同时关注 geom_ratio 与 scatter_factor_1sigma。\n');
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
%   优先按完整名称匹配，例如 {'logX_c:magnitude_c','magnitude_c:logX_c'}
%
% contains_terms:
%   如果 exact_names 找不到，则找同时包含这些词的项。

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


function [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme)
% 安全读取 AIC / BIC / LogLikelihood

    aic_val = NaN;
    bic_val = NaN;
    loglik_val = NaN;

    try
        crit = lme.ModelCriterion;

        if istable(crit)
            aic_val = crit{1, 'AIC'};
            bic_val = crit{1, 'BIC'};
            loglik_val = crit{1, 'LogLikelihood'};
            return;
        end

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
        try
            loglik_val = lme.LogLikelihood;
        catch
        end
    end
end


function [r2_cond, r2_fixed] = pseudo_r2_lme(lme, Tfit, response_name)
% 简单 pseudo-R2
%
% conditional:
%   使用包含随机效应的 fitted 值
%
% fixed:
%   使用仅固定效应预测

    y = Tfit.(response_name);
    sst = sum((y - mean(y, 'omitnan')).^2);

    r2_cond = NaN;
    r2_fixed = NaN;

    if sst <= 0
        return;
    end

    try
        yhat_cond = fitted(lme);
        sse_cond = sum((y - yhat_cond).^2);
        r2_cond = 1 - sse_cond / sst;
    catch
    end

    try
        yhat_fixed = predict(lme, Tfit, 'Conditional', false);
        sse_fixed = sum((y - yhat_fixed).^2);
        r2_fixed = 1 - sse_fixed / sst;
    catch
    end
end