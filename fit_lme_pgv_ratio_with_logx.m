function ratio_result = fit_lme_pgv_ratio_with_logx(metric_all, varargin)
% FIT_LME_PGV_RATIO_WITH_LOGX
% 对高频/低频 PGV 比值做混合效应模型分析。
%
% 模型：
%   log_ratio ~ magnitude_c + logX + group_m + (1 + logX | event_id)
%
% 其中：
%   log_ratio = log10(pgv_T0_020_h) - log10(pgv_Tx_xxx_h)
%   logX      = log10(pgv_Tx_xxx_h)
%
% 重点看：
%   beta_magnitude_c:
%       控制 logX 后，震级是否影响高频/低频比值
%
%   beta_logX:
%       高频/低频比值是否随低频 PGV 水平变化

    %% ---------- 参数 ----------
    p = inputParser;
    addRequired(p, 'metric_all', @istable);

    addParameter(p, 'TargetPeriods', [1 2 5], @isnumeric);
    addParameter(p, 'HighFreqField', 'pgv_T0_020_h', @(x)ischar(x) || isstring(x));
    addParameter(p, 'UseEventLevelMagnitudeMean', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'FitMethod', 'ML', @(x)ischar(x) || isstring(x));
    addParameter(p, 'MakePlots', true, @(x)islogical(x) || isnumeric(x));

    parse(p, metric_all, varargin{:});

    target_periods = p.Results.TargetPeriods;
    high_field = char(p.Results.HighFreqField);
    use_event_level_mag_mean = logical(p.Results.UseEventLevelMagnitudeMean);
    fit_method = char(p.Results.FitMethod);
    make_plots = logical(p.Results.MakePlots);

    %% ---------- 基本检查 ----------
    required_vars = {'group_m', 'event_n', 'event_code', ...
                     'station_code', 'magnitude', high_field};

    for i = 1:numel(required_vars)
        if ~ismember(required_vars{i}, metric_all.Properties.VariableNames)
            error('metric_all 缺少字段：%s', required_vars{i});
        end
    end

    %% ---------- 预处理 ----------
    Tbase = metric_all;

    Tbase.event_id = categorical( ...
        "G" + string(Tbase.group_m) + ...
        "_E" + string(Tbase.event_n) + ...
        "_" + string(Tbase.event_code));

    Tbase.group_m = categorical(Tbase.group_m);
    Tbase.station_code = categorical(string(Tbase.station_code));

    % 震级中心化：默认按事件级均值
    if use_event_level_mag_mean
        event_mag_table = unique(Tbase(:, {'event_id', 'magnitude'}), 'rows');
        mag_mean = mean(event_mag_table.magnitude, 'omitnan');
    else
        mag_mean = mean(Tbase.magnitude, 'omitnan');
    end

    Tbase.magnitude_c = Tbase.magnitude - mag_mean;

    formula = 'log_ratio ~ magnitude_c + logX + group_m + (1 + logX | event_id)';

    fprintf('\n============================================================\n');
    fprintf('开始 ratio + logX 混合效应模型\n');
    fprintf('Formula: %s\n', formula);
    fprintf('High-frequency field: %s\n', high_field);
    fprintf('Magnitude center = %.4f\n', mag_mean);
    fprintf('Target periods: %s s\n', mat2str(target_periods));
    fprintf('FitMethod: %s\n', fit_method);
    fprintf('============================================================\n');

    %% ---------- 容器 ----------
    nT = numel(target_periods);

    models = cell(nT, 1);
    fit_tables = cell(nT, 1);
    coef_tables = cell(nT, 1);
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

        high_raw = Tfit.(high_field);
        low_raw  = Tfit.(low_field);

        mask = isfinite(high_raw) & isfinite(low_raw) & ...
               high_raw > 0 & low_raw > 0 & ...
               isfinite(Tfit.magnitude_c);

        Tfit = Tfit(mask, :);

        Tfit.log_high = log10(Tfit.(high_field));
        Tfit.log_low  = log10(Tfit.(low_field));

        % 响应变量与自变量
        Tfit.log_ratio = Tfit.log_high - Tfit.log_low;
        Tfit.logX = Tfit.log_low;

        if height(Tfit) < 20
            warning('T=%.3fs 有效样本太少，跳过。N=%d', T0, height(Tfit));
            continue;
        end

        n_events = numel(unique(Tfit.event_id));
        n_groups = numel(unique(Tfit.group_m));
        n_stations = numel(unique(Tfit.station_code));

        fprintf('\n--- 拟合 ratio + logX, T = %.3f s ---\n', T0);
        fprintf('N rows = %d, events = %d, groups = %d, stations = %d\n', ...
            height(Tfit), n_events, n_groups, n_stations);

        %% ---------- 拟合 LME ----------
        lme = fitlme(Tfit, formula, 'FitMethod', fit_method);

        models{k} = lme;
        fit_tables{k} = Tfit;
        coef_tables{k} = lme.Coefficients;

        coef = lme.Coefficients;

        beta_mag  = get_coef_flexible(coef, {'magnitude_c'});
        beta_logX = get_coef_flexible(coef, {'logX'});

        [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme);
        [r2_cond, r2_fixed] = pseudo_r2_lme(lme, Tfit, 'log_ratio');

        ratio_mean = mean(Tfit.log_ratio, 'omitnan');
        ratio_std  = std(Tfit.log_ratio, 'omitnan');

        new_row = table( ...
            T0, ...
            string(low_field), ...
            height(Tfit), ...
            n_events, ...
            n_groups, ...
            n_stations, ...
            mag_mean, ...
            ratio_mean, ...
            ratio_std, ...
            beta_mag.Estimate, ...
            beta_mag.SE, ...
            beta_mag.pValue, ...
            beta_logX.Estimate, ...
            beta_logX.SE, ...
            beta_logX.pValue, ...
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
                'mean_log_ratio', ...
                'std_log_ratio', ...
                'beta_magnitude_c', ...
                'SE_magnitude_c', ...
                'p_magnitude_c', ...
                'beta_logX', ...
                'SE_logX', ...
                'p_logX', ...
                'pseudo_R2_conditional', ...
                'pseudo_R2_fixed', ...
                'AIC', ...
                'BIC', ...
                'LogLikelihood'});

        summary = [summary; new_row]; %#ok<AGROW>

        disp(lme);

        fprintf('关键项 magnitude_c = %.6f, p = %.4g\n', ...
            beta_mag.Estimate, beta_mag.pValue);

        fprintf('关键项 logX = %.6f, p = %.4g\n', ...
            beta_logX.Estimate, beta_logX.pValue);

        %% ---------- 可选图 ----------
        if make_plots
            % 残差图
            figure('Color', 'w', 'Name', sprintf('Ratio + logX residuals T=%.3fs', T0));
            plotResiduals(lme, 'fitted');
            title(sprintf('Ratio + logX LME residuals, T = %.3f s', T0));
            grid on;

            % log_ratio vs logX，按 group 着色
            figure('Color', 'w', 'Name', sprintf('log ratio vs logX T=%.3fs', T0));
            hold on; box on; grid on;

            g = categories(Tfit.group_m);
            cmap = lines(numel(g));

            for ig = 1:numel(g)
                mask_g = Tfit.group_m == g{ig};

                scatter(Tfit.logX(mask_g), Tfit.log_ratio(mask_g), ...
                    12, cmap(ig, :), 'filled', ...
                    'MarkerFaceAlpha', 0.35, ...
                    'MarkerEdgeAlpha', 0.35, ...
                    'DisplayName', ['G' char(g{ig})]);
            end

            xlabel(['log_{10}(' strrep(low_field, '_', '\_') ')']);
            ylabel(sprintf('log_{10}(%s / %s)', ...
                strrep(high_field, '_', '\_'), strrep(low_field, '_', '\_')));

            title(sprintf('log ratio vs logX, T = %.3f s', T0));
            legend('Location', 'bestoutside');
            yline(0, 'k--', 'LineWidth', 1);

            hold off;
        end
    end

    %% ---------- 输出 ----------
    ratio_result = struct();
    ratio_result.models = models;
    ratio_result.fit_tables = fit_tables;
    ratio_result.coef_tables = coef_tables;
    ratio_result.summary = summary;
    ratio_result.formula = formula;
    ratio_result.high_field = high_field;
    ratio_result.target_periods = target_periods;
    ratio_result.magnitude_center = mag_mean;
    ratio_result.fit_method = fit_method;

    fprintf('\n============================================================\n');
    fprintf('ratio + logX 混合效应模型拟合完成。\n');
    fprintf('重点看 beta_magnitude_c / p_magnitude_c 以及 beta_logX / p_logX。\n');
    fprintf('============================================================\n');

    disp(summary);
end


function fn = pgv_field_name(T0)
    T_str = sprintf('T%.3f', T0);
    T_str = strrep(T_str, '.', '_');
    fn = ['pgv_' T_str '_h'];
end


function out = get_coef_flexible(coef_table, exact_names, contains_terms)
    if nargin < 3
        contains_terms = {};
    end

    names = string(coef_table.Name);
    idx = [];

    for i = 1:numel(exact_names)
        idx = find(names == string(exact_names{i}), 1);
        if ~isempty(idx)
            break;
        end
    end

    if isempty(idx) && ~isempty(contains_terms)
        mask = true(size(names));
        for i = 1:numel(contains_terms)
            mask = mask & contains(names, string(contains_terms{i}));
        end
        idx = find(mask, 1);
    end

    if isempty(idx)
        warning('没有找到系数项：%s', strjoin(string(exact_names), ' / '));
        out.Estimate = NaN;
        out.SE = NaN;
        out.pValue = NaN;
    else
        out.Estimate = coef_table.Estimate(idx);
        out.SE = coef_table.SE(idx);
        out.pValue = coef_table.pValue(idx);
    end
end


function [aic_val, bic_val, loglik_val] = get_model_criterion_safe(lme)
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