%% run_analysis_metric_G1.m
clear; close all; clc;

load('D:\FAST\metric_all_G5.mat', 'metric_all');
%%
result_pgv = analysis_metric_table(metric_all, ...
    'Groups', [], ...
    'AnalysisType', 'PGV', ...
    'BasisSuffix', 'T0_020_h', ...
    'MakePlots', true, ...
    'ColorByEvent', true);

summary_pgv = result_pgv.summary;

%%
metric_all=[g1;g2;g3;g4;g5;g6;g7;g8];

%% 混合效应模型 basic 震级影响

lme_result = fit_lme_pgv_basic(metric_all, ...
    'TargetPeriods', [1 2 5], ...
    'BasisField', 'pgv_T0_020_h', ...
    'MakePlots', true);

%% 混合效应模型  +事件随机斜率

lme_rs_result = fit_lme_pgv_event_random_slope(metric_all, ...
    'TargetPeriods', [1 2 5], ...
    'BasisField', 'pgv_T0_020_h', ...
    'MakePlots', true, ...
    'FitMethod', 'ML');

%% 从事件随机斜率 LME 模型中提取每个事件的条件斜率，并画 slope vs magnitude

event_slope_table = extract_event_random_slopes(lme_rs_result, ...
    'MakePlot', true);

%%

figure;
hold on; box on; grid on;

periods = unique(event_slope_table.T_sec, 'stable');
cmap = lines(numel(periods));
markers = {'o', 's', '^'};

stat_text = strings(numel(periods), 1);

for k = 1:numel(periods)
    T0 = periods(k);
    mask_T = event_slope_table.T_sec == T0;

    x = event_slope_table.magnitude(mask_T);
    y = event_slope_table.random_slope_logX(mask_T);

    mask = isfinite(x) & isfinite(y);

    scatter(x(mask), y(mask), 48, markers{k}, ...
        'MarkerEdgeColor', cmap(k,:), ...
        'LineWidth', 1.1, ...
        'DisplayName', sprintf('T = %.0f s random slope', T0));

    mdl = fitlm(x(mask), y(mask));

    xfit = linspace(min(x(mask)), max(x(mask)), 100)';
    yfit = predict(mdl, xfit);

    plot(xfit, yfit, '-', ...
        'Color', cmap(k,:), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('T = %.0f s fit', T0));

    [rval, ~] = corr(x(mask), y(mask), 'Type', 'Pearson');
    R2 = mdl.Rsquared.Ordinary;
    pval = mdl.Coefficients.pValue(2);

    stat_text(k) = sprintf('T = %.0f s: r=%.3f, R^2=%.3f, p=%.3g', ...
        T0, rval, R2, pval);
end

yline(0, 'k--', 'LineWidth', 1);

xlabel('Magnitude');
ylabel('Random slope of logX');
title('Residual event random slope vs Magnitude');

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
%%  ratio 严谨版

ratio_result = fit_lme_pgv_ratio_basic(metric_all, ...
    'TargetPeriods', [1 2 5], ...
    'BroadBandField', 'pgv_T0_020_h', ...
    'FitMethod', 'ML', ...
    'MakePlots', true);
%%
for k = 1:3
    Tfit = ratio_result.fit_tables{k};
    lme = ratio_result.models{k};

    y = Tfit.log_ratio;
    yhat_fixed = predict(lme, Tfit, 'Conditional', false);
    yhat_cond = fitted(lme);

    rmse_fixed = sqrt(mean((y - yhat_fixed).^2, 'omitnan'));
    rmse_cond = sqrt(mean((y - yhat_cond).^2, 'omitnan'));

    fprintf('T%d: RMSE fixed = %.4f, factor = %.3f; RMSE cond = %.4f, factor = %.3f\n', ...
        k, rmse_fixed, 10^rmse_fixed, rmse_cond, 10^rmse_cond);
end


%% ratio 考虑振幅

ratio_logx_result = fit_lme_pgv_ratio_with_logx(metric_all, ...
    'TargetPeriods', [1 2 5], ...
    'HighFreqField', 'pgv_T0_020_h', ...
    'FitMethod', 'ML', ...
    'MakePlots', true);
%% 区分地理组

figure;
tiledlayout(1, numel(periods), 'TileSpacing', 'compact');

for k = 1:numel(periods)
    T0 = periods(k);
    mask_T = event_slope_table.T_sec == T0;

    nexttile;
    boxchart(categorical(event_slope_table.group_m(mask_T)), ...
             event_slope_table.random_slope_logX(mask_T));

    yline(0, 'k--');
    grid on; box on;
    xlabel('Group');
    ylabel('Random slope');
    title(sprintf('T = %.0f s', T0));
end

sgtitle('Event random slope by group');