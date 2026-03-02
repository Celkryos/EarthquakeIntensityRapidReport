%% 批量地震动参数相关性分析
% 功能：
% 1. 支持切换 PGA / PGV 分析模式
% 2. 以最短周期（通常为 T=0.02s）为基准，分析与其他周期的相关性
% 3. 绘制所有周期的回归散点图 (Figure 1)
% 4.绘制 R^2 随周期变化的趋势图 (Figure 2)
%
% 依赖变量：stations (结构体)

clearvars -except stations adata; 
clc; close all;

%% ================= 配置区域 =================

% 1. 分析类型选择: 'PGA' 或 'PGV'
analysis_type = 'PGA'; 

% 2. 基准周期的后缀（通常选最短周期 T=0.02s 作为基准）
%    脚本会自动拼接前缀，如 'pga_' + basis_suffix
basis_suffix  = 'T0_020_h'; 

% 3. 是否使用 Log10 坐标 (推荐 true)
use_log = true;

% ===========================================

% 根据配置生成字段名模式
switch upper(analysis_type)
    case 'PGA'
        basis_field    = ['pga_' basis_suffix];
        target_pattern = 'pga_T.*_h';
    case 'PGV'
        basis_field    = ['pgv_' basis_suffix];
        target_pattern = 'pgv_T.*_h';
    otherwise
        error('analysis_type 必须是 ''PGA'' 或 ''PGV''');
end

fprintf('--- 开始批量分析 (%s 模式) ---\n', analysis_type);
fprintf('基准变量 (X): %s\n', basis_field);

%% 1. 数据提取与准备
names = fieldnames(stations);
is_valid_sta = false(numel(names), 1);
all_vals_basis = nan(numel(names), 1);

% 提取基准数据 (X)
for i = 1:numel(names)
    s = stations.(names{i});
    if isfield(s, 'is_valid') && s.is_valid && isfield(s, basis_field)
        val = s.(basis_field);
        if ~isempty(val) && (~use_log || val > 0)
            all_vals_basis(i) = val;
            is_valid_sta(i) = true;
        end
    end
end

if sum(is_valid_sta) < 5
    error('有效数据点太少 (<5)，请检查 basis_field 设置是否正确。');
end

% 查找所有目标字段 (Y)
sample_sta = stations.(names{find(is_valid_sta, 1)});
all_fields = fieldnames(sample_sta);
match_idx  = ~cellfun(@isempty, regexp(all_fields, target_pattern));
target_fields = all_fields(match_idx);

% 排除基准字段自身
target_fields = setdiff(target_fields, {basis_field});

% 解析周期并排序
T_vals = zeros(numel(target_fields), 1);
for k = 1:numel(target_fields)
    t_str = regexp(target_fields{k}, 'T(\d+)_(\d+)', 'tokens');
    if ~isempty(t_str)
        sec = str2double(t_str{1}{1});
        dec = str2double(t_str{1}{2});
        T_vals(k) = sec + dec/1000;
    else
        T_vals(k) = 999; % 解析失败排到最后
    end
end
[T_vals, sort_idx] = sort(T_vals);
target_fields = target_fields(sort_idx);

fprintf('共找到 %d 个目标字段，开始拟合...\n', numel(target_fields));

%% 2. 循环拟合与绘制散点图 (Figure 1)
n_targets = numel(target_fields);
colors = parula(n_targets);   % 

% 记录 R2 数据用于 Figure 2
R2_history = zeros(n_targets, 1);
T_history  = T_vals;

% 创建散点图窗口
fig1 = figure('Name', ['Regression Analysis - ' analysis_type], 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.1 0.6 0.8]);

% 自动计算子图行列
n_cols = ceil(sqrt(n_targets));
n_rows = ceil(n_targets / n_cols);

for k = 1:n_targets
    y_field = target_fields{k};
    col = colors(k, :);
    
    % 提取 Y 数据
    vals_y = nan(size(all_vals_basis));
    mask   = is_valid_sta;
    
    for i = 1:numel(names)
        if ~mask(i), continue; end
        s = stations.(names{i});
        if isfield(s, y_field)
            v = s.(y_field);
            if ~isempty(v) && (~use_log || v > 0)
                vals_y(i) = v;
            else
                mask(i) = false; 
            end
        else
            mask(i) = false;
        end
    end
    
    X = all_vals_basis(mask);
    Y = vals_y(mask);
    
    if isempty(X)
        R2_history(k) = NaN;
        continue; 
    end
    
    % 坐标转换
    if use_log
        XX = log10(X);
        YY = log10(Y);
        label_prefix = 'log_{10}';
    else
        XX = X;
        YY = Y;
        label_prefix = '';
    end
    
    % 线性回归
    mdl = fitlm(XX, YY);
    b = mdl.Coefficients.Estimate(1); 
    k_slope = mdl.Coefficients.Estimate(2); 
    r2 = mdl.Rsquared.Ordinary;
    
    % 记录 R2
    R2_history(k) = r2;
    
    % --- 绘图 (Subplot) ---
    figure(fig1); % 确保画在 fig1 上
    subplot(n_rows, n_cols, k);
    hold on; box on; grid on;
    
    % 散点与拟合线
    scatter(XX, YY, 12, col, 'filled', 'MarkerFaceAlpha', 0.5);
    
    x_min = min(XX); x_max = max(XX);
    x_range = x_max - x_min;
    x_fit = linspace(x_min, x_max, 100);
    y_fit = k_slope * x_fit + b;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5);
    
    % 标题 (T=xxx s)
    title_str = sprintf('T = %.3f s', T_vals(k));
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    % 坐标标签
    if k > (n_rows-1)*n_cols % 最后一行显示X轴标签
        xlabel([label_prefix ' (' strrep(basis_field, '_', '\_') ')'], 'FontSize', 8);
    end
    if mod(k, n_cols) == 1 % 第一列显示Y轴标签
        ylabel([label_prefix ' (Target)'], 'FontSize', 8);
    end
    
    % 显示公式与 R2 (精确到小数点后4位)
    y_range = max(YY) - min(YY);
    text_x = min(XX) + 0.05 * x_range;
    text_y = max(YY) - 0.15 * y_range;
    
    eq_str = sprintf('y=%.2fx%+.2f\nR^2=%.4f', k_slope, b, r2);
    text(text_x, text_y, eq_str, 'BackgroundColor', 'w', 'EdgeColor', 'none', ...
        'FontSize', 8, 'Margin', 1);
    
    axis tight;
    xlim([min(XX)-0.1*x_range, max(XX)+0.1*x_range]);
    ylim([min(YY)-0.1*y_range, max(YY)+0.1*y_range]);
    hold off;
end
sgtitle(fig1, ['Correlation Analysis: ' analysis_type ' (Base: ' strrep(basis_field, '_', '\_') ')'], ...
    'FontSize', 14, 'FontWeight', 'bold');

%% 3. 绘制 R^2 随周期变化趋势图 (Figure 2 - 新窗口)
% 过滤掉无效点 (T=999 或 NaN)
valid_idx = (T_vals < 900) & ~isnan(R2_history);
T_plot = T_vals(valid_idx);
R2_plot = R2_history(valid_idx);

if ~isempty(T_plot)
    fig2 = figure('Name', ['R2 vs Period - ' analysis_type], 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.66 0.1 0.3 0.4]); % 放在屏幕右侧
    
    plot(T_plot, R2_plot, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'b', 'Color', 'b');
    grid on; box on;
    
    xlabel('Period T (s)');
    ylabel('Coefficient of Determination (R^2)');
    title(['R^2 vs. Period (' analysis_type ')']);
    
    % 标记每个点的数值 (精确到4位)
    for i = 1:length(T_plot)
        text(T_plot(i), R2_plot(i), sprintf(' %.4f', R2_plot(i)), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
            'FontSize', 9, 'Rotation', 45);
    end
    
    % 适当调整 Y 轴范围，留出文字空间
    ylim([min(R2_plot)-0.05, 1.05]);
    xlim([0, max(T_plot)*1.1]);
    
    fprintf('--- 绘图完成: 查看 Figure 1 (散点) 和 Figure 2 (R2趋势) ---\n');
else
    warning('没有有效的 T 和 R2 数据用于绘制趋势图。');
end