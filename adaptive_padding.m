function [x_final, total_length, n_prepended] = adaptive_padding(x, varargin)
% ADAPTIVE_PADDING - 对信号进行自适应的“镜像-渐隐-零填充”处理
%
% Syntax: [x_final, total_length, n_prepended] = adaptive_padding(x, p)
%
% Description:
%   本函数根据信号两端的“干净”程度，自动调整镜像填充的长度和
%   渐隐窗口的参数，以最优方式为后续的FFT做准备。
%   该策略旨在保护原始信号的幅度信息不受窗口函数影响。
%
% Inputs:
%   x - 输入的列向量信号
%
% Optional Name-Value Pair Arguments:
%   'BoundaryPercent' - 用于计算边界能量的窗口占总长的百分比 (默认: 5)
%   'MaxMirrorPercent' - 允许的最大镜像长度占总长的百分比 (默认: 30)
%   'BERThreshold' - 判断边界是否“干净”的能量比阈值 (默认: 0.05)
%
% Outputs:
%   x_final     - 经过完整处理后，可以直接用于FFT的信号
%   total_length- x_final 的总长度 (通常是2的幂次)
%   n_prepended - 在原始信号前端增加的总点数(镜像+零填充)，用于后续截取

    % --- 1. 解析输入参数 ---
    p = inputParser;
    addRequired(p, 'x', @isnumeric);
    addParameter(p, 'BoundaryPercent', 5, @isnumeric);
    addParameter(p, 'MaxMirrorPercent', 30, @isnumeric);
    addParameter(p, 'BERThreshold', 0.05, @isnumeric);
    parse(p, x, varargin{:});

    x = x(:); % 确保是列向量
    npts = length(x);
    boundary_percent = p.Results.BoundaryPercent / 100;
    max_mirror_percent = p.Results.MaxMirrorPercent / 100;
    ber_threshold = p.Results.BERThreshold;

    % --- 2. 计算自适应指标 (BER) ---
    n_boundary = floor(npts * boundary_percent);
    if n_boundary < 1
        n_boundary = 1;
    end
    
    x_start_window = x(1:n_boundary);
    x_end_window = x(end-n_boundary+1:end);
    
    rms_start = sqrt(mean(x_start_window.^2));
    rms_end = sqrt(mean(x_end_window.^2));
    A_max = max(abs(x));
    
    if A_max == 0
        ber_metric = 0;
    else
        ber_metric = max(rms_start, rms_end) / A_max;
    end
disp(ber_metric);
    % --- 3. 根据BER自适应计算镜像长度 ---
    % ber_metric/ber_threshold 的比值被限制在[0, 1]区间
    % 当BER很小时，镜像比例接近0；当BER很大时，镜像比例接近最大值
    mirror_ratio = max_mirror_percent * min(ber_metric / ber_threshold, 1.0);

disp(mirror_ratio);
    
    

    n_mirror = floor(npts * mirror_ratio);
    
    % --- 4. 执行镜像填充 (偶填充，确保信号连续性) ---
    if n_mirror > 0
        % 注意镜像的索引，避免复制边界点自身
        mirrored_part_start = flipud(x(2 : n_mirror + 1));
        mirrored_part_end = flipud(x(end - n_mirror : end - 1));
        x_mirrored = [mirrored_part_start; x; mirrored_part_end];
    else
        x_mirrored = x;
    end
    
    % --- 5. 执行渐隐加窗 (Tapering) ---
    % 关键：窗口的渐隐区只作用于镜像段，平台区(值为1)覆盖整个原始信号
    if n_mirror > 0
        n_mirrored_total = length(x_mirrored);
        win = ones(n_mirrored_total, 1);
        
        % 创建一个长度为 2*n_mirror 的全余弦窗 
        taper_win = tukeywin(n_mirror * 2, 1);
        
        % 将窗口的前半部分应用到信号前端的镜像段
        win(1:n_mirror) = taper_win(1:n_mirror);
        
        % 将窗口的后半部分应用到信号后端的镜像段
        win(end-n_mirror+1:end) = taper_win(n_mirror+1:end);
        
        x_tapered = x_mirrored .* win;
    else
        % 如果没有镜像，可以施加一个非常温和的窗，例如1%
        win = tukeywin(npts, 0.01 * 2);
        x_tapered = x .* win;
    end
    
    % --- 6. 执行最终零填充 (Zero Padding) ---
    n_tapered = length(x_tapered);
    total_length = 2^nextpow2(n_tapered); % 找到下一个2的幂次
    
    % 将信号居中放置在零填充数组中
    n_pre_pad = floor((total_length - n_tapered) / 2);
    n_post_pad = total_length - n_tapered - n_pre_pad;
    
    x_final = [zeros(n_pre_pad, 1); x_tapered; zeros(n_post_pad, 1)];
    
    % 计算在原始信号前端增加的总点数，用于后续截取
    n_prepended = n_mirror + n_pre_pad;

end