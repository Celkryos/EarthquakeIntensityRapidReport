function [x_final, total_length, n_prepended] = taper_zeropad(x, varargin)
% TAPER_ZEROPAD - 对信号进行直接加窗并零填充（日本PPT推荐流程）
%
% Syntax: [x_final, total_length, n_prepended] = taper_zeropad(x, varargin)
%
% Description:
%   本函数对输入信号直接施加一个余弦渐隐窗（Tukey window），然后进行零填充
%
% Inputs:
%   x - 输入的列向量信号
%
% Optional Name-Value Pair Arguments:
%   'TaperPercent' - 渐隐窗占总信号长度的百分比 (默认: 10, 即两端各5%)
%   'PadFactor'    - 填充因子，总长度至少为 npts*PadFactor (默认: 2)
%
% Outputs:
%   x_final     - 经过完整处理后，可以直接用于FFT的信号
%   total_length- x_final 的总长度 (通常补到2的幂次)
%   n_prepended - 在原始信号前端增加的零填充点数，用于后续截取

    % --- 1. 解析输入参数 ---
    p = inputParser;
    addRequired(p, 'x', @isnumeric);
    addParameter(p, 'TaperPercent', 10, @isnumeric); % 10% total taper is a common choice
    addParameter(p, 'PadFactor', 2, @isnumeric);
    parse(p, x, varargin{:});

    x = x(:); % 确保是列向量
    npts = length(x);
    taper_percent = p.Results.TaperPercent / 100;
    pad_factor = p.Results.PadFactor;

    % --- 2. 直接在原始信号上施加渐隐窗 ---
    % tukeywin的第二个参数是渐隐部分占总长的比例
    win = tukeywin(npts, taper_percent);
    x_tapered = x .* win;
    
    % --- 3. 执行零填充 (Zero Padding) ---
    % 确保填充后的总长度足够长，以容纳滤波器的瞬态响应
    % 然后再扩展到下一个2的幂次
    min_length = npts * pad_factor;
    total_length = 2^nextpow2(min_length);
    
    % 只在尾部零填充（避免前端平移），保持信号索引一致性
    n_post_pad = total_length - npts;
    
    x_final = [x_tapered; zeros(n_post_pad, 1)];
    
    % 前端无填充，n_prepended为0，后续处理无需修正索引
    n_prepended = 0;

end