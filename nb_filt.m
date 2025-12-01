function [acc_untrimmed, trim_info, freq, tf_bandpass] = nb_filt(data, fch, T0)
% NB_FILT - 频域零相位带通滤波（纯被调用函数，不修改输入 data）
%
% 输入：
%   data.sampling_freq_hz  - 采样频率 (Hz)
%   data.acceleration_jjh  - 加速度时程 (gal)，列向量或行向量
%   data.time              - (可选) 时间轴 (s)，仅用于可视化
%   fch                    - 高截止频率 (Hz)，对应你说的“最低频”
%   T0                     - "最高频率对应周期"(s)，低截止频率 fc_low = 1/T0
%
% 输出：
%   acc_untrimmed - 带 pad 的滤波后加速度时程 (gal)
%   trim_info     - 结构体，包含：
%                       .n_prepended      - 前面镜像/渐隐填充的样点数
%                       .original_length  - 原始信号长度 N
%                       .totallength      - FFT 用的总长度
%                       .dt               - 采样间隔
%                       .fc_low           - 低截止频率
%                       .fc_high          - 高截止频率
%   freq          - 频率轴 (0 ~ f_Nyq)，长度 Nnyq
%   tf_bandpass   - 对应的带通滤波器幅频 (与 freq 同长)

    %------------------ 基本参数 ------------------%
    if ~isfield(data, 'sampling_freq_hz')
        error('nb_filt:MissingField', 'data.sampling_freq_hz 缺失。');
    end
    if ~isfield(data, 'acceleration_jjh')
        error('nb_filt:MissingField', 'data.acceleration_jjh 缺失。');
    end

    fs   = data.sampling_freq_hz;
    xacc = data.acceleration_jjh(:);    % 列向量
    N    = length(xacc);

    keshihua = 0;   % 可视化开关

    % 时间轴（仅用于画图）
    if isfield(data,'time') && ~isempty(data.time)
        t = data.time(:);
        if length(t) ~= N
            error('data.time 长度与 acceleration_jjh 不一致。');
        end
    else
        t = (0:N-1)'/fs;
    end

    % 截止频率
    if T0 <= 0
        error('nb_filt:InvalidT0', 'T0 必须为正。');
    end
    fc_low  = 1 / T0;   % 低截止频率
    fc_high = fch;      % 高截止频率

    dt = 1/fs;

    %------------------ 自适应填充 (镜像+渐隐+零填充) ------------------%
    moshi = 0;    % 1 用 adaptive_padding，0 用 taper_zeropad

    if moshi
        [xpad, totallength, tbegin] = adaptive_padding(xacc);
    else
        [xpad, totallength, tbegin] = taper_zeropad(xacc, 'TaperPercent', 10);
    end

    %------------------ 频域滤波参数 ------------------%
    df   = 1/(totallength*dt);          % 频率分辨率
    Nnyq = totallength/2 + 1;           % 奈奎斯特频率对应点数
    freq = ((1:Nnyq) - 1)' * df;        % [0, f_Nyq]

    jie    = 4;
    n_low  = jie;
    n_high = jie;

    % 低/高通传递函数（幅度）
    tf_low  = sqrt(1.0 ./ (1.0 + (freq / fc_low ).^(2 * n_low )));
    tf_high = sqrt((freq / fc_high).^(2 * n_high) ./ (1.0 + (freq / fc_high).^(2 * n_high)));
    tf_high(1) = 0;                     % 直流分量置零
    tf_bandpass = tf_low .* tf_high;

    %------------------ 频域滤波 ------------------%
    % FFT
    X_freq = fft(xpad);

    % 只在 0~f_Nyq 区间乘以滤波器
    X_freq_filtered = zeros(size(X_freq));
    X_freq_filtered(1:Nnyq) = X_freq(1:Nnyq) .* tf_bandpass;

    % 构造共轭对称
    X_freq_filtered(1)     = real(X_freq_filtered(1));
    X_freq_filtered(Nnyq)  = real(X_freq_filtered(Nnyq));
    X_freq_filtered(totallength+2-(2:Nnyq-1)) = conj(X_freq_filtered(2:Nnyq-1));

    % IFFT 回时域
    y_acc_pad = real(ifft(X_freq_filtered));

    % 输出：保留 pad 段，不裁剪
    acc_untrimmed = y_acc_pad(:);

    % 填写 trim_info（只统计信息，不直接裁剪）
    trim_info = struct();
    trim_info.n_prepended     = tbegin;
    trim_info.original_length = N;
    trim_info.totallength     = totallength;
    trim_info.dt              = dt;
    trim_info.fc_low          = fc_low;
    trim_info.fc_high         = fc_high;

    %------------------ 可视化（可选） ------------------%
    if keshihua
        % 画原始时程和“裁剪后”的滤波结果
        idx_start = tbegin + 1;
        idx_end   = tbegin + N;
        y_acc = y_acc_pad(idx_start:idx_end);

        figure;
        subplot(2,1,1);
        plot(t, xacc, 'k-'); grid on;
        title('原始 acc'); xlabel('t (s)'); ylabel('gal');

        subplot(2,1,2);
        plot(t, y_acc, 'b-'); grid on;
        ttl = sprintf('频域滤波(acc): fclow=%.4f Hz (T0=%.3f s), fchigh=%.4f Hz, Butter-%d', ...
            fc_low, T0, fc_high, n_low);
        title(ttl); xlabel('t (s)'); ylabel('gal');
    end

end
