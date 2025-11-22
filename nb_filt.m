function data = nb_filt(data, T0)%------------------ 参数 ------------------%
% 用频域滤波对加速度滤波（零相位）
% 截止频率 fc = 1/T0。结果直接存回 data，字段名 nb_acc_T0
%
% 输入：
%   data.sampling_freq_hz   - 采样频率 (Hz)
%   data.acceleration_jjh   - 加速度时程 (gal)
%   data.time               - (可选) 时间轴 (s)
%   T0                      - "最高频率对应周期"(s)，fc = 1/T0
%
% 输出：
%   data.nb_acc_T0  - 滤波后的加速度（gal）
    
    %------------------ 参数 ------------------%
    fs   = data.sampling_freq_hz;
    xacc = data.acceleration_jjh(:);
    N    = length(xacc);
    keshihua = 0; % 可视化开关
    
    % 时间轴
    if ~isfield(data,'time') || isempty(data.time)
        data.time = (0:N-1)'/fs;
    else
        data.time = data.time(:);
    end
    if length(data.time) ~= N
        error('data.time 长度与 acceleration_jjh 不一致。');
    end

    % 截止频率
    if T0 <= 0, error('T0 必须为正。'); end

    fc_low = 1 / T0;                   % 低通截止频率 
    fc_high = 0.05;                    % 高通截止频率

    %fc_high = zishiying();            % 根据信噪比自适应调整最低频

    dt = 1/fs;                         % 时间步长

    %------------------ 频域滤波参数 ------------------%
    % 滤波器阶数 当信号较短时需考虑自适应调整 当前数据集（9.0）可固定6阶
    if T0==10 
        jie=2;
    elseif T0==5
        jie=3;
    elseif T0==2
        jie=4;
    else
        jie=6;
    end

    jie=4;

    
    n_low = jie;   % 低通滤波器阶数
    n_high = jie;  % 高通滤波器阶数
    
    %------------------ 准备时域信号 (去趋势 & 自适应填充) ------------------%
    % 1. 去均值和线性趋势（这一步已分进单独的初步预处理函数）
    % x_d = detrend(xacc, 1);
    
    % 2. 自适应 镜像填充→渐隐→零填充
    % 确定最终用于FFT的信号 xpad 和它的长度 totallength

    moshi=1;    % 内置变量 1采用我的处理函数 0采用经典流程函数
    
    if moshi
        [xpad, totallength, tbegin] = adaptive_padding(xacc);
    else
        [xpad, totallength, tbegin] = taper_zeropad(xacc, 'TaperPercent', 10);
    end


    % 3. 设计频域滤波器
    df = 1/(totallength*dt);           % 频率分辨率
    Nnyq = totallength/2 + 1;          % 奈奎斯特频率对应点数
    freq = ((1:Nnyq) - 1)' * df;       % 频率向量 (建议转为列向量以匹配tf_bandpass')
    %disp(df);

    % 4. 滤波器传递函数
    tf_low = sqrt(1.0 ./ (1.0 + (freq / fc_low).^(2 * n_low)));
    tf_high = sqrt((freq / fc_high).^(2 * n_high) ./ (1.0 + (freq / fc_high).^(2 * n_high)));
    tf_high(1) = 0; % 明确将直流分量(freq(1)=0)设置为0
    tf_bandpass = tf_low .* tf_high;

    %------------------ 频域滤波 ------------------%
    % 5. FFT变换到频域
    X_freq = fft(xpad);
    data.fft=X_freq;
    % 6. 应用滤波器 
    X_freq_filtered = zeros(size(X_freq));
    X_freq_filtered(1:Nnyq) = X_freq(1:Nnyq) .* tf_bandpass; 
    
    
    data.lvboqi=tf_bandpass;


    % 7. 构造共轭对称 
    X_freq_filtered(1) = real(X_freq_filtered(1));
    X_freq_filtered(Nnyq) = real(X_freq_filtered(Nnyq));
    X_freq_filtered(totallength+2-(2:Nnyq-1)) = conj(X_freq_filtered(2:Nnyq-1));
    
    % 8. IFFT变换回时域
    y_acc_pad = real(ifft(X_freq_filtered));

    data.acc_untrimmed = y_acc_pad;

    % 9. 提取有效信号部分 (使用 adaptive_padding 返回的 tbegin)
    y_acc = y_acc_pad(tbegin+1 : tbegin+N);

    % 写回 data（字段名带 T0，替换小数点为下划线）
    fname_acc = strrep(sprintf('nb_acc_%g', T0), '.', '_');
    data.(fname_acc) = y_acc;



    data.trim_info.n_prepended = tbegin;
    data.trim_info.original_length = N;

    %------------------ 可视化 ------------------%
    if keshihua
        figure;
        subplot(2,1,1);
        plot(data.time, xacc, 'k-'); grid on;
        title('原始 acc'); xlabel('t (s)'); ylabel('gal');

        subplot(2,1,2);
        plot(data.time, y_acc, 'b-'); grid on;
        ttl = sprintf('频域滤波(acc): fclow=%.4f Hz (T0=%.3f s), fchigh=%.4f Hz, Butter-%d', fc_low, T0, fc_high, n_low);
        title(ttl); xlabel('t (s)'); ylabel('gal');

    end
end