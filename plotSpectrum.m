function plotSpectrum(y, Fs, fmin)
% plotSpectrum: 绘制信号 y 的单边振幅谱，可选高通去 DC。
%
% y    : 输入信号
% Fs   : 采样频率 (Hz)
% fmin : (可选) 高通截止频率 (Hz)，默认 0 (不滤波)

    if nargin < 3
        fmin = 0;   % 默认不做高通
    end
y = detrend(y, 0);
    % --- 可选高通滤波，去除 DC 与低频漂移 ---
    if fmin > 0
        Wn = fmin / (Fs/2);                % 归一化截止频率
        hpFilt = designfilt('highpassiir', ...
                            'FilterOrder', 4, ...
                            'HalfPowerFrequency', Wn, ...
                            'DesignMethod', 'butter');
        y = filtfilt(hpFilt, y);           % 零相位高通滤波
    end

    % --- FFT 计算 ---
    L = length(y); 
    Y = fft(y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(L/2))/L;

    % --- 绘图 ---
    plot(f, P1)
    if fmin > 0
        title(sprintf('单边振幅谱 (高通 %.2f Hz)', fmin));
    else
        title('单边振幅谱 (未高通)');
    end
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    grid on
end
