function [rms_noise, rms_signal, ratio] = demo_check_truncated(acc, fs, noise_win_sec)
% DEMO_CHECK_TRUNCATED  Slide3 演示用：
% 输入一条加速度序列，自动计算噪声/强震 RMS，并画图标出窗口
%
% acc           加速度序列（列/行向量均可）
% fs            采样率 (Hz)
% noise_win_sec 噪声窗口长度，单位秒（可选，默认 5 s）
%
% 输出：
%   rms_noise   噪声窗口 RMS
%   rms_signal  强震窗口 RMS
%   ratio       rms_noise / rms_signal

if nargin < 3
    noise_win_sec = 5;   % 默认取前 5 秒作为噪声窗
end

acc = acc(:);                 % 转列向量
n   = length(acc);
t   = (0:n-1)'/fs;

%% 1. 计算 Arias 强度并找 5% / 95% 能量点
AI = cumsum(acc.^2) / fs;     % 这里少了常数，对“相对能量”足够
if AI(end) == 0
    AI(:) = 0;
else
    AI = AI / AI(end);        % 归一化到 0–1
end

idx5  = find(AI >= 0.05, 1, 'first');
idx95 = find(AI >= 0.95, 1, 'first');

% 防御：如果记录很短或几乎没能量，给个兜底
if isempty(idx5),  idx5  = round(0.1*n); end
if isempty(idx95), idx95 = round(0.9*n); end

%% 2. 定义噪声 & 强震窗口索引

% ---- 噪声窗：从开头往后 noise_win_sec 秒 ----
noise_len  = round(noise_win_sec * fs);
noise_start = 1;
% 尽量不让噪声窗踩到 5% 能量点之前的样本
noise_end   = min([noise_len, idx5-1, n]);
if noise_end < noise_start
    % 极端情况：5% 点就在最开始，那就只能强行取从头到 noise_len
    noise_end = min(noise_len, n);
end

noise_idx  = noise_start:noise_end;

% ---- 强震窗：5%–95% AI ----
signal_idx = idx5:idx95;

%% 3. 计算 RMS
rms_noise  = sqrt(mean(acc(noise_idx).^2));
rms_signal = sqrt(mean(acc(signal_idx).^2));
ratio      = rms_noise / max(rms_signal, eps);   % 防止除 0

%% 4. 画图：上面时程 + 窗口，下方 Arias 强度
figure('Color','w');

% (a) 加速度时程
subplot(2,1,1);
plot(t, acc, 'k'); hold on;
yl = ylim;

% 噪声窗（前 5 秒）
patch([t(noise_start) t(noise_end) t(noise_end) t(noise_start)], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.8 0.8 1.0], 'FaceAlpha',0.3, 'EdgeColor','none');   % 淡蓝

% 强震窗（5%–95%）
patch([t(idx5) t(idx95) t(idx95) t(idx5)], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [1.0 0.8 0.8], 'FaceAlpha',0.3, 'EdgeColor','none');   % 淡红

plot(t, acc, 'k');   % 再画一遍保证波形在最前面
xline(t(idx5), '--r', '5% AI');
xline(t(idx95),'--r', '95% AI');

title('加速度时程：噪声窗（前 5 s）与强震窗（5%–95% AI）');
xlabel('Time (s)');
ylabel('Acceleration');

legend({'加速度',...
        '噪声窗',...
        '强震窗',...
        '5% AI','95% AI'}, ...
        'Location','best');

text(0.01,0.95, sprintf('RMS_{noise}=%.3g\nRMS_{signal}=%.3g\nratio=%.3g', ...
     rms_noise, rms_signal, ratio), ...
     'Units','normalized', 'VerticalAlignment','top', ...
     'BackgroundColor','w','Margin',5);

% (b) Arias 强度
subplot(2,1,2);
plot(t, AI, 'k','LineWidth',1.2); hold on;
yline(0.05,'--r','5%');
yline(0.95,'--r','95%');
xline(t(idx5), '--r');
xline(t(idx95),'--r');

xlabel('Time (s)');
ylabel('Normalized AI');
title('归一化 Arias 强度');
grid on;

end
