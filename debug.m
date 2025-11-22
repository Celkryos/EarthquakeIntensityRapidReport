%%  流程概述
%   读+配对
%   整体性预处理（日本+美国那个流程，含：
%   自适应加窗及首尾pad填充（镜像填充→加窗→零填充）
%   保留填充段的前提下积分
%   多项式回扣
%   统一剪裁
%   合成时程

%% 批量读
%adata=batch_read_earthquake_data('D:\FAST\IWT0211103111509');
adata=batch_read_earthquake_data('D:\FAST\20250602035200knt');
%% 初步矫正
adata=chubu(adata);


%% 台站配对
[stations, unpaired] = pair_earthquake_stations(adata);

%% 参数设置

T0=1;
dt=0.01;
%%
%n=56;
for n=1:length(adata)
% 1.滤波
adata{n, 1} = nb_filt(adata{n, 1}, T0); 

% 2.积分
adata{n, 1} = acc2vel(adata{n, 1});

% 3.基线校正 
t_untrimmed = (0 : length(adata{n, 1}.disp_untrimmed)-1)' * (1/adata{n,1}.sampling_freq_hz);
t_normalized = (t_untrimmed - t_untrimmed(1)) / (t_untrimmed(end) - t_untrimmed(1));    %时间轴归一化，防止数值过小时报错亏秩

[corrected_disp_untrimmed,trend, ~] = polyfit_constrained_zerobase(...
    t_normalized, adata{n, 1}.disp_untrimmed, 6);

% (可选) 从校正后的位移反推速度和加速度
temp=diff([0; trend]) / dt;     %保长度

quduoxiang_vel_untrimmed = adata{n, 1}.vel_untrimmed-temp;
quduoxiang_acc_untrimmed = adata{n, 1}.acc_untrimmed-diff([0; temp]) / dt;

% 保存校正后的未裁剪位移
adata{n, 1}.disp_corrected_untrimmed = corrected_disp_untrimmed;


% 4.统一裁剪

info = adata{n, 1}.trim_info;
start_idx = info.n_prepended + 1;
end_idx   = info.n_prepended + info.original_length;

% 执行裁剪
adata{n, 1}.acceleration_gal_processed = adata{n, 1}.acc_untrimmed(start_idx:end_idx);
adata{n, 1}.velocity_cms = adata{n, 1}.vel_untrimmed(start_idx:end_idx);
adata{n, 1}.displacement_cm = adata{n, 1}.disp_corrected_untrimmed(start_idx:end_idx);

% 对裁剪后的位移再做一次detrend，移除因裁剪边界不为零而引入的最终漂移（但牵涉到永久位移时，这一步应谨慎）
%adata{n, 1}.displacement_cm = detrend(adata{n, 1}.displacement_cm, 1);

disp(max(abs(adata{n, 1}.acceleration_gal_processed))); %检查，显示处理后的pga


%plot(adata{n, 1}.time, adata{n, 1}.velocity_cms)
%disp(max(abs(adata{n, 1}.velocity_cms))); %检查，显示处理后的pgv
end

%% 批量合成时程
stations = hcsc(adata, stations);

%%
%plotSpectrum(x, fs);

