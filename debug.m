%%  流程概述
%   读+配对
%   剔除截断记录
%   加窗→零填充
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

%% 剔除截断记录
[adata, stations, rejected_info] = check_truncated_records(adata, stations);

%% 多频段滤波
T_list  = [5, 2, 1, 0.5];   % 
fc_high = 0.1;
adata = multiband_filter_controller(adata, stations, T_list, fc_high);
%% 积
acc_fields = {'acc_T5_000s','acc_T2_000s','acc_T1_000s','acc_T0_500s'};
for i = 1:numel(adata)
    if ~(adata{i}.is_valid), continue; end
    adata{i} = acc2vel(adata{i}, acc_fields);
end
%%
% 4.统一裁剪

info = adata{n, 1}.trim_info;
start_idx = info.n_prepended + 1;
end_idx   = info.n_prepended + info.original_length;

% 执行裁剪
adata{n, 1}.acceleration_gal_processed = adata{n, 1}.acc_untrimmed(start_idx:end_idx);
adata{n, 1}.velocity_cms = adata{n, 1}.vel_untrimmed(start_idx:end_idx);
adata{n, 1}.displacement_cm = adata{n, 1}.disp_corrected_untrimmed(start_idx:end_idx);

disp(max(abs(adata{n, 1}.acceleration_gal_processed))); %检查，显示处理后的pga


%plot(adata{n, 1}.time, adata{n, 1}.velocity_cms)
%disp(max(abs(adata{n, 1}.velocity_cms))); %检查，显示处理后的pgv
end

%% 批量合成时程
stations = hcsc(adata, stations);

%%
%plotSpectrum(x, fs);

