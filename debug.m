%%  流程概述
%   读+配对
%   剔除截断记录
%   加窗→零填充
%   保留填充段的前提下积分
%   多项式回扣
%   统一剪裁
%   合成时程
clear; close all; clc;
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
T_list  = [5, 2, 1, 0.5,0.1];   % 周期
fc_high = [0.1,0.1,0.1,0.1,0.1];    % Hz
adata = multiband_filter_controller(adata, stations, T_list, fc_high);
%% 积
acc_fields = {'acc_T5_000s','acc_T2_000s','acc_T1_000s','acc_T0_500s','acc_T0_100s'};
for i = 1:numel(adata)
    if ~(adata{i}.is_valid), continue; end
    adata{i} = acc2vel(adata{i}, acc_fields);
end
%% 合成时程（水平 / 三分量）
stations = hcsc(adata, stations, T_list);
%% 计算国标烈度
stations = compute_intensity_gb(stations);
%% 插值并输出等高线图
% 1. 构造输入数组
%ziduan='pgv_T0_100';kong=true;
ziduan='pgv_T0_500_h';kong=true;
%ziduan='I_gb';kong=false;
names = fieldnames(stations);
lon_sta = [];
lat_sta = [];
val_sta = [];

for i = 1:numel(names)
    s = stations.(names{i});
    % 只要 is_valid的台站
    if ~isfield(s, 'is_valid') || ~s.is_valid
        continue;   % 跳
    end

    % 需存在相应字段
    if ~isfield(s, ziduan)
        continue;   % 跳
    end

    % 纳入插值列表
    lon_sta(end+1,1) = s.station_long;
    lat_sta(end+1,1) = s.station_lat;
    val_sta(end+1,1) = s.(ziduan);
end

n_valid = numel(val_sta); % 合法台站数

%

% 插值

[LonG, LatG, ZZ, metaZ] = interp2d_field( ...
    lon_sta, lat_sta, val_sta, ...
    'UseLog', kong, ...
    'GridStep', 0.02, ...
    'NumLevels', 12);
%
% 画图
figure;
plot_contour_field(LonG, LatG, ZZ, metaZ, lon_sta, lat_sta);




%%
%plotSpectrum(x, fs);
% 测试另一个函数
[Xq, Yq, Vq] = interp_seismic_field(lon_sta, lat_sta, val_sta, 'linear');
