function [adata, paired_stations, rejected_info] = check_truncated_records(adata, paired_stations, varargin)
% CHECK_TRUNCATED_RECORDS - 检查并剔除截断记录，并为保留记录添加 lag, dist, noise_rms 参数
%
% 输入:
%   adata           - 包含地震数据的 cell 数组
%   paired_stations - 配对好的台站结构体
%   varargin        - 可选参数对:
%       'Vp'            - P波平均波速 (km/s), 默认 5.0
%       'CheckWindow'   - 检查起始噪声的时间窗长度 (s)
%       'NoiseThreshold'- 起始RMS/PGA 的阈值，超过此值且时间晚到则判为截断。
%       'timelag'       - 允许的最大滞后时间阈值 (s)
%
% 输出:
%   adata           - 更新后的数据 (包含 .lag, .dist_km, .noise_rms, .is_valid)
%   paired_stations - 更新后的配对结构体 (包含上述字段)
%   rejected_info   - 被剔除台站的详细信息列表

    % --- 参数解析 ---
    p = inputParser;
    addRequired(p, 'adata');
    addRequired(p, 'paired_stations');
    addParameter(p, 'Vp', 5.0, @isnumeric);           % km/s
    addParameter(p, 'CheckWindow', 2.0, @isnumeric);  % 秒
    addParameter(p, 'NoiseThreshold', 0.15, @isnumeric); % 阈值
    addParameter(p, 'timelag', -2.0, @isnumeric); % 时间滞后容忍阈值
    parse(p, adata, paired_stations, varargin{:});
    
    Vp = p.Results.Vp;
    win_len = p.Results.CheckWindow;
    noise_ratio_thres = p.Results.NoiseThreshold;
    lag_thres = p.Results.timelag;
    
    station_names = fieldnames(paired_stations);
    rejected_info = struct('station', {}, 'reason', {});
    
    fprintf('--- 开始进行截断记录筛查 (Vp=%.1f km/s, LagThres=%.1fs) ---\n', Vp, lag_thres);
    
    count_rejected = 0;
    
    for i = 1:length(station_names)
        name = station_names{i};
        
        % 获取 EW 和 NS 的索引
        idx_ew = paired_stations.(name).ew;
        idx_ns = paired_stations.(name).ns;
        idx_ud = paired_stations.(name).ud;
        % 提取数据 (以EW为例进行计算)
        data = adata{idx_ud}; 
        
        try
            % 1. 时间解析与校正
            fmt = 'yyyy/MM/dd HH:mm:ss';
            t_origin = datetime(data.origin_time, 'InputFormat', fmt);
            t_record_header = datetime(data.record_time, 'InputFormat', fmt);
            
            % 扣除 15s Pre-trigger delay
            t_start_real = t_record_header - seconds(15);
            
            % 2. 计算震中距和理论到时
            dist_km = epi_station_distance(data.latitude, data.longitude, ...
                                           data.station_lat, data.station_long, 'km');
            
            travel_time = dist_km / Vp;
            t_arrival_theo = t_origin + seconds(travel_time);
            
            % 计算滞后时间 (秒)
            % 正值 = 记录开始晚于理论P波到时 (可能丢失初至)
            time_lag = seconds(t_start_real - t_arrival_theo);
            
            % 3. 计算起始段能量 (RMS)
            % 确保使用已进行初步基线校正的数据 acceleration_jjh
            fs = data.sampling_freq_hz;
            n_samples_check = round(win_len * fs);
            if n_samples_check > length(data.acceleration_jjh)
                n_samples_check = length(data.acceleration_jjh);
            end
            
            acc_segment = detrend(data.acceleration_jjh(1:n_samples_check), 'constant');
            start_rms = rms(acc_segment);
            
            full_pga = max(abs(detrend(data.acceleration_jjh, 'constant')));
            if full_pga == 0, full_pga = eps; end
            
            noise_ratio = start_rms / full_pga;
            
            % 4. 综合判别逻辑
            is_time_late  = time_lag > lag_thres;
            is_high_noise = noise_ratio > noise_ratio_thres;
            highhigh_noise= noise_ratio>(noise_ratio_thres+0.05);
            
            if (is_time_late && is_high_noise)||highhigh_noise
                % --- 剔除分支 ---
                reason = sprintf('疑似截断记录: 滞后 %.2fs (>%.1f), 起始RMS占比 %.1f%%', ...
                    time_lag, lag_thres, noise_ratio*100);
                
                count_rejected = count_rejected + 1;
                rejected_info(count_rejected).station = name;
                rejected_info(count_rejected).reason  = reason;
                
                fprintf('剔除台站 %s: %s (Dist: %.1f km)\n', name, reason, dist_km);
                
                % 标记为无效
                adata{idx_ew}.is_valid = false;
                adata{idx_ns}.is_valid = false;
                adata{idx_ud}.is_valid = false;
                
                paired_stations.(name).is_valid = false;
                
            else
                % --- 保留分支：将计算出的参数写回结构体 ---
                
                % 1. 更新 paired_stations
                paired_stations.(name).lag       = time_lag;
                paired_stations.(name).dist_km   = dist_km;
                paired_stations.(name).noise_rms = start_rms; % 新增：写入初始噪声RMS
                paired_stations.(name).is_valid  = true;
                
                % 2. 更新 adata (EW分量)
                adata{idx_ew}.lag       = time_lag;
                adata{idx_ew}.dist_km   = dist_km;       % 新增：写入距离
                adata{idx_ew}.noise_rms = start_rms;     % 新增：写入初始噪声RMS
                adata{idx_ew}.is_valid  = true;
                
                % 3. 更新 adata (NS分量) - 赋予相同的值
                adata{idx_ns}.lag       = time_lag;
                adata{idx_ns}.dist_km   = dist_km;
                adata{idx_ns}.noise_rms = start_rms; 
                adata{idx_ns}.is_valid  = true;

                % 4. 更新 adata (UD分量) - 赋予相同的值
                adata{idx_ud}.lag       = time_lag;
                adata{idx_ud}.dist_km   = dist_km;
                adata{idx_ud}.noise_rms = start_rms; 
                adata{idx_ud}.is_valid  = true;
            end
            
        catch ME
            warning('处理台站 %s 时出错: %s', name, ME.message);
        end
    end
    
    % 统计有效台站数量
    names = fieldnames(paired_stations);
    valid_mask = false(size(names));
    for i = 1:numel(names)
        name = names{i};
        if isfield(paired_stations.(name), 'is_valid') && paired_stations.(name).is_valid
            valid_mask(i) = true;
        end
    end
    n_valid = sum(valid_mask);

    fprintf('--- 筛查完成，共剔除 %d 个台站，剩余 %d 个有效台站 ---\n', ...
        count_rejected, n_valid);
end