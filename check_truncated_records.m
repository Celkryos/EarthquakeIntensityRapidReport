function [adata, paired_stations, rejected_info] = check_truncated_records(adata, paired_stations, varargin)
% CHECK_TRUNCATED_RECORDS - 仅基于信噪比(SNR)筛查并剔除低质量记录
%
% 输入:
%   adata           - 包含地震数据的 cell 数组
%   paired_stations - 配对好的台站结构体
%   varargin        - 可选参数对:
%       'CheckWindow'    - 噪声时间窗长度 (s)，默认 2.0
%       'NoiseThreshold' - 噪声RMS/强震RMS 的阈值，超过则剔除（默认 0.2）
%       'Vp'             - (已弃用) 仅为兼容旧调用保留，不再参与判据
%       'timelag'        - (已弃用) 仅为兼容旧调用保留，不再参与判据
%
% 输出:
%   adata           - 更新后的数据 (包含 .dist_km, .noise_rms, .noise_ratio, .is_valid)
%   paired_stations - 更新后的配对结构体 (包含上述字段)
%   rejected_info   - 被剔除台站的详细信息列表

    % --- 参数解析 ---
    p = inputParser;
    addRequired(p, 'adata');
    addRequired(p, 'paired_stations');
    addParameter(p, 'CheckWindow', 2.0, @isnumeric);  % 秒
    addParameter(p, 'NoiseThreshold', 0.3, @isnumeric); % 阈值
    % 兼容旧接口：保留但不使用
    addParameter(p, 'Vp', 4.5, @isnumeric);
    addParameter(p, 'timelag', 55.0, @isnumeric);
    parse(p, adata, paired_stations, varargin{:});

    win_len = p.Results.CheckWindow;
    noise_ratio_thres = p.Results.NoiseThreshold;
    
    station_names = fieldnames(paired_stations);
    rejected_info = struct('station', {}, 'reason', {});
    
    fprintf('--- 开始筛查 (SNR筛查, NoiseRatioThres=%.3g, NoiseWin=%.1fs) ---\n', noise_ratio_thres, win_len);
    
    count_rejected = 0;
    
    for i = 1:length(station_names)
        name = station_names{i};
        
        % 获取索引（EW/NS 必需，UD 可选）
        idx_ew = paired_stations.(name).ew;
        idx_ns = paired_stations.(name).ns;
        if isfield(paired_stations.(name), 'ud')
            idx_ud = paired_stations.(name).ud;
        else
            idx_ud = [];
        end

        % 参考头部信息：优先用 EW，否则用 NS
        data_ref = adata{idx_ew};
        if isempty(data_ref)
            data_ref = adata{idx_ns};
        end
        
        try
            % 1) 计算震中距（可选：便于日志/后续分析；与筛查判据无关）
            dist_km = NaN;
            if isfield(data_ref, 'latitude') && isfield(data_ref, 'longitude') && ...
               isfield(data_ref, 'station_lat') && isfield(data_ref, 'station_long')
                dist_km = epi_station_distance(data_ref.latitude, data_ref.longitude, ...
                    data_ref.station_lat, data_ref.station_long, 'km');
            end


            %% 3 & 4. 计算 Arias 强度并确定噪声窗口

            fs = data_ref.sampling_freq_hz;
            dt = 1 / fs;

            % 仅使用水平分量：EW/NS 合成幅值来做 SNR/窗选取
            if isfield(adata{idx_ew}, 'acceleration_jjh')
                acc_ew0 = adata{idx_ew}.acceleration_jjh(:);
            else
                acc_ew0 = adata{idx_ew}.acceleration_gal(:);
            end
            if isfield(adata{idx_ns}, 'acceleration_jjh')
                acc_ns0 = adata{idx_ns}.acceleration_jjh(:);
            else
                acc_ns0 = adata{idx_ns}.acceleration_gal(:);
            end

            n_total = min(numel(acc_ew0), numel(acc_ns0));
            acc_ew0 = acc_ew0(1:n_total);
            acc_ns0 = acc_ns0(1:n_total);

            acc_full= sqrt(acc_ew0.^2 + acc_ns0.^2);
            n_total  = length(acc_full);

            % 先计算 Arias 强度，获取 idx_5 / idx_95
            ai     = cumsum(acc_full.^2) * dt;
            ai_max = ai(end);

            if ai_max < eps
                % 极端情况：信号全 0 或能量极小
                idx_5  = n_total;
                idx_95 = n_total;
            else
                ai_norm = ai / ai_max;

                idx_5  = find(ai_norm >= 0.05, 1, 'first');
                idx_95 = find(ai_norm >= 0.95, 1, 'first');

                if isempty(idx_5),  idx_5  = 1;      end
                if isempty(idx_95), idx_95 = n_total; end
            end

            % 噪声窗口长度：取 win_len 对应点数 与 idx_5 的较小值
            % 防止强震早到污染噪声窗口
            n_noise_samples = round(win_len * fs);
            n_noise_samples = min([n_noise_samples, idx_5 - 1, n_total]);
            n_noise_samples = max(n_noise_samples, 1);  % 至少 1 个样本

            acc_noise_seg = acc_full(1:n_noise_samples);
            start_rms     = rms(acc_noise_seg);

            % 强震段 RMS（5%–95% 能量窗口）
            if idx_95 - idx_5 < 2
                strong_rms = rms(acc_full);              % 窗口太短就退回全段
            else
                strong_rms = rms(acc_full(idx_5:idx_95));
            end

            %% 5. 噪声比例
            strong_rms  = max(strong_rms, eps);
            noise_ratio = start_rms / strong_rms;

            % 判别：仅基于噪声RMS / 强震RMS
            is_high_noise = noise_ratio > noise_ratio_thres;

            if is_high_noise
                % --- 剔除分支 ---
                reason = sprintf('低信噪比: 噪声RMS占比 %.1f%% (>%.1f%%)', ...
                    noise_ratio*100, noise_ratio_thres*100);
                
                count_rejected = count_rejected + 1;
                rejected_info(count_rejected).station = name;
                rejected_info(count_rejected).reason  = reason;
                
                fprintf('剔除台站 %s: %s (Dist: %.1f km)\n', name, reason, dist_km);
                
                % 标记为无效
                adata{idx_ew}.is_valid = false;
                adata{idx_ns}.is_valid = false;
                if ~isempty(idx_ud) && idx_ud <= numel(adata) && ~isempty(adata{idx_ud})
                    adata{idx_ud}.is_valid = false;
                end
                
                paired_stations.(name).is_valid = false;
                
            else
                % --- 保留分支：将计算出的参数写回结构体 ---
                
                % 1. 更新 paired_stations
                paired_stations.(name).dist_km   = dist_km;
                paired_stations.(name).noise_rms = start_rms; % 新增：写入初始噪声RMS
                paired_stations.(name).noise_ratio = noise_ratio;
                paired_stations.(name).is_valid  = true;
                
                % 2. 更新 adata (EW分量)
                adata{idx_ew}.dist_km   = dist_km;       % 新增：写入距离
                adata{idx_ew}.noise_rms = start_rms;     % 新增：写入初始噪声RMS
                adata{idx_ew}.noise_ratio = noise_ratio;
                adata{idx_ew}.is_valid  = true;
                
                % 3. 更新 adata (NS分量) - 赋予相同的值
                adata{idx_ns}.dist_km   = dist_km;
                adata{idx_ns}.noise_rms = start_rms; 
                adata{idx_ns}.noise_ratio = noise_ratio;
                adata{idx_ns}.is_valid  = true;

                % 4. 更新 adata (UD分量) - 可选
                if ~isempty(idx_ud) && idx_ud <= numel(adata) && ~isempty(adata{idx_ud})
                    adata{idx_ud}.dist_km   = dist_km;
                    adata{idx_ud}.noise_rms = start_rms;
                    adata{idx_ud}.noise_ratio = noise_ratio;
                    adata{idx_ud}.is_valid  = true;
                end
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