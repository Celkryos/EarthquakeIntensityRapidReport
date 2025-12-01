function paired_stations = calc_gb_intensity(paired_stations, adata)
% CALC_GB_INTENSITY - 根据 GB/T 17742-2020 计算国标仪器地震烈度
%
% 输入:
%   paired_stations - 包含 .ew, .ns, .ud 索引的配对结构体
%   adata           - 包含原始数据的 cell 数组
%
% 输出:
%   paired_stations - 更新后的结构体，新增以下字段:
%       .gb_intensity  - 最终仪器烈度值 (保留一位小数)
%       .gb_pga        - 用于计算的合成 PGA (m/s^2)
%       .gb_pgv        - 用于计算的合成 PGV (m/s)
%       .I_A           - 加速度烈度计算值
%       .I_V           - 速度烈度计算值

    station_names = fieldnames(paired_stations);
    fprintf('--- 计算国标仪器烈度 (GB/T 17742-2020) ---\n');
    
    % 滤波器参数：国标要求 0.1 Hz ~ 10 Hz 带通 [GB/T 17742-2020 A.5]
    gb_fch = 0.1;       % 高通截止频率 0.1 Hz
    gb_T0  = 0.1;       % 低通截止周期 0.1s => 频率 10 Hz
    
    for i = 1:length(station_names)
        name = station_names{i};
        st_info = paired_stations.(name);
        
        try
            % 1. 检查数据有效性
            if isfield(st_info, 'is_valid') && ~st_info.is_valid
                continue; % 跳过无效台站
            end
            
            % 2. 检查是否有三个分量
            if ~isfield(st_info, 'ew') || ~isfield(st_info, 'ns') || ~isfield(st_info, 'ud')
                warning('台站 %s 缺少三分量数据，无法严格按国标计算，跳过。', name);
                continue;
            end
            
            % 提取三个分量的原始数据
            d_ew = adata{st_info.ew};
            d_ns = adata{st_info.ns};
            d_ud = adata{st_info.ud};
            
            % 3. 滤波与积分 (0.1-10Hz)
            % 注意：nb_filt 会重写 acc_untrimmed，所以必须现算现用
            
            % --- EW ---
            d_ew = nb_filt(d_ew, gb_fch, gb_T0);
            d_ew = acc2vel(d_ew); 
            
            % --- NS ---
            d_ns = nb_filt(d_ns, gb_fch, gb_T0);
            d_ns = acc2vel(d_ns);
            
            % --- UD ---
            d_ud = nb_filt(d_ud, gb_fch, gb_T0);
            d_ud = acc2vel(d_ud);
            
            % 4. 统一裁剪 (取交集长度)
            % 利用 trim_info 对齐三个分量
            len_ew = length(d_ew.velocity_cms); % acc2vel后velocity_cms由acc_untrimmed截取而来吗？
            % 检查 acc2vel 的输出，它生成 vel_untrimmed。
            % 为了安全，我们这里直接使用 untrimmed 数据并通过 trim_info 对齐
            % 实际上 nb_filt 的 padding 逻辑应该是一致的，只要原始长度差不多
            
            % 这里为了严谨，我们取三个分量 untrimmed 结果的最小长度
            % (通常 K-NET 同一台站记录长度是一样的)
            n_pts = min([length(d_ew.acc_untrimmed), length(d_ns.acc_untrimmed), length(d_ud.acc_untrimmed)]);
            
            % 提取前 n_pts 点 (去除尾部可能的差异)
            acc_ew = d_ew.acc_untrimmed(1:n_pts);
            acc_ns = d_ns.acc_untrimmed(1:n_pts);
            acc_ud = d_ud.acc_untrimmed(1:n_pts);
            
            vel_ew = d_ew.vel_untrimmed(1:n_pts);
            vel_ns = d_ns.vel_untrimmed(1:n_pts);
            vel_ud = d_ud.vel_untrimmed(1:n_pts);
            
            % 5. 三分量合成 [GB/T 17742-2020 A.6]
            acc_syn = sqrt(acc_ew.^2 + acc_ns.^2 + acc_ud.^2);
            vel_syn = sqrt(vel_ew.^2 + vel_ns.^2 + vel_ud.^2);
            
            % 6. 计算 PGA 和 PGV [GB/T 17742-2020 A.7]
            % 注意单位转换！输入数据通常是 gal 和 cm/s
            % 国标要求: Acceleration (m/s^2), Velocity (m/s)
            
            PGA_gal = max(acc_syn);
            PGV_cms = max(vel_syn);
            
            PGA_SI = PGA_gal / 100.0; % m/s^2
            PGV_SI = PGV_cms / 100.0; % m/s
            
            % 7. 计算仪器烈度 [GB/T 17742-2020 A.8]
            
            % 公式 A.5: Ia = 3.17 * log10(PGA) + 6.59
            if PGA_SI > 0
                I_A = 3.17 * log10(PGA_SI) + 6.59;
            else
                I_A = 0;
            end
            
            % 公式 A.6: Iv = 3.00 * log10(PGV) + 9.77
            if PGV_SI > 0
                I_V = 3.00 * log10(PGV_SI) + 9.77;
            else
                I_V = 0;
            end
            
            % 8. 判定逻辑 [GB/T 17742-2020 公式 A.7]
            % 如果 Ia >= 6.0 或 Iv >= 6.0，则 I = Iv
            if I_A >= 6.0 || I_V >= 6.0
                I_inst = I_V;
            else
                % 否则取平均
                I_inst = (I_A + I_V) / 2.0;
            end
            
            % 9. 边界处理 [A.8.3]
            if I_inst < 1.0
                I_inst = 1.0;
            elseif I_inst > 12.0
                I_inst = 12.0;
            end
            
            % 10. 结果存回
            % 结果可取小数点后一位 (标准建议)，这里我们存原始值，显示时再约
            paired_stations.(name).gb_intensity = I_inst;
            paired_stations.(name).gb_pga = PGA_SI; % m/s^2
            paired_stations.(name).gb_pgv = PGV_SI; % m/s
            paired_stations.(name).I_A = I_A;
            paired_stations.(name).I_V = I_V;
            
            % fprintf('台站 %s: 烈度 %.1f (Ia=%.2f, Iv=%.2f)\n', name, I_inst, I_A, I_V);
            
        catch ME
            warning('计算台站 %s 国标烈度时出错: %s', name, ME.message);
        end
    end
    
    fprintf('--- 国标烈度计算完成 ---\n');
end