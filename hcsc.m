function paired_stations = hcsc(adata, paired_stations)
% 根据 sqrt(EW.^2 + NS.^2) 公式合成水平向时程。合成后的结果将作为新字段添加到输入的 paired_stations 结构体中。
%
% 输入:
%   adata           - 包含所有台站完整处理结果的 n*1 cell数组。
%   paired_stations - 由 pair_earthquake_stations 生成的配对信息结构体。
%
% 输出:
%   paired_stations - 更新后的结构体。每个台站字段下会新增以下字段：
%                     .acc_h_gal   (合成后的水平向加速度)
%                     .vel_h_cms   (合成后的水平向速度)
%                     .time        (与分量一致的时间轴)

    % 获取所有已配-对的台站名
    station_names = fieldnames(paired_stations);
    
    if isempty(station_names)
        warning('输入的 paired_stations 结构体为空，没有可合成的台站。');
        return;
    end
    
    fprintf('--- 开始批量合成 %d 个台站的水平分量 ---\n', length(station_names));
    
    for i = 1:length(station_names)
        station_name = station_names{i};
        
        try
            % 1. 获取该台站EW和NS分量在adata中的索引
            ew_index = paired_stations.(station_name).ew;
            ns_index = paired_stations.(station_name).ns;
            
            % 2. 提取处理后的时程数据
            %    我们假设处理后的数据存储在 'acceleration_gal_processed' 和 'velocity_cms' 字段中
            if ~isfield(adata{ew_index}, 'acceleration_gal_processed') || ...
               ~isfield(adata{ns_index}, 'acceleration_gal_processed')
                warning('跳过台站 %s: 缺少 "acceleration_gal_processed" 字段。', station_name);
                continue;
            end
            
            acc_ew = adata{ew_index}.acceleration_gal_processed;
            acc_ns = adata{ns_index}.acceleration_gal_processed;
            
            vel_ew = adata{ew_index}.velocity_cms;
            vel_ns = adata{ns_index}.velocity_cms;
            
            % 3. 检查数据长度是否一致 (一个重要的安全检查)
            if length(acc_ew) ~= length(acc_ns)
                warning('跳过台站 %s: EW 和 NS 分量加速度长度不一致。', station_name);
                continue;
            end
            if length(vel_ew) ~= length(vel_ns)
                warning('跳过台站 %s: EW 和 NS 分量速度长度不一致。', station_name);
                continue;
            end
            
            % 4. 执行时程合成
            acc_h = sqrt(acc_ew.^2 + acc_ns.^2);
            vel_h = sqrt(vel_ew.^2 + vel_ns.^2);
            
            % 5. 将合成结果和时间轴存回 paired_stations 结构体
            paired_stations.(station_name).acc_h_gal = acc_h;
            paired_stations.(station_name).vel_h_cms = vel_h;

            % 6. 记录pga和pgv
            paired_stations.(station_name).pga=max(abs(acc_h));     %其实abs可以去掉，因为按说不会有负值
            paired_stations.(station_name).pgv=max(abs(vel_h));
            % 随便从EW或NS分量里拿一个时间轴即可，因为它们是一样的
            paired_stations.(station_name).time = adata{ew_index}.time;
            
            fprintf('台站 %s 合成成功。\n', station_name);
            
        catch ME
            fprintf('处理台站 %s 时发生错误: %s\n', station_name, ME.message);
        end
    end
    
    fprintf('--- 批量合成完成 ---\n');
    
end