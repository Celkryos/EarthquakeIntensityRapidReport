function [paired_stations, unpaired_info] = pair_earthquake_stations(adata)
% PAIR_EARTHQUAKE_STATIONS - 将地震动数据的EW和NS分量进行配对
%
%   接收一个由 batch_read_earthquake_data 生成的cell数组 adata，
%   并根据台站代码 (station_code) 将EW和NS两个水平分量进行配对。
%   函数会优先尝试快速匹配(i vs i+n/2)，失败后会进行全局搜索。:
% 输入：adata
%
% 输出:
%   paired_stations - 一个结构体，其字段名为台站代码。每个字段下包含 .ew 和 .ns 两个属性，值为该分量在adata 中的原始索引。
%   unpaired_info   - 一个结构体数组，记录了所有未能成功配对的台站信息，包括台站代码、原始索引和分量方向。

    % 初始化
    nn = length(adata);
    if nn == 0 || mod(nn, 2) ~= 0
        warning('输入数据为空或数量为奇数，可能无法正确配对。');
    end
    n_half = floor(nn / 2);
    
    paired_stations = struct();
    unpaired_list = []; % 用于临时存储未配对的台站信息
    
    % 创建一个“已处理”标记数组，避免重复匹配
    processed_flags = false(nn, 1);

    % 优先尝试快速匹配 (i vs i+n/2)
    fprintf('--- 开始快速匹配 ---\n');
    for i = 1:n_half
        % 获取当前EW分量的信息
        ew_candidate_idx = i;
        ew_station_code = adata{ew_candidate_idx}.station_code;
        
        % 获取对应的NS分量候选信息
        ns_candidate_idx = i + n_half;
        ns_station_code = adata{ns_candidate_idx}.station_code;
        
        % 检查台站代码是否匹配
        if strcmp(ew_station_code, ns_station_code)
            % 匹配成功！
            station_code = ew_station_code;
            % 确保字段名是合法的MATLAB变量名
            station_field_name = matlab.lang.makeValidName(station_code);
            paired_stations.(station_field_name).ew = ew_candidate_idx;
            paired_stations.(station_field_name).ns = ns_candidate_idx;
            
            % 标记这两个台站为已处理
            processed_flags(ew_candidate_idx) = true;
            processed_flags(ns_candidate_idx) = true;
            
            fprintf('快速匹配成功: %s (EW: %d, NS: %d)\n', ...
                station_code, ew_candidate_idx, ns_candidate_idx);
        end
    end
    
    % --- 3. 对所有未处理的台站进行全局搜索匹配 ---
    unprocessed_indices = find(~processed_flags);
    if ~isempty(unprocessed_indices)
        fprintf('--- 快速匹配完成，对剩余 %d 个记录进行全局搜索 ---\n', length(unprocessed_indices));
        
        for i_idx = 1:length(unprocessed_indices)
            current_idx = unprocessed_indices(i_idx);
            
            % 如果当前记录在本次循环中已经被配对，则跳过
            if processed_flags(current_idx)
                continue;
            end
            
            current_code = adata{current_idx}.station_code;
            current_dir = adata{current_idx}.direction;
            
            found_partner = false;
            
            % 在所有其他未处理的记录中寻找配对
            for j_idx = i_idx + 1 : length(unprocessed_indices)
                partner_idx = unprocessed_indices(j_idx);
                
                % 如果伙伴已被处理，则跳过
                if processed_flags(partner_idx)
                    continue;
                end
                
                partner_code = adata{partner_idx}.station_code;
                partner_dir = adata{partner_idx}.direction;
                
                % 判断条件：台站代码相同，且方向不同
                if strcmp(current_code, partner_code) && ~strcmp(current_dir, partner_dir)
                    % 配对成功！
                    station_code = current_code;
                    station_field_name = matlab.lang.makeValidName(station_code);
                    
                    % 判断哪个是EW，哪个是NS
                    if contains(upper(current_dir), 'E-W')
                        paired_stations.(station_field_name).ew = current_idx;
                        paired_stations.(station_field_name).ns = partner_idx;
                    else
                        paired_stations.(station_field_name).ew = partner_idx;
                        paired_stations.(station_field_name).ns = current_idx;
                    end
                    
                    % 标记两者都已处理
                    processed_flags(current_idx) = true;
                    processed_flags(partner_idx) = true;
                    
                    fprintf('全局匹配成功: %s (Indices: %d, %d)\n', ...
                        station_code, current_idx, partner_idx);
                        
                    found_partner = true;
                    break; % 已找到伙伴，跳出内层循环
                end
            end
            
            % 如果遍历完所有都找不到伙伴，则记录为未配对
            if ~found_partner
                unpaired_list(end+1).station_code = current_code;
                unpaired_list(end).index = current_idx;
                unpaired_list(end).direction = current_dir;
            end
        end
    end
    
    % --- 4. 整理并输出最终结果 ---
    unpaired_info = unpaired_list;
    num_paired = length(fieldnames(paired_stations));
    num_unpaired = length(unpaired_info);
    
    fprintf('--- 配对完成 ---\n');
    fprintf('成功配对台站: %d 个\n', num_paired);
    fprintf('未能配对记录: %d 个\n', num_unpaired);
    
    if num_unpaired > 0
        fprintf('未能配对的台站列表:\n');
        for k = 1:num_unpaired
            fprintf('  - Code: %s, Index: %d, Dir: %s\n', ...
                unpaired_info(k).station_code, unpaired_info(k).index, unpaired_info(k).direction);
        end
    end
end