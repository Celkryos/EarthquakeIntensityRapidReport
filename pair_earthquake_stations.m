function [paired_stations, unpaired_info] = pair_earthquake_stations(adata)
% PAIR_EARTHQUAKE_STATIONS - 将地震动数据的EW、NS和UD分量进行配对
%
%   接收一个由 batch_read_earthquake_data 生成的cell数组 adata，
%   并根据台站代码 (station_code) 将三个分量进行配对。
%   函数会优先尝试快速匹配(假设数据是按块读取的)，失败后会进行全局搜索。
%
% 输入：adata (包含 .EW, .NS, .UD 数据)
%
% 输出:
%   paired_stations - 结构体，字段名为台站代码。包含 .ew, .ns, .ud 三个属性(索引)。
%   unpaired_info   - 结构体数组，记录未能凑齐3个分量的台站信息。

    % 初始化
    nn = length(adata);
    if nn == 0 
        warning('输入数据为空。');
    end
    
    % 检查是否能被3整除（只是个警告，不影响后续处理）
    if mod(nn, 3) ~= 0
        warning('输入数据总数 (%d) 不是3的倍数，可能存在缺失分量的台站。', nn);
    end
    
    % 假设数据是按 [所有EW; 所有NS; 所有UD] 排列的，计算步长
    n_third = floor(nn / 3);
    
    paired_stations = struct();
    unpaired_list = []; 
    
    % 创建一个“已处理”标记数组
    processed_flags = false(nn, 1);

    % --- 1. 快速匹配 (i, i+n/3, i+2n/3) ---
    fprintf('--- 开始快速匹配 (三分量) ---\n');
    % 只有当数据量足够且大概率是对齐时才尝试
    if n_third > 0
        for i = 1:n_third
            % 候选索引
            idx1 = i;
            idx2 = i + n_third;
            idx3 = i + 2 * n_third;
            
            % 越界保护
            if idx3 > nn, continue; end
            
            code1 = adata{idx1}.station_code;
            code2 = adata{idx2}.station_code;
            code3 = adata{idx3}.station_code;
            
            % 检查三个位置的台站代码是否一致
            if strcmp(code1, code2) && strcmp(code1, code3)
                
                % 匹配成功，开始分配方向
                indices = [idx1, idx2, idx3];
                dirs = {adata{idx1}.direction, adata{idx2}.direction, adata{idx3}.direction};
                
                % 临时结构体用于存放当前组
                temp_s = struct();
                
                for k = 1:3
                    d_str = upper(dirs{k});
                    if contains(d_str, 'E') % E-W
                        temp_s.ew = indices(k);
                    elseif contains(d_str, 'N') % N-S
                        temp_s.ns = indices(k);
                    elseif contains(d_str, 'U') % U-D
                        temp_s.ud = indices(k);
                    end
                end
                
                % 只有当三个方向都齐了才算成功
                if isfield(temp_s, 'ew') && isfield(temp_s, 'ns') && isfield(temp_s, 'ud')
                    station_field_name = matlab.lang.makeValidName(code1);
                    paired_stations.(station_field_name) = temp_s;
                    
                    % 标记已处理
                    processed_flags([idx1, idx2, idx3]) = true;
                    
                    % 仅打印部分日志以免刷屏
                    % fprintf('快速匹配成功: %s\n', code1);
                end
            end
        end
    end
    
    % --- 2. 对剩余数据进行全局搜索 (分组匹配) ---
    unprocessed_indices = find(~processed_flags);
    if ~isempty(unprocessed_indices)
        fprintf('--- 快速匹配完成，对剩余 %d 个记录进行全局搜索 ---\n', length(unprocessed_indices));
        
        % 为了效率，先把剩余数据的 idx 分组归类到 Map 中
        % Key: StationCode, Value: [idx1, idx2, ...]
        station_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        for k = 1:length(unprocessed_indices)
            curr_idx = unprocessed_indices(k);
            code = adata{curr_idx}.station_code;
            
            if isKey(station_map, code)
                station_map(code) = [station_map(code), curr_idx];
            else
                station_map(code) = curr_idx;
            end
        end
        
        % 遍历 Map 中的每个台站，检查分量是否齐全
        keys = station_map.keys;
        for k = 1:length(keys)
            s_code = keys{k};
            idxs = station_map(s_code);
            
            temp_s = struct();
            
            % 遍历该台站下找到的所有索引
            for m = 1:length(idxs)
                idx = idxs(m);
                d_str = upper(adata{idx}.direction);
                
                if contains(d_str, 'E')
                    temp_s.ew = idx;
                elseif contains(d_str, 'N')
                    temp_s.ns = idx;
                elseif contains(d_str, 'U')
                    temp_s.ud = idx;
                end
            end
            
            % 检查是否凑齐三兄弟
            if isfield(temp_s, 'ew') && isfield(temp_s, 'ns') && isfield(temp_s, 'ud')
                % 配对成功
                station_field_name = matlab.lang.makeValidName(s_code);
                paired_stations.(station_field_name) = temp_s;
                
                processed_flags(idxs) = true; % 标记(虽然不需要了，但保持逻辑一致)
                fprintf('全局匹配成功: %s (EW:%d, NS:%d, UD:%d)\n', ...
                    s_code, temp_s.ew, temp_s.ns, temp_s.ud);
            else
                % 未凑齐 (可能是只有2个分量，或者有重复分量)
                % 将这些索引加入未配对列表
                for m = 1:length(idxs)
                    idx = idxs(m);
                    unpaired_list(end+1).station_code = s_code;
                    unpaired_list(end).index = idx;
                    unpaired_list(end).direction = adata{idx}.direction;
                end
            end
        end
    end
    
    % --- 3. 整理并输出最终结果 ---
    unpaired_info = unpaired_list;
    num_paired = length(fieldnames(paired_stations));
    num_unpaired = length(unpaired_info);
    
    fprintf('--- 配对完成 ---\n');
    fprintf('成功配对台站 (3分量齐全): %d 个\n', num_paired);
    fprintf('未能配对/分量缺失记录: %d 个\n', num_unpaired);
    
    if num_unpaired > 0
        fprintf('未能配对的台站列表 (前10个):\n');
        for k = 1:min(10, num_unpaired)
            fprintf('  - Code: %s, Index: %d, Dir: %s\n', ...
                unpaired_info(k).station_code, unpaired_info(k).index, unpaired_info(k).direction);
        end
        if num_unpaired > 10
            fprintf('  ... (剩余 %d 个略)\n', num_unpaired - 10);
        end
    end
end