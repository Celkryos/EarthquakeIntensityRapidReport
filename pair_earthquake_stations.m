function [paired_stations, unpaired_info] = pair_earthquake_stations(adata)
% PAIR_EARTHQUAKE_STATIONS - 将地震动数据按台站代码配对（至少 EW+NS），UD 可选
%
% 输入:
%   adata - 由 batch_read_earthquake_data 生成的 cell 数组
%
% 输出:
%   paired_stations - 结构体，字段名为台站代码。包含：
%       .ew, .ns               EW/NS 分量在 adata 中的索引（必需）
%       .ud                    UD 分量索引（可选，若存在则给出）
%       .station_lat/.station_long  台站经纬度（优先取 EW，否则 NS）
%   unpaired_info - 结构体数组，记录未能凑齐 EW+NS 的台站信息
%       .station_code, .has_ew, .has_ns, .has_ud, .indices

    nn = numel(adata);
    paired_stations = struct();
    unpaired_info = struct('station_code', {}, 'has_ew', {}, 'has_ns', {}, 'has_ud', {}, 'indices', {});

    if nn == 0
        warning('输入数据为空。');
        return;
    end

    % --- 建立 station_code -> indices 的映射 ---
    station_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:nn
        if isempty(adata{i}) || ~isfield(adata{i}, 'station_code')
            continue;
        end
        code = char(adata{i}.station_code);
        if isKey(station_map, code)
            station_map(code) = [station_map(code), i];
        else
            station_map(code) = i;
        end
    end

    keys = station_map.keys;
    fprintf('--- 开始台站分量配对（要求 EW+NS，UD 可选）---\n');

    for k = 1:numel(keys)
        s_code = keys{k};
        idxs = station_map(s_code);

        temp_s = struct();
        for m = 1:numel(idxs)
            idx = idxs(m);
            if isempty(adata{idx}) || ~isfield(adata{idx}, 'direction')
                continue;
            end
            d_str = upper(char(adata{idx}.direction));
            if contains(d_str, 'E')
                temp_s.ew = idx;
            elseif contains(d_str, 'N')
                temp_s.ns = idx;
            elseif contains(d_str, 'U')
                temp_s.ud = idx;
            end
        end

        has_ew = isfield(temp_s, 'ew');
        has_ns = isfield(temp_s, 'ns');
        has_ud = isfield(temp_s, 'ud');

        if has_ew && has_ns
            station_field_name = matlab.lang.makeValidName(s_code);

            % 台站经纬度：优先取 EW，否则取 NS
            if has_ew && isfield(adata{temp_s.ew}, 'station_lat')
                temp_s.station_lat  = adata{temp_s.ew}.station_lat;
                temp_s.station_long = adata{temp_s.ew}.station_long;
            else
                temp_s.station_lat  = adata{temp_s.ns}.station_lat;
                temp_s.station_long = adata{temp_s.ns}.station_long;
            end

            paired_stations.(station_field_name) = temp_s;
            if has_ud
                fprintf('配对成功: %s (EW:%d, NS:%d, UD:%d)\n', s_code, temp_s.ew, temp_s.ns, temp_s.ud);
            else
                fprintf('配对成功: %s (EW:%d, NS:%d, UD:缺失)\n', s_code, temp_s.ew, temp_s.ns);
            end
        else
            % 未能凑齐 EW+NS → 记录
            rec.station_code = s_code;
            rec.has_ew = has_ew;
            rec.has_ns = has_ns;
            rec.has_ud = has_ud;
            rec.indices = idxs;
            unpaired_info(end+1) = rec; %#ok<AGROW>
        end
    end

    fprintf('--- 配对完成 ---\n');
    fprintf('成功配对台站 (EW+NS): %d 个\n', numel(fieldnames(paired_stations)));
    fprintf('未能配对 (EW/NS 缺失): %d 个\n', numel(unpaired_info));

    if ~isempty(unpaired_info)
        fprintf('未能配对的台站列表 (前10个):\n');
        for k = 1:min(10, numel(unpaired_info))
            r = unpaired_info(k);
            fprintf('  - Code: %s, hasEW:%d, hasNS:%d, hasUD:%d, nIdx:%d\n', ...
                r.station_code, r.has_ew, r.has_ns, r.has_ud, numel(r.indices));
        end
        if numel(unpaired_info) > 10
            fprintf('  ... (剩余 %d 个略)\n', numel(unpaired_info) - 10);
        end
    end
end
