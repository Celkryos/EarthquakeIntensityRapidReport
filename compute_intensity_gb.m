function stations = compute_intensity_gb(stations, varargin)
% 按 GB/T 17742-2020 附录A 计算仪器地震烈度 I_gb
% 输入：
%   stations  - hcsc 之后的台站结构体，每个台站至少包含：
%               .is_valid      (check_truncated_records 生成)
%               .pga_T0_100    0.1s 三分量合成 PGA
%               .pgv_T0_100    0.1s 三分量合成 PGV
% 可选参数（Name-Value）：
%   'PGAField' - PGA 字段名，默认 'pga_T0_100'
%   'PGVField' - PGV 字段名，默认 'pgv_T0_100'
% 输出：
%   stations.(sta) 中新增字段：
%       .I_A_gb   - 按 A.5 计算的 I_A
%       .I_V_gb   - 按 A.6 计算的 I_V
%       .I_gb     - 仪器地震烈度 I_1（保留一位小数，截断到[1.0, 12.0]）

    p = inputParser;
    addParameter(p, 'PGAField', 'pga_T0_100', @(s)ischar(s) || isstring(s));
    addParameter(p, 'PGVField', 'pgv_T0_100', @(s)ischar(s) || isstring(s));
    parse(p, varargin{:});

    pga_field = char(p.Results.PGAField);
    pgv_field = char(p.Results.PGVField);

    names = fieldnames(stations);
    n_sta = numel(names);

    fprintf('--- 按国标计算仪器烈度 I_gb（使用 %s / %s）---\n', ...
        pga_field, pgv_field);
    for i = 1:n_sta
        sname = names{i};
        s = stations.(sname);
        % 1) 只对 is_valid 的台站计算
        if isfield(s, 'is_valid') && ~s.is_valid
            continue;
        end
        % 2) 必须有0.1s字段
        if ~isfield(s, pga_field) || ~isfield(s, pgv_field)
            fprintf('  台站 %s 缺少 %s 或 %s，跳过。\n', ...
                sname, pga_field, pgv_field);
            continue;
        end
        PGA = s.(pga_field);
        PGV = s.(pgv_field);

        
        
        % 3) 值必须是正值
        if ~isfinite(PGA) || ~isfinite(PGV) || PGA <= 0 || PGV <= 0
            fprintf('  台站 %s 的 PGA/PGV 非正或非法，设 I_gb = NaN。\n', sname);
            stations.(sname).I_A_gb = NaN;
            stations.(sname).I_V_gb = NaN;
            stations.(sname).I_gb   = NaN;
            continue;
        end
        % 附录A公式   记得单位转换 cms/cms2到ms/ms2
        I_A = 3.17 * log10(PGA*0.01) + 6.59;
        I_V = 3.00 * log10(PGV*0.01) + 9.77;
        if (I_A >= 6.0) && (I_V >= 6.0)
            I_1 = I_V;
        else
            I_1 = (I_A + I_V) / 2;
        end
        % 保留一位小数
        I_1 = round(I_1 * 10) / 10;

        % 限制在 [1.0, 12.0]
        if I_1 < 1.0
            I_1 = 1.0;
        elseif I_1 > 12.0
            I_1 = 12.0;
        end

        % 写回 stations
        stations.(sname).I_A_gb = I_A;
        stations.(sname).I_V_gb = I_V;
        stations.(sname).I_gb   = I_1;
    end

    fprintf('--- 仪器烈度计算完成 ---\n');
end
