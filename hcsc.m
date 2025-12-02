function paired_stations = hcsc(adata, paired_stations, T_list)
% HCSC - 按周期合成水平/三分量时程，并根据 trim_info 裁剪到原始长度。
% 前提：
%   - adata{idx} 中已存在各周期的滤波 + 积分结果：
%       acc_T5_000s, vel_T5_000s, acc_T0_100s, vel_T0_100s, ...
%   - adata{idx}.trim_info 中有：
%       .n_prepended, .original_length
%   - paired_stations.(sta) 里有：
%       .ew, .ns, .ud (索引)，以及 .is_valid (由 check_truncated_records 生成)
%
% 输出：
%   在 paired_stations 下新增字段（举例 T0=5, 0.1）：
%       .sc_T5_000_acc_h   水平合成加速度
%       .pga_T5_000_h      水平 PGA
%       .sc_T5_000_v_h     水平合成速度
%       .pgv_T5_000_h      水平 PGV
%       .sc_T0_100_acc     三分量合成加速度
%       .pga_T0_100        三分量 PGA
%       .sc_T0_100_v       三分量合成速度
%       .pgv_T0_100        三分量 PGV
%   均为根据 trim_info 裁剪为原始长度后的时程。

    sta_names = fieldnames(paired_stations);
    fprintf('--- 开始按周期合成时程 (共 %d 个台站) ---\n', numel(sta_names));

    for iSta = 1:numel(sta_names)
        sta_name = sta_names{iSta};
        sta_info = paired_stations.(sta_name);

        % 跳过无效台站
        if isfield(sta_info, 'is_valid') && ~sta_info.is_valid
            fprintf('跳过无效台站 %s\n', sta_name);
            continue;
        end

        % 必须至少有 EW/NS
        if ~isfield(sta_info, 'ew') || ~isfield(sta_info, 'ns')
            fprintf('台站 %s 缺少 EW/NS 索引，跳过。\n', sta_name);
            continue;
        end

        idx_ew = sta_info.ew;
        idx_ns = sta_info.ns;
        idx_ud = sta_info.ud;

        if isempty(adata{idx_ew}) || isempty(adata{idx_ns})
            fprintf('台站 %s 的 EW/NS 记录为空，跳过。\n', sta_name);
            continue;
        end

        % ---------- 获取 trim_info ----------
        if ~isfield(adata{idx_ew}, 'trim_info')
            warning('台站 %s 缺少 trim_info，跳过合成。', sta_name);
            continue;
        end
        trim = adata{idx_ew}.trim_info;

        n_pre   = trim.n_prepended;
        N_orig  = trim.original_length;
        start_i = n_pre + 1;
        end_i   = n_pre + N_orig;

        % 台站存一个 time，沿用 adata 中分量
        paired_stations.(sta_name).time = adata{idx_ew}.time(:);

        fprintf('台站 %s：合成与裁剪 [%d, %d] 样点。\n', ...
            sta_name, start_i, end_i);

        % ---------- 逐周期处理 ----------
        for kT = 1:numel(T_list)
            T0 = T_list(kT);

            % 构造加速度 / 速度字段名，例如 acc_T5_000s
            T_str_full  = sprintf('T%.3fs', T0);            % "T5.000s"
            acc_field   = matlab.lang.makeValidName(['acc_' T_str_full]);
            vel_field   = regexprep(acc_field, '^acc', 'vel', 1);

            % 生成台站级别的 sc 字段名前缀，例如 "sc_T5_000"
            %   先去掉末尾的 's'，再把小数点变成下划线
            T_token     = regexprep(T_str_full, 's$', '');  % "T5.000"
            T_short     = strrep(T_token, '.', '_');        % "T5_000"
            sc_prefix   = ['sc_' T_short];

            is_T01 = abs(T0 - 0.1) < 1e-6;

            % 检查字段是否存在
            if ~isfield(adata{idx_ew}, acc_field) || ~isfield(adata{idx_ns}, acc_field) ...
                    || ~isfield(adata{idx_ew}, vel_field) || ~isfield(adata{idx_ns}, vel_field)
                fprintf('  台站 %s: 周期 T=%.3fs 缺少 EW/NS 字段 %s / %s，跳过。\n', ...
                    sta_name, T0, acc_field, vel_field);
                continue;
            end

            % ---------- 读取 & 裁剪 EW/NS ----------
            acc_ew_full = adata{idx_ew}.(acc_field)(:);
            acc_ns_full = adata{idx_ns}.(acc_field)(:);
            vel_ew_full = adata{idx_ew}.(vel_field)(:);
            vel_ns_full = adata{idx_ns}.(vel_field)(:);

            if end_i <= start_i
                warning('台站 %s: T=%.3fs 有效长度 <= 0，跳过。', sta_name, T0);
                continue;
            end

            acc_ew = acc_ew_full(start_i:end_i);
            acc_ns = acc_ns_full(start_i:end_i);
            vel_ew = vel_ew_full(start_i:end_i);
            vel_ns = vel_ns_full(start_i:end_i);

            % ---------- 水平合成 ----------
            acc_h = sqrt(acc_ew.^2 + acc_ns.^2);
            vel_h = sqrt(vel_ew.^2 + vel_ns.^2);

            if is_T01
                % ===== T0 = 0.1 s: 三分量合成 =====
                if isempty(idx_ud) || isempty(adata{idx_ud}) ...
                        || ~isfield(adata{idx_ud}, acc_field) ...
                        || ~isfield(adata{idx_ud}, vel_field)
                    fprintf('  台站 %s: T=0.1 s 缺少 UD 分量，按只有水平处理。\n', sta_name);
                    % 退化成水平，只存 _acc_h/_v_h 以及对应的水平 PGA/PGV
                    field_acc_h = matlab.lang.makeValidName([sc_prefix '_acc_h']);
                    field_v_h   = matlab.lang.makeValidName([sc_prefix '_v_h']);
                    paired_stations.(sta_name).(field_acc_h) = acc_h;
                    paired_stations.(sta_name).(field_v_h)   = vel_h;

                    % 水平 PGA/PGV（退化情况）
                    pga_field_h = matlab.lang.makeValidName(['pga_' T_short '_h']);
                    pgv_field_h = matlab.lang.makeValidName(['pgv_' T_short '_h']);
                    paired_stations.(sta_name).(pga_field_h) = max(abs(acc_h));
                    paired_stations.(sta_name).(pgv_field_h) = max(abs(vel_h));
                else
                    % 读取并裁剪 UD
                    acc_ud_full = adata{idx_ud}.(acc_field)(:);
                    vel_ud_full = adata{idx_ud}.(vel_field)(:);

                    acc_ud = acc_ud_full(start_i:end_i);
                    vel_ud = vel_ud_full(start_i:end_i);

                    % 三分量合成
                    acc_3c = sqrt(acc_ew.^2 + acc_ns.^2 + acc_ud.^2);
                    vel_3c = sqrt(vel_ew.^2 + vel_ns.^2 + vel_ud.^2);

                    field_acc = matlab.lang.makeValidName([sc_prefix '_acc']);
                    field_v   = matlab.lang.makeValidName([sc_prefix '_v']);

                    paired_stations.(sta_name).(field_acc) = acc_3c;
                    paired_stations.(sta_name).(field_v)   = vel_3c;

                    % 三分量 PGA/PGV，命名无 _h
                    pga_field = matlab.lang.makeValidName(['pga_' T_short]);
                    pgv_field = matlab.lang.makeValidName(['pgv_' T_short]);
                    paired_stations.(sta_name).(pga_field) = max(abs(acc_3c));
                    paired_stations.(sta_name).(pgv_field) = max(abs(vel_3c));
                end
            else
                % ===== T0 ≠ 0.1 s: 只存水平 =====
                field_acc_h = matlab.lang.makeValidName([sc_prefix '_acc_h']);
                field_v_h   = matlab.lang.makeValidName([sc_prefix '_v_h']);

                paired_stations.(sta_name).(field_acc_h) = acc_h;
                paired_stations.(sta_name).(field_v_h)   = vel_h;

                % 水平 PGA/PGV
                pga_field_h = matlab.lang.makeValidName(['pga_' T_short '_h']);
                pgv_field_h = matlab.lang.makeValidName(['pgv_' T_short '_h']);
                paired_stations.(sta_name).(pga_field_h) = max(abs(acc_h));
                paired_stations.(sta_name).(pgv_field_h) = max(abs(vel_h));
            end
        end % for kT
    end % for iSta

    fprintf('--- 合成 & 裁剪 + PGA/PGV 提取完成 ---\n');
end
