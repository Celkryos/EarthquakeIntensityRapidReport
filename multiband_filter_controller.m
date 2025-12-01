function adata = multiband_filter_controller(adata, paired_stations, T_list, fc_high)
% MULTIBAND_FILTER_CONTROLLER
% 做多频段 nb_filt 滤波
%   - 对 T0 ~= 0.1 s：只处理水平分量 (EW/NS)
%   - 对 T0 == 0.1 s：处理三分量 (EW/NS/UD)
% 保存“带 pad 的 acc_Txxx”，暂不裁剪。

    if isempty(T_list)
        warning('T_list 为空，multiband_filter_controller 未执行任何滤波。');
        return;
    end

    station_names = fieldnames(paired_stations);
    fprintf('--- 多频段滤波开始，共 %d 个台站，%d 个周期 ---\n', ...
        numel(station_names), numel(T_list));

    for iSta = 1:numel(station_names)
        sta_name = station_names{iSta};
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
            fprintf('台站 %s 在 adata 中记录为空，跳过。\n', sta_name);
            continue;
        end

        % 构造要处理的分量列表
        comp_indices = [idx_ew, idx_ns,idx_ud];
        comp_labels  = {'ew', 'ns','ud'};        
        fprintf('处理台站 %s\n', sta_name);

        for ic = 1:numel(comp_indices)
            comp_idx = comp_indices(ic);
            comp_lab = comp_labels{ic};

            data_orig = adata{comp_idx};
            if isempty(data_orig)
                continue;
            end

            for iT = 1:numel(T_list)
                T0 = T_list(iT);
                % T0 == 0.1 s：三分量都处理
                % T0 ~= 0.1 s：只处理 EW/NS，跳过 UD
                is_T01 = abs(T0 - 0.1) < 1e-6;
                if ~is_T01 && strcmp(comp_lab, 'ud')
                    % 非 0.1 s 且是 UD 分量 → 跳过
                    continue;
                end

                % 调用纯函数 nb_filt
                [acc_untrim, trim_info] = nb_filt(data_orig, fc_high, T0);

                % 生成字段名，例如 acc_T5_000s, acc_T0_100s
                T_str   = sprintf('T%.3fs', T0);                   % "T5.000s"
                fieldnm = matlab.lang.makeValidName(['acc_' T_str]);
                adata{comp_idx}.(fieldnm) = acc_untrim;

                % trim_info 只存一份（各 T0 共用同一 padding）
                if ~isfield(adata{comp_idx}, 'trim_info')
                    adata{comp_idx}.trim_info = trim_info;
                end

                if ~isfield(adata{comp_idx}, 'filtered_T_list')
                    adata{comp_idx}.filtered_T_list = T_list(:).';
                end

                fprintf('   [%s] %s -> %s 完成，T0=%.3fs, fc_high=%.3f Hz\n', ...
                    sta_name, comp_lab, fieldnm, T0, fc_high);
            end
        end
    end

    fprintf('--- 多频段滤波完成（保留 pad 段） ---\n');
end
