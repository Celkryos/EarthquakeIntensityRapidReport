function T_all = append_metric_table(T_all, T_new)
% APPEND_METRIC_TABLE
% 安全纵向合并两个 metric table。
% 如果变量名不完全一致，会自动补齐缺失变量。

    if isempty(T_new) || height(T_new) == 0
        return;
    end

    if isempty(T_all) || height(T_all) == 0
        T_all = T_new;
        return;
    end

    vars_all = T_all.Properties.VariableNames;
    vars_new = T_new.Properties.VariableNames;

    all_vars = unique([vars_all, vars_new], 'stable');

    for i = 1:numel(all_vars)
        v = all_vars{i};

        if ~ismember(v, vars_all)
            T_all.(v) = make_missing_column(T_all, T_new, v);
        end

        if ~ismember(v, vars_new)
            T_new.(v) = make_missing_column(T_new, T_all, v);
        end
    end

    T_all = T_all(:, all_vars);
    T_new = T_new(:, all_vars);

    T_all = [T_all; T_new];
end


function col = make_missing_column(T_target, T_ref, var_name)
% 根据参考表中的变量类型，为目标表生成缺失列。

    n = height(T_target);

    if ismember(var_name, T_ref.Properties.VariableNames)
        ref_col = T_ref.(var_name);

        if isstring(ref_col)
            col = strings(n, 1);
        elseif iscellstr(ref_col)
            col = repmat({''}, n, 1);
        elseif isnumeric(ref_col)
            col = NaN(n, 1);
        elseif islogical(ref_col)
            col = false(n, 1);
        else
            col = repmat(missing, n, 1);
        end
    else
        col = NaN(n, 1);
    end
end