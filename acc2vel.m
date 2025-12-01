function data = acc2vel(data, acc_field_names)
% ACC2VEL - 对一个或多个加速度时程字段积分，生成对应的速度和位移字段。
%
% --- 采样间隔 ---
    if ~isfield(data, 'sampling_freq_hz')
        error('acc2vel:MissingField', 'data.sampling_freq_hz 缺失。');
    end
    dt = 1 / data.sampling_freq_hz;

    % --- 解析 acc_field_names 参数 ---
    if nargin < 2 || isempty(acc_field_names)
    % 默认使用 acc_untrimmed
        acc_field_names = {'acc_untrimmed'};
    elseif ischar(acc_field_names)
        acc_field_names = {acc_field_names};
    elseif ~iscell(acc_field_names)
        error('acc2vel:InvalidInput', 'acc_field_names 必须是字符串或字符串 cell。');
    end

    % --- 对每一个加速度字段进行积分 ---
    for k = 1:numel(acc_field_names)
        acc_field = acc_field_names{k};

        if ~isfield(data, acc_field)
            warning('acc2vel:MissingAccField', ...
                '字段 %s 在 data 中不存在，跳过。', acc_field);
            continue;
        end

        acc = data.(acc_field)(:);  % 列向量

        % 1. 积分得到未裁剪的速度
        vel_untrimmed = cumtrapz(acc) * dt;

        % 2. 积分得到未裁剪的位移
        % 对速度做一次 detrend(,1) 去掉积累漂移
        disp_untrimmed = cumtrapz(detrend(vel_untrimmed, 1)) * dt;

        % 3. 生成速度/位移字段名
        %    把前缀 acc 换成 vel/disp
        vel_field  = regexprep(acc_field, '^acc', 'vel', 1);
        disp_field = regexprep(acc_field, '^acc', 'disp', 1);

       
        % 4. 写回 data
        data.(vel_field)  = vel_untrimmed;
        data.(disp_field) = disp_untrimmed;
    end

end
