function data = acc2vel(data)

    dt = 1 / data.sampling_freq_hz;
    
    % 确保使用的是未裁剪的加速度
    if ~isfield(data, 'acc_untrimmed')
        error('错误: acc2vel 需要由 nb_filt 处理后的 data.acc_untrimmed 字段。');
    end

    % 1. 积分得到未裁剪的速度
    % 滤波后的加速度均值理论上已接近0，detrend(..., 0) 可以省略
    vel_untrimmed = cumtrapz(data.acc_untrimmed) * dt;
    
    % 2. 积分得到未裁剪的位移
    % 这里需要对速度进行detrend(...,0)或detrend(...,1)，以移除累积误差
    disp_untrimmed = cumtrapz(detrend(vel_untrimmed, 1)) * dt;
    
    % 3. 保存未裁剪的速度和位移结果
    data.vel_untrimmed = vel_untrimmed;
    data.disp_untrimmed = disp_untrimmed;
    
    
end