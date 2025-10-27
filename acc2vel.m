function data = acc2vel(data)
%%%%%%%%%%%%%%%%%%待修改%%%%%%%%%%%%%%%%%%%%%
    % 积 得 v 把加速度 (gal) 积分成速度 (cm/s)
    dt = 1 / data.sampling_freq_hz;
    % gal = cm/s^2 直接积分得到 cm/s
    data.velocity_cms = cumtrapz(detrend(data.nb_acc_0_02,1)) * dt;
    data.d_cm=cumtrapz(detrend(data.velocity_cms,1)) * dt;
    data.d_cm=detrend(data.d_cm,1);
end
