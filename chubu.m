function jdata = chubu(data)
    for i = 1:length(data)
        data{i, 1}.acceleration_jjh=detrend(data{i, 1}.acceleration_gal,1); % jjh=基线矫正后
    end
    jdata=data;
end