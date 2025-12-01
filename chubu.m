function jdata = chubu(data)
    for i = 1:length(data)
        if isempty(data{i,1}), continue; end
        if isfield(data{i,1}, 'is_valid') && ~data{i,1}.is_valid
            continue;
        end
        data{i, 1}.acceleration_jjh=detrend(data{i, 1}.acceleration_gal,1); % jjh=基线矫正后
    end
    jdata=data;
end