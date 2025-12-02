%ziduan='pgv_T0_100';kong=true;
%ziduan='pga_T0_500_h';kong=true;
ziduan='I_gb';kong=false;
names = fieldnames(stations);
lon_sta = [];
lat_sta = [];
val_sta = [];

for i = 1:numel(names)
    s = stations.(names{i});
    % 只要 is_valid的台站
    if ~isfield(s, 'is_valid') || ~s.is_valid
        continue;   % 跳
    end
    % 需存在相应字段
    if ~isfield(s, ziduan)
        continue;   % 跳
    end

    % 纳入插值列表
    lon_sta(end+1,1) = s.station_long;
    lat_sta(end+1,1) = s.station_lat;
    val_sta(end+1,1) = s.(ziduan);
end
In=val_sta;
%%
pinduan = {'pga_T0_100','pga_T0_500_h', 'pga_T1_000_h', 'pga_T2_000_h', 'pga_T5_000_h'};
%pinduan = {'pgv_T0_100','pgv_T0_500_h', 'pgv_T1_000_h', 'pgv_T2_000_h', 'pgv_T5_000_h'};

R2=[];
for p = 1:length(pinduan)
    ziduan = pinduan{p}; 
    names = fieldnames(stations);
    lon_sta = [];
    lat_sta = [];
    val_sta = [];

    for i = 1:numel(names)
        s = stations.(names{i});
        if ~isfield(s, 'is_valid') || ~s.is_valid
            continue;
        end
        if ~isfield(s, ziduan)
            continue; 
        end
        lon_sta(end+1,1) = s.station_long;
        lat_sta(end+1,1) = s.station_lat;
        val_sta(end+1,1) = s.(ziduan);
    end
    pgv = val_sta;
    mdl = fitlm(log10(pgv), In);
    R2(p) = mdl. Rsquared. Ordinary;

    %figure;
    %plotResiduals(mdl, 'fitted');  % 残差 vs 拟合值
    %title(sprintf('R^2 = %.3f', mdl.Rsquared.Ordinary));
end

