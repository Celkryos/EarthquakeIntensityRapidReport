%ziduan='pgv_T0_100';kong=true;
%ziduan='pga_T0_500_h';kong=true;
ziduana='pgv_T0_020_h';  %【【【】】】
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
    if ~isfield(s, ziduana)
        continue;   % 跳
    end

    % 纳入列表
    lon_sta(end+1,1) = s.station_long;
    lat_sta(end+1,1) = s.station_lat;
    val_sta(end+1,1) = s.(ziduana); %基准值
end
jichu=val_sta;

%%
pinduana = {'pga_T0_100_h','pga_T0_500_h', 'pga_T1_000_h', 'pga_T2_000_h', 'pga_T5_000_h'};
pinduanv = {'pgv_T0_100_h','pgv_T0_500_h', 'pgv_T1_000_h', 'pgv_T2_000_h', 'pgv_T5_000_h'};

%R2=[];
%for p = 1:length(pinduana)
p=5;
    ziduana = pinduana{p}; 
    ziduanv=pinduanv{p};
    names = fieldnames(stations);
    lon_sta = [];
    lat_sta = [];
    val_staa = [];
    val_stav=[];

    for i = 1:numel(names)
        s = stations.(names{i});
        if ~isfield(s, 'is_valid') || ~s.is_valid
            continue;
        end
        if ~isfield(s, ziduana)
            continue; 
        end
        lon_sta(end+1,1) = s.station_long;
        lat_sta(end+1,1) = s.station_lat;
        val_staa(end+1,1) = s.(ziduana);
        val_stav(end+1,1) = s.(ziduanv);
    end
   

    mdl = fitlm(log10(val_stav), log10(jichu)); %【【【】】】
    R2 = mdl. Rsquared. Ordinary;

    
    figure;
    plotResiduals(mdl, 'fitted');  % 残差 vs 拟合值
    title(sprintf('R^2 = %.3f', mdl.Rsquared.Ordinary));
    figure;
    scatter(log10(jichu), log10(val_stav), 18, 'filled'); grid on;
    xlabel('jichu = log10 pgv\_T0\_020\_h');
    ylabel(ziduanv);  % 比如 pgv_T1_000_h
    title('val\_stav vs jichu log-log下');

    %figure;
%plotResiduals(mdl, 'probability');
%title('Residual Q-Q plot (log-log model)');

%plotDiagnostics(mdl,'cookd')

%end
%figure;
%plot(R2);
