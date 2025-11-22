function d = epi_station_distance(eq_lat, eq_lon, sta_lat, sta_lon, unit)
% 计算震源与台站之间的球面大圆距离
% 输入（标量或同尺寸矩阵）：
%   eq_lat, eq_lon   - 震源纬度、经度（单位：度）
%   sta_lat, sta_lon - 台站纬度、经度（单位：度）
%   unit             - 输出单位: 'km' 或 'm'，默认 'km'
% 输出：
%   d                - 震源与台站之间的大圆距离（单位：km 或 m）

    if nargin < 5 || isempty(unit)
        unit = 'km';
    end

    % 地球平均半径
    R_km = 6371.0;   % km

    switch lower(unit)
        case 'km'
            R = R_km;
        case 'm'
            R = R_km * 1000;
        otherwise
            error('unit 只能是 ''km'' 或 ''m''。');
    end

    % 角度 -> 弧度
    phi1 = deg2rad(eq_lat);
    phi2 = deg2rad(sta_lat);
    dphi = deg2rad(sta_lat - eq_lat);
    dlambda = deg2rad(sta_lon - eq_lon);

    % Haversine 公式
    a = sin(dphi/2).^2 + cos(phi1) .* cos(phi2) .* sin(dlambda/2).^2;
    c = 2 .* atan2(sqrt(a), sqrt(1 - a));

    % 大圆距离
    d = R .* c;
end
