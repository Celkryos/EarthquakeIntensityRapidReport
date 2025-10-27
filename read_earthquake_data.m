function earthquake_data = read_earthquake_data(filename)
% 读取地震数据文件 (.EW 或 .NS 格式) 
% 输入: filename - 文件名 (如 'MYZ0132408081643.EW')
% 输出: earthquake_data - 包含所有信息的结构体

    % 检查文件是否存在
    if ~exist(filename, 'file')
        error('文件不存在: %s', filename);
    end
    
    % 打开文件
    fid = fopen(filename, 'r', 'n', 'UTF-8');  % 指定编码
    if fid == -1
        error('无法打开文件: %s', filename);
    end
    
    try
        % 初始化结构体
        earthquake_data = struct();
        earthquake_data.filename = filename;
        
        % 读取所有行
        all_lines = {};
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line)
                all_lines{end+1} = line;
            end
        end
        fclose(fid);
        
        % 找到 Memo. 行的位置
        memo_line_idx = 0;
        for i = 1:length(all_lines)
            if contains(all_lines{i}, 'Memo.')
                memo_line_idx = i;
                break;
            end
        end
        
        if memo_line_idx == 0
            error('未找到 Memo. 行，文件格式可能有误');
        end
        
        % 解析头部信息（Memo.行之前的所有行）
        for i = 1:memo_line_idx
            line = all_lines{i};
            earthquake_data = parse_header_line(line, earthquake_data);
        end
        
        % 读取加速度数据（Memo.行之后的所有行）
        data_lines = all_lines((memo_line_idx+1):end);
        acceleration_data = parse_acceleration_data(data_lines);
        
        earthquake_data.acceleration_raw = acceleration_data;
        earthquake_data.num_samples = length(acceleration_data);
        
        % 转换为实际加速度值 (gal)
        if isfield(earthquake_data, 'scale_factor') && ~isnan(earthquake_data.scale_factor)
            earthquake_data.acceleration_gal = acceleration_data * earthquake_data.scale_factor;
        else
            earthquake_data.acceleration_gal = acceleration_data;
            warning('未找到有效的scale_factor，使用原始数据');
        end
        
        % 生成时间序列
        if isfield(earthquake_data, 'sampling_freq_hz') && ~isnan(earthquake_data.sampling_freq_hz)
            dt = 1 / earthquake_data.sampling_freq_hz;
            earthquake_data.time = (0:dt:(length(acceleration_data)-1)*dt)';
        else
            earthquake_data.time = (1:length(acceleration_data))';
            warning('未找到采样频率，使用样本序号作为时间轴');
        end
        
        % 计算统计信息
        if isfield(earthquake_data, 'acceleration_gal')
            earthquake_data.computed_max_acc = max(abs(earthquake_data.acceleration_gal));
            earthquake_data.computed_min_acc = min(earthquake_data.acceleration_gal);
            earthquake_data.computed_mean_acc = mean(earthquake_data.acceleration_gal);
            earthquake_data.computed_std_acc = std(earthquake_data.acceleration_gal);
        end
        
    catch ME
        if fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end
end

function earthquake_data = parse_header_line(line, earthquake_data)
% 解析单行头部信息
    
    % 去除首尾空格
    line = strtrim(line);
    
    if isempty(line)
        return;
    end
    
    try
        if contains(line, 'Origin Time')
            % Origin Time       发震时刻
            tokens = regexp(line, 'Origin Time\s+(.+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.origin_time = strtrim(tokens{1}{1});
            end
            
        elseif contains(line, 'Lat.') && ~contains(line, 'Station')
            % Lat.              震源纬度
            tokens = regexp(line, 'Lat\.\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.latitude = str2double(tokens{1}{1});
            else
                earthquake_data.latitude = NaN;
            end
            
        elseif contains(line, 'Long.') && ~contains(line, 'Station')
            % Long.             震源经度
            tokens = regexp(line, 'Long\.\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.longitude = str2double(tokens{1}{1});
            else
                earthquake_data.longitude = NaN;
            end
            
        elseif contains(line, 'Depth. (km)')
            % Depth. (km)       震源深度
            tokens = regexp(line, 'Depth\.\s*\(km\)\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.depth_km = str2double(tokens{1}{1});
            else
                earthquake_data.depth_km = NaN;
            end
            
        elseif contains(line, 'Mag.') && ~contains(line, 'Max')
            % Mag.              震级
            tokens = regexp(line, 'Mag\.\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.magnitude = str2double(tokens{1}{1});
            else
                earthquake_data.magnitude = NaN;
            end
            
        elseif contains(line, 'Station Code')
            % Station Code      台站代码
            tokens = regexp(line, 'Station Code\s+(\S+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.station_code = tokens{1}{1};
            end
            
        elseif contains(line, 'Station Lat.')
            % Station Lat.      台站纬度
            tokens = regexp(line, 'Station Lat\.\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.station_lat = str2double(tokens{1}{1});
            else
                earthquake_data.station_lat = NaN;
            end
            
        elseif contains(line, 'Station Long.')
            % Station Long.     台站经度
            tokens = regexp(line, 'Station Long\.\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.station_long = str2double(tokens{1}{1});
            else
                earthquake_data.station_long = NaN;
            end
            
        elseif contains(line, 'Station Height(m)')
            % Station Height(m)  台站海拔 
            tokens = regexp(line, 'Station Height\(m\)\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.station_height_m = str2double(tokens{1}{1});
            else
                earthquake_data.station_height_m = NaN;
            end
            
        elseif contains(line, 'Record Time')
            % Record Time       本记录起始时间
            % 注意！！↓ 来自https://www.kyoshin.bosai.go.jp/kyoshin/man/knetform_en.html
            % (Line 10th) Recording start time
            % This time includes the effect of 15-second trigger delay in the data logger. 
            % The true start time of the time series is obtained by subtracting 15 seconds from this "Recording start time."
            tokens = regexp(line, 'Record Time\s+(.+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.record_time = strtrim(tokens{1}{1});
            end
            
        elseif contains(line, 'Sampling Freq(Hz)')
            % Sampling Freq(Hz)  采样率 通常为100Hz
            tokens = regexp(line, 'Sampling Freq\(Hz\)\s+([0-9.]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.sampling_freq_hz = str2double(tokens{1}{1});
            else
                earthquake_data.sampling_freq_hz = NaN;
            end
            
        elseif contains(line, 'Duration Time(s)')
            % Duration Time(s)  记录持续时长
            tokens = regexp(line, 'Duration Time\(s\)\s+([0-9.]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.duration_time_s = str2double(tokens{1}{1});
            else
                earthquake_data.duration_time_s = NaN;
            end
            
        elseif contains(line, 'Dir.') && ~contains(line, 'Direction')
            % Dir.              东西or南北
            tokens = regexp(line, 'Dir\.\s+(\S+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.direction = tokens{1}{1};
            end
            
        elseif contains(line, 'Scale Factor')
            % Scale Factor      aaaa(gal)/bbbbbbb
            tokens = regexp(line, 'Scale Factor\s+(.+)', 'tokens');
            if ~isempty(tokens)
                scale_str = strtrim(tokens{1}{1});
                earthquake_data.scale_factor_raw = scale_str;
                
                % 解析数值
                scale_tokens = regexp(scale_str, '(\d+)\(gal\)/(\d+)', 'tokens');
                if ~isempty(scale_tokens)
                    earthquake_data.scale_factor_gal = str2double(scale_tokens{1}{1});
                    earthquake_data.scale_factor_denominator = str2double(scale_tokens{1}{2});
                    earthquake_data.scale_factor = earthquake_data.scale_factor_gal / earthquake_data.scale_factor_denominator;
                else
                    earthquake_data.scale_factor = NaN;
                end
            end
            
        elseif contains(line, 'Max. Acc. (gal)')
            % Max. Acc. (gal)   最大acc
            tokens = regexp(line, 'Max\. Acc\. \(gal\)\s+([0-9.-]+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.max_acc_gal = str2double(tokens{1}{1});
            else
                earthquake_data.max_acc_gal = NaN;
            end
            
        elseif contains(line, 'Last Correction')
            % Last Correction   最近一次仪器时间校准的时间点 通常可忽略
            tokens = regexp(line, 'Last Correction\s+(.+)', 'tokens');
            if ~isempty(tokens)
                earthquake_data.last_correction = strtrim(tokens{1}{1});
            end
            
        elseif contains(line, 'Memo.')
            earthquake_data.memo = '';
        end
        
    catch
        % 如果某行解析失败，跳过该行
        warning('解析失败的行: %s', line);
    end
end

function acceleration_data = parse_acceleration_data(data_lines)
% 解析加速度数据行
    
    % 合并所有数据行
    data_text = strjoin(data_lines, ' ');
    
    % 使用正则表达式提取所有数字（包括负数）
    numbers_str = regexp(data_text, '-?\d+', 'match');
    
    % 转换为数值
    acceleration_data = str2double(numbers_str);
    
    % 移除NaN值
    acceleration_data = acceleration_data(~isnan(acceleration_data));
    
    % 转换为列向量
    acceleration_data = acceleration_data(:);
end