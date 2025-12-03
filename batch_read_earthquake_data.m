function all_data = batch_read_earthquake_data(folder_path)
% 批量读取文件夹中的所有地震数据文件
% 输入: folder_path - 文件夹路径
% 输出: all_data - 包含所有文件数据的cell数组

    if nargin < 1
        folder_path = '.'; % 当前文件夹
    end
    % 保留 .EW.. 或 .EW2..  但排除 .EW1..
    all_ew = dir(fullfile(folder_path, '*.EW*'));
    ew_files = all_ew(~cellfun(@isempty, regexp({all_ew.name}, '\.EW([2]?)$')));
    all_ns = dir(fullfile(folder_path, '*.NS*'));
    ns_files = all_ns(~cellfun(@isempty, regexp({all_ns.name}, '\.NS([2]?)$')));

    all_ud = dir(fullfile(folder_path, '*.UD*'));
    ud_files = all_ud(~cellfun(@isempty, regexp({all_ud.name}, '\.UD([2]?)$')));


    % 查找所有 .EW 和 .NS 文件（knet 用，kiknet加入后看上边↑）
    %ew_files = dir(fullfile(folder_path, '*.EW*'));
    %ns_files = dir(fullfile(folder_path, '*.NS*'));
    %ud_files = dir(fullfile(folder_path, '*.UD*'));
    
    all_files = [ew_files; ns_files;ud_files];
    
    if isempty(all_files)
        warning('未找到任何 .EW 或 .NS 文件');
        all_data = {};
        return;
    end
    
    fprintf('找到 %d 个文件，开始读取...\n', length(all_files));
    
    all_data = cell(length(all_files), 1);
    
    for i = 1:length(all_files)
        filename = fullfile(folder_path, all_files(i).name);
        fprintf('读取文件 %d/%d: %s\n', i, length(all_files), all_files(i).name);
        
        try
            all_data{i} = read_earthquake_data(filename);
            fprintf('  成功读取 %d 个数据点\n', all_data{i}.num_samples);
        catch ME
            fprintf('  读取失败: %s\n', ME.message);
            all_data{i} = [];
        end
    end
    
    % 移除失败的读取
    all_data = all_data(~cellfun(@isempty, all_data));
    
    fprintf('成功读取 %d 个文件\n', length(all_data));
end