function [adata, event_info] = read_earthquake_event(event_dirs, m, n)
% READ_EARTHQUAKE_EVENT
% 按组 m、组内事件 n 读取一个地震事件。
% kik 和 knet 会依次读取，并合并到同一个 adata 里。
%
% 用法：
%   [adata, event_info] = read_earthquake_event(event_dirs, 1, 3);

    if isempty(event_dirs)
        error('event_dirs 为空，请先运行 build_earthquake_read_dirs。');
    end

    group_list = [event_dirs.m];
    iG = find(group_list == m, 1);

    if isempty(iG)
        error('没有找到第 %d 组。', m);
    end

    events = event_dirs(iG).events;

    if n < 1 || n > numel(events)
        error('第 %d 组中没有第 %d 个事件。该组共有 %d 个事件。', ...
              m, n, numel(events));
    end

    event_info = events(n);
    adata = {};

    fprintf('--- 读取事件：组 m=%d，事件 n=%d，%s ---\n', ...
            m, n, event_info.event_folder);

    for k = 1:numel(event_info.networks)
        net = event_info.networks{k};
        folder_path = event_info.data_dirs{k};

        if ~isfolder(folder_path)
            fprintf('跳过不存在目录 [%s]：%s\n', net, folder_path);
            continue;
        end

        fprintf('读取 [%s]：%s\n', net, folder_path);

        tmp = batch_read_earthquake_data(folder_path);

        % 给每条记录补充来源信息，后面调试时方便；不影响现有流程
        for i = 1:numel(tmp)
            if isempty(tmp{i}), continue; end
            tmp{i}.source_network = net;
            tmp{i}.event_code = event_info.event_code;
            tmp{i}.group_m = m;
            tmp{i}.event_n = n;
        end

        adata = [adata; tmp(:)]; %#ok<AGROW>
    end

    if isempty(adata)
        warning('事件 m=%d, n=%d 没有读到任何数据。', m, n);
    else
        fprintf('事件读取完成：共 %d 条分量记录。\n', numel(adata));
    end
end