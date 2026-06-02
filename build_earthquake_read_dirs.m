function event_dirs = build_earthquake_read_dirs(root_dir, varargin)
% BUILD_EARTHQUAKE_READ_DIRS
% 构造两层事件读取目录：
%   第 1 层：组别 m，例如 root_dir\1
%   第 2 层：组内第 n 个事件，例如 20220122010830_ascii
%
% 每个事件下自动构造：
%   ...\事件编号\kik\ascii
%   ...\事件编号\knet\ascii
%
% 输出 event_dirs:
%   event_dirs(iGroup).m
%   event_dirs(iGroup).events(n).event_folder
%   event_dirs(iGroup).events(n).event_code
%   event_dirs(iGroup).events(n).data_dirs
%   event_dirs(iGroup).events(n).network_exists

    p = inputParser;
    addRequired(p, 'root_dir', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Networks', {'kik', 'knet'}, @(x)iscell(x) || isstring(x));
    parse(p, root_dir, varargin{:});

    root_dir = char(root_dir);
    networks = cellstr(string(p.Results.Networks));

    if ~isfolder(root_dir)
        error('根目录不存在：%s', root_dir);
    end

    % ---------- 找组目录：只接受纯数字文件夹，如 1,2,3 ----------
    group_raw = dir(root_dir);
    group_raw = group_raw([group_raw.isdir]);

    group_names = {group_raw.name};
    is_group = ~ismember(group_names, {'.', '..'}) & ...
               ~cellfun(@isempty, regexp(group_names, '^\d+$', 'once'));

    group_raw = group_raw(is_group);

    if isempty(group_raw)
        warning('根目录下没有找到数字组目录：%s', root_dir);
        event_dirs = struct('m', {}, 'group_folder', {}, 'group_path', {}, 'events', {});
        return;
    end

    group_nums = str2double({group_raw.name});
    [group_nums, order] = sort(group_nums);
    group_raw = group_raw(order);

    event_dirs = struct('m', {}, 'group_folder', {}, 'group_path', {}, 'events', {});

    % ---------- 逐组扫描事件 ----------
    for iG = 1:numel(group_raw)
        m = group_nums(iG);
        group_folder = group_raw(iG).name;
        group_path = fullfile(root_dir, group_folder);

        ev_raw = dir(fullfile(group_path, '*_ascii'));
        ev_raw = ev_raw([ev_raw.isdir]);

        if isempty(ev_raw)
            events = struct('m', {}, 'n', {}, 'event_folder', {}, ...
                            'event_code', {}, 'event_base', {}, ...
                            'networks', {}, 'data_dirs', {}, ...
                            'network_exists', {}, 'read_dirs', {});
        else
            [~, ev_order] = sort({ev_raw.name});
            ev_raw = ev_raw(ev_order);

            events = struct('m', {}, 'n', {}, 'event_folder', {}, ...
                            'event_code', {}, 'event_base', {}, ...
                            'networks', {}, 'data_dirs', {}, ...
                            'network_exists', {}, 'read_dirs', {});

            for n = 1:numel(ev_raw)
                event_folder = ev_raw(n).name;                 % eg. 20220122010830_ascii
                event_code   = regexprep(event_folder, '_ascii$', '');  % eg. 20220122010830

                event_base = fullfile(group_path, event_folder, event_code);

                data_dirs = cell(1, numel(networks));
                network_exists = false(1, numel(networks));

                for k = 1:numel(networks)
                    data_dirs{k} = fullfile(event_base, networks{k}, 'ascii');
                    network_exists(k) = isfolder(data_dirs{k});
                end

                events(n).m = m;
                events(n).n = n;
                events(n).event_folder = event_folder;
                events(n).event_code = event_code;
                events(n).event_base = event_base;
                events(n).networks = networks;
                events(n).data_dirs = data_dirs;
                events(n).network_exists = network_exists;
                events(n).read_dirs = data_dirs(network_exists);
            end
        end

        event_dirs(iG).m = m;
        event_dirs(iG).group_folder = group_folder;
        event_dirs(iG).group_path = group_path;
        event_dirs(iG).events = events;
    end

    n_events = sum(arrayfun(@(g)numel(g.events), event_dirs));
    fprintf('已构建目录索引：%d 个组，%d 个事件。\n', numel(event_dirs), n_events);
end