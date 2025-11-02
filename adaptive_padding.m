function [x_final, total_length, n_prepended, n_mirror] = adaptive_padding(x, varargin)
% 对信号进行自适应的“镜像-渐隐-零填充”处理
% 只在镜像段渐隐；零填充只加在尾部；无镜像时不乘窗
    p = inputParser;
    addRequired(p, 'x', @isnumeric);
    addParameter(p, 'BoundaryPercent', 5, @isnumeric);      %计算边界能量的窗口占总长的百分比
    addParameter(p, 'MaxMirrorPercent', 50, @isnumeric);    %允许的最大镜像长度占总长的百分比
    addParameter(p, 'BERThreshold', 0.05, @isnumeric);      %判断边界是否“干净”的能量比阈值
    addParameter(p, 'MinMirrorSec', 0, @isnumeric);         % 可选参数 保证最短镜像（点数由外部fs换算）
    addParameter(p, 'FIRLen', 0, @isnumeric);               % 可选参数 若做频域FIR，填冲激响应长度
    parse(p, x, varargin{:});

    x = x(:);
    npts = length(x);
    bp   = max(1, floor(npts * (p.Results.BoundaryPercent/100)));

    x_start = x(1:bp);
    x_end   = x(end-bp+1:end);
    Amax    = max(abs(x));
    if Amax==0, ber=0; else, ber = max(rms(x_start), rms(x_end))/Amax; end

    mirror_ratio = (p.Results.MaxMirrorPercent/100) * min(ber / p.Results.BERThreshold, 1.0);
    n_mirror = floor(npts * mirror_ratio);

disp(mirror_ratio);

    % 安全夹逼
    n_mirror = min(n_mirror, max(npts-2,0));
    if p.Results.MinMirrorSec>0
        % 可在外面用 fs 换算为点数传进来
        n_mirror = max(n_mirror, round(p.Results.MinMirrorSec)); 
        n_mirror = min(n_mirror, max(npts-2,0));
    end

    % 偶对称镜像（不包含端点）
    if n_mirror>0
        left  = flipud(x(2:n_mirror+1));
        right = flipud(x(end-n_mirror:end-1));
        x_m = [left; x; right];

        % 只在镜像段淡入/淡出，拼接点处=1
        win = ones(size(x_m));
        tw  = tukeywin(2*n_mirror, 1);    % 全余弦
        win(1:n_mirror) = tw(1:n_mirror);
        win(end-n_mirror+1:end) = tw(n_mirror+1:end);
        x_t = x_m .* win;
    else
        x_t = x; % 无镜像，不乘窗，避免压掉起始峰值
    end

    % 计算目标长度
    grow = length(x_t);
    Lh   = max(0, p.Results.FIRLen - 1);
    total_length = 2^nextpow2(grow + Lh);

    % 只在尾部零填充（避免前端平移）
    n_post = total_length - grow;
    x_final = [x_t; zeros(n_post,1)];

    % 返回前端需要丢弃的点数 即镜像的长度 后续IFFT剪掉
    n_prepended = n_mirror;
end
