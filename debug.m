%% 流程概述
% 读+配对
% 整体性预处理（日本+美国那个流程，含：
%   自适应加窗及首尾pad填充（镜像填充→加窗→零填充）
%   多项式回扣

%% 批量读
adata=batch_read_earthquake_data('D:\FAST\20250602035200knt');
%adata=batch_read_earthquake_data('D:\FAST\data');



%%  窄带滤波
for i = 1:numel(adata)
    adata{i, 1} = nb_filt(adata{i, 1},2);
end
%data = nb_filt(data, );   % 会在 data 里新增 nb_acc_5, nb_v_5


%% 积分得v
for i = 1:numel(adata)
    adata{i, 1} = acc2vel(adata{i, 1});
end
%% 
x=adata{1, 1}.acceleration_gal;
%x=adata{1, 1}.velocity_cms;
%x=adata{1, 1}.nb_acc_2;
%x=adata{1, 1}.nb_v_1;

fs=adata{1, 1}.sampling_freq_hz;

%%
plotSpectrum(x, fs);

