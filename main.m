clear; close all; clc;

root_dir = 'D:\FAST\earthquakedir';
event_dirs = build_earthquake_read_dirs(root_dir);

T_list  = [5, 2, 1, 0.5, 0.1, 0.02];
fc_high = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

batch_result = process_selected_events( ...
    event_dirs, 7, T_list, fc_high, ...
    'SaveEachEvent', true, ...
    'SaveFolder', 'D:\FAST\processed_metric_tables', ...
    'KeepFullResult', false);

metric_all = batch_result.metric_all;

height(metric_all)
unique(metric_all.event_code)
metric_all.Properties.VariableNames'