function [y_corrected, trend, coeffs] = polyfit_constrained_zerobase(t, y, order)
% POLYFIT_CONSTRAINED_ZEROBASE - 拟合一个N阶多项式，强制其0阶和1阶系数为0。
%
% Inputs:
%   t     - 时间向量 (列)
%   y     - 需要校正的数据向量 (列)
%   order - 多项式的阶数 (6)
%
% Outputs:
%   y_corrected - 校正后的数据 (y - trend)
%   trend       - 拟合出的多项式漂移趋势
%   coeffs      - 拟合出的系数 [c_order, c_{order-1}, ..., c2]

    % 确保输入是列向量
    t_col = t(:);
    y_col = y(:);
        

    % 构建设计矩阵 A
    % 要拟合 y = c_n*t^n + ... + c_2*t^2 所以矩阵的每一列是 t 的某个幂次
    num_coeffs = order - 1; % 只求解 c2 到 c_order 这几个系数
    A = zeros(length(t_col), num_coeffs);
    
    for i = 1:num_coeffs
        % A的第一列是 t^order, 第二列是 t^(order-1), ..., 最后一列是 t^2
        A(:, i) = t_col .^ (order - i + 1);
    end

    % --- 求解线性方程组 Ac = y 的最小二乘解 ---
    coeffs = A \ y_col;

    % --- 计算拟合出的漂移趋势 ---
    trend = A * coeffs;
    
    % --- 从原始数据中减去漂移 ---
    y_corrected = y_col - trend;
end