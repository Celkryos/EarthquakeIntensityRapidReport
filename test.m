

%%
%x=adata{n, 1}.acceleration_gal;
%x=detrend(x,1);
%[x_final, total_length, n_prepended] = adaptive_padding(x);

%%
n=1;
adata{n, 1} = nb_filt(adata{n, 1},0.02);

%%
adata{n, 1} = acc2vel(adata{n, 1});
%%
[corrected_d, fitted_trend, coeffs] = polyfit_constrained_zerobase(adata{n, 1}.time, adata{n, 1}.d_cm,6);


