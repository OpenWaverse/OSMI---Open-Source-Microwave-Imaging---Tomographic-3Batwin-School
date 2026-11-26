function [nrmse_re, nrmse_im, nrmse_avg] = compute_nrmse(A, B)
% A: recovered complex contrast vector [Ncell x 1]
% B: ground truth complex contrast vector [Ncell x 1]
%
% Formula:
%   NRMSE(q) = sqrt( mean_j[ (A_j(q) - B_j(q))^2 ] ) / Delta(q)
%   Delta(q) = min( max(A(q)) - min(A(q)), max(B(q)) - min(B(q)) )

% Real component (q = real)
A_re     = real(A(:));
B_re     = real(B(:));
delta_re = min(max(A_re) - min(A_re), max(B_re) - min(B_re));
rmse_re  = sqrt(mean((A_re - B_re).^2));
nrmse_re = rmse_re / max(delta_re, eps);

% Imaginary component (q = imag)
A_im     = imag(A(:));
B_im     = imag(B(:));
delta_im = min(max(A_im) - min(A_im), max(B_im) - min(B_im));
rmse_im  = sqrt(mean((A_im - B_im).^2));
nrmse_im = rmse_im / max(delta_im, eps);

% Mean across both components
nrmse_avg = 0.5 * (nrmse_re + nrmse_im);
end