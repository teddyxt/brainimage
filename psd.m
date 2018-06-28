function [psd] = (U,V,fmin,fmax,dsrate)
% Function calculate the psd of the video after pca(dimensionality
% reduciton): fn(t)=Un*V=sum[Un(i)Vt(i)], i=1:ncomps;
% fn(w)=sum[Un(i)V(w)], i=1:ncomps; psd = fn(w)*fn(w)',w=freq band of
% interest, fn(w)': complex conjugate of fn(w)

% video fs after downsampling: 40/dsrate = 40/5 = 8Hz

%Un: U of the video; Vt: V of the video; fmin, fmax: frequency band of
%interest; dsrate:downsampling rate
Fs = 40/dsrate;       % Sampling rate
T = 1/Fs;             % Sampling period       
L = size(V,2);        % Length of signal
t = (0:L-1)*T;        % Time vector

Un = reshape(U,[],size(U,3));
V_w = fft(V);
fn_w = sum(Un*V_w,2);



