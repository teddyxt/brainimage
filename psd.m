function [psd_i] = psd(U,V,fmin,fmax,dsrate)
% Function calculate the psd of the video after pca(dimensionality
% reduciton): fn(t)=Un*V=sum[Un(i)Vt(i)], i=1:ncomps;
% fn(w)=sum[Un(i)V(w)], i=1:ncomps; psd = fn(w)*fn(w)',w=freq band of
% interest, fn(w)': complex conjugate of fn(w)


%Un: U of the video; Vt: V of the video; fmin, fmax: frequency band of
%interest; dsrate:downsampling rate

Fs = 40/dsrate;       % Sampling rate   
                      % video fs after downsampling: 40/dsrate = 40/5 = 8Hz
L = size(V,2);        % Length of signal
% f = Fs*(1/L:1/L:1); % frequency vector : L number of points make up Fs
frange = ceil((fmin/Fs)/(1/L):(fmax/Fs)/(1/L)); % bandwidth of interest

% Un = reshape(U,[],size(U,3));
% V_w = fft(V);
% fn(w)=sum[Un(i)V(w)]
fn_w = reshape(U,[],size(U,3))*fft(V); %fft(V) apply FT to every column of V
psd_w = fn_w.* conj(fn_w); % power spectrum in 3D
psd_i = reshape(sum(psd_w(:,frange),2),size(U,1),size(U,2));% psd of the frequency range of interest for each pixel
% save the psd as mat.file
psd_w = reshape(psd_w,size(U,1),size(U,2),[]);
fileID = fopen('video_psd.mat','w');
fwrite(fileID, psd_w, 'single');
fclose(fileID);
clear psd_w fn_w




