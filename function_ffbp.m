function PA_Image = function_ffbp(Signal,POINT,COS)

% figure;
% subplot(211);
% plot(Signal(1:1024,65));
fL = 1e3;
fH = 10e6;
fs = 25.51e6;
width           = 2^nextpow2(size(Signal,1));
bw              = round((fL)*(width/2)/(fs/2));
wc              = round((fH)*(width/2)/(fs/2));
proj_fft        = fft(Signal,width);
filter          = ones(1,width);

filter(1:bw+2) = 0;           % filter(2:10)与filter(2040:2048)为低频部分filter(1)为常值分量
filter(width-bw:width) = 0;
filter = repmat(filter,128,1);
proj_filtered = proj_fft.*filter';
Signal          = ifft(proj_filtered);
% subplot(212);
% plot(Signal(1:1024,65));

N = length(POINT(:,1));
PA_Image = zeros(N,128);
for i = 1:128
    for t = 1:N
        PA_Image(t,i) = -Signal(POINT(t,i),i);
    end
end
PA_Image = PA_Image.* COS;
PA_Image = sum(PA_Image,2);

end