%% ffbp.m
% 一代快速反投影，将坐标参数的计算剥离
% 沿着每一列遍历矩阵速度更快

%% parameter used
load('rf3.mat');
Signal = ZEG1(:,:,3);

% system parameter
fs          = 25.51e6;                      % sampling frequency [samples/sec]
c           = 1540;                         % speed of sound [m/s]
Ts          = 1/fs;                         % sampling interval 
Ds          = c*Ts;                         % distance unit for one sampling interval
delay       = 0*Ts;
sensorNum   = 128;

% ROI parameter
center_x    =  0;
center_y    =  0;
Width       =  0.035;
Higth       =  0.035;
left_top_x  = -Width/2;
left_top_y  =  0;
Npx         =  256;
Npy         =  256;
N           =  Npx*Npy;
ddx         =  Width/Npx;
ddy         =  Higth/Npy;
ddl         =  sqrt(ddx^2+ddy^2);
dds         =  0.038/sensorNum;

%% 滤波
fL = 200e3;
fH = 10e6;
width           = 2^nextpow2(size(Signal,1));
bw              = round((fL)*(width/2)/(fs/2));
wc              = round((fH)*(width/2)/(fs/2));
proj_fft        = fft(Signal,width);    
filter          = ones(1,width);

filter(1:bw+2) = 0;           % filter(2:10)与filter(2040:2048)为低频部分filter(1)为常值分量
filter(width-bw:width) = 0;
filter = repmat(filter,sensorNum,1);
proj_filtered = proj_fft.*filter';

Signal          = ifft(proj_filtered);
Signal          = Signal(1:Nline,:);

%% ROI坐标
rij      = zeros(N,2,'double'); 
ROIx     = linspace(left_top_x,-left_top_x,Npx);
ROIx     = repmat(ROIx,Npy,1);
rij(:,1) = reshape(ROIx,N,1);
ROIy     = linspace(left_top_y,-Higth,Npy)';
ROIy     = repmat(ROIy,1,Npx);
rij(:,2) = reshape(ROIy,N,1);

%% 传感器坐标
pitch   = 0.0003;
rm      = zeros(sensorNum,2,'double');    
rm(:,1) = linspace(-0.019,0.019,sensorNum);
rm(:,2) = left_top_y+ddy;

%% 网格文件生成
D = zeros(N,sensorNum);
COS = zeros(N,sensorNum);

for i = 1:128
    for t = 1:N
        D(t,i) = sqrt((rm(i,1)-rij(t,1))^2 + (rm(i,2)-rij(t,2))^2);
        COS(t,i)   = (left_top_y + ddy - rij(t,2))/D(t,i);
    end
end
POINT = ceil(D./Ds);

%% 快速反投影
tic;
P = zeros(N,128);
for i = 1:128
    for t = 1:N
        P(t,i) = -Signal(POINT(t,i),i);
    end
end
P = P.* COS;
P = sum(P,2);
toc;

recon_bp        = reshape(P, Npy, Npx);
recon_nor1       = recon_bp/max(recon_bp(:));
img_index_x = linspace(left_top_x, -left_top_x, Npx); 
img_index_y = linspace(0, -Higth, Npy);
figure;
hold on; imagesc(img_index_x*1e3, img_index_y*1e3, recon_nor1,[0.05,1] );
axis image;
axis off;
colormap(hot);
% colorbar;
xlabel('mm');
ylabel('depth(mm)');
