%% ffbp.m
% 一代快速反投影，将坐标参数的计算剥离
% 沿着每一列遍历矩阵速度更快
%% ffbp2.m
% 二代快速反投影，将坐标参数计算封装成函数

%% parameter used
load('rf3.mat');
Signal = ZEG1(:,:,3);
Npx = 512;
Npy = 512;
[POINT,COS,N,img_index_x,img_index_y] = Grid_para(0.035,0.035,Npx,Npy);

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
recon_nor       = recon_bp/max(recon_bp(:));
figure;
hold on; imagesc(img_index_x*1e3, img_index_y*1e3, recon_nor,[0.05,1] );
axis image;
axis off;
colormap(hot);
% colorbar;
xlabel('mm');
ylabel('depth(mm)');
