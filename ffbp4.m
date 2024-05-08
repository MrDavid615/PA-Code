% ffbp.m
% 一代快速反投影，将坐标参数的计算剥离
% 沿着每一列遍历矩阵速度更快
%% ffbp2.m
% 二代快速反投影，将坐标参数计算封装成函数
%% ffbp3.m
% 三代快速反投影，将反投影封装为函数
%% ffbp4.m
% 四代快速反投影，将反投影封装为C语言函数

%%
load('rf3.mat');
Signal = ZEG1(:,:,3);
Npx = 512;
Npy = 512;
[POINT,COS,N,img_index_x,img_index_y] = Grid_para(0.035,0.035,Npx,Npy);

%% 快速反投影
tic;
POINT = int32(POINT);
% P = FFBP_fun_C(N,Signal',POINT',COS');
% P = FFBP_fun_C_2(N,Signal',POINT,COS);
toc;

recon_bp        = reshape(P(1,:), Npy, Npx);
recon_nor       = recon_bp/max(recon_bp(:));
figure;
hold on; imagesc(img_index_x*1e3, img_index_y*1e3, recon_nor,[0,1] );
axis image;
axis off;
colormap(hot);
% colorbar;
xlabel('mm');
ylabel('depth(mm)');
