%% ffbp.m
% 一代快速反投影，将坐标参数的计算剥离
% 沿着每一列遍历矩阵速度更快
%% ffbp2.m
% 二代快速反投影，将坐标参数计算封装成函数
%% ffbp3.m
% 三代快速反投影，将反投影封装为函数

%%
% load('rf3.mat');
% Signal = Filtered_ZEG;
Signal = ZEG1(:,:,1);
Signal(10:4096,:) = Signal(1:4087,:);
Npx = 512;
Npy = 512;
[POINT,COS,N,img_index_x,img_index_y] = Grid_para(0.038,0.038,Npx,Npy,1504);

%% 快速反投影
tic;
P = function_ffbp(Signal,POINT,COS);
toc;

recon_bp        = reshape(P, Npy, Npx);
recon_nor       = recon_bp/max(recon_bp(:));
recon_tem       = abs(hilbert(recon_nor));
recon_tem       = recon_tem/max(recon_tem(:));
% L = recon_nor<0;
% recon_nor(L) = 0.01;
% L = recon_nor==0;
% recon_nor(L) = 0.01;
G = figure;
% hold on; imagesc(img_index_x*1e3, img_index_y*1e3, log10(recon_tem),[-1,0] );
% hold on; imagesc(img_index_x*1e3, img_index_y*1e3, log10(recon_nor),[-1,0] );
hold on; imagesc(img_index_x*1e3, img_index_y*1e3, recon_tem,[0,1] );
% hold on; imagesc(img_index_x*1e3, img_index_y*1e3, recon_nor,[0.05,1] );
axis image;
% axis off;
colormap(jet);
colorbar;
xlabel('mm');
ylabel('depth(mm)');
