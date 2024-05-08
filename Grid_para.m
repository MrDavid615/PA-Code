function [POINT,COS,N,img_index_x,img_index_y] = Grid_para(Width,Higth,Npx,Npy,vs)

fs          = 25.51e6;                      % sampling frequency [samples/sec]
c           = vs;                         % speed of sound [m/s]
Ts          = 1/fs;                         % sampling interval 
Ds          = c*Ts;                         % distance unit for one sampling interval
sensorNum   = 128;

left_top_x  = -Width/2;
left_top_y  =  0;
N           =  Npx*Npy;
ddx         =  Width/Npx;
ddy         =  Higth/Npy;
ddl         =  sqrt(ddx^2+ddy^2);
dds         =  0.038/sensorNum;

rij      = zeros(N,2,'double'); 
ROIx     = linspace(left_top_x,-left_top_x,Npx);
ROIx     = repmat(ROIx,Npy,1);
rij(:,1) = reshape(ROIx,N,1);
ROIy     = linspace(left_top_y,-Higth,Npy)';
ROIy     = repmat(ROIy,1,Npx);
rij(:,2) = reshape(ROIy,N,1);

pitch   = 0.0003;
rm      = zeros(sensorNum,2,'double');    
rm(:,1) = linspace(-0.019,0.019,sensorNum);
rm(:,2) = left_top_y+ddy;

D = zeros(N,sensorNum);
COS = zeros(N,sensorNum);

for i = 1:128
    for t = 1:N
        D(t,i) = sqrt((rm(i,1)-rij(t,1))^2 + (rm(i,2)-rij(t,2))^2);
        COS(t,i)   = (left_top_y + ddy - rij(t,2))/D(t,i);
    end
end
POINT = ceil(D./Ds);

img_index_x = linspace(left_top_x, -left_top_x, Npx); 
img_index_y = linspace(0, -Higth, Npy);

end