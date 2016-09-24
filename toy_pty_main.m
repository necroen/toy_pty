% implementation of the ptychography algorithm
% this is a very simple toy program
% it will make two dir: det_mod and recon
% det_mod is the modulus of diffraction pattern get by the detector at
% different position.
% recon is the reconstruct image during the iteration.
clear all;clc;close all;

sampleMod = rgb2gray(imread('lena.jpg')); % sample image
sampleMod = uint8(MatMap(sampleMod,0,255));

sample = sampleMod;
%===========
% if the sample include phase information, it could write as:
% sample = sampleMod.*exp(i*samplePhase);
% the size of sample phase matrix should be same as sampleMod and range from -pi to
% pi
%===========

figure;
imshow(sample);
sample = MatMap(sample,0,1);  
%%
scanDim = 3;
scanNum = scanDim*scanDim;

Npic = size(sample,1);

startPos = 64+2;
endPos = 192-2; % this two args will control the overlap rate

step = (endPos-startPos)/(scanDim-1);   % 两个方向步长一致
[xx yy] = meshgrid(startPos:step:endPos);
xx = round(xx(:));
yy = round(yy(:));  % 注意这些数要是整数  center of the support at different postion

% m = 2*(min(xx)-2);  % pinhole(support)的直径  
m = 100; % the radius of support
%%
maskDB = zeros(Npic,Npic,scanNum);
maskAll = zeros(Npic,Npic);
for i = 1:scanNum;
    maskDB(:,:,i) = circle_mask(Npic,m,xx(i),yy(i));
    maskAll = maskAll + maskDB(:,:,i);
end

% maskAll = (maskAll>0);
figure;
imagesc(maskAll);
axis square;
clear endPos i startPos step;
%%
mkdir('notCenter');
caijianDB = zeros(size(maskDB));
for i = 1:scanNum;
    tmp = maskDB(:,:,i);
    caijianDB(:,:,i) = tmp.*sample;
    imwrite(caijianDB(:,:,i),strcat(pwd,'\notCenter\',num2str(i),'.jpg'));
end
clear i tmp;
%%
% mkdir('Center');
caijianCenterDB = zeros(size(caijianDB));
for i = 1:scanNum;
    tmp = caijianDB(:,:,i);
    caijianCenterDB(:,:,i) = circshift(tmp,[Npic/2 - xx(i),Npic/2 - yy(i)]);
%     imwrite(caijianCenterDB(:,:,i),strcat(pwd,'\Center\',num2str(i),'.jpg'));
end
clear i tmp caijianDB;
%% 镶嵌
N = 512;                         % 衍射图片的尺寸
det_mod = zeros(N,N,scanNum);
L = Npic/2 -1;

for i = 1:scanNum
    det_mod(N/2-L:N/2+1+L , N/2-L:N/2+1+L,i) = caijianCenterDB(:,:,i);
end
clear i;
%% 得到 det_mod
for i = 1:scanNum
    det_mod(:,:,i) = abs(fftshift(fft2(det_mod(:,:,i))));
end
clear i;
%% 将det_mod 适当处理后输出
mkdir('det_mod');

for i = 1:scanNum
    tmp = det_mod(:,:,i);
    tmp = log(1+tmp);
    tmp = uint8(MatMap(tmp,0,255));
    imwrite(tmp,strcat(pwd,'\det_mod\',num2str(i),'.jpg'));
end
clear tmp i;
%%
support1 = circle_mask(N,m,N/2,N/2);
%% 开始用 HIO算法 做相位恢复
recon = rand(size(sample));

itnum = 50;   % max iteration number
figure;
%%
mkdir('recon');
for k = 1:itnum
    disp(strcat('=========================',num2str(k)));
    for i = 1:scanNum
        disp(strcat('========scanNum:',num2str(i)));
        % 先从recon中截取一部分出来放到一个 N x N 的矩阵中
        gInit = zeros(N,N);
        tmp = recon.*maskDB(:,:,i);  % 这是一个 Npic x Npic 大小的矩阵
        tmp = circshift(tmp,[Npic/2 - xx(i),Npic/2 - yy(i)]);
        gInit(N/2-L:N/2+1+L , N/2-L:N/2+1+L) = tmp; % 放到 N x N 的矩阵中
        % 做 HIO 更新幅值得到 g1
        g1 = HIO2D_ver1(det_mod(:,:,i),support1,1,0.7,gInit); %
        % 将 g1 裁剪后 再平移成和maskDB一致的位置 再回填到 recon 中
        g1 = g1(N/2-L:N/2+1+L , N/2-L:N/2+1+L);
        g1 = circshift(g1,[xx(i) - Npic/2,yy(i) - Npic/2]);
        recon(find(maskDB(:,:,i)==1)) = g1(find(maskDB(:,:,i)==1));
    end
    % 记录 err
    % 显示 recon
    imshow(recon);
    title(strcat('第',num2str(k),'次迭代结果'));
    pause(0.1);
    imwrite(recon,strcat(pwd,'\recon\',num2str(k),'.jpg'));
end
clear tmp k i g1 gInit;