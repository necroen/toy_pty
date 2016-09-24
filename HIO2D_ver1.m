function [gEnd err] = HIO2D_ver1(S,support,itnum,alpha,gInit)
% 这是标准版本的 二维 HIO 算法

g = gInit;
err = [];

for k=1:itnum
    G  = fftshift(fft2(g));
    
    errEnd = sum(sum( (abs(G) - S).^2 ));
    err = [err errEnd];
    
    G2 = S .* (G ./ abs(G));
    g2 = ifft2(ifftshift(G2));
    
    % 对 g2 旁边那些相位噪音做去除
    g2 = g2.*support;     % 这句非常重要，不然恢复不出来，幅值乱飘，相位一团糟
    %
    
    g  = (g2>=0).* g2 + (g2<0).*(g - alpha * g2);
%     disp(strcat('HIO:',num2str(k)));
end

gEnd = g;