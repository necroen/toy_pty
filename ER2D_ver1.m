function [gEnd err] = ER2D_ver1(S,support,itnum,gInit)
% ���Ǳ�׼�汾�� ��ά ER �㷨

g = gInit;
err = [];

for k=1:itnum
    G  = fftshift(fft2(g));
    
    errEnd = sum(sum( (abs(G) - S).^2 ));
    err = [err errEnd];
    
    G2 = S .* (G ./ abs(G));
    g2 = ifft2(ifftshift(G2));
    
    % �� g2 �Ա���Щ��λ������ȥ��
    g2 = g2.*support;     % ���ǳ���Ҫ����Ȼ�ָ�����������ֵ��Ʈ����λһ����
    %
    
    g  = (g2>=0).*g2;
%    disp(strcat('ER:',num2str(k)));
end

gEnd = g;