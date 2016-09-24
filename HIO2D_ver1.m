function [gEnd err] = HIO2D_ver1(S,support,itnum,alpha,gInit)
% ���Ǳ�׼�汾�� ��ά HIO �㷨

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
    
    g  = (g2>=0).* g2 + (g2<0).*(g - alpha * g2);
%     disp(strcat('HIO:',num2str(k)));
end

gEnd = g;