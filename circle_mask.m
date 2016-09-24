function maski = circle_mask(N,m,x,y)
% �ĸ����������������������������ƶ�
% ����һ�� N x N �ߴ��������� (x,y) �� ����һ���뾶Ϊ m ��Բ
% (x,y)ΪԲ�ĵ�
% ����ͼƬ����ϵ�����ϵ���Ϊ x ����������Ϊ y
% m ����Ϊż��
% 
% �����������Բ�ľ���

mOKflag = 1;
if mod(m,2)==1
    disp('ֱ��m����ż����');
    mOKflag = 0;
end

r = m/2;

xyOKflag = 1;
if (x>=r)&(x<=(N-r))&(y>=r)&(y<=(N-r))
%     disp('x y ok');
    xyOKflag = 1;
else
    disp('x y error.��������Բ���ᴥ�ߣ�');
    xyOKflag = 0;
end

if (xyOKflag == 1)&(mOKflag == 1)
    % ��ʼ��Բ �������ƶ���ָ��λ����
    [xx yy] = meshgrid(-N/2:N/2-1);
    z = sqrt(xx.^2 + yy.^2);
    clear xx yy xyOKflag;
    z = (z<r);       % ����z<=r ʹ�� sum(sum(z)) ���ӽ��� round(pi*r*r)
    z = circshift(z,[x-N/2,y-N/2]);
else
    disp('PIEmask something is wrong!');
    z = [];
end

maski = double(z);