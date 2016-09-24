function Mout = MatMap(M,ymin,ymax)
% ��������Ԫ�صķ�Χת���� [ymin ymax] ��Χ�ڣ�����double����

data = M(:);
data = double(data);
mapdata = (ymax - ymin)*((data - min(data))/(max(data) - min(data)))...
    + ymin;

Mout = reshape(mapdata,size(M));