function [Fx,Fy]=F2calcu(deltsita,ps,PT,DZ,R,L,n,m)
%% __________________input variables_______________________
% deltsita ����ڵ���� ps  ���ѹ��������ϵ��
% PT  ��ѹ��            DZ  ˮ������
% R   ��а뾶          L   ��п��
% n,m   �ܡ�����ڵ���
%% __________________output variables_______________________
% Fx,Fy ���X��Y��������� 
%% __________________formal functions_______________________
dbstop if error
%   ���µ���
PT=PT*ps; %�ָ�����������
sitas=(0:n*DZ-1)*deltsita; %����
for i=1:m
    A(:,i)=sin(sitas');
    B(:,i)=cos(sitas');
end
PA=PT.*A;PB=PT.*B;
Fx=-(sum(sum(PB))*deltsita*R*L/m);
Fy=(sum(sum(PA))*deltsita*R*L/m);
end