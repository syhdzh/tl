function [Fx,Fy]=F2calcu(deltsita,ps,PT,DZ,R,L,n,m)
%% __________________input variables_______________________
% deltsita 周向节点距离 ps  入口压力无量纲系数
% PT  动压力            DZ  水槽数量
% R   轴承半径          L   轴承宽度
% n,m   周、轴向节点数
%% __________________output variables_______________________
% Fx,Fy 轴承X、Y方向承载力 
%% __________________formal functions_______________________
dbstop if error
%   更新迭代
PT=PT*ps; %恢复成有量纲量
sitas=(0:n*DZ-1)*deltsita; %周向
for i=1:m
    A(:,i)=sin(sitas');
    B(:,i)=cos(sitas');
end
PA=PT.*A;PB=PT.*B;
Fx=-(sum(sum(PB))*deltsita*R*L/m);
Fy=(sum(sum(PA))*deltsita*R*L/m);
end