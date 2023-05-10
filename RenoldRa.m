function Qcz=RenoldRa(DH,DZ,Deep,e)
%% _________input variables_________________
% D   轴颈直径        c   半径间隙 单位为m
% L   轴承长度        l   内衬厚2
% nd0 水的粘度 Pa.s   w   自转转速 r/min--rad/s
% E1  内衬弹性模量    v1  橡胶泊松比
% E2  轴弹性模量      v2  轴泊松比
% v0  当量泊松比      W   轴承外载荷
% DZ  水槽数量        DH  水槽角度
% ps  进口压力无量纲参数 DZ  水槽数量
% M   轴承重量        g   重力加速度
% Ra1 轴承表面粗糙度  Ra2 轴径表面粗糙度

%% _________output variable________
% PT  轴承水膜        HT  轴承水膜厚度分布
% Fx  轴承油膜承载力  Fy  轴承油膜y方向油膜承载力
% Fp  由液膜压力流产生的摩擦力
% Fc  由液膜剪切流产生的摩擦力
% Fmc 总摩擦力
% fmc 总摩擦系数

%% ___________testing data_________________
dbstop if error
clc;close all;
ps=1;
D=0.320; c=0.0005; L=0.561;l=0.0145;nd0=0.0008994;
E1=440*10^6; v1=0.45; E2=207*10^9;v2=0.3; M=11.5; g=9.816;
Ra1=0.8*10^(-6); Rendax1=10*Ra1;Renday1=10*Ra1;
Ra2=2.4*10^(-6);Rendax2=10*Ra2;Renday2=10*Ra2;
Cmd=5*10^(-6)/(Ra1^2);Cbeta=10^4*Ra1; %粗糙度计算模型
%% __________formal function_________________
R=D/2;E=2*((1-v1^2)/E1+(1-v2^2)/E2)^(-1);
m=360; n=100; %轴、周向等分数；
deltL=L/m; %轴向等分间距
deltsita=(2*pi-DZ*DH)/DZ/n; %周向角度等分大小
DDL=0:L/m:L;
Dr=deltsita;Dz=deltL;
Dc=Dccal(m,n,Dr,Dz);
Rc1=Cucao(Ra1,Rendax1,Renday1,n,m,D,L,DZ);%轴承内表面随机粗糙峰分布
Rc2=Cucao(Ra2,Rendax2,Renday2,n,m,D,L,DZ);%轴径表面随机粗糙峰分布
save DATA.mat c Cbeta Cmd D DDL deltL deltsita DZ E E1 E2 g l L m M n nd0 ps R Ra1 Ra2 Rendax1 Rendax2 Renday1 Renday2 v1 v2 DH Dr Dz Rc1 Rc2
Qx=-2*1000;Qy=0;
%% 主轴变形计算程序
w=200*2*pi/60;fai0=0;
VV=6.*nd0.*w.*R.^2./(c.^2); %速度无量纲参数fai0=0; %初始偏心量和偏位角
JHcs=[DH*R,Deep];ks=4;err=1;
% while err>0.005
[PT,HT,WY,Fx,Fy,~]=DD2(w,VV,e,fai0,JHcs,ks,Dc,Rc1,Rc2);
fai=atan(-Fx/Fy);
% err=abs((fai-fai0)/fai0);
% fai0=fai;
% end
Qcz=Fx*sin(abs(fai))+Fy*cos(abs(fai));
% [Fp,Fc,Fmc,fmc]=ucalcu(HT,PT,c,deltsita,deltL,R,w,Qx,Qy,M,g,nd0,VV);
%膜厚绘图
% PHplot(m,L,sitas,PT,HT)
end
