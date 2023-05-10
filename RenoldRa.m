function Qcz=RenoldRa(DH,DZ,Deep,e)
%% _________input variables_________________
% D   �ᾱֱ��        c   �뾶��϶ ��λΪm
% L   ��г���        l   �ڳĺ�2
% nd0 ˮ��ճ�� Pa.s   w   ��תת�� r/min--rad/s
% E1  �ڳĵ���ģ��    v1  �𽺲��ɱ�
% E2  �ᵯ��ģ��      v2  �Ჴ�ɱ�
% v0  �������ɱ�      W   ������غ�
% DZ  ˮ������        DH  ˮ�۽Ƕ�
% ps  ����ѹ�������ٲ��� DZ  ˮ������
% M   �������        g   �������ٶ�
% Ra1 ��б���ֲڶ�  Ra2 �ᾶ����ֲڶ�

%% _________output variable________
% PT  ���ˮĤ        HT  ���ˮĤ��ȷֲ�
% Fx  �����Ĥ������  Fy  �����Ĥy������Ĥ������
% Fp  ��ҺĤѹ����������Ħ����
% Fc  ��ҺĤ������������Ħ����
% Fmc ��Ħ����
% fmc ��Ħ��ϵ��

%% ___________testing data_________________
dbstop if error
clc;close all;
ps=1;
D=0.320; c=0.0005; L=0.561;l=0.0145;nd0=0.0008994;
E1=440*10^6; v1=0.45; E2=207*10^9;v2=0.3; M=11.5; g=9.816;
Ra1=0.8*10^(-6); Rendax1=10*Ra1;Renday1=10*Ra1;
Ra2=2.4*10^(-6);Rendax2=10*Ra2;Renday2=10*Ra2;
Cmd=5*10^(-6)/(Ra1^2);Cbeta=10^4*Ra1; %�ֲڶȼ���ģ��
%% __________formal function_________________
R=D/2;E=2*((1-v1^2)/E1+(1-v2^2)/E2)^(-1);
m=360; n=100; %�ᡢ����ȷ�����
deltL=L/m; %����ȷּ��
deltsita=(2*pi-DZ*DH)/DZ/n; %����Ƕȵȷִ�С
DDL=0:L/m:L;
Dr=deltsita;Dz=deltL;
Dc=Dccal(m,n,Dr,Dz);
Rc1=Cucao(Ra1,Rendax1,Renday1,n,m,D,L,DZ);%����ڱ�������ֲڷ�ֲ�
Rc2=Cucao(Ra2,Rendax2,Renday2,n,m,D,L,DZ);%�ᾶ��������ֲڷ�ֲ�
save DATA.mat c Cbeta Cmd D DDL deltL deltsita DZ E E1 E2 g l L m M n nd0 ps R Ra1 Ra2 Rendax1 Rendax2 Renday1 Renday2 v1 v2 DH Dr Dz Rc1 Rc2
Qx=-2*1000;Qy=0;
%% ������μ������
w=200*2*pi/60;fai0=0;
VV=6.*nd0.*w.*R.^2./(c.^2); %�ٶ������ٲ���fai0=0; %��ʼƫ������ƫλ��
JHcs=[DH*R,Deep];ks=4;err=1;
% while err>0.005
[PT,HT,WY,Fx,Fy,~]=DD2(w,VV,e,fai0,JHcs,ks,Dc,Rc1,Rc2);
fai=atan(-Fx/Fy);
% err=abs((fai-fai0)/fai0);
% fai0=fai;
% end
Qcz=Fx*sin(abs(fai))+Fy*cos(abs(fai));
% [Fp,Fc,Fmc,fmc]=ucalcu(HT,PT,c,deltsita,deltL,R,w,Qx,Qy,M,g,nd0,VV);
%Ĥ���ͼ
% PHplot(m,L,sitas,PT,HT)
end
