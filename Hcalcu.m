function [H,Hz,Wy,jd]=Hcalcu(csj,deltsita,fai,e,c,P,ps,Dc,Dh,n,m,E,Rc1,Rc2)
% ������㴦ˮĤ���
%% __________________input variables_______________________
% csj �Ӵ��������Ų��� deltsita ����ڵ����
% c   ��м�϶           P    �Ӵ�ѹ��
% ps  ��ˮѹ��������ϵ�� Dh   ˮ�ۼ��β���
% e   ���ƫ����         fai  ƫλ��
% Rc1 ����ڴֲڱ������ Rc2  ����ֲڱ������
% n   �Ӵ�������ڵ���   m    �Ӵ���������ڵ���
% Dc  ���ξ���ϵ��       DDBX ���������
%% __________________output variables_______________________
% H   ��ѹĤ��           Hz   ѹ��С��0������
% WY  �ĵ����           jd   �Ƕȱ���
%% __________________formal function_______________________
dbstop if error
H=zeros(n+1,m+1);PP=zeros(n+1,m+1);Hz=zeros(n+1,m+1);
Wy=WyCalcu(P,m,n,E,ps,Dc);
KK=zeros(n+1,m+1);
Dh=[Dh;Dh(1)]/c;jd=zeros(n+1,m+1);
for i=1:n+1 %����
    for j=1:m+1 %����
        jd(i,j)=csj+(i-1)*deltsita-fai;
        PP(i,j)=0;
        KK(i,j)=Wy(i,j)/c+Dh(i);
        if KK(i,j)>max(Dh(i))
            KK(i,j)=max(Dh(i));
        end
        H(i,j)=1+e*(cos(jd(i,j)))+Wy(i,j)/c+Rc1(i,j)/c+Rc2(i,j)/c+PP(i,j)+KK(i,j);%��ʼ���Ա�����
        if (0<H(i,j))&&(H(i,j)<=20*10^(-4))
            Hz(i,j)=1;
        elseif H(i,j)<=0
            Hz(i,j)=1;
            H(i,j)=0;
        end
    end
end
havera=c*mean(H(:));
if isnan(havera)
    disp Ĥ��ΪNaN������
end
end
%         lz=(j-1)*L/m;
%         ez(j)=sqrt(lz^2*DDBX(j)^2-2*lz*eZJ(j)*DDBX(j)*cos(fai)+eZJ(j)^2);
