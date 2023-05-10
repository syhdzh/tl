function [H,Hz,Wy,jd]=Hcalcu(csj,deltsita,fai,e,c,P,ps,Dc,Dh,n,m,E,Rc1,Rc2)
% 计算各点处水膜厚度
%% __________________input variables_______________________
% csj 接触区域的序号参数 deltsita 周向节点距离
% c   轴承间隙           P    接触压力
% ps  进水压力无量纲系数 Dh   水槽几何参数
% e   轴承偏心率         fai  偏位角
% Rc1 轴承内粗糙表面参数 Rc2  主轴粗糙表面参数
% n   接触区周向节点数   m    接触区域轴向节点数
% Dc  变形矩阵系数       DDBX 主轴变形量
%% __________________output variables_______________________
% H   动压膜厚           Hz   压力小于0的区域
% WY  衬垫变形           jd   角度变形
%% __________________formal function_______________________
dbstop if error
H=zeros(n+1,m+1);PP=zeros(n+1,m+1);Hz=zeros(n+1,m+1);
Wy=WyCalcu(P,m,n,E,ps,Dc);
KK=zeros(n+1,m+1);
Dh=[Dh;Dh(1)]/c;jd=zeros(n+1,m+1);
for i=1:n+1 %周向
    for j=1:m+1 %轴向
        jd(i,j)=csj+(i-1)*deltsita-fai;
        PP(i,j)=0;
        KK(i,j)=Wy(i,j)/c+Dh(i);
        if KK(i,j)>max(Dh(i))
            KK(i,j)=max(Dh(i));
        end
        H(i,j)=1+e*(cos(jd(i,j)))+Wy(i,j)/c+Rc1(i,j)/c+Rc2(i,j)/c+PP(i,j)+KK(i,j);%初始弹性变形量
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
    disp 膜厚为NaN，出错
end
end
%         lz=(j-1)*L/m;
%         ez(j)=sqrt(lz^2*DDBX(j)^2-2*lz*eZJ(j)*DDBX(j)*cos(fai)+eZJ(j)^2);
