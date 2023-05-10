function [PT,HT,WY,Fx,Fy,fai]=DD2(w,VV,e,fai,JHcs,ks,Dc,Rc1,Rc2)
% 主轴变形量>轴承间隙时，压力迭代计算
%% __________________input variables_______________________
% WYz 主轴变形量        w    主轴转速
% VV  速度无量纲参数    DDBX 主轴变形量
% e   轴承偏心率        fai  偏位角         
% ks  水槽几何参数系数  JHcs 水槽几何参数
% Dc   变形矩阵系数
%% __________________output variables_______________________
% PT  压力              HT   膜厚
% XT  轴承中心X方向位移 YT   轴承中心Y方向位移
% WY  衬垫变形          Fx、Fy  油膜所产生的X、Y方向的力
% fai 新的偏位角        WYz    新的变形量
% jx、jy  边界接触的区域节点号
%% __________________formal function_______________________
load DATA.mat 
xs=1;PT=[];HT=[];WY=[];PPJ=[];
for dk=1:DZ
    csj=(dk-1)*(2*pi/DZ);ccs=csj:deltsita:csj+(n-1)*deltsita;%参数角
    Dh=DhCalcu(R,ccs,DZ,JHcs,ks); %求解水槽几何参数
    Ar=0;Hr=0;ER=zeros(10000,2);ERR3=5.0e-3; k=1; %迭代初始值
    PK=zeros(n+1,m+1);
    while k>0
        P=PK;%下次迭代赋值
        %% 膜厚求解
        [H,Hz,Wy,jd]=Hcalcu(csj,deltsita,fai,e,c,P,ps,Dc,Dh,n,m,E,Rc1,Rc2);
        havera=c*mean(H(:));
        if isnan(havera)
            disp 膜厚为NaN，出错
            return
        end
        %% 压力求解与迭代误差计算
        PK2=P2calcu(H,Hz,deltsita,deltL,P,n,m,VV,c,e,w,R,jd,ps);
        PK=PK+(PK2-PK)*xs;
        [sumt,~,sum2]=SumErr(PK,P,n,m);
        if (sum2<10^(-5))&&(k>10)
            break;
        end
        if sumt<=ERR3
            break;
        end
        if k>10000
            break;
        end
        ER(k,1)=sumt;ER(k,2)=sum2;k=k+1;
    end
    PPK=PK(1:n,1:m);HK=H(1:n,1:m);
    PT=[PT;PPK];HT=[HT;HK];WY=[WY;Wy];AR(dk)=Ar;HR(dk)=Hr;
end
if norm(PT)==0
    disp 压力计算值为空
    PT=zeros(n*DZ,m);
end
%% 求解载荷
[Fx,Fy]=F2calcu(deltsita,ps,PT,DZ,R,L,n,m);
if sum(sum(PPJ))~=0
    disp 膜厚小于粗糙度，进入混合润滑状态
end
end