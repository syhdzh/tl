function Rc=Cucao(Ra,Rendax,Renday,n,m,D,L,ZD)
%%  计算  三维粗糙表面的高度模型
%% input character 
% Ra   高度起伏均方根,粗糙度           n   x方向采样点数
% Rendax 横向相干长度为5um             Rendax 纵向相干长度为5um
% m    y方向采样点数                   D   轴承直径
% L    轴承宽度                        ZD  沟槽数量
%% output character 
% Rc 一个沟槽计算区域内的随机粗糙表面
%% 初始参数设置
% Ra=0.5*10^-7;Rendax=5*10^-7;Renday=5*10^-7;n=50;m=100;D=0.236;L=0.2;ZD=14;
Lx=pi*D/ZD; Ly=L; ux=-n/2:n/2;   uy=-m/2:m/2;%X方向、Y方向节点
Ux=2*pi*ux/Lx; Uy=2*pi*uy/Ly;%X、Y方向离散波数
[vx,vy]=meshgrid(ux,uy);
deltax=Lx/n;   deltay=Ly/m;%X、Y方向采样间隔
Vx=vx*deltax;  Vy=vy*deltay;%X、Y方向
%% 计算孔径函数
P=zeros(n+1,m+1);
for i=1:n+1
    for j=1:m+1        
        P(i,j)=(sum(sum(Ra^2*exp(-(Vx.^2/(Rendax^2)+Vy.^2/(Renday^2))).*exp(1i*(Ux(i)*Vx+Uy(j)*Vy))*deltax*deltay,1),2));%上光所
    end
end
% figure(1);
% surf(real(P),'Edgecolor','none');title('孔径函数');
%% 计算复高度分布
X0=ux*deltax;Y0=uy*deltay; [Ux0,Uy0]=meshgrid(Uy,Ux);
yita=(randn(n+1,m+1)+1i*randn(n+1,m+1))/sqrt(2);
hc=zeros(n+1,m+1);
for i=1:n+1
    for j=1:m+1
        k1=exp(-1i*(Ux0*X0(i)+Uy0*Y0(j)));
        k2=sqrt(P).*yita;K2=k2.*k1;
        hc(i,j)=sqrt(2)*pi*sum(sum(K2))/sqrt(Lx*Ly)/(sqrt(2*pi));
    end
end
Rc=real(hc);
end
% figure(2);
% surf(real(hc),'Edgecolor','none');
% colormap jet
% colorbar
% view(2)
% title('二维随机粗糙表面');