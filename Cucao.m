function Rc=Cucao(Ra,Rendax,Renday,n,m,D,L,ZD)
%%  ����  ��ά�ֲڱ���ĸ߶�ģ��
%% input character 
% Ra   �߶����������,�ֲڶ�           n   x�����������
% Rendax ������ɳ���Ϊ5um             Rendax ������ɳ���Ϊ5um
% m    y�����������                   D   ���ֱ��
% L    ��п��                        ZD  ��������
%% output character 
% Rc һ�����ۼ��������ڵ�����ֲڱ���
%% ��ʼ��������
% Ra=0.5*10^-7;Rendax=5*10^-7;Renday=5*10^-7;n=50;m=100;D=0.236;L=0.2;ZD=14;
Lx=pi*D/ZD; Ly=L; ux=-n/2:n/2;   uy=-m/2:m/2;%X����Y����ڵ�
Ux=2*pi*ux/Lx; Uy=2*pi*uy/Ly;%X��Y������ɢ����
[vx,vy]=meshgrid(ux,uy);
deltax=Lx/n;   deltay=Ly/m;%X��Y����������
Vx=vx*deltax;  Vy=vy*deltay;%X��Y����
%% ����׾�����
P=zeros(n+1,m+1);
for i=1:n+1
    for j=1:m+1        
        P(i,j)=(sum(sum(Ra^2*exp(-(Vx.^2/(Rendax^2)+Vy.^2/(Renday^2))).*exp(1i*(Ux(i)*Vx+Uy(j)*Vy))*deltax*deltay,1),2));%�Ϲ���
    end
end
% figure(1);
% surf(real(P),'Edgecolor','none');title('�׾�����');
%% ���㸴�߶ȷֲ�
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
% title('��ά����ֲڱ���');