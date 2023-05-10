function [PT,HT,WY,Fx,Fy,fai]=DD2(w,VV,e,fai,JHcs,ks,Dc,Rc1,Rc2)
% ���������>��м�϶ʱ��ѹ����������
%% __________________input variables_______________________
% WYz ���������        w    ����ת��
% VV  �ٶ������ٲ���    DDBX ���������
% e   ���ƫ����        fai  ƫλ��         
% ks  ˮ�ۼ��β���ϵ��  JHcs ˮ�ۼ��β���
% Dc   ���ξ���ϵ��
%% __________________output variables_______________________
% PT  ѹ��              HT   Ĥ��
% XT  �������X����λ�� YT   �������Y����λ��
% WY  �ĵ����          Fx��Fy  ��Ĥ��������X��Y�������
% fai �µ�ƫλ��        WYz    �µı�����
% jx��jy  �߽�Ӵ�������ڵ��
%% __________________formal function_______________________
load DATA.mat 
xs=1;PT=[];HT=[];WY=[];PPJ=[];
for dk=1:DZ
    csj=(dk-1)*(2*pi/DZ);ccs=csj:deltsita:csj+(n-1)*deltsita;%������
    Dh=DhCalcu(R,ccs,DZ,JHcs,ks); %���ˮ�ۼ��β���
    Ar=0;Hr=0;ER=zeros(10000,2);ERR3=5.0e-3; k=1; %������ʼֵ
    PK=zeros(n+1,m+1);
    while k>0
        P=PK;%�´ε�����ֵ
        %% Ĥ�����
        [H,Hz,Wy,jd]=Hcalcu(csj,deltsita,fai,e,c,P,ps,Dc,Dh,n,m,E,Rc1,Rc2);
        havera=c*mean(H(:));
        if isnan(havera)
            disp Ĥ��ΪNaN������
            return
        end
        %% ѹ����������������
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
    disp ѹ������ֵΪ��
    PT=zeros(n*DZ,m);
end
%% ����غ�
[Fx,Fy]=F2calcu(deltsita,ps,PT,DZ,R,L,n,m);
if sum(sum(PPJ))~=0
    disp Ĥ��С�ڴֲڶȣ���������״̬
end
end