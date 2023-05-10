function Wy=WyCalcu(P,m,n,E,ps,Dc)
P=P*ps;
Wy=zeros(n+1,m+1);
for i=1:n+1
    for j=1:m+1
        dc=Dc{i,j};
        Wy(i,j)=sum(sum(dc.*P))*2/(pi*E);
    end
end
end


