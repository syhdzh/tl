function ZXplot(XT,YT,XvT,YvT,XaT,YaT)
%%�����Ĺ켣ͼ
load
scatter(XT,YT);xlabel('���Ĺ켣x����(m)'); ylabel('���Ĺ켣y����(m)')
saveas(gcf, '���Ĺ켣λ������.jpg'); 
scatter(XvT,YvT);xlabel('���Ĺ켣x�ٶ�(m/s)'); ylabel('���Ĺ켣y�ٶ�(m/s)')
saveas(gcf, '���Ĺ켣�ٶ�����.jpg');
scatter(XaT,YaT);xlabel('���Ĺ켣x���ٶ�(m2/s)'); ylabel('���Ĺ켣y���ٶ�(m2/s)')
saveas(gcf, '���Ĺ켣���ٶ�����.jpg');
TTT=(1:2001)/400;
plot(TTT,XT,TTT,YT);xlim([0 5]);xlabel('ʱ��t(s)'); ylabel('���Ĺ켣x,yλ��(m)');legend('xλ��','yλ��')
saveas(gcf, '����x��y����λ��.jpg');
end