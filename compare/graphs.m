clear,clc
% A=readsgy('../compare/fu.sgy');
A=readsgy('../downs/model_data_down2.sgy');
A=A';
%% sgy graph
figure
% A=loog(A);%AGC by 
m=max(abs(A(:)))/1.5;%AGC
imagesc(A,[-m,m])
mycolor=load('color.csv');
colormap(mycolor)
axis off
set(gca,'position',[0,0,1,1]);
set(gcf,'unit','centimeters','position',[10,10,4.8,6]);


%% Dip graph
% dip=load('../pos.csv');dip=dip';
% angle_r=80;
% 
% figure;
% imagesc(dip,[-angle_r,angle_r]);colormap turbo;axis off;
% set(gcf,'unit','centimeters','position',[10,10,4.8,6]);set(gca,'position',[0,0,1,1]);

% figure;
% imagesc(ddip*90/pi,[-angle_r,angle_r]);colormap turbo;axis off;
% set(gcf,'unit','centimeters','position',[20,20,6.88,8.62]);set(gca,'position',[0,0,1,1]);