clc,clear;
filename='../data/field_data.sgy';
[D,hdr]=readsgy(filename);
[Nx,Nt]=size(D);
%% 
Fk=fft2(D,Nx,Nt);
FkD=Fk;
FkU=Fk;

att=0.03;% Attenuation factor
%% either:deivide by quadrant
for i=1:Nx
    for j=1:Nt
        if (i>Nx/2&&j>Nt/2)||(i<Nx/2&&j<Nt/2)
            FkD(i,j)=FkD(i,j)*att;
        % elseif (i<Nx/2&&j>Nt/2)||(i>Nx/2&&j<Nt/2)
        else
            FkU(i,j)=FkD(i,j)*att;
        end
    end
end
%% either:manual
% maskU=mkpolygon(loog(abs(FkU)),2);%draw the polygon
% maskD=mkpolygon(loog(abs(FkD)),2);
% FkU=FkU.*maskU;
% FkD=FkD.*maskD;
%%
DD=ifft2(FkD);
DU=ifft2(FkU);

DD=real(DD);
DU=real(DU);
%% paint

subplot(2,2,1)
imagesc(D')
title('raw')
subplot(2,2,2)
imagesc(loog(abs(Fk')))
title('f-k')
subplot(2,2,3)
imagesc(DU')
title('Up')
subplot(2,2,4)
imagesc(DD')
title('Down')

writesgy(DU,'fu.sgy',hdr);
writesgy(DD,'fd.sgy',hdr);
writesgy(DD+DU-D,'lefted.sgy',hdr);

function [mask,points]=mkpolygon(X,n)
% press 'enter' to end the drawing after one polygon is obtained.
if nargin==1, n=1; end
[x,y]=size(X);
fig1=figure('Name','Press enter to end drawing'); 
imagesc(X);
title('To stop the drawing, press Enter');
hold on;

mask=zeros(x,y,n,'logical');
for i=1:n
    points=[];
    pre_pt=[];
    while 1 
        pt = ginput(1);
        if isequal(pt,[]), break, end
        scatter(pt(1),pt(2),'ro'); % 标记点
        if ~isempty(pre_pt)
            line([pre_pt(1),pt(1)],[pre_pt(2),pt(2)]);%标记线
        end
        points=[points;pt];
        pre_pt=pt;
    end
    mask(:,:,i)=poly2mask(points(:,1),points(:,2),size(X,1),size(X,2));
end
close(fig1)

mask=any(mask,3);
mask=squeeze(mask);
% points(end+1,:)=points(1,:);
% imagesc(X);hold on;line(points(:,1),points(:,2));
end
