clc,clear;
tic;
filename='../data/model_data_n.sgy';
[Xt,hdr]=readsgy(filename);

dt=2;
[Nx,Nt]=size(Xt);

p=31;%
first=load('../data/model_data_fb.txt');%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xt1=zeros(Nx,Nt);
D1=zeros(Nx,Nt);
U1=zeros(Nx,Nt);
for i=1:Nx
    for j=1:Nt
        if(j-(first(i)-first(1))>0&&j-(first(i)-first(1))<Nt)
            Xt1(i,j-(first(i)-first(1)))=Xt(i,j);
        end
    end
end% flatten

D=medfilt2(Xt1,[p,1]);
U=Xt1-D;

for i=1:Nx
    for j=1:Nt
        if(j+(first(i)-first(1))>0&&j+(first(i)-first(1))<Nt)
            D1(i,j+(first(i)-first(1)))=D(i,j);
            U1(i,j+(first(i)-first(1)))=U(i,j);
        end
    end
end
%% paint
subplot(2,2,1)
imagesc(Xt1')
title('Flatten')
subplot(2,2,2)
imagesc(D1')
title('Down')
subplot(2,2,3)
imagesc(U1')
title('Up')
subplot(2,2,4)
imagesc(Xt')
title('raw')

writesgy(U1,'fu.sgy',hdr);
writesgy(D1,'fd.sgy',hdr);
% writesgy(DD+DU-D,'lefted.sgy',hdr);
toc