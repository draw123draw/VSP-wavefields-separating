clc,clear;
tic;
filename='../data/model_data_n.sgy';
[Xt,hdr]=readsgy(filename);Xt=Xt';
[Nt,Nx]=size(Xt);
Xt=Xt';

Lp=100;%Lp is the length of the slope, step satifies (-1/step:step:+1/step), where x=1/step satisfies 2*x*x=Lp
step=sqrt(2/Lp);
Lt=Nt;
Tp=zeros(Lp,Lt);
for i=1:Lp
    p=step*i-1/step;
    for j=1:Lt
        tao=j;
        for k=1:Nx
            t=tao+p*k;
            if(t>1&&t<Nt)
                t=round(t);
                Tp(i,j)=Tp(i,j)+Xt(k,t);
            end
        end
    end
end
Xt0=zeros(Nx,Nt);
Xtu=zeros(Nx,Nt);
Xtd=zeros(Nx,Nt);
% Tp(abs(Tp)<mean(Tp))=0;
Tpu=Tp;
Tpd=Tp;
Tpu(Lp/2-4:Lp,:)=0;
Tpd(1:Lp/2+5,:)=0;
% Tpu(Lp/2:Lp,:)=0;
% Tpd(1:Lp/2,:)=0;
for i=1:Nx
    for j=1:Nt
        for k=1:Lp
            tao=j-i*(step*k-1/step);
            if(tao>1&&tao<Nt)
                tao=round(tao);
                Xtu(i,j)=Xtu(i,j)+Tpu(k,tao);
                Xtd(i,j)=Xtd(i,j)+Tpd(k,tao);
            end
        end
    end
end

%% firstbreak
first=load('../data/model_data_fb.txt');
for i=1:Nx
    for j=1:Nt
        if(j<first(i))
            Xtu(i,j)=0;
            Xtd(i,j)=0;
        end
    end
end

subplot(2,2,1)
imagesc(Xt')
title('raw')
subplot(2,2,2)
temp_t=1:Nt;
temp_p=-Lp/10:Lp/10;
imagesc(temp_p,temp_t,loog(Tp'))
title('Tp')
subplot(2,2,3)
imagesc(Xtu')
title('Up')
subplot(2,2,4)
imagesc(Xtd')
title('Down')

writesgy(Xtu,'fu.sgy',hdr);
writesgy(Xtd,'fd.sgy',hdr);
toc
writesgy(Xtd+Xtu-Xt,'lefted.sgy',hdr);