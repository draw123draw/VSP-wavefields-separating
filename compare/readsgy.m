function [Xt,header]=readsgy(inputfile)

filesize=dir(inputfile);
fp=fopen(inputfile,'rb');
if fp==-1
    disp("ERROR: file not found");
end
X=fread(fp,3600);
if X(3225)+X(3226)==5%pc
    if X(3226)==0
        Nt=X(3221)+X(3222)*256;
%         sprinf('big endian,小字节在前\n');
        pcs='ieee-le';
    else
        Nt=X(3221)*256+X(3222);
%         sprintf('little endian,大字节在前\n');
        pcs='ieee-be';
    end
    Nx=(filesize.bytes-3600)/(240+4*Nt);
    header=zeros(240*Nx+3600,1);
    header(1:3600)=X;
    Xt=zeros(Nx,Nt);
    for i=1:Nx
        header(3600+(i-1)*240+1:3600+i*240)=fread(fp,240);
        Xt(i,:)=fread(fp,Nt,'float',pcs);
    end
elseif X(3225)+X(3226)==1 || X(3225)+X(3226)==2%sun
    if X(3225)==0
        Nt=X(3221)*256+X(3222);
        pcs='ieee-be';
%         sprintf('little endian,大字节在前\n');
    else
        Nt=X(3221)+X(3222)*256;
        pcs='ieee-le';
%         sprintf('big endian,小字节在前\n');
    end
    Nx=(filesize.bytes-3600)/(240+4*Nt);
    header=zeros(240*Nx+3600,1);
    header(1:3600)=X;
    Xt=zeros(Nx,Nt,'uint32');
    for i=1:Nx
        header(3600+(i-1)*240+1:3600+i*240)=fread(fp,240);
        Xt(i,:)=fread(fp,Nt,'uint32',pcs);
    end
    if X(3225)+X(3226)==1
        Xt=u32_ibm_2double(Xt);
    else
        tmp=single(Xt);
        for i=1:Nx
            for j=1:Nt
                if Xt(i,j)>2^31
                    tmp(i,j)=-single(4294967296-Xt(i,j));
                end
            end
        end
        Xt=tmp;
        clear tmp;
    end
else
    disp(["format unknown,byte(3225,3226)=",num2str(X(3225)),num2str(X(3226))]);
end
fclose(fp);
return
end
function data2=u32_ibm_2double(data)%摘抄自S4M文件下的Geophysics_3.0文件夹其中的ibm2double
data2=(1-2*double(bitget(data,32))).*16.^ ...
  (double(bitshift(bitand(data,uint32(hex2dec('7f000000'))),-24))-64) .* ...
  (double(bitand(data,uint32(hex2dec('00ffffff'))))/2^24);
end