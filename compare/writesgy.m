function writesgy(X,outputfile,header)
[Nx,Nt]=size(X);
fp=fopen(outputfile,'wb');
if nargin<3
    header=uint8(zeros(3600,1));
    theader=uint8(zeros(240,1));
    dt=1000;
    fmt=5;%此处可修改默认格式,1=sun,5=pc
    if fmt==5
        header(3213)=mod(Nx,256);
        header(3214)=round((Nx-mod(Nx,256))/256);
        header(3217)=mod(dt,256);
        header(3218)=round((dt-mod(dt,256))/256);
        header(3219)=header(3217);
        header(3220)=header(3218);
        header(3221)=mod(Nt,256);
        header(3222)=round((Nt-mod(Nt,256))/256);
        header(3223)=header(3221);
        header(3224)=header(3222);
        header(3225)=5;
        header(3255)=1;

        theader(29)=1;
        theader(35)=1;
        theader(115)=mod(Nt,256);
        theader(116)=round((Nt-mod(Nt,256))/256);
        theader(117)=mod(dt,256);
        theader(118)=round((dt-mod(dt,256))/256);
        
        fwrite(fp,header);
        for i=1:Nx
            theader(1:4)=typecast(int32(i),'uint8');
            theader(5:8)=typecast(int32(1),'uint8');
            fwrite(fp,theader);
            fwrite(fp,X(i,:),'float');
        end
    elseif fmt==1
        header(3214)=mod(Nx,256);
        header(3213)=round((Nx-mod(Nx,256))/256);
        header(3218)=mod(dt,256);
        header(3217)=round((dt-mod(dt,256))/256);
        header(3219)=header(3217);
        header(3220)=header(3218);
        header(3222)=mod(Nt,256);
        header(3221)=round((Nt-mod(Nt,256))/256);
        header(3223)=header(3221);
        header(3224)=header(3222);
        header(3226)=1;
        header(3256)=1;

        theader(30)=1;
        theader(36)=1;
        theader(115)=round((Nt-mod(Nt,256))/256);
        theader(116)=mod(Nt,256);
        theader(117)=round((dt-mod(dt,256))/256);
        theader(118)=mod(dt,256);
        
        fwrite(fp,header);
        for i=1:Nx
            theader(1:4)=typecast(swapbytes(int32(i)),'uint8');
            theader(5:8)=typecast(swapbytes(int32(1)),'uint8');
            fwrite(fp,theader);
            % for j=1:Nt
            %     fwrite(fp,float2ibm(X(i,j)),'ubit8');
            % end
            fwrite(fp,swapbytes(ieee2ibm(X(i,:))),'uint32');
        end
    end
else
    fwrite(fp,header(1:3600),'ubit8');
    if header(3225)==5&&header(3226)==0%pc
        for i=1:Nx
            fwrite(fp,header(3601+(i-1)*240:3600+i*240),'ubit8');
            fwrite(fp,X(i,:),'float');
        end
    elseif header(3225)==0&&header(3226)==5
        for i=1:Nx
            fwrite(fp,header(3601+(i-1)*240:3600+i*240),'ubit8');
            fwrite(fp,X(i,:),'float','ieee-be');
        end
    elseif header(3225)+header(3226)==1%sun
        for i=1:Nx
            fwrite(fp,header(3601+(i-1)*240:3600+i*240),'ubit8');
            fwrite(fp,swapbytes(ieee2ibm(X(i,:))),'uint32');
            % for j=1:Nt
            %     fwrite(fp,float2ibm(X(i,j)),'ubit8');
            % end
        end
    end
end
fclose(fp);
fclose all;
return