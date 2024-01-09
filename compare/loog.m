function y=loog(x)%取对数以增强图片图片内容
if(ismatrix(x))
    [row,vol]=size(x);
    y=zeros(row,vol);
    for i=1:row
        for j=1:vol
            if(x(i,j)>=0)
                y(i,j)=log(abs(1+x(i,j)));%这里不取abs的话会莫名其妙出复数，随之出错
            else
                y(i,j)=-log(abs(1-x(i,j)));
            end
        end
    end
end
end