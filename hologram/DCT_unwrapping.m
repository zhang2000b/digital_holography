%% 利用DCT进行相位解包裹
%% image为相位，n为迭代次数，一般为2到3
function [fai]=DCT_unwrapping(image,n)
old_fai=image;
for i=1:n
    [M,N]=size(image);
    y=linspace(-M/2,M/2,M);
    x=linspace(-N/2,N/2,N);
    [xx,yy]=meshgrid(x,y);
    step_1=idct2(dct2((cos(old_fai).*idct2((xx.^2+yy.^2).*dct2(sin(old_fai)))))./(xx.^2+yy.^2));
    step_2=idct2(dct2((sin(old_fai).*idct2((xx.^2+yy.^2).*dct2(cos(old_fai)))))./(xx.^2+yy.^2));
    fai_e=step_1-step_2;
    fai=old_fai+2*pi*round((fai_e-old_fai)/2/pi);
    old_fai=fai;
end
% fai=abs(fai);
end