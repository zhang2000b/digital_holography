%% 利用FFT进行相位解包裹
%% image为相位，n为迭代次数，一般为2到3
function [fai]=FFT_unwrapping(image,n)
% linhas=size(image,1);
% colunas=size(image,2);
% 
% image=padarray(image,[abs(linhas-colunas),abs(linhas-colunas)],'replicate','post');%Faz padarray se a imagem nao for quadrada
% if mod(size(image,1),2)~=0
%     image=padarray(image,[1,0],'replicate','post');%Faz padarray se as dimensoes nao forem par
%     %在第一个维度上添加一行
% end
% if mod(size(image,2),2)~=0
%     image=padarray(image,[0,1],'replicate','post');%Faz padarray se as dimensoes nao forem par
%     %在第二个维度上添加一行
% end
[M,N]=size(image);
old_fai=image;
for i=1:n
    
    y=linspace(-M/2,M/2,M);
    x=linspace(-N/2,N/2,N);
    [xx,yy]=meshgrid(x,y);
    step_1=ifft2(fft2((cos(old_fai).*ifft2((xx.^2+yy.^2).*fft2(sin(old_fai)))))./(xx.^2+yy.^2));
    step_2=ifft2(fft2((sin(old_fai).*ifft2((xx.^2+yy.^2).*fft2(cos(old_fai)))))./(xx.^2+yy.^2));
    fai_e=step_1-step_2;
    fai=old_fai+2*pi*round((fai_e-old_fai)/2/pi);
    old_fai=fai;
end
fai=abs(fai);
end