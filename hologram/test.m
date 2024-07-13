clc;clear
R=10; 
%Set grid number
n=100; 
theta = (-n:2:n)/n*pi;
phi = ([0,0:2:n])'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(end) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(end) = 0;
x = R*cosphi*cos(theta);
y = R*cosphi*sintheta;
z = R*sin(phi)*ones(1,n+1);
figure();surf(x,y,z);
title('半球');

fai=exp(1i*z);
re=real(fai);
im=imag(fai);
ang=atan2(im,re);

image=ang;

linhas=size(image,1);
colunas=size(image,2);

image=padarray(image,[abs(linhas-colunas),abs(linhas-colunas)],'replicate','post');%Faz padarray se a imagem nao for quadrada
if mod(size(image,1),2)~=0
    image=padarray(image,[1,0],'replicate','post');%Faz padarray se as dimensoes nao forem par
    %在第一个维度上添加一行
end
if mod(size(image,2),2)~=0
    image=padarray(image,[0,1],'replicate','post');%Faz padarray se as dimensoes nao forem par
    %在第二个维度上添加一行
end
[M,N]=size(image);
old_fai=image;
for i=1:2
    
    y=linspace(-M/2,M/2,M);
    x=linspace(-N/2,N/2,N);
    [xx,yy]=meshgrid(x,y);
    step_1=ifft2(fft2((cos(old_fai).*ifft2((xx.^2+yy.^2).*fft2(sin(old_fai)))))./(xx.^2+yy.^2));
    step_2=ifft2(fft2((sin(old_fai).*ifft2((xx.^2+yy.^2).*fft2(cos(old_fai)))))./(xx.^2+yy.^2));
    fai_e=step_1-step_2;
    fai=old_fai+2*pi*round((fai_e-old_fai)/2/pi);
    old_fai=fai;
end
% fai=abs(fai);

figure();surf(x,y,fai);title("解包裹半球")