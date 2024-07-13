%Hemisphere
clear,clc,close all
%Sphere radius
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
figure();mesh(ang);
% % agg=Phase_unwrapping_volkovt(ang);
% agg=MCF(ang);
% figure();surf(x,y,agg);title("解包裹半球")


% image=zeros(200);
% image(10:190,10:190)=10;
% [M,N]=size(image);
% x=linspace(-M/2,M/2,M);
% y=linspace(-N/2,N/2,N);
% [xx,yy]=meshgrid(x,y);
% mesh(x,y,image)
% 
% z=image;
% fai=exp(1i*z);
% re=real(fai);
% im=imag(fai);
% ang=atan2(im,re);
% % agg=Phase_unwrapping_volkovt(ang);
% agg=FFT_unwrapping(ang,2);
% figure();surf(x,y,agg);title("解包裹半球")