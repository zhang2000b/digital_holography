clc;clear;
image=double(imread("hologram_xidian.bmp"));
[M,N]=size(image);
lambda=0.6328e-3;           %波长
k=2*pi/lambda;              %波矢量
pix=0.00465;                %ccd像素大小
L=N*pix;                    %全息图大小
z0=1000;                    %物体到ccd的距离
L0=lambda*z0*N/L;           %重建平面大小
%%    1-fft重建
x=linspace(-L/2,L/2,N);
y=linspace(-L/2,L/2,N);
[xx,yy]=meshgrid(x,y);
Fr=exp(1i*k/2/z0*(xx.^2+yy.^2));
f=image.*Fr;
Uf=fftshift(fft2(f));       %振幅
X=linspace(-L0/2,L0/2,N);
Y=linspace(-L0/2,L0/2,N);
[XX,YY]=meshgrid(X,Y);
phase=exp(1i*k*z0)/(1i*lambda*z0).*exp(1i*k/2/z0*(XX.^2+YY.^2));
Uf=Uf.*phase;               %物平面光场
Gmax=max(max(abs(Uf)));Gmin=min(min(abs(Uf)));
figure();imshow(abs(Uf),[Gmin Gmax/100]);colormap("gray");title("1-fft重建")
%%    图像截取
UU=zeros(N,N);
p=getrect;
UU(N/2-p(4)/2:N/2+p(4)/2,N/2-p(3)/2:N/2+p(3)/2) = Uf(p(2):p(2)+p(4),p(1):p(1)+p(3)); 
Ns=max(p(3),p(4));
% figure();imshow(abs(UU),[Gmin Gmax/100]);colormap("gray");title('截取图像')

%%    虚拟全息
z1=lambda*Ns/L/L*z0^2;      %虚拟全息图距离1-fft再现像面距离
L1=lambda*z1*N/L0;          %虚拟全息图宽度
x=linspace(-L0/2,L0/2,N);
y=linspace(-L0/2,L0/2,N);
[xx,yy]=meshgrid(x,y);
Fr=exp(-1i*k/2/z1*(xx.^2+yy.^2));
f1=UU.*Fr;                  %振幅
Ut=fftshift(fft2(f1));
X=linspace(-L1/2,L1/2,N);
Y=linspace(-L1/2,L1/2,N);
[XX,YY]=meshgrid(X,Y);
phase=exp(-1i*k*z1)/(-1i*lambda*z1).*exp(-1i*k/2/z1*(XX.^2+YY.^2));
Ut=Ut.*phase;
% Gmax1=max(max(abs(Ut)));Gmin1=min(min(abs(Ut)));
% figure();imshow(abs(Ut),[Gmin1 Gmax1]);colormap("gray");title("虚拟全息图")    %虚拟平面全息图
%%   角谱重建
Us=fft2(Ut);
fx=lambda*linspace(-N/2/L1,N/2/L1,N);
fy=lambda*linspace(-N/2/L1,N/2/L1,N);
[Fx,Fy]=meshgrid(fy,fx);
H=exp(1i*k*z1*sqrt(1-Fx.^2-Fy.^2));
result=ifft2(Us.*H);
I=result.*conj(result);
Gmax2=max(max(abs(I)));Gmin2=min(min(abs(I)));
I=imrotate(I,180);
figure();imshow(abs(I),[Gmin2 Gmax2/10]);colormap("gray");title("VDH4FFT重建");
% 
real_image=real(result);
imag_image=imag(result);
ang=atan2(imag_image,real_image);
figure();imshow(abs(ang));title("未解包裹相位图")
ang=Phase_unwrapping_volkovt(ang);
figure();mesh(ang);colorbar;title("解包裹相位图")

%    第二幅

% 
% image=double(imread("hologram_0.bmp"));
% [M,N]=size(image);
% lambda=0.6328e-3;           %波长
% k=2*pi/lambda;              %波矢量
% pix=0.00465;                %ccd像素大小
% L=N*pix;                    %全息图大小
% z0=1000;                    %物体到ccd的距离
% L0=lambda*z0*N/L;           %重建平面大小
% %%    1-fft重建
% x=linspace(-L/2,L/2,N);
% y=linspace(-L/2,L/2,N);
% [xx,yy]=meshgrid(x,y);
% Fr=exp(1i*k/2/z0*(xx.^2+yy.^2));
% f=image.*Fr;
% Uf=fftshift(fft2(f));       %振幅
% X=linspace(-L0/2,L0/2,N);
% Y=linspace(-L0/2,L0/2,N);
% [XX,YY]=meshgrid(X,Y);
% phase=exp(1i*k*z0)/(1i*lambda*z0).*exp(1i*k/2/z0*(XX.^2+YY.^2));
% Uf=Uf.*phase;               %物平面光场
% Gmax=max(max(abs(Uf)));Gmin=min(min(abs(Uf)));
% % figure();imshow(abs(Uf),[Gmin Gmax/100]);colormap("gray");title("1-fft重建")
% %%    图像截取
% UU=zeros(N,N);
% UU(N/2-p(4)/2:N/2+p(4)/2,N/2-p(3)/2:N/2+p(3)/2) = Uf(p(2):p(2)+p(4),p(1):p(1)+p(3)); 
% Ns=max(p(3),p(4));
% % figure();imshow(abs(UU),[Gmin Gmax/100]);colormap("gray");title('截取图像')
% 
% %%    虚拟全息
% z1=lambda*Ns/L/L*z0^2;      %虚拟全息图距离1-fft再现像面距离
% L1=lambda*z1*N/L0;          %虚拟全息图宽度
% x=linspace(-L0/2,L0/2,N);
% y=linspace(-L0/2,L0/2,N);
% [xx,yy]=meshgrid(x,y);
% Fr=exp(-1i*k/2/z1*(xx.^2+yy.^2));
% f1=UU.*Fr;                  %振幅
% Ut=fftshift(fft2(f1));
% X=linspace(-L1/2,L1/2,N);
% Y=linspace(-L1/2,L1/2,N);
% [XX,YY]=meshgrid(X,Y);
% phase=exp(-1i*k*z1)/(-1i*lambda*z1).*exp(-1i*k/2/z1*(XX.^2+YY.^2));
% Ut=Ut.*phase;
% Gmax1=max(max(abs(Ut)));Gmin1=min(min(abs(Ut)));
% % figure();imshow(abs(Ut),[Gmin1 Gmax1]);colormap("gray");title("虚拟全息图")    %虚拟平面全息图
% %%   角谱重建
% Us=fft2(Ut);
% fx=lambda*linspace(-N/2/L1,N/2/L1,N);
% fy=lambda*linspace(-N/2/L1,N/2/L1,N);
% [Fx,Fy]=meshgrid(fy,fx);
% H=exp(1i*k*z1*sqrt(1-Fx.^2-Fy.^2));
% result1=ifft2(Us.*H);
% I=result1.*conj(result1);
% Gmax2=max(max(abs(I)));Gmin2=min(min(abs(I)));
% I=imrotate(I,180);
% % figure();imshow(abs(I),[Gmin2 Gmax2/10]);colormap("gray");title("VDH4FFT重建");
% 
% 
% 
% 
% %%   相位解包裹
% real_image=real(result);
% imag_image=imag(result);
% ang=atan2(imag_image,real_image);
% angg=Phase_unwrapping_volkovt(ang);
% figure();mesh(angg);colorbar;title("解包裹相位图1");
% 
% 
% real_image1=real(result1);
% imag_image1=imag(result1);
% ang1=atan2(imag_image1,real_image1);
% angg1=Phase_unwrapping_volkovt(ang1);
% figure();mesh(angg1);colorbar;title("解包裹相位图2");
% % ang=ang-ang1;
% % angg=Phase_unwrapping_volkovt(ang);
% % 
% % 
% % figure();mesh(angg);colorbar;title("解包裹相位图");
% 
% height_1=max(max(angg));
% height_2=max(max(angg1));
% height=height_1-height_2;