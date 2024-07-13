%% 角谱衍射法
clc;clear;
image=double(imread("images.png"));
N=512;                      %像素大小
z0=1000;                    %物体到ccd的距离
pix=0.00465;                %ccd像素大小
L=N*pix;                    %ccd的大小
lambda=0.632e-3;            %波长
k=2*pi/lambda;              %波矢量
[M0,N0]=size(image);
X1=imresize(image,N/N0);
[M1,N1]=size(X1);
X=zeros(N,N);
X(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2)=X1;
L0=sqrt(lambda*z0*N);
Uf=fftshift(fft2(X));
x=lambda*linspace(-N/L0/2,N/L0/2,N);x1=linspace(-L0/2,L0/2,N);
y=lambda*linspace(-N/L0/2,N/L0/2,N);y1=linspace(-L0/2,L0/2,N);
[xx,yy]=meshgrid(x,y);[xx1,yy1]=meshgrid(x1,y1);
count=input("1-基尔霍夫传递函数，2-角谱传递函数，3-菲涅尔传递函数,4-瑞丽传递函数，5-菲涅尔傅里叶形式:");
switch count
    case 1
        H1=exp(1i*k*sqrt(z0^2+xx1.^2+yy1.^2)).*(sqrt(z0^2+xx1.^2+yy1.^2)+z0)...
            ./(1i*2*lambda.*(z0^2+xx1.^2+yy1.^2));
        H1=fftshift(fft2(H1));
    case 2
        H1=exp(1i*k*z0*sqrt(1-xx.^2-yy.^2));
    case 3
        H1=exp(1i*k*z0.*(1-lambda^2/2.*(xx.^2+yy.^2)));
    case 4
        H1=z0*exp(1i*k.*sqrt(z0^2+xx1.^2+yy1.^2))./(1i*lambda.*sqrt(z0^2+xx1.^2+yy1.^2));
        H1=fftshift(fft2(H1));
    case 5
        H1=exp(1i*k*z0).*exp(1i*k/2/z0.*(xx1.^2+yy1.^2))/(1i*lambda*z0);
        H1=fftshift(fft2(H1));
end
UU=ifft2(Uf.*H1);

%%   生成参考光
X=linspace(-L/2,L/2,N);             %ccd平面
Y=linspace(-L/2,L/2,N);
[XX,YY]=meshgrid(X,Y);
ang=pi/20;                          %离轴角度有一定限制，ang越大说明离轴角度越小，ang越小，离轴角度越大。
Qx=cos(ang);
Qy=Qx;
AR=max(max(abs(UU)));         %获取最大振幅
R=AR*exp(1i*k*(XX.*Qx+YY.*Qy)); %参考光  
E=R+UU;
II=E.*conj(E);                       %全息图
II=uint8(255*II./(max(max(II))));
% imwrite(II,"hologram.bmp")
figure();imshow(II,[]);colormap("gray");title("全息图");
%%   再现
% Er=exp(1i*k*(XX.*Qx+YY.*Qy));
Er=1;
z=z0;
U1=fftshift(fft2(Er.*double(II)));
x=lambda*linspace(-N/L/2,N/L/2,N);x1=linspace(-L/2,L/2,N);
y=lambda*linspace(-N/L/2,N/L/2,N);y1=linspace(-L/2,L/2,N);
[xx,yy]=meshgrid(x,y);
switch count
    case 1
        H2=exp(1i*k*sqrt(z^2+xx1.^2+yy1.^2)).*(sqrt(z^2+xx1.^2+yy1.^2)+z)...
            ./(1i*2*lambda.*(z^2+xx1.^2+yy1.^2));
        H2=fftshift(fft2(H2));
    case 2
        H2=exp(1i*k*z*sqrt(1-xx.^2-yy.^2));
    case 3
        H2=exp(1i*k*z.*(1-lambda^2/2.*(xx.^2+yy.^2)));
    case 4
        H2=z*exp(1i*k.*sqrt(z^2+xx1.^2+yy1.^2))./(1i*lambda.*sqrt(z^2+xx1.^2+yy1.^2));
        H2=fftshift(fft2(H2));
    case 5
        H2=exp(1i*k*z).*exp(1i*k/2/z.*(xx1.^2+yy1.^2))/(1i*lambda*z);
        H2=fftshift(fft2(H2));
end
Uv=ifft2(U1.*H2);
III=Uv.*conj(Uv);
Gmax=max(max(III));
Gmin=min(min(III));
figure();imshow(III,[Gmin Gmax/2]);colormap("gray");title("角谱再现");