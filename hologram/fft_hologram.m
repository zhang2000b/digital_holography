clc;clear;
%%   数据初始化(单位为mm）
image=(imread("images.png"));
[M0,N0]=size(image);
n=min(M0,N0);
N=1024;                         %形成全息图的取样数
image=imresize(image,N/n/4);
[M1,N1]=size(image);
new_image=zeros(N,N);
new_image(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2)=image; %图像的扩充
new_image=double(new_image);
lambda=0.6328e-3;           %波长
k=2*pi/lambda;              %波矢量
pix=0.00465;                %ccd像素大小
L=N*pix;                    %ccd的大小
z0=1000;                    %衍射距离
L0=40;           %物平面的大小
LL=lambda*N*z0/L;           %重建平面大小
%%    增加相位
d=1e-3; n=1.4;            %厚度  折射率
T=1-((n-1)/(n+1))^2;       %透射率
fai=k*(n-1)*d;
% fai=rand(N,N)*2*pi;
ang_image=T*new_image.*exp(1i.*fai);
% ang_image=new_image;
figure();imshow(ang_image,[]);colormap("gray");title("物平面图像");

%%    衍射计算
x=linspace(-L0/2,L0/2,N);
y=linspace(-L0/2,L0/2,N);

[xx,yy]=meshgrid(x,y);              %物平面
F=exp(1i*k/2/z0*(xx.^2+yy.^2));
F_image=fftshift(fft2(ang_image.*F));

X=linspace(-L/2,L/2,N);             %ccd平面
Y=linspace(-L/2,L/2,N);
[XX,YY]=meshgrid(X,Y);
phase=exp(1i*k*z0)/(1i*lambda*z0)*exp(1i*k/2/z0*(XX.^2+YY.^2));
FF_image=F_image.*phase;            %到达ccd平面的物光振幅

%%   生成参考光
ang=pi/20;                          %离轴角度有一定限制，ang越大说明离轴角度越小，ang越小，离轴角度越大。
Qx=cos(ang);
Qy=Qx;
AR=max(max(abs(FF_image)));         %获取最大振幅
R=AR*exp(1i*k*(XX.*Qx+YY.*Qy)); %参考光
E=R+FF_image;
I=E.*conj(E);                       %全息图
I=uint8(255*I./(max(max(I))));
imwrite(I,"hologram_xidian.bmp")
% figure();imshow(I,[]);colormap("gray");title("全息图");

%%  全息图再现
% EE=exp(1i*k*(XX.*Qx+YY.*Qy));   %再现光
EE=1;
z=1000;             %重建距离
Fr=exp(1i*k/2/z*(XX.^2+YY.^2));
Fr_image=fftshift(fft2(double(I).*EE.*Fr));

x1=linspace(-LL/2,LL/2,N);      %再现平面
y1=linspace(-LL/2,LL/2,N);
[xx1,yy1]=meshgrid(x1,y1);
phase=exp(1i*k*z)/(1i*lambda*z)*exp(1i*k/2/z*(xx1.^2+yy1.^2));
re_image=Fr_image.*phase;
II=re_image.*conj(re_image);
Max=max(max(abs(II)));
Min=min(min(abs(II)));
figure();imshow(II,[Min Max/1000]);colormap("gray");title("菲涅尔重建");




