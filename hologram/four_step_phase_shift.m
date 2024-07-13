clc;clear;
%%   数据初始化(单位为mm）
image=imread("images.png");
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
d=1e-5; n=1.4;            %厚度  折射率
T=1-((n-1)/(n+1))^2;       %透射率
fai=k*(n-1)*d;
% a=rand(N,N)*2*pi;
ang_image=T*new_image.*exp(1i.*fai);
% figure();imshow(ang_image,[]);colormap("gray");title("物平面图像");

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

%%   生成参考光                  相移为0
fai1=0;
AR=mean(mean(abs(FF_image)));         %获取最大振幅
R=AR*exp(1i*fai1); %参考光
E=R+FF_image;
I1=E.*conj(E);                       %全息图
Gmax=max(max(I1));Gmin=min(min(I1));
subplot(221);imshow(I1,[Gmin Gmax/10]);colormap("gray");title("同轴0相移全息图");

%%   生成参考光                  相移为pi/2
fai1=pi/2;
AR=mean(mean(abs(FF_image)));         %获取最大振幅
R=AR*exp(1i*fai1); %参考光
E=R+FF_image;
I2=E.*conj(E);                       %全息图
Gmax=max(max(I2));Gmin=min(min(I2));
subplot(222);imshow(I2,[Gmin Gmax/10]);colormap("gray");title("同轴pi/2相移全息图");
% imwrite(I2,"tongzhou.bmp");
%%   生成参考光                  相移为3pi/2
fai1=3*pi/2;
AR=mean(mean(abs(FF_image)));         %获取平均振幅
R=AR*exp(1i*fai1); %参考光
E=R+FF_image;
I3=E.*conj(E);                       %全息图
Gmax=max(max(I3));Gmin=min(min(I3));
subplot(223);imshow(I3,[Gmin Gmax/10]);colormap("gray");title("同轴3pi/2相移全息图");

%%   生成参考光                  相移为2pi
fai1=pi*2;
AR=mean(mean(abs(FF_image)));         %获取最大振幅
R=AR*exp(1i*fai1); %参考光
E=R+FF_image;
I4=E.*conj(E);                       %全息图
Gmax=max(max(I4));Gmin=min(min(I4));
subplot(224);imshow(I4,[Gmin Gmax/10]);colormap("gray");title("同轴2pi相移全息图");

%%  计算物光场振幅
I1=double(I1);
I2=double(I2);
I3=double(I3);
I4=double(I4);
U=(I4-I2)+1i*(I1-I3)/4/AR;
% U=sqrt((I1-I3).^2+(I4-I2).^2)/4/AR;

%%  1-fft重建
z=1000;             %重建距离
Fr=exp(1i*k/2/z*(XX.^2+YY.^2));
Fr_image=(fft2(U.*Fr));

x1=linspace(-LL/2,LL/2,N);      %再现平面
y1=linspace(-LL/2,LL/2,N);
[xx1,yy1]=meshgrid(x1,y1);
phase=exp(1i*k*z)/(1i*lambda*z)*exp(1i*k/2/z*(xx1.^2+yy1.^2));
re_image=Fr_image.*phase;
II=re_image.*conj(re_image);
Max=max(max(abs(II)));
Min=min(min(abs(II)));
figure();imshow(II,[Min Max/10]);colormap("gray");title("四步相移重建图");

%%  相位
ang=atan((I4-I2)/(I1-I3));
deg=Phase_unwrapping_volkovt(ang);
% figure();mesh(xx1,yy1,deg);