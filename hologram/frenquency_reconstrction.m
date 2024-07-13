clc;

clear;

%读入全息图

Hol2=imread('hologram_xidian.bmp');

% Hol2=double(Hol (:,:,1)); figure(1);s

% imshow(Hol2/255); 
Hol2=double(Hol2);
[M,N]=size(Hol2); 

%再现参数
pix=0.00465;                %ccd像素大小
L=N*pix;                    %ccd的大小
lambda=0.6328e-3;%单位mm
k=2*pi/lambda;

DDist=0.05; %衍射距离起始值
LL=sqrt(lambda*DDist*N);           %重建平面大小
dz=0.05; %步长

Repeats=10; %步数

%设置照明光
X=linspace(-L/2,L/2,N);             %ccd平面
Y=linspace(-L/2,L/2,N);
[XX,YY]=meshgrid(X,Y);
[x,y] = meshgrid(-(M/2-1):M/2,-(N/2-1):N/2);

x1=linspace(-LL/2,LL/2,N);      %再现平面
y1=linspace(-LL/2,LL/2,N);
[xx1,yy1]=meshgrid(x1,y1);

choice=menu('Indicate the central points of the +1 diffraction orders','Yes','No');

if choice==1

    Hol2F=fftshift(fft2((Hol2))); 

    figure(1); 

    [x1,y1,a1]=impixel(255*abs(Hol2F)/max(max(abs(Hol2F)))); %归一化为255

    close(1);

else

     x1=0; y1=0;

end

ThetaX=asin((x1-M/2)*lambda/(M*pix));

ThetaY=asin((y1-N/2)*lambda/(N*pix));

Ref=exp(-1i*2*pi*(pix*x*sin(ThetaX)+pix*y*sin(ThetaY))/lambda);% 照明光

%-------------全息再现----------------

Hol2=Hol2.*Ref; % 加照明光
T_grad=zeros(1,Repeats);
var_grad=zeros(1,Repeats);
distance=zeros(1,Repeats);

for mm=1:Repeats

    Image=Diffract(Hol2,lambda,pix,DDist);

    Image_A=abs(Image)/max(max(abs(Image)));

    Max=max(max(abs(Image_A)));
    Min=min(min(abs(Image_A)));
    figure();imshow(Image_A,[Min Max]);colormap("gray");title("角谱重建图");
    distance(mm)=DDist;
    DDist=DDist+dz;

    pause(1);%
    T_grad(mm)=TenenGrad(Image_A);
    var_grad(mm)=var(Image_A);
end
T_grad=1-T_grad/(max(T_grad));
var_grad=1-var_grad/(max(var_grad));
figure();
plot(distance,T_grad,distance,var_grad)
title("重建距离对重建像清晰度影响")
xlabel("重建距离");ylabel("归一化评价函数");legend("能量梯度函数","方差");


function DiffractZ=Diffract(Object, Lamda, dx, Zd)



    Object=double(Object);

    [M,N]=size(Object);

    dy=dx;



    du=1/(M*dx);dv=1/(N*dy); %空间频谱面的抽样间隔

    [u,v]=meshgrid(-du*(M/2-1):du:du*M/2,-dv*(N/2-1):dv:dv*N/2); % 空间频率坐标

    TransFunct=exp(1i*2.0*pi*Zd*((1/Lamda)^2-u.^2-v.^2).^0.5);% 传递函数

    Object_F=fftshift(fft2(fftshift(Object))); % 输入物体的傅立叶变换

    DiffractZ=fftshift(ifft2(fftshift(Object_F.*TransFunct))); % 计算衍射光场



end


