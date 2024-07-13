%参考光法

clc;

clear;

%读入全息图

Hol=imread('hologram_xidian.bmp');

Hol2=double(Hol (:,:,1)); figure(1);

% imshow(Hol2/255); 

[M,N]=size(Hol2); 

%再现参数
pix=0.00465;                %ccd像素大小
L=N*pix;                    %ccd的大小
Lamda=0.6328e-3;%单位mm
k=2*pi/Lamda;

DDist=1000; %衍射距离起始值
LL=Lamda*N*DDist/L;           %重建平面大小
dz=100; %步长

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

    Hol2F=fftshift(fft2(Hol2)); 

    figure(1); 

    [x1,y1,a1]=impixel(255*abs(Hol2F)/max(max(abs(Hol2F)))); %归一化为255

    close(1);

else

     x1=0; y1=0;

end

ThetaX=asin((x1-M/2)*Lamda/(M*pix));

ThetaY=asin((y1-N/2)*Lamda/(N*pix));

Ref=exp(-1i*2*pi*(pix*x*sin(ThetaX)+pix*y*sin(ThetaY))/Lamda);% 照明光

%-------------全息再现----------------

Hol2=Hol2.*Ref; % 加照明光
T_grad=zeros(1,Repeats);
var_grad=zeros(1,Repeats);
distance=zeros(1,Repeats);
for mm=1:Repeats
   
    Fr=exp(1i*k/2/DDist*(XX.^2+YY.^2));
    Fr_image=fftshift(fft2(Hol2.*Fr));
    phase=exp(1i*k*DDist)/(1i*Lamda*DDist)*exp(1i*k/2/DDist*(xx1.^2+yy1.^2));
    re_image=Fr_image.*phase;
    II=re_image.*conj(re_image);
    Max=max(max(abs(II)));
    Min=min(min(abs(II)));
    figure();imshow(II,[Min Max/1000]);colormap("gray");title("与参考光平行入射重建图");
    distance(mm)=DDist;
    DDist=DDist+dz;

    pause(1);%
    III=II(M/2-200:M/2+200,N/2-200:N/2+200);
    III(III>=Max/100)=255;III(III<=Min)=0;
    III=III*255/max(max(III));
    % figure();imshow(abs(III));
    T_grad(mm)=TenenGrad(III);
    var_grad(mm)=var(III);
end
T_grad=1-T_grad/(max(T_grad));
var_grad=1-var_grad/(max(var_grad));
figure();
plot(distance,T_grad,distance,var_grad)
title("重建距离对重建像清晰度影响")
xlabel("重建距离");ylabel("归一化评价函数");legend("能量梯度函数","方差");











% 
% dx=5;
% 
% dy=dx; 

% DDist=800; %衍射距离起始值
% 
% dz=10.0; %步长
% 
% Repeats=10; %步数
% 
% %设置照明光
% 
% [x,y] = meshgrid(-(M/2-1):M/2,-(N/2-1):N/2);
% 
% choice=menu('Indicate the central points of the +1 diffraction orders','Yes','No');
% 
% if choice==1
% 
%     Hol2F=fftshift(fft2(fftshift(Hol2))); 
% 
%     figure(1); 
% 
%     [x1,y1,a1]=impixel(100*abs(Hol2F)/max(max(abs(Hol2F)))); 
% 
%     close(1);
% 
% else
% 
%      x1=0; y1=0;
% 
% end
% 
% ThetaX=asin((x1-M/2)*Lamda/(M*dx));
% 
% ThetaY=asin((y1-N/2)*Lamda/(N*dy));
% 
% Ref=exp(-1i*2*pi*(dx*x*sin(ThetaX)+dy*y*sin(ThetaY))/Lamda);% 照明光
% 
% %-------------全息再现----------------
% 
% Hol2=Hol2.*Ref; % 加照明光
% 
% for mm=1:Repeats
% 
%     Image=Diffract(Hol2,Lamda,dx,DDist);
% 
%     Image_A=abs(Image)/max(max(abs(Image)));
% 
%     figure(2);
% 
%     imshow(Image_A,[0 1]);
% 
%     DDist=DDist+dz;
% 
%     pause(1);%
% 
% end
% 


% function DiffractZ=Diffract(Object, Lamda, dx, Zd)
% 
% 
% 
%     Object=double(Object);
% 
%     [M,N]=size(Object);
% 
%     dy=dx;
% 
% 
% 
%     du=1/(M*dx);dv=1/(N*dy); %空间频谱面的抽样间隔
% 
%     [u,v]=meshgrid(-du*(M/2-1):du:du*M/2,-dv*(N/2-1):dv:dv*N/2); % 空间频率坐标
% 
%     TransFunct=exp(1i*2.0*pi*Zd*((1/Lamda)^2-u.^2-v.^2).^0.5);% 传递函数
% 
%     Object_F=fftshift(fft2(fftshift(Object))); % 输入物体的傅立叶变换
% 
%     DiffractZ=fftshift(ifft2(fftshift(Object_F.*TransFunct))); % 计算衍射光场
% 
% 
% 
% end



%空间滤波法


%读入全息图
% 
% Hol=imread('test.bmp');
% 
% Hol2=double(Hol (:,:,1)); figure(1);
% 
% imshow(Hol1/255); 
% 
% [M,N]=size(Hol2); 
% 
% %再现参数
% 
% Lamda=0.6328;%单位um
% 
% dx=5.0;
% 
% dy=dx; 
% 
% DDist=300.0e+3; %衍射距离起始值
% 
% dz=10.0e+3; %步长
% 
% Repeats=10; %步数
% 
% %设置照明光
% 
% [x,y] = meshgrid(-(M/2-1):M/2,-(N/2-1):N/2);
% 
% choice=menu('Indicate the central points of the +1 diffraction orders','Yes','No');
% 
% if choice==1
% 
%     Hol2F=fftshift(fft2(fftshift(Hol2))); 
% 
%     figure(1); 
% 
%     [x1,y1,a1]=impixel(100*abs(Hol2F)/max(max(abs(Hol2F)))); 
% 
%     close(1);
% 
% else
% 
%      x1=0; y1=0;
% 
% end
% 
% Hol2F=circshift(Hol2F,[-y1+M/2 -x1+M/2]);
% 
% Hol2F=cyl(x,y,abs(x1-M/2)).*Hol2F;
% 
% Hol2=ifftshift(ifft2(ifftshift(Hol2F)));
% 
% %-------------全息再现----------------
% 
% for mm=1:Repeats
% 
%     Image=Diffract(Hol2,Lamda,dx,DDist);
% 
% Image_A=abs(Image)/max(max(abs(Image)));
% 
% figure(2);
% 
%     imshow(Image_A,[0 1]);
% 
%     DDist=DDist+dz;
% 
%     pause(1);%
% 
% end
% 
% 
% 
% function DiffractZ=Diffract(Object, Lamda, dx, Zd)
% 
% 
% 
    % Object=double(Object);
    % 
    % [M,N]=size(Object);
    % 
    % dy=dx;
% 
% 
% 
%     du=1/(M*dx);dv=1/(N*dy); %空间频谱面的抽样间隔
% 
%     [u,v]=meshgrid(-du*(M/2-1):du:du*M/2,-dv*(N/2-1):dv:dv*N/2); % 空间频率坐标
% 
%     TransFunct=exp(1i*2.0*pi*Zd*((1/Lamda)^2-u.^2-v.^2).^0.5);% 传递函数
% 
%     Object_F=fftshift(fft2(fftshift(Object))); % 输入物体的傅立叶变换
% 
%     DiffractZ=fftshift(ifft2(fftshift(Object_F.*TransFunct))); % 计算衍射光场
% 
% 
% 
% end
% 
% function z = cyl(x,y,d)
% 
% 
% 
%     r = sqrt(x.*x+y.*y);
% 
%     z = zeros(size(r));
% 
%     z(r<d/2)=1.0;
% 
%     z(r==d/2)=0.5;
% end