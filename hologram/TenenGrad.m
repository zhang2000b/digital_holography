%% 能量梯度函数
function f=TenenGrad(image)   % 值越大，图像越锐利
[m,n]=size(image);
t=0;
for i=1:m-1
    for j=1:n-1
        grad=(image(i+1,j)-image(i,j))^2+(image(i,j+1)-image(i,j))^2+t;
        t=grad;
    end
end
f=t;
end