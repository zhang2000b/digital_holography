%%  平方梯度函数
function value=squared_gradient(image)
[m,n]=size(image);
count=0;
for i=1:m
    for j=1:n-1
        value=count+sum(sum((image(i,j+1)-image(i,j))^2));
    end
end

end
 