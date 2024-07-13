%% 方差评价函数
function var_value=var(image)
average=mean(mean(image));
var_value=sum(sum((image-average).^2));  
end