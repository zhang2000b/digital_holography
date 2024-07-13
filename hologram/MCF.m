%%  MCF相位解包裹
function [output]=MCF(phase)
[M,N]=size(phase);
dx=zeros(M,N);dy=zeros(M,N);
for i=1:M-1
    for j=1:N-1
        dx(i,j)=phase(i+1,j)-phase(i,j);
        dy(i,j)=phase(i,j+1)-phase(i,j);
    end
end
output=zeros(M,N);
for i=1:M-1
    for j=1:N-1
        if(dx(i,j)+dy(i+1,j)-dx(i,j+1)-dy(i,j)==0)
            output(i,j)=sum(sum(dx(1:i,1:j)))+sum(sum(dy(1:i,1:j)));
        end
    end
end
end

