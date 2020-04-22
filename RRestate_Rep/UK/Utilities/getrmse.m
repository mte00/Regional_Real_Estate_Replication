function [ out ] = getrmse( resid,horz )

out=[];
for i=1:length(horz)
    temp=sqrt(mean(resid(1:horz(i),:).^2,1));
    out=[out temp];
end


end

