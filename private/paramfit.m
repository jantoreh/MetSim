% Define fitting function
function val = paramfit(x,p)
    val = p(:,1)+p(:,2).*x.^p(:,3) + p(:,4).*exp(p(:,5).*(x+p(:,6)).^p(:,7));
end

