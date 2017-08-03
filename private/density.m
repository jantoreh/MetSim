% Define probability density functions
function out = density(type,x,p) 
    % Input: Type, x-values, parameters
    switch type
        case 'w2'
            out = wblpdf(x,p(1),p(2));
        case 'w3'
            out = wblpdf(x-p(3),p(1),p(2));
        case 'n'
            out = normpdf(x,p(1),p(2));
        case 'n2'
            obj = gmdistribution((p(3:4)),permute(p(5:6),[3,2,1]),p(1:2)');
            out = pdf(obj,x');
        case 'ln'
            out = lognpdf(x,p(1),p(2));
        case 'vm'
            out = exp(p(2).*cosd(x-p(1)))./(2*pi*besseli(0,p(2)))*pi/180;
        case 'vm2'
            f1 = p(1)*exp(p(5).*cosd(x-p(3)))./(2*pi*besseli(0,p(5)))*pi/180;
            f2 = p(2)*exp(p(6).*cosd(x-p(4)))./(2*pi*besseli(0,p(6)))*pi/180;
            out = f1+f2;
        case'vm3'
            f1 = p(1)*exp(p(7).*cosd(x-p(4)))./(2*pi*besseli(0,p(7)))*pi/180;
            f2 = p(2)*exp(p(8).*cosd(x-p(5)))./(2*pi*besseli(0,p(8)))*pi/180;
            f3 = p(3)*exp(p(9).*cosd(x-p(6)))./(2*pi*besseli(0,p(9)))*pi/180;
            out = f1+f2+f3;
    end
end