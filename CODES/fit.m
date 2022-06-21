function [Plg,ygrid]=fit(Xn,Yn)
    if size(Xn,1)==1;Xn=Xn';end
    if size(Yn,1)==1;Yn=Yn';end

    modelFun =  @(p,x) p(1).*(1-(0.5*exp(-(x-Xn(1))./p(2))+0.5*(exp(-(x-Xn(1))./p(3)))))+ Yn(1);
    startingVals = [(max(Yn)-Yn(1)) 0.1 0.2];

%     modelFun =  @(p,x) p(1)+(1-(1./(1+((x/p(2)).^p(3)))));
%     startingVals = [(max(Yn)-Yn(1)) 0.1 3];

    RNCF = @(p) norm(Yn - modelFun(p,Xn));
    beta = patternsearch(RNCF,startingVals);
    Plg= nlinfit(Xn,Yn,modelFun,beta);
%     [Plg,~,~,CovB] = nlinfit(Xn,Yn,modelFun,beta);
%     error=sqrt(diag(CovB));
%     error= error(1);
    ygrid=(modelFun(Plg,Xn));
end
