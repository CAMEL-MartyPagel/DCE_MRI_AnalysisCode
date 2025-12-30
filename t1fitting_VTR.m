function [M0, T1map,c] = t1fitting_VTR(S, TRs)
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'...
        ,'Display','off','MaxFunEvals',100000,'MaxIter',100000,'TolFun',...
        1e-8,'TolX',1e-8);
    lb = [];
    ub = [];
    modelfun = @(p,x) (p(1).*(1-exp(-x./p(2)))+p(3));
    x0 = [400 1000 0];
    p = lsqcurvefit(modelfun,x0,TRs,S,lb,ub,options);
    M0 = p(1);
    T1map = p(2);
    c=p(3);
    
    %plot
%     times = linspace(TRs(1),TRs(end));
%     figure;
%     plot(TRs,S,'ko',times,modelfun(p,times),'b-');
end

function [M0, T1map,c] = t1fitting_VFA(S, TR, FAs)
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'...
        ,'Display','off','MaxFunEvals',100000,'MaxIter',100000,'TolFun',...
        1e-8,'TolX',1e-8);
    lb = [];
    ub = [];
    modelfun = @(p,x) (p(1).*((1-exp(-TR./p(2)))/(1-cos(x)*exp(-TR./p(2))))+p(3));
    x0 = [400 1000 0];
    p = lsqcurvefit(modelfun,x0,TRs,S,lb,ub,options);
    M0 = p(1);
    T1map = p(2);
    c=p(3);
    
    %plot
%     times = linspace(TRs(1),TRs(end));
%     figure;
%     plot(TRs,S,'ko',times,modelfun(p,times),'b-');
end
