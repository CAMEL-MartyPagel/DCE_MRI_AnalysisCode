function pars = fitdcemri(varargin)
%  DCE_MRI_FIT performs fitting to estimate DCE MRI parameters for the following models:
%
%    1) Linear Reference Region Model
%           pars = FITDCEMRI(toi,rr,time,'lsq')
%
%    2) Non-Linear Reference Region Model
%           pars = FITDCEMRI(toi,rr,time,x0,lb,ub,'RRM')
%
%    3) Linear Tofts Model
%           pars = FITDCEMRI(Ctoi,Cp,time)
%
%    4) Non-linear Tofts Model
%           pars = FITDCEMRI(Ctoi,Cp,time,x0,lb,ub,'Tofts')
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Linear Reference Region Model (LRRM):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   pars = fitdcemri(toi,rr,time,'fit')
%    /::INPUTS::\
%       toi: a column vector of either delta R1 values or Contrast Agent(CA)
%           Concentration over time for the Tissue of Interest
%       rr: a column vector of either delta R1 values or CA Concentration
%           over time for the Reference Region
%       time: a column vector of time in minutes corresponding to toi and
%           rr
%       'fit': a string designating the desired fitting algorithm to be
%           used. Can be:
%                   'robust': uses robustfit
%                   'lsq': uses linear least squared fitting (mldivide)
%                   'nonneg': uses linear least squared fitting with a
%                       non-negative constraint
%                   'lasso': uses lasso
%    /::OUTPUT::\
%       pars: a vector with the following components
%           pars(1)=  Relative Ktrans (RKtrans)
%           pars(2)=  kep for the tissue of interest, kepTOI
%           pars(3)=  kep for the reference region,   kepRR
%           pars(4)=  rsquare for the fitting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Non-Linear Reference Region Model (NRRM):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   pars = fitdcemri(toi,rr,time,x0,lb,ub,'RRM')
%    /::INPUTS::\
%       toi: a column vector of either delta R1 values or Contrast Agent(CA)
%           Concentration over time for the Tissue of Interest
%       rr: a column vector of either delta R1 values or CA Concentration
%           over time for the Reference Region
%       time: a column vector of time in minutes corresponding to toi and
%           rr
%       x0: a vector of best guesses for parameters from the non-linear 
%           least squared fitting (Ex:[RKtrans_guess,kepTOI_guess,kepRR_guess])
%       lb: vector of lowerbounds for the parameters (same form as above)
%       ub: vector of upperbounds for the parameters (same form as above)
%       'RRM': This string indicates that you want the non-linear Reference
%           Region Model (versus the non-linear Tofts Model)
%    /::OUTPUT::\
%       pars: a vector with the following components
%           pars(1)=  Relative Ktrans (RKtrans)
%           pars(2)=  kep for the tissue of interest, kepTOI
%           pars(3)=  kep for the reference region,   kepRR
%           pars(4)=  rsquare for the fitting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Linear Toft's Model (LTM):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   pars = fitdcemri(Ctoi,Cp,time)
%    /::INPUTS::\
%       Ctoi: a column vector of Contrast Agent(CA) Concentration over time 
%           for the Tissue of Interest.
%           NOTE: Must be concentration, not delta R1 (deltaR1 = r1*Conc
%           where r1 is the relaxivity of the CA)
%       Cp: an Arterial Input Function (AIF) for the tissue in question in
%           the form of a column vector equal in size to toi
%       time: a column vector of time in minutes corresponding to toi and
%           Cp
%    /::OUTPUT::\
%       pars: a vector with the following components
%           pars(1)=  Ktrans
%           pars(2)=  kep for the tissue of interest, kepTOI
%           pars(3)=  rsquare for the fitting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Non-Linear Toft's Model (NLTM):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   pars = fitdcemri(Ctoi,Cp,time,x0,lb,ub,'Tofts')
%    /::INPUTS::\
%       Ctoi: a column vector of Contrast Agent(CA) Concentration over time 
%           for the Tissue of Interest.
%           NOTE: Must be concentration, not delta R1 (deltaR1 = r1*Conc
%           where r1 is the relaxivity of the CA)
%       Cp: an Arterial Input Function (AIF) for the tissue in question in
%           the form of a column vector equal in size to toi
%       time: a column vector of time in minutes corresponding to toi and
%           Cp
%       x0: a vector of best guesses for parameters from the non-linear 
%           least squared fitting (Ex:[Ktrans_guess,kepTOI_guess])
%       lb: vector of lowerbounds for the parameters (same form as above)
%       ub: vector of upperbounds for the parameters (same form as above)
%       'Tofts': This string indicates that you want the non-linear Tofts
%           Model (versus the non-linear RRM)
%    /::OUTPUT::\
%       pars: a vector with the following components
%           pars(1)=  Ktrans
%           pars(2)=  kep for the tissue of interest, kepTOI
%           pars(3)=  rsquare for the fitting
% Authors:
% Joey DeGrandchamp                 Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jdegrandchamp@email.arizona.edu   cardenaj@email.arizona.edu
%
%                       www.cardenaslab.org/resources
% v2.0 09/24/2015


% Check that the number of arguments is reasonable
narginchk(3,7);


if nargin == 4 %% LRRM
    toi = varargin{1};
    rr = varargin{2};
    time = varargin{3};
    fit = varargin{4};
    
    % Make sure all vectors are column, not row
    if size(toi,1) == 1
        toi = toi';
    end
    if size(rr,1) == 1
        rr = rr';
    end
    if size(time,1) == 1
        time = time';
    end
    
    % Make sure data vectors are same size
    if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(rr,1) || size(rr,1) ~= size(time,1)
       error('dce_mri_fit:toi,rr,and time vectors must all be same size'); 
    end
    
    % Set up LRRM (see references)
    y=toi;
        x1=rr;      
            x2=cumtrapz(time,rr); 
                x3=cumtrapz(time,toi);
                    X=[x1,x2,-x3];
    
    % Switch to users desired fitting algorithm
    switch fit
        case {'robust','Robust','ROBUST'}
            
            % Estimate parameters using robust method
            B = robustfit(X,y);
            B = B(2:end);
            
        case {'lsq','LSQ','least-squared'} 
            
            % Estimate parameters using LLSQ method (mldivide)
            B=X\y;
            
        case {'nonneg','Nonneg','NONNEG','non-neg'}
            
            % Estimate parameters using LLSQ method with non-negative constraint
            B=lsqnonneg(X,y);
            
        case {'lasso','Lasso','LASSO'}
            
            % Estimate parameters using lasso
            [b,FitInfo]=lasso(X,y);
            [~,ii]=min(FitInfo.MSE);    
            B=b(:,ii);
            
        otherwise
            
            error('dce_mri_fit:Invalid fit option for the LRRM');
            
    end
    
    %  Store each parameter of output pars (see reference) ##
    pars(1,1)=B(1);             %RKtrans;
    pars(2,1)=B(3);             %kepTOI;
    pars(3,1)=(B(2))*B(1)^-1;   %kepRR;

    % Calculate predicted curve for the TOI

    prediction=X*B; %#ok<*NBRAK>

    % Calculate the the RMSE for the predicted curve and allocate to pars

    pars(4,1)=rsquare(y,prediction);

elseif nargin == 5 || nargin == 6
    error('dce_mri_fit:Invalid number of arguments');

elseif nargin == 7 %% Non-linear fits (RRM/Tofts)
    
    model = varargin{7};
    
    switch model
        case {'RRM','NLRRM','NRRM','rrm','nlrrm','nrrm'}
    
            toi = varargin{1};
            rr = varargin{2};
            time = varargin{3};
            x0 = varargin{4};
            lb = varargin{5};
            ub = varargin{6};

            % Make sure all vectors are column, not row
            if size(toi,1) == 1
                toi = toi';
            end
            if size(rr,1) == 1
                rr = rr';
            end
            if size(time,1) == 1
                time = time';
            end

            % Make sure data vectors are same size
            if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(rr,1) || size(rr,1) ~= size(time,1)
               error('dce_mri_fit:toi,rr,and time vectors must all be same size'); 
            end

            % Set options for non-linear LSQ fitting
            S=optimset; S.Algorithm='trust-region-reflective'; S.Display='off';
            S.TolFun=1e-3; S.TolX=1e-3;
            S.MaxIter=1000;

            % Perform the fitting (the function "rrm" is in the nested functions below)
            [B,~,~] = lsqcurvefit(@rrm,x0,[time,rr],toi,lb,ub,S);

            % Store each parameter (see reference)
            pars(1,1)=B(1);    %RKtrans
            pars(2,1)=B(2);    %kepTOI
            pars(3,1)=B(3);    %kepRR

            pred=rrm(B,[time,rr]);
            pars(4)=rsquare(toi,pred);
       
        case {'Tofts','tofts','NLTM','nltm'}
            
            toi = varargin{1};
            Cp = varargin{2};
            time = varargin{3};
            x0 = varargin{4};
            lb = varargin{5};
            ub = varargin{6};
            
            % Make sure all vectors are column, not row
            if size(toi,1) == 1
                toi = toi';
            end
            if size(Cp,1) == 1
                Cp = Cp';
            end
            if size(time,1) == 1
                time = time';
            end

            % Make sure data vectors are same size
            if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(Cp,1) || size(Cp,1) ~= size(time,1)
               error('dce_mri_fit:toi,Cp,and time vectors must all be same size'); 
            end
            
            % Set options for non-linear LSQ fitting
            S=optimset; S.Algorithm='trust-region-reflective'; S.Display='off';
            S.TolFun=1e-3; S.TolX=1e-3;
            S.MaxIter=1000;

            % Perform the fitting (the function "rrm" is in the nested functions below)
            [B,~,~] = lsqcurvefit(@tofts,x0,[time,Cp],toi,lb,ub,S);

            % Store each parameter (see reference)
            pars(1,1)=B(1);    %Ktrans
            pars(2,1)=B(2);    %kepTOI

            pred=tofts(B,[time,Cp]);
            pars(3,1)=rsquare(toi,pred);
            
        otherwise
            
            error('dce_mri_fit:Must choose ''Tofts'' or ''RRM'' for last option');
            
    end

elseif nargin ==3 %% LTM
    toi = varargin{1};
    Cp = varargin{2};
    time = varargin{3};
    
    % Make sure all vectors are column, not row
    if size(toi,1) == 1
        toi = toi';
    end
    if size(Cp,1) == 1
        Cp = Cp';
    end
    if size(time,1) == 1
        time = time';
    end
    
    % Make sure data vectors are same size
    if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(Cp,1) || size(Cp,1) ~= size(time,1)
       error('dce_mri_fit:toi,Cp,and time vectors must all be same size'); 
    end
    
        y=toi;      
         x1=cumtrapz(time,Cp); 
             x2=cumtrapz(time,toi);
                 X=[x1,-x2,ones(size(y))];

    % Estimate parameters using a LLSQ with non-negative constraint

     [B,~,~]=lsqnonneg(X,y);
    % Calculate each parameter of output pars (see reference) ##

    pars(1,1)=B(1);             %Ktrans;
    pars(2,1)=B(2);             %kep;

    % Calculate predicted curve for the TOI

    prediction=X*B; %#ok<*NBRAK>

    % Calculate the the r-square for the predicted curve and allocated to pars

    pars(3,1)=rsquare(y,prediction);
end

end



%% ####### Internal Functions ##########


function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));

end

function C_toi = rrm(beta,X)
    
time=X(:,1); 
C_rr=X(:,2);

RKtrans= beta(1);  %ktrans_toi / ktrans_rr ;   RKtrans
kepRR=   beta(3);  %ktrans_rr  / ve_rr;    KepRR
kepTOI=  beta(2);  %ktrans_toi  / ve_toi;  KepTOI

R4= RKtrans*(kepRR-kepTOI);

C_toi=RKtrans.*C_rr + R4.*convolution(C_rr,exp(-kepTOI.*time));

end

function C_toi = tofts(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
    
    Ktrans = beta(1);
    kepTOI = beta(2);
    
    C_toi = Ktrans*convolution(Cp,exp(-kepTOI.*time));

end

function c = convolution(a, b, shape)
%CONVOLUTION MODIFIED BY JULIO CARDENAS, MAY OF 2011.
%   SAME THAN CONV BUT RETURN A VECTOR FROM a(1) to a(end), not the central
%   section as described for the usual convolution function.
%  

if ~isvector(a) || ~isvector(b)
  error(message('MATLAB:conv:AorBNotVector'));
end

if nargin < 3
    shape = 'full';
end

if ~ischar(shape)
  error(message('MATLAB:conv:unknownShapeParameter'));
end

% compute as if both inputs are column vectors
[rows,~]=size(a);
c = conv2(a(:),b(:),shape);
c=c(1:rows);

% restore orientation
if shape(1) == 'f'
    if length(a) > length(b)
        if size(a,1) == 1 %row vector
            c = c.';
        end
    else
        if size(b,1) == 1 %row vector
            c = c.';
        end
    end
else
    if size(a,1) == 1 %row vector
        c = c.';
    end
end

end