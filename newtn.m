function [x,J,iflag] = newtn(x0, F, tol,vararin)
%[x,J,iflag] = newt(x0,F,tol);
% Simple function that applies Newton's method to find the root of

%  F(x) = 0
% 
% Input: 
% x0: initial iterate
% F: function for which we seek a zero
% tol: convergence criteria ||F(x)||<tol
%
% Output:
% x: steady state solution
% J: Jacobian = [partial F/partial x]
% iflag: 0 ==> Newton's method converged to the desired
%                      tolerance
%        1 ==> Newton's method did not converge 
    iprint = 0;
    MAXIT = 40;
    x = x0;
    if (nargin==4)
        [F0,iJ] = F(x);
    else
        [F0, J] = F(x);
    end
    iflag = 0; itno = 0;
    while (((norm(F0) > tol) & (itno<MAXIT)) )
        if (nargin==4)            
            dx = -iJ(F0);
            x = x + dx;
            [F0, iJ] = F(x);
        else
            warning('off');
            dx = -J\F0;
            x = x + dx;
            [F0, J] = F(x);
            warning('on');
        end
        itno = itno+1;
        if (iprint)
            fprintf('%i: norm(F0) = %e norm(F0) = %e \n',itno,norm(F0),norm(F0(1:end-1)));
        end
    end
    if (itno>=MAXIT)
        iflag = 1;
        fprintf('Warning Newton''s Method did not converge.\n')
    end
    if (nargin==4)
        J = iJ;
    end
end