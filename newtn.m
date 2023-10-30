function [x,J,iflag] = newtn(x0, F, tol)
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
    iprint = 1;
    MAXIT = 50;
    x = x0;
    if (nargin==4)
        [F0,iJ] = F(x);
    else
        [F0, J] = F(x); 
    end
    % keyboard
    iflag = 0; itno = 0;
    while (((norm(F0) > tol) && (itno<MAXIT)) )
        if (nargin==4)            
            dx = -iJ(F0);
            x = x + dx;
            [F0, iJ] = F(x);
        else
            %[P,R,C] = equilibrate(J);
            %B = R*P*J*C;
            %dx = -C*(B\(R*P*F0));
            %x = x+dx;            
            % warning('off');
            dx = -J\F0;
            x = x + dx;
            %warning('on');
            [F0, J] = F(x);
            % keyboard
        end
        itno = itno+1;
        if (iprint)
            fprintf('%i: norm(F0) = %e norm(F0) = %e \n',itno,norm(F0),norm(F0(1:end-1)));
        end
    end
    if (itno>=MAXIT)
        if norm(F0(1:end-1)) < tol*1e1
            iflag = 2;
            fprintf('Warniing Newton''s Method did not converge.\n ')
            fprintf('At max iteration, value was one order of magnitude more than tolerance...\n')
            fprintf('Consider keeping. \n')
        else
            iflag = 1;
            fprintf('Warning Newton''s Method did not converge.\n')
        end
    end
    if (nargin==4)
        J = iJ;
    end
end
