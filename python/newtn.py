import numpy as np

def newtn(x0, F, tol, iJ = None):
# [x,J,iflag] = newt(x0,F,tol);
#  Simple function that applies Newton's method to find the root of

#   F(x) = 0
#  
#  Input: 
#  x0: initial iterate
#  F: function for which we seek a zero
#  tol: convergence criteria ||F(x)||<tol
# 
#  Output:
#  x: steady state solution
#  J: Jacobian = [partial F/partial x]
#  iflag: 0 ==> Newton's method converged to the desired
#                       tolerance
#         1 ==> Newton's method did not converge 
    iprint = 0
    MAXIT = 50
    x = x0
    if iJ is not None:
        F0, _ = F(x)
    else:
        F0, J = F(x)
    iflag = 0
    itno = 0

    while np.linalg.norm(F0) > tol and itno < MAXIT:
        if iJ is not None:
            dx = -iJ(F0)
            x = x + dx
            F0, _ = F(x)
        else:
            dx = -np.linalg.solve(J, F0)
            x = x + dx
            F0, J = F(x)
        
        itno += 1
        if iprint:
            print(f'{itno}: norm(F0) = {np.linalg.norm(F0):e} norm(F0[:-1]) = {np.linalg.norm(F0[:-1]):e}')
    
    if itno >= MAXIT:
        if np.linalg.norm(F0[:-1]) < tol * 1e1:
            iflag = 2
            print("Warning: Newton's Method did not converge.")
            print("But value was only 10x more than tolerance...")
            print("so it was close, recommend keeping.")
        else:
            iflag = 1
            print("Warning: Newton's Method did not converge.")
    
    if iJ is not None:
        J = iJ

    return x, J, iflag

    