import numpy as np

def limp(z, y, w, obs, sys, opt):
    #
    # Negative log probability for the CO2 system a.k.a. log improbability i.e., limp!
    #
    # INPUT:
    # y := measured components, non-measured components set to NaN
    # w := measurement precisions, same size and order as y, non-measured components set to anything, they are ignored
    # gpK := first order derivatives of pK wrt T, S, P
    # ggpK := second order derivatives of pK wrt T, S, P

    # OUTPUT:
    #
    # f := limp
    # g := grad f w.r.t. x
    # h := hessian of f w.r.t. x

    p = sys['p']
    q = sys['q']
    M = sys['M']
    K = sys['K']
    nrk = K.shape[0]
    nTP = len(sys['tp'])
    nv = M.shape[1]
    x = z[:nv]  # thermodynamic system state, PP*x is modeled measured vars

    # fill ypK, gypK, ggypK with associated calculated pK, gpK, and ggpK values
    # update the pK values based on the new estimate of (T, P)
    y, gy, ggy = update_y(y, x, obs, sys, opt)
    # Make a vector of measured quantities
    id = np.where(~np.isnan(y))[0]
    y = y[id]
    gy = gy[id]
    ggy = ggy[id]

    # Make a precision matrix
    W = np.diag(w[id])

    # Build a matrix that picks out the measured components of x
    I = np.eye(nv)  # for chain rule
    PP = I[id, :]  # picking/pick out the measured ones
    e = PP @ x - y  # calculated - measured (minus)
    ge = PP - gy

    nlam = M.shape[0] + K.shape[0]
    lam = z[nv:]  # Lagrange multipliers 

    # Constraint equations
    c = np.concatenate((M @ q(x), K @ x))
    f = 0.5 * e.T @ W @ e + lam.T @ c  # limp, method of Lagrange multipliers    
    # -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    dcdx = np.vstack([M @ np.diag(sys['dqdx'](x)), K])

    # iTSP = [[sys.tp(:).iT], sys.isal, [sys.tp(:).iP]];
    g = np.concatenate((e.T @ W @ ge + lam.T @ dcdx, c.T))

    ddq = np.diag(sys['d2qdx2'](x))  # q"

    nr, nc = M.shape
    gg = np.zeros((nc, nc))
    for row in range(nr):
        gg += lam[row] * np.diag(M[row, :]) @ ddq

    tmp = np.zeros((len(x), len(x)))
    eW = e.T @ W
    for jj in range(ggy.shape[0]):
        tmp += eW[jj] * ggy[jj, :, :]

    H = np.block([
        [ge.T @ W @ ge - tmp + gg, dcdx.T],  # derivatives wrt lambdas
        [dcdx, np.zeros((nlam, nlam))]  # derivatives wrt variables
    ])

    g = g.flatten()  # make sure g is returned as a column vector
    return g, H, f
