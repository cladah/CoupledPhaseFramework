def eulerbackwardstep(f,modelvar,mesh,ti,dt,yi,bi):
    from scipy.optimize import fsolve
    import numpy as np
    error = 1
    #while error >1E-6:
    t0 = ti
    y0 = yi
    tp = t0 + dt
    yp = y0 + dt * f(t0, y0)

    yp = fsolve(backward_euler_residual, yp, args=(f, t0, y0, tp))

    return yp

def backward_euler_residual(yp,f,t0,y0,tp):
    value = yp - y0 - (tp - t0) * f(tp, yp);
    return value

def rhs():
    b = yi
def pde(y,t):
    fsolve(lambda dydt : dydt - rhs(y,dydt, t),np.zeros())