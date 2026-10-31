"""Independent cross-check of the RCAM-Flight-Dynamics MATLAB model.

The aerodynamic and moment equations below were re-implemented in Python
from a separate open implementation (github.com/reggiejohanes/RCAM,
RCAM_script_version/RCAM_script_dynamics.m), without reusing any code from
this repository. The script trims the aircraft in straight and level flight,
linearizes it by central differences and prints the trim and the four
lateral-directional and four longitudinal eigenvalues, for comparison with
the README results table.

Requires numpy and scipy.  Run:  python rcam_crosscheck.py
"""
import numpy as np
from scipy.optimize import minimize

m, cbar, lt, S, St = 120000.0, 6.6, 24.8, 260.0, 64.0
Ib = m * np.array([[40.07, 0, -2.0923], [0, 64.0, 0], [-2.0923, 0, 99.92]])
invIb = np.linalg.inv(Ib)
rac = np.array([0.12 * cbar, 0.0, 0.0])
r1, r2 = np.array([0.0, -7.94, -1.9]), np.array([0.0, 7.94, -1.9])
rho, g, depsda = 1.225, 9.81, 0.25
aL0, asw = np.deg2rad(-11.5), np.deg2rad(14.5)


def f(X, U, rcg):
    u, v, w, p, q, r, phi, th, _ = X
    uA, uE, uR, T1, T2 = U
    V, W = np.array([u, v, w]), np.array([p, q, r])
    Va = np.linalg.norm(V)
    al, be, Q = np.arctan2(w, u), np.arcsin(v / Va), 0.5 * rho * Va**2
    CLwb = 5.5 * (al - aL0) if al <= asw else -768.5*al**3 + 609.2*al**2 - 155.2*al + 15.212
    eps = depsda * (al - aL0)
    CL = CLwb + 3.1 * (St / S) * (al - eps + uE + 1.3 * q * lt / Va)
    CD = 0.13 + 0.07 * (5.5 * al + 0.654)**2
    CY = -1.6 * be + 0.24 * uR
    ca, sa = np.cos(al), np.sin(al)
    FA = np.array([[ca, 0, -sa], [0, 1, 0], [sa, 0, ca]]) @ (np.array([-CD, CY, -CL]) * Q * S)
    eta = np.array([-1.4 * be,
                    -0.59 - (3.1 * St * lt / (S * cbar)) * (al - eps),
                    (1 - al * 180 / (15 * np.pi)) * be])
    dCMdx = (cbar / Va) * np.array([[-11, 0, 5],
                                    [0, -4.03 * St * lt**2 / (S * cbar**2), 0],
                                    [1.7, 0, -11.5]])
    dCMdu = np.array([[-0.6, 0, 0.22],
                      [0, -3.1 * St * lt / (S * cbar), 0],
                      [0, 0, -0.63]])
    MA = (eta + dCMdx @ W + dCMdu @ np.array([uA, uE, uR])) * Q * S * cbar
    MA = MA + np.cross(FA, rcg - rac)
    F1, F2 = np.array([T1 * m * g, 0, 0]), np.array([T2 * m * g, 0, 0])
    mw = lambda ra: np.array([rcg[0] - ra[0], ra[1] - rcg[1], rcg[2] - ra[2]])
    ME = np.cross(mw(r1), F1) + np.cross(mw(r2), F2)
    Fg = m * g * np.array([-np.sin(th), np.cos(th) * np.sin(phi), np.cos(th) * np.cos(phi)])
    vdot = (Fg + F1 + F2 + FA) / m - np.cross(W, V)
    wdot = invIb @ (MA + ME - np.cross(W, Ib @ W))
    H = np.array([[1, np.sin(phi) * np.tan(th), np.cos(phi) * np.tan(th)],
                  [0, np.cos(phi), -np.sin(phi)],
                  [0, np.sin(phi) / np.cos(th), np.cos(phi) / np.cos(th)]])
    return np.concatenate([vdot, wdot, H @ W])


def trim(Va, rcg):
    def build(z):
        al, uE, T = z
        return (np.array([Va * np.cos(al), 0, Va * np.sin(al), 0, 0, 0, 0, al, 0]),
                np.array([0, uE, 0, T, T]))
    res = minimize(lambda z: np.sum(f(*build(z), rcg)[:6]**2), [0.05, -0.15, 0.08],
                   method="Nelder-Mead",
                   options=dict(xatol=1e-12, fatol=1e-14, maxiter=40000, maxfev=40000))
    return (*build(res.x), res.x, res.fun)


def jacobian(X, U, rcg, h=1e-6):
    A = np.zeros((9, 9))
    for j in range(9):
        d = np.zeros(9); d[j] = h
        A[:, j] = (f(X + d, U, rcg) - f(X - d, U, rcg)) / (2 * h)
    return A


def describe(lam):
    out = []
    for l in sorted(lam, key=lambda z: (abs(z.imag) < 1e-6, z.real)):
        if abs(l.imag) > 1e-6:
            if l.imag > 0:
                out.append(f"{l.real:+.3f} ± {l.imag:.3f}j (wn {abs(l):.3f}, zeta {-l.real/abs(l):.3f})")
        else:
            out.append(f"{l.real:+.3f}")
    return "; ".join(out)


if __name__ == "__main__":
    LAT, LON = [1, 3, 5, 6], [0, 2, 4, 7]   # (v,p,r,phi), (u,w,q,theta)
    for zcg in (0.10, 0.0):
        rcg = np.array([0.23 * cbar, 0.0, zcg * cbar])
        for Va in (85.0, 80.0):
            X, U, z, res = trim(Va, rcg)
            A = jacobian(X, U, rcg)
            print(f"Zcg = {zcg:.2f} cbar, Va = {Va:.0f} m/s: alpha {np.rad2deg(z[0]):.2f} deg, "
                  f"uE {np.rad2deg(z[1]):.2f} deg, uTH {np.rad2deg(z[2]):.2f} deg (residual {res:.1e})")
            print("  lateral     :", describe(np.linalg.eigvals(A[np.ix_(LAT, LAT)])))
            print("  longitudinal:", describe(np.linalg.eigvals(A[np.ix_(LON, LON)])))
