import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 2D Laplace Potential Flow using BEM
#
# Domain:
#     0 <= x <= 1
#     0 <= y <= 1
#
# Boundary condition:
#     phi = x
#
# Exact solution:
#     phi = x
#     velocity = (1, 0)
#
# Constant boundary elements
# ============================================================

N_SIDE = 16
N_GAUSS = 10


# ============================================================
# 1. Boundary mesh
# ============================================================

def make_square_boundary(n_side):

    pts = []

    # bottom
    # (0,0) -> (1,0)
    for k in range(n_side):
        pts.append([k / n_side, 0.0])

    # right
    # (1,0) -> (1,1)
    for k in range(n_side):
        pts.append([1.0, k / n_side])

    # top
    # (1,1) -> (0,1)
    for k in range(n_side):
        pts.append([1.0 - k / n_side, 1.0])

    # left
    # (0,1) -> (0,0)
    for k in range(n_side):
        pts.append([0.0, 1.0 - k / n_side])

    pts = np.asarray(pts, dtype=float)

    # close boundary
    closed = np.vstack([pts, pts[0]])

    p1 = closed[:-1]
    p2 = closed[1:]

    # element center
    mid = 0.5 * (p1 + p2)

    d = p2 - p1

    # element length
    L = np.linalg.norm(d, axis=1)

    # CCW boundary
    # outward normal = right-hand normal
    normal = np.column_stack(
        [d[:, 1], -d[:, 0]]
    ) / L[:, None]

    return p1, p2, mid, L, normal


p1, p2, mid, L, normal = make_square_boundary(N_SIDE)

NE = len(L)


# ============================================================
# 2. Gauss integration
# ============================================================

gp, gw = np.polynomial.legendre.leggauss(N_GAUSS)


# ============================================================
# 3. Fundamental solution
#
# G = -(1/2pi) ln(r)
#
# dG/dn =
# -(1/2pi) ((x-xi) dot n) / r^2
# ============================================================

def integrate_element(source, a, b, nvec):

    length = np.linalg.norm(b - a)

    IG = 0.0
    IH = 0.0

    for z, w in zip(gp, gw):

        # natural coordinate
        x = 0.5 * (
            (1 - z) * a
            +
            (1 + z) * b
        )

        J = length / 2.0

        rvec = x - source

        r = np.linalg.norm(rvec)

        # Green function
        G = -(1.0 / (2 * np.pi)) * np.log(r)

        # normal derivative of Green function
dGdn = -(1.0 / (2 * np.pi)) * (
    np.dot(rvec, nvec) / (r * r)
)
        IG += w * G * J

        IH += w * dGdn * J

    return IG, IH


# ============================================================
# 4. Assemble H and G matrices
# ============================================================

H = np.zeros((NE, NE))

Gmat = np.zeros((NE, NE))


for i in range(NE):

    # collocation/source point
    source = mid[i]

    for j in range(NE):

        # -----------------------------------------
        # singular self element
        # -----------------------------------------
        if i == j:

            # CPV / jump term
            H[i, j] = 0.5

            # analytic logarithmic singular integral
            Gmat[i, j] = (
                L[j]
                / (2 * np.pi)
                * (
                    1.0
                    - np.log(L[j] / 2.0)
                )
            )

        # -----------------------------------------
        # regular element
        # -----------------------------------------
        else:

            IG, IH = integrate_element(
                source,
                p1[j],
                p2[j],
                normal[j]
            )

            Gmat[i, j] = IG

            H[i, j] = IH


# ============================================================
# 5. Boundary condition
#
# phi = x
# ============================================================

phi_boundary = mid[:, 0]


# BEM equation:
#
# H phi = G q
#
# unknown q

q_boundary = np.linalg.solve(
    Gmat,
    H @ phi_boundary
)


# ============================================================
# 6. Compare boundary flux with exact solution
#
# Exact:
# phi = x
#
# grad(phi) = (1,0)
#
# q = grad(phi) dot n
#   = n_x
# ============================================================

q_exact = normal[:, 0]

error_q = np.sqrt(
    np.mean(
        (q_boundary - q_exact)**2
    )
)

print(
    "Boundary q RMS error =",
    error_q
)


# ============================================================
# 7. Interior potential
#
# phi(xi)
# =
# integral [
# G q
# -
# phi dG/dn
# ] dGamma
# ============================================================

def interior_phi(xi):

    value = 0.0

    for j in range(NE):

        IG, IH = integrate_element(
            xi,
            p1[j],
            p2[j],
            normal[j]
        )

        value += (
            IG * q_boundary[j]
            -
            IH * phi_boundary[j]
        )

    return value


# ============================================================
# 8. Interior grid
# ============================================================

xs = np.linspace(
    0.05,
    0.95,
    25
)

ys = np.linspace(
    0.05,
    0.95,
    25
)

X, Y = np.meshgrid(xs, ys)

Phi = np.zeros_like(X)


for iy in range(len(ys)):

    for ix in range(len(xs)):

        point = np.array(
            [
                X[iy, ix],
                Y[iy, ix]
            ]
        )

        Phi[iy, ix] = interior_phi(point)


# ============================================================
# 9. Velocity
#
# u = grad(phi)
# ============================================================

dphi_dy, dphi_dx = np.gradient(
    Phi,
    ys,
    xs
)

U = dphi_dx

V = dphi_dy


# ============================================================
# 10. Exact solution comparison
# ============================================================

phi_exact = X

error_phi = np.max(
    np.abs(
        Phi - phi_exact
    )
)

print(
    "Max phi error =",
    error_phi
)

print(
    "Mean U =",
    np.mean(U)
)

print(
    "Mean V =",
    np.mean(V)
)


# ============================================================
# 11. Plot
# ============================================================

plt.figure(
    figsize=(8, 6)
)

contour = plt.contourf(
    X,
    Y,
    Phi,
    levels=20
)

plt.colorbar(
    contour,
    label="Potential phi"
)

# velocity vectors
plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2]
)

# boundary
plt.plot(
    np.r_[p1[:, 0], p1[0, 0]],
    np.r_[p1[:, 1], p1[0, 1]],
    "k-"
)

plt.xlabel("x")
plt.ylabel("y")

plt.title(
    "2D Laplace Potential Flow by BEM"
)

plt.axis("equal")

plt.tight_layout()

plt.show()
