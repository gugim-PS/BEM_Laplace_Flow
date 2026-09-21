2D Laplace Potential Flow Using the Boundary Element Method (BEM)

1. Overview

This example solves a simple two-dimensional potential-flow problem using the Boundary Element Method (BEM).

The computational domain is the unit square

$$0 \le x \le 1, \qquad 0 \le y \le 1.$$

The prescribed velocity potential on the boundary is

$$\phi(x,y)=x.$$

The exact solution inside the domain is also

$$\phi(x,y)=x.$$

Therefore, the velocity field is

$$\mathbf{u}=\nabla\phi,$$

or

$$u=\frac{\partial\phi}{\partial x}=1, \qquad v=\frac{\partial\phi}{\partial y}=0.$$

Thus, the exact flow is a uniform flow from left to right:

$$\mathbf{u}=(1,0).$$

The purpose of this program is to demonstrate the basic BEM procedure:

Generate boundary elements.

Define the Laplace fundamental solution.

Evaluate boundary integrals using Gauss quadrature.

Assemble the BEM matrices $H$ and $G$.

Solve for the unknown normal derivative $q=\partial\phi/\partial n$.

Reconstruct the potential inside the domain.

Obtain the internal velocity field from $\nabla\phi$.

2. Governing Equation

For incompressible and irrotational potential flow,

$$\mathbf{u}=\nabla\phi,$$

where $\phi$ is the velocity potential.

For incompressible flow,

$$\nabla\cdot\mathbf{u}=0.$$

Substituting $\mathbf{u}=\nabla\phi$ gives

$$\nabla\cdot(\nabla\phi)=0,$$

and therefore

$$\boxed{\nabla^2\phi=0}.$$

In two dimensions,

$$\boxed{\frac{\partial^2\phi}{\partial x^2}+\frac{\partial^2\phi}{\partial y^2}=0}.$$

This is the Laplace equation.

3. Boundary Integral Equation

BEM transforms the differential equation in the domain into an integral equation defined only on the boundary.

For the two-dimensional Laplace equation,

\int_{\Gamma}
G(\mathbf{x},\boldsymbol{\xi})
q(\mathbf{x})
,d\Gamma.
$$

Here,

$$q=\frac{\partial\phi}{\partial n}.$$

The symbols are:

$\boldsymbol{\xi}$: source or collocation point,

$\mathbf{x}$: integration point on the boundary,

$\Gamma$: boundary of the domain,

$G$: fundamental solution,

$n$: outward unit normal direction.

For a smooth two-dimensional boundary,

$$\boxed{c(\boldsymbol{\xi})=\frac{1}{2}}.$$

4. Fundamental Solution

For the 2D Laplace equation, the fundamental solution is

$$\boxed{G(\mathbf{x},\boldsymbol{\xi})=-\frac{1}{2\pi}\ln r},$$

where

$$r=\left|\mathbf{x}-\boldsymbol{\xi}\right|.$$

The normal derivative of the fundamental solution is

-\frac{1}{2\pi}
\frac{(\mathbf{x}-\boldsymbol{\xi})\cdot\mathbf{n}}{r^2}
}.
$$

These two kernels are used to construct the BEM matrices.

5. Boundary Discretization

Unlike FEM, BEM does not mesh the entire two-dimensional domain. Only the boundary is discretized.

For the square domain, each side is divided into N_SIDE constant boundary elements.

N_SIDE = 16

Since the square has four sides,

$$N_E=4N_{\mathrm{SIDE}}.$$

With $N_{\mathrm{SIDE}}=16$,

$$\boxed{N_E=64}.$$

The function

make_square_boundary(n_side)

creates the four boundary sides in counterclockwise order.

For every element, the program computes:

first endpoint p1,

second endpoint p2,

element midpoint mid,

element length L,

outward unit normal normal.

For one element,

$$\mathbf{d}=\mathbf{x}_2-\mathbf{x}_1,$$

and

$$L_e=|\mathbf{d}|.$$

Because the boundary is ordered counterclockwise, the outward normal is

\frac{1}{L_e}
\begin{bmatrix}
d_y\
-d_x
\end{bmatrix}
}.
$$

6. Constant Boundary Elements

This program uses constant boundary elements.

Within one element,

$$\phi(\mathbf{x})\approx\phi_j,$$

and

$$q(\mathbf{x})\approx q_j.$$

Here, $\phi_j$ and $q_j$ are represented by the values at the element midpoint.

After dividing the boundary into $N_E$ elements, the continuous boundary integral equation becomes

\sum_{j=1}^{N_E}
q_j
\int_{\Gamma_j}
G_i
,d\Gamma.
$$

Define the matrix coefficient

c_i\delta_{ij}
+
\int_{\Gamma_j}
\frac{\partial G_i}{\partial n}
,d\Gamma
},
$$

and

\int_{\Gamma_j}
G_i
,d\Gamma
}.
$$

Then the discretized BEM equation becomes

$$\boxed{H\Phi=GQ}.$$

The boundary-potential vector is

$$
\Phi=
\begin{bmatrix}
\phi_1\
\phi_2\
\vdots\
\phi_{N_E}
\end{bmatrix},
$$

and the boundary-flux vector is

$$
Q=
\begin{bmatrix}
q_1\
q_2\
\vdots\
q_{N_E}
\end{bmatrix}.
$$

7. Gauss Quadrature

For regular boundary elements, the required integrals are evaluated numerically using Gauss-Legendre quadrature.

The standard formula is

$$
\boxed{
\int_{-1}^{1} f(z),dz
\approx
\sum_{k=1}^{N_G} w_k f(z_k)
}.
$$

Here,

$z_k$: Gauss points,

$w_k$: Gauss weights,

$N_G$: number of Gauss points.

The code generates them with

gp, gw = np.polynomial.legendre.leggauss(N_GAUSS)

For example,

N_GAUSS = 10

uses ten Gauss points for each regular element integral.

8. Mapping the Gauss Point to the Physical Element

Gauss quadrature is defined on the standard interval

$$-1\le z\le1.$$

A physical boundary element extends from point $\mathbf{a}$ to $\mathbf{b}$.

The coordinate mapping is

\frac12\left[(1-z)\mathbf{a}+(1+z)\mathbf{b}\right]
}.
$$

In the code,

x = 0.5 * (
    (1 - z) * a
    +
    (1 + z) * b
)

The Jacobian of this transformation is

$$\boxed{J=\frac{L_e}{2}},$$

so

$$d\Gamma=J,dz.$$

Therefore,

\int_{-1}^{1} f(\mathbf{x}(z))J,dz,
$$

and Gauss quadrature gives

$$
\boxed{
\int_{\Gamma_e} f,d\Gamma
\approx
\sum_{k=1}^{N_G} w_k f(z_k)J
}.
$$

9. Element Integration Function

The function

integrate_element(source, a, b, nvec)

calculates

$$I_G=\int_{\Gamma_e}G,d\Gamma,$$

and

$$I_H=\int_{\Gamma_e}\frac{\partial G}{\partial n},d\Gamma.$$

For each Gauss point,

rvec = x - source
r = np.linalg.norm(rvec)

computes

$$\mathbf{r}=\mathbf{x}-\boldsymbol{\xi},$$

and

$$r=|\mathbf{r}|.$$

The fundamental solution is evaluated as

G = -(1.0 / (2 * np.pi)) * np.log(r)

which corresponds to

$$G=-\frac{1}{2\pi}\ln r.$$

The normal derivative is evaluated as

dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)

which corresponds to

-\frac{1}{2\pi}
\frac{\mathbf{r}\cdot\mathbf{n}}{r^2}.
$$

Finally,

IG += w * G * J
IH += w * dGdn * J

implements

$$I_G\approx\sum_k w_k G_k J,$$

and

$$
I_H
\approx
\sum_k
w_k
\left(\frac{\partial G}{\partial n}\right)_k
J.
$$

10. Singular Self Element

A special case occurs when the source point lies on the same boundary element being integrated.

Then

$$r\rightarrow0,$$

and the Green function becomes singular.

For the $H$ matrix, the Cauchy Principal Value treatment produces the jump term

$$\boxed{H_{ii}=\frac12}$$

for a smooth boundary.

The code uses

H[i, j] = 0.5

when

i == j

For the logarithmic singularity in the $G$ matrix, the self-element integral for a constant straight element can be evaluated analytically:

\frac{L_i}{2\pi}
\left[1-\ln\left(\frac{L_i}{2}\right)\right]
}.
$$

The code implements this as

Gmat[i, j] = (
    L[j]
    / (2 * np.pi)
    * (
        1.0
        - np.log(L[j] / 2.0)
    )
)

This avoids applying ordinary Gauss quadrature directly at the singularity.

11. Assembly of the Global BEM Matrices

The program loops over all source elements $i$ and integration elements $j$.

for i in range(NE):
    source = mid[i]

    for j in range(NE):

Each pair $(i,j)$ generates one matrix coefficient.

Thus,

$$
H=
\begin{bmatrix}
H_{11}&H_{12}&\cdots&H_{1N}\
H_{21}&H_{22}&\cdots&H_{2N}\
\vdots&\vdots&\ddots&\vdots\
H_{N1}&H_{N2}&\cdots&H_{NN}
\end{bmatrix},
$$

and similarly for $G$.

Unlike a typical FEM stiffness matrix, BEM matrices are generally dense, because every boundary element influences every collocation point through the Green function.

12. Boundary Condition

The prescribed boundary condition is

$$\boxed{\phi=x}.$$

Since constant boundary elements are used, the prescribed value is evaluated at each element midpoint:

phi_boundary = mid[:, 0]

Therefore,

$$\phi_j=x_j$$

for every boundary element.

The unknown is

$$q_j=\frac{\partial\phi}{\partial n}.$$

The BEM equation is

$$H\Phi=GQ,$$

so formally,

$$\boxed{Q=G^{-1}H\Phi}.$$

The program solves this system without explicitly computing the inverse:

q_boundary = np.linalg.solve(
    Gmat,
    H @ phi_boundary
)

Using np.linalg.solve is numerically preferable to explicitly forming $G^{-1}$.

13. Exact Boundary Flux

Since

$$\phi=x,$$

the exact gradient is

$$
\nabla\phi=
\begin{bmatrix}
1\
0
\end{bmatrix}.
$$

The normal derivative is

$$q=\nabla\phi\cdot\mathbf{n}.$$

Therefore,

$$\boxed{q=n_x}.$$

This is implemented as

q_exact = normal[:, 0]

Physically:

left wall: $q=-1$,

right wall: $q=+1$,

bottom wall: $q=0$,

top wall: $q=0$.

The RMS boundary error is

\sqrt{
\frac{1}{N_E}
\sum_{j=1}^{N_E}
\left(q_j-q_j^{\mathrm{exact}}\right)^2
}
}.
$$

It is computed by

error_q = np.sqrt(
    np.mean(
        (q_boundary - q_exact)**2
    )
)

14. Interior Potential Reconstruction

An important feature of BEM is that the interior does not need to be meshed.

Once the boundary values $\phi$ and $q$ are known, the potential at any interior point $\boldsymbol{\xi}$ can be reconstructed from the boundary solution:

\phi\frac{\partial G}{\partial n}
\right]
,d\Gamma
}.
$$

For an interior point, the coefficient multiplying $\phi(\boldsymbol{\xi})$ is 1.

The code evaluates the integral as

value += (
    IG * q_boundary[j]
    -
    IH * phi_boundary[j]
)

for every boundary element.

This is why BEM can display an internal potential field even though only the boundary was discretized.

The internal points are evaluation points, not additional BEM unknowns.

15. Interior Grid

For visualization, the program creates a set of interior field points:

xs = np.linspace(0.05, 0.95, 25)
ys = np.linspace(0.05, 0.95, 25)

The endpoints are kept slightly away from the boundary to reduce near-singular integration errors in this simple example.

A regular Cartesian grid is formed:

X, Y = np.meshgrid(xs, ys)

For each point,

Phi[iy, ix] = interior_phi(point)

evaluates the boundary integral representation.

The resulting array Phi contains

$$\phi(x_i,y_j)$$

throughout the interior.

16. Velocity Field

For potential flow,

$$\boxed{\mathbf{u}=\nabla\phi}.$$

Therefore,

$$u=\frac{\partial\phi}{\partial x},$$

and

$$v=\frac{\partial\phi}{\partial y}.$$

The code obtains these derivatives numerically:

dphi_dy, dphi_dx = np.gradient(
    Phi,
    ys,
    xs
)

and then

U = dphi_dx
V = dphi_dy

Hence,

$$\boxed{U\approx\frac{\partial\phi}{\partial x}},$$

and

$$\boxed{V\approx\frac{\partial\phi}{\partial y}}.$$

For this example,

$$U\approx1, \qquad V\approx0$$

throughout the domain.

17. Exact Interior Solution

The exact potential is

$$\phi^{\mathrm{exact}}=x.$$

In the code,

phi_exact = X

The maximum absolute potential error is

\max
\left|
\phi^{\mathrm{BEM}}-\phi^{\mathrm{exact}}
\right|
}.
$$

It is calculated with

error_phi = np.max(
    np.abs(
        Phi - phi_exact
    )
)

The program also reports

np.mean(U)
np.mean(V)

which should approach

$$\boxed{\overline{U}\rightarrow1, \qquad \overline{V}\rightarrow0}$$

as the discretization becomes sufficiently accurate.

18. Visualization

The potential field is drawn using

plt.contourf(
    X,
    Y,
    Phi,
    levels=20
)

This gives a contour map of $\phi(x,y)$.

The velocity field is drawn using

plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2]
)

The expected result is a set of nearly horizontal arrows pointing in the positive $x$-direction:

→ → → → → → →
→ → → → → → →
→ → → → → → →
→ → → → → → →

19. Complete Computational Flow

The complete BEM calculation can be summarized as follows.

Step 1: Governing equation

$$\nabla^2\phi=0.$$

Step 2: Boundary integral equation

\int_{\Gamma}
Gq
,d\Gamma.
$$

Step 3: Boundary discretization

$$\Gamma\rightarrow\sum_e\Gamma_e.$$

Step 4: Gauss quadrature

$$
\int_{\Gamma_e}(\cdots),d\Gamma
\approx
\sum_k w_k(\cdots)_kJ.
$$

Step 5: Matrix assembly

$$\boxed{H\Phi=GQ}.$$

Step 6: Solve boundary unknowns

$$\boxed{Q=G^{-1}H\Phi}.$$

Step 7: Interior reconstruction

\phi\frac{\partial G}{\partial n}
\right]
,d\Gamma.
$$

Step 8: Velocity field

$$\boxed{\mathbf{u}=\nabla\phi}.$$

20. Why This Example Is Useful

This problem has an exact solution, so it is useful for checking whether the BEM implementation is correct.

It demonstrates the essential BEM ideas:

only the boundary is discretized,

the Green function transfers boundary information into the domain,

Gauss quadrature evaluates regular boundary integrals,

singular self-integrals require special treatment,

BEM matrices are generally dense,

interior values are reconstructed after the boundary solution is known.

21. Important Numerical Notes

21.1 Python line continuation

A line-continuation backslash with trailing spaces can cause a syntax error. A safer style is to use parentheses:

dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)

21.2 Near-singular integration

The present example treats exactly singular self-elements analytically.

However, if an interior evaluation point lies extremely close to a boundary element, ordinary Gauss quadrature may lose accuracy even though the integral is not mathematically singular.

More advanced BEM programs may use:

adaptive quadrature,

element subdivision,

special near-singular integration,

singularity subtraction.

21.3 Higher-order elements

The current program uses constant elements because they are easy to understand.

More accurate BEM formulations may use linear or quadratic boundary elements:

$$\phi(\xi)=\sum_jN_j(\xi)\phi_j,$$

and

$$q(\xi)=\sum_jN_j(\xi)q_j,$$

where $N_j$ are shape functions.

22. Requirements

The script requires:

numpy
matplotlib

Install them with

pip install numpy matplotlib

23. Running the Program

Save the Python code as, for example,

bem_laplace_flow.py

and run

python bem_laplace_flow.py

The program prints:

boundary-flux RMS error,

maximum potential error,

mean $u$-velocity,

mean $v$-velocity.

It also displays the potential contour and velocity vectors.

24. Expected Physical Result

Because

$$\phi=x,$$

the velocity should be

$$\boxed{\mathbf{u}=(1,0)}.$$

Therefore:

potential contours should be approximately vertical,

velocity arrows should point horizontally to the right,

$U$ should be approximately 1,

$V$ should be approximately 0.

25. Main Concept

The key idea is

$$
\boxed{
\text{Solve only on the boundary}
\quad\Longrightarrow\quad
\text{reconstruct the solution anywhere inside}
}.
$$

This is the defining computational idea behind the Boundary Element Method for Laplace-type problems.

26. Complete Python Example

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


def make_square_boundary(n_side):
    pts = []

    # bottom: (0,0) -> (1,0)
    for k in range(n_side):
        pts.append([k / n_side, 0.0])

    # right: (1,0) -> (1,1)
    for k in range(n_side):
        pts.append([1.0, k / n_side])

    # top: (1,1) -> (0,1)
    for k in range(n_side):
        pts.append([1.0 - k / n_side, 1.0])

    # left: (0,1) -> (0,0)
    for k in range(n_side):
        pts.append([0.0, 1.0 - k / n_side])

    pts = np.asarray(pts, dtype=float)
    closed = np.vstack([pts, pts[0]])

    p1 = closed[:-1]
    p2 = closed[1:]

    mid = 0.5 * (p1 + p2)
    d = p2 - p1
    L = np.linalg.norm(d, axis=1)

    # CCW boundary: right-hand normal points outward
    normal = np.column_stack([d[:, 1], -d[:, 0]]) / L[:, None]

    return p1, p2, mid, L, normal


p1, p2, mid, L, normal = make_square_boundary(N_SIDE)
NE = len(L)

# Gauss-Legendre integration points and weights
gp, gw = np.polynomial.legendre.leggauss(N_GAUSS)


def integrate_element(source, a, b, nvec):
    length = np.linalg.norm(b - a)

    IG = 0.0
    IH = 0.0

    for z, w in zip(gp, gw):
        x = 0.5 * ((1 - z) * a + (1 + z) * b)
        J = length / 2.0

        rvec = x - source
        r = np.linalg.norm(rvec)

        G = -(1.0 / (2 * np.pi)) * np.log(r)

        dGdn = (
            -(1.0 / (2 * np.pi))
            * np.dot(rvec, nvec)
            / (r * r)
        )

        IG += w * G * J
        IH += w * dGdn * J

    return IG, IH


# Assemble BEM matrices
H = np.zeros((NE, NE))
Gmat = np.zeros((NE, NE))

for i in range(NE):
    source = mid[i]

    for j in range(NE):
        if i == j:
            H[i, j] = 0.5

            Gmat[i, j] = (
                L[j]
                / (2 * np.pi)
                * (1.0 - np.log(L[j] / 2.0))
            )
        else:
            IG, IH = integrate_element(
                source,
                p1[j],
                p2[j],
                normal[j],
            )

            Gmat[i, j] = IG
            H[i, j] = IH


# Boundary condition: phi = x
phi_boundary = mid[:, 0]

# Solve H phi = G q
q_boundary = np.linalg.solve(
    Gmat,
    H @ phi_boundary,
)

# Exact boundary flux
q_exact = normal[:, 0]

error_q = np.sqrt(
    np.mean((q_boundary - q_exact) ** 2)
)

print("Boundary q RMS error =", error_q)


def interior_phi(xi):
    value = 0.0

    for j in range(NE):
        IG, IH = integrate_element(
            xi,
            p1[j],
            p2[j],
            normal[j],
        )

        value += (
            IG * q_boundary[j]
            - IH * phi_boundary[j]
        )

    return value


# Interior evaluation grid
xs = np.linspace(0.05, 0.95, 25)
ys = np.linspace(0.05, 0.95, 25)

X, Y = np.meshgrid(xs, ys)
Phi = np.zeros_like(X)

for iy in range(len(ys)):
    for ix in range(len(xs)):
        point = np.array([X[iy, ix], Y[iy, ix]])
        Phi[iy, ix] = interior_phi(point)


# Velocity = grad(phi)
dphi_dy, dphi_dx = np.gradient(Phi, ys, xs)
U = dphi_dx
V = dphi_dy

# Exact solution comparison
phi_exact = X
error_phi = np.max(np.abs(Phi - phi_exact))

print("Max phi error =", error_phi)
print("Mean U =", np.mean(U))
print("Mean V =", np.mean(V))


# Plot
plt.figure(figsize=(8, 6))

contour = plt.contourf(
    X,
    Y,
    Phi,
    levels=20,
)

plt.colorbar(contour, label="Potential phi")

plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2],
)

plt.plot(
    np.r_[p1[:, 0], p1[0, 0]],
    np.r_[p1[:, 1], p1[0, 1]],
    "k-",
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("2D Laplace Potential Flow by BEM")
plt.axis("equal")
plt.tight_layout()
plt.show()
