2D Laplace Potential Flow Using the Boundary Element Method (BEM)

GitHub-safe version: all equations are written in plain Unicode text.
No LaTeX math delimiters are used, so the equations will not break if GitHub math rendering is unavailable.

1. Overview

This example solves a simple two-dimensional potential-flow problem using the Boundary Element Method (BEM).

The computational domain is the unit square:

0 ≤ x ≤ 1
0 ≤ y ≤ 1

The prescribed velocity potential on the boundary is:

φ(x, y) = x

The exact solution inside the domain is also:

φ(x, y) = x

The velocity field is obtained from the gradient of the potential:

u⃗ = ∇φ

Therefore,

u = ∂φ/∂x = 1
v = ∂φ/∂y = 0

Hence, the exact velocity field is:

u⃗ = (1, 0)

The purpose of this program is to demonstrate the basic BEM procedure:

Generate boundary elements.

Define the Laplace fundamental solution.

Evaluate boundary integrals using Gauss quadrature.

Assemble the BEM matrices H and G.

Solve for the unknown normal derivative q = ∂φ/∂n.

Reconstruct the potential inside the domain.

Obtain the internal velocity field from ∇φ.

2. Governing Equation

For incompressible and irrotational potential flow,

u⃗ = ∇φ

where φ is the velocity potential.

For incompressible flow,

∇ · u⃗ = 0

Substituting u⃗ = ∇φ gives:

∇ · (∇φ) = 0

Therefore,

∇²φ = 0

In two dimensions,

∂²φ/∂x² + ∂²φ/∂y² = 0

This is the Laplace equation.

3. Boundary Integral Equation

BEM converts the differential equation inside the domain into an integral equation defined only on the boundary.

For the two-dimensional Laplace equation:

c(ξ) φ(ξ)
+ ∫Γ φ(x) [∂G(x, ξ)/∂n] dΓ
= ∫Γ G(x, ξ) q(x) dΓ

where:

q = ∂φ/∂n

The symbols are:

ξ: source or collocation point

x: integration point on the boundary

Γ: boundary of the domain

G: fundamental solution

n: outward unit normal direction

q: normal derivative of the velocity potential

For a smooth two-dimensional boundary:

c(ξ) = 1/2

4. Fundamental Solution

For the 2D Laplace equation, the fundamental solution is:

G(x, ξ) = -(1 / 2π) ln(r)

where:

r = |x - ξ|

The normal derivative of the fundamental solution is:

∂G/∂n = -(1 / 2π) [(x - ξ) · n] / r²

These two kernels are used to construct the BEM matrices H and G.

5. Boundary Discretization

Unlike FEM, BEM does not mesh the entire two-dimensional domain.

Only the boundary is discretized.

For the square domain, each side is divided into N_SIDE constant boundary elements.

N_SIDE = 16

Since the square has four sides:

N_E = 4 N_SIDE

For:

N_SIDE = 16

the total number of boundary elements is:

N_E = 64

The function:

make_square_boundary(n_side)

creates the four sides in counterclockwise order.

For every boundary element, the program computes:

first endpoint p1

second endpoint p2

element midpoint mid

element length L

outward unit normal normal

For one element:

d = x₂ - x₁

and:

L_e = |d|

Because the boundary is ordered counterclockwise, the outward normal is:

n = (1 / L_e) (d_y, -d_x)ᵀ

This is the same quantity implemented in the code by:

normal = np.column_stack(
    [d[:, 1], -d[:, 0]]
) / L[:, None]

6. Constant Boundary Elements

This program uses constant boundary elements.

Within one element:

φ(x) ≈ φ_j

and:

q(x) ≈ q_j

The values φ_j and q_j are represented by values at the element midpoint.

After dividing the boundary into N_E elements, the boundary integral equation becomes:

c_i φ_i
+ Σ(j=1 to N_E) [ φ_j ∫Γ_j (∂G_i/∂n) dΓ ]
= Σ(j=1 to N_E) [ q_j ∫Γ_j G_i dΓ ]

Define:

H_ij = c_i δ_ij + ∫Γ_j (∂G_i/∂n) dΓ

and:

G_ij = ∫Γ_j G_i dΓ

Then the discretized BEM equation is:

H Φ = G Q

where:

Φ = (φ₁, φ₂, ..., φ_N_E)ᵀ
Q = (q₁, q₂, ..., q_N_E)ᵀ

This is the main algebraic system used by the BEM code.

7. Gauss Quadrature

For regular boundary elements, the required integrals are evaluated numerically using Gauss-Legendre quadrature.

The standard formula is:

∫[-1 to 1] f(z) dz ≈ Σ(k=1 to N_G) w_k f(z_k)

where:

z_k: Gauss point

w_k: Gauss weight

N_G: number of Gauss points

The code generates the Gauss points and weights with:

gp, gw = np.polynomial.legendre.leggauss(N_GAUSS)

For example:

N_GAUSS = 10

uses ten Gauss points for each regular boundary-element integral.

8. Mapping a Gauss Point to the Physical Element

Gauss quadrature is defined on the standard interval:

-1 ≤ z ≤ 1

A physical boundary element extends from point a to point b.

The mapping is:

x(z) = 0.5 [(1 - z)a + (1 + z)b]

In the code:

x = 0.5 * (
    (1 - z) * a
    +
    (1 + z) * b
)

The Jacobian is:

J = L_e / 2

Therefore:

dΓ = J dz

The physical boundary integral becomes:

∫Γ_e f(x) dΓ = ∫[-1 to 1] f(x(z)) J dz

Using Gauss quadrature:

∫Γ_e f dΓ ≈ Σ(k=1 to N_G) w_k f(z_k) J

9. Element Integration Function

The function:

integrate_element(source, a, b, nvec)

calculates two quantities:

I_G = ∫Γ_e G dΓ

and:

I_H = ∫Γ_e (∂G/∂n) dΓ

At each Gauss point, the code evaluates:

rvec = x - source
r = np.linalg.norm(rvec)

which corresponds to:

r⃗ = x - ξ
r = |r⃗|

The Green function is:

G = -(1.0 / (2 * np.pi)) * np.log(r)

which represents:

G = -(1 / 2π) ln(r)

The normal derivative is calculated as:

dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)

which represents:

∂G/∂n = -(1 / 2π) (r⃗ · n) / r²

Finally:

IG += w * G * J
IH += w * dGdn * J

implements:

I_G ≈ Σ w_k G_k J
I_H ≈ Σ w_k (∂G/∂n)_k J

10. Singular Self Element

A special case occurs when the source point lies on the same boundary element that is being integrated.

Then:

r → 0

and the Green function becomes singular.

For the H matrix, the Cauchy Principal Value treatment produces the jump term:

H_ii = 1/2

for a smooth boundary.

The code uses:

H[i, j] = 0.5

when:

i == j

For the logarithmic singularity in the G matrix, the self-element integral for a constant straight element can be evaluated analytically:

G_ii = [L_i / (2π)] [1 - ln(L_i / 2)]

The code implements this as:

Gmat[i, j] = (
    L[j]
    / (2 * np.pi)
    * (
        1.0
        - np.log(L[j] / 2.0)
    )
)

This avoids applying ordinary Gauss quadrature directly at the singular point.

11. Assembly of the Global BEM Matrices

The program loops over all source elements i and all integration elements j.

for i in range(NE):
    source = mid[i]

    for j in range(NE):

Each pair (i, j) generates one matrix coefficient.

The matrices can be written symbolically as:

H = [H_ij]
G = [G_ij]

with dimensions:

N_E × N_E

Unlike a typical FEM stiffness matrix, BEM matrices are generally dense.

This is because every boundary element influences every collocation point through the Green function.

12. Boundary Condition

The prescribed boundary condition is:

φ = x

Since constant boundary elements are used, the prescribed value is evaluated at each element midpoint:

phi_boundary = mid[:, 0]

Therefore:

φ_j = x_j

for every boundary element.

The unknown quantity is:

q_j = ∂φ/∂n

The BEM equation is:

H Φ = G Q

Formally:

Q = G⁻¹ H Φ

However, the code does not explicitly compute the inverse.

Instead, it solves the linear system directly:

q_boundary = np.linalg.solve(
    Gmat,
    H @ phi_boundary
)

Using np.linalg.solve is numerically preferable to explicitly forming G⁻¹.

13. Exact Boundary Flux

Since:

φ = x

the exact gradient is:

∇φ = (1, 0)ᵀ

The normal derivative is:

q = ∇φ · n

Therefore:

q = n_x

This is implemented as:

q_exact = normal[:, 0]

Physically:

left wall: q = -1

right wall: q = +1

bottom wall: q = 0

top wall: q = 0

The RMS boundary error is:

e_q = sqrt[ (1 / N_E) Σ(j=1 to N_E) (q_j - q_j,exact)² ]

The code calculates it as:

error_q = np.sqrt(
    np.mean(
        (q_boundary - q_exact)**2
    )
)

14. Interior Potential Reconstruction

An important feature of BEM is that the interior does not need to be meshed.

Once the boundary values φ and q are known, the potential at any interior point ξ can be reconstructed from the boundary solution:

φ(ξ) = ∫Γ [ G q - φ (∂G/∂n) ] dΓ

For an interior field point, the coefficient multiplying φ(ξ) is 1.

The code evaluates the boundary contribution as:

value += (
    IG * q_boundary[j]
    -
    IH * phi_boundary[j]
)

for every boundary element.

Therefore, the interior grid contains evaluation points only.

They are not additional unknown BEM nodes.

15. Interior Grid

For visualization, the program creates a set of interior field points:

xs = np.linspace(0.05, 0.95, 25)
ys = np.linspace(0.05, 0.95, 25)

The endpoints are kept slightly away from the boundary to reduce near-singular integration error in this simple example.

A regular Cartesian grid is formed:

X, Y = np.meshgrid(xs, ys)

At each point:

Phi[iy, ix] = interior_phi(point)

evaluates the boundary integral representation.

The array Phi therefore contains:

φ(x_i, y_j)

throughout the interior.

16. Velocity Field

For potential flow:

u⃗ = ∇φ

Therefore:

u = ∂φ/∂x
v = ∂φ/∂y

The code obtains these derivatives numerically:

dphi_dy, dphi_dx = np.gradient(
    Phi,
    ys,
    xs
)

and then:

U = dphi_dx
V = dphi_dy

Hence:

U ≈ ∂φ/∂x
V ≈ ∂φ/∂y

For this example, the expected values are:

U ≈ 1
V ≈ 0

throughout the domain.

17. Exact Interior Solution

The exact potential is:

φ_exact = x

In the code:

phi_exact = X

The maximum absolute potential error is:

e_φ = max |φ_BEM - φ_exact|

It is calculated as:

error_phi = np.max(
    np.abs(
        Phi - phi_exact
    )
)

The program also reports:

np.mean(U)
np.mean(V)

which should approach:

mean(U) → 1
mean(V) → 0

as the discretization becomes sufficiently accurate.

18. Visualization

The potential field is drawn using:

plt.contourf(
    X,
    Y,
    Phi,
    levels=20
)

This gives a contour map of the potential φ(x, y).

The velocity field is drawn using:

plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2]
)

The expected velocity vectors are nearly horizontal and point in the positive x direction:

→ → → → → → →
→ → → → → → →
→ → → → → → →
→ → → → → → →

19. Complete Computational Flow

The complete BEM calculation can be summarized as:

1. Governing equation

   ∇²φ = 0

          ↓

2. Boundary integral equation

   c φ + ∫Γ φ (∂G/∂n) dΓ = ∫Γ G q dΓ

          ↓

3. Boundary discretization

   Γ → Γ₁ + Γ₂ + ... + Γ_N

          ↓

4. Gauss quadrature

   ∫Γ_e (...) dΓ ≈ Σ w_k (...)_k J

          ↓

5. Matrix assembly

   H Φ = G Q

          ↓

6. Solve boundary unknowns

   Q = G⁻¹ H Φ

          ↓

7. Interior reconstruction

   φ(ξ) = ∫Γ [Gq - φ(∂G/∂n)] dΓ

          ↓

8. Velocity field

   u⃗ = ∇φ

20. Why This Example Is Useful

This problem has an exact solution, so it is useful for checking whether the BEM implementation is correct.

It demonstrates the essential BEM ideas:

only the boundary is discretized

the Green function transfers boundary information into the domain

Gauss quadrature evaluates regular boundary integrals

singular self-integrals require special treatment

BEM matrices are generally dense

interior values are reconstructed after the boundary solution is known

21. Important Numerical Notes

21.1 Python line continuation

Avoid using a backslash for line continuation when it may contain trailing spaces.

Instead of a fragile form, use parentheses:

dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)

21.2 Near-singular integration

The present example treats exactly singular self-elements analytically.

However, if an interior evaluation point lies extremely close to a boundary element, ordinary Gauss quadrature can lose accuracy even though the integral is not mathematically singular.

More advanced BEM programs may use:

adaptive quadrature

element subdivision

special near-singular integration

singularity subtraction

21.3 Constant elements

The current program uses constant elements because they are easy to understand.

More accurate BEM formulations may use linear or quadratic elements.

For example:

φ(ξ) = Σ_j N_j(ξ) φ_j
q(ξ) = Σ_j N_j(ξ) q_j

where N_j is the corresponding shape function.

22. Requirements

The script requires:

numpy
matplotlib

Install them with:

pip install numpy matplotlib

23. Running the Program

Save the Python code as:

bem_laplace_flow.py

and run:

python bem_laplace_flow.py

The program prints:

boundary-flux RMS error

maximum potential error

mean x-direction velocity

mean y-direction velocity

It also displays the potential contour and velocity vectors.

24. Expected Physical Result

Because:

φ = x

the velocity is:

u⃗ = (1, 0)

Therefore:

potential contours should be approximately vertical

velocity arrows should point horizontally to the right

U should be approximately 1

V should be approximately 0

25. Main Concept

The most important idea of this example is:

Solve only on the boundary
        ↓
Reconstruct the solution anywhere inside

This is the defining computational idea behind the Boundary Element Method for Laplace-type problems.

26. Original Python Program

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
    for k in range(n_side):
        pts.append([k / n_side, 0.0])

    # right
    for k in range(n_side):
        pts.append([1.0, k / n_side])

    # top
    for k in range(n_side):
        pts.append([1.0 - k / n_side, 1.0])

    # left
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

    # CCW boundary:
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
        dGdn = (
            -(1.0 / (2 * np.pi))
            * np.dot(rvec, nvec)
            / (r * r)
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

        # singular self element
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

        # regular element
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
# 5. Boundary condition: phi = x
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
# q = grad(phi) dot n = n_x
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

plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2]
)

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
