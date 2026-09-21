# 2D Laplace Potential Flow Using the Boundary Element Method (BEM)

## 1. Overview

This example solves a simple two-dimensional potential-flow problem with the **Boundary Element Method (BEM)**.

The computational domain is a unit square:

\[
0 \le x \le 1, \qquad 0 \le y \le 1
\]

The prescribed velocity potential on the boundary is

\[
\phi(x,y)=x
\]

The exact solution inside the domain is also

\[
\phi(x,y)=x
\]

Therefore, the velocity field is

\[
\mathbf{u}=\nabla\phi
\]

or

\[
u=\frac{\partial\phi}{\partial x}=1,
\qquad
v=\frac{\partial\phi}{\partial y}=0
\]

Thus, the exact flow is a uniform flow from left to right:

\[
\boxed{\mathbf{u}=(1,0)}
\]

The purpose of this program is not to solve a difficult flow problem, but to demonstrate the basic BEM procedure:

1. Generate boundary elements.
2. Define the Laplace fundamental solution.
3. Evaluate boundary integrals using Gauss quadrature.
4. Assemble the BEM matrices \(H\) and \(G\).
5. Solve for the unknown normal derivative \(q=\partial\phi/\partial n\).
6. Reconstruct the potential inside the domain.
7. Obtain the internal velocity field from \(\nabla\phi\).

---

# 2. Governing Equation

For incompressible and irrotational potential flow,

\[
\mathbf{u}=\nabla\phi
\]

where \(\phi\) is the velocity potential.

For incompressible flow,

\[
\nabla\cdot\mathbf{u}=0
\]

Substituting

\[
\mathbf{u}=\nabla\phi
\]

gives

\[
\nabla\cdot(\nabla\phi)=0
\]

and therefore

\[
\boxed{\nabla^2\phi=0}
\]

In two dimensions,

\[
\boxed{
\frac{\partial^2\phi}{\partial x^2}
+
\frac{\partial^2\phi}{\partial y^2}
=0
}
\]

This is the **Laplace equation**.

---

# 3. Boundary Integral Equation

BEM transforms the differential equation in the domain into an integral equation defined only on the boundary.

For the two-dimensional Laplace equation,

\[
c(\boldsymbol{\xi})\phi(\boldsymbol{\xi})
+
\int_\Gamma
\phi(\mathbf{x})
\frac{\partial G(\mathbf{x},\boldsymbol{\xi})}{\partial n}
\,d\Gamma
=
\int_\Gamma
G(\mathbf{x},\boldsymbol{\xi})
q(\mathbf{x})
\,d\Gamma
\]

where

\[
q=\frac{\partial\phi}{\partial n}
\]

and:

- \(\boldsymbol{\xi}\): source or collocation point,
- \(\mathbf{x}\): integration point on the boundary,
- \(\Gamma\): boundary of the domain,
- \(G\): fundamental solution,
- \(n\): outward unit normal vector.

For a smooth two-dimensional boundary,

\[
\boxed{c(\boldsymbol{\xi})=\frac12}
\]

---

# 4. Fundamental Solution

For the 2D Laplace equation, the fundamental solution is

\[
\boxed{
G(\mathbf{x},\boldsymbol{\xi})
=
-\frac{1}{2\pi}\ln r
}
\]

where

\[
r=
\left|
\mathbf{x}-\boldsymbol{\xi}
\right|
\]

The normal derivative of the fundamental solution is

\[
\boxed{
\frac{\partial G}{\partial n}
=
-\frac{1}{2\pi}
\frac{
(\mathbf{x}-\boldsymbol{\xi})\cdot\mathbf{n}
}{
r^2
}
}
\]

These two kernels are used to construct the BEM matrices.

---

# 5. Boundary Discretization

Unlike FEM, BEM does not mesh the entire two-dimensional domain.

Only the boundary is discretized.

For the square domain, each side is divided into `N_SIDE` constant boundary elements.

```python
N_SIDE = 16
```

Since the square has four sides,

\[
N_E=4N_{\text{SIDE}}
\]

Therefore, with

\[
N_{\text{SIDE}}=16
\]

the total number of boundary elements is

\[
\boxed{N_E=64}
\]

The function

```python
make_square_boundary(n_side)
```

creates the four boundary sides in counterclockwise order.

For every element it computes:

- first endpoint `p1`,
- second endpoint `p2`,
- element midpoint `mid`,
- element length `L`,
- outward unit normal `normal`.

For one element,

\[
\mathbf{d}
=
\mathbf{x}_2-\mathbf{x}_1
\]

and

\[
L_e=|\mathbf{d}|
\]

Since the boundary is ordered counterclockwise, the outward normal is

\[
\boxed{
\mathbf{n}
=
\frac{1}{L_e}
\begin{bmatrix}
d_y\\
-d_x
\end{bmatrix}
}
\]

---

# 6. Constant Boundary Elements

This program uses **constant boundary elements**.

This means that within one element,

\[
\phi(\mathbf{x})\approx\phi_j
\]

and

\[
q(\mathbf{x})\approx q_j
\]

where \(\phi_j\) and \(q_j\) are represented by values at the element midpoint.

The continuous boundary integral equation is therefore discretized as

\[
c_i\phi_i
+
\sum_{j=1}^{N_E}
\phi_j
\int_{\Gamma_j}
\frac{\partial G_i}{\partial n}
\,d\Gamma
=
\sum_{j=1}^{N_E}
q_j
\int_{\Gamma_j}
G_i
\,d\Gamma
\]

Define

\[
H_{ij}
=
c_i\delta_{ij}
+
\int_{\Gamma_j}
\frac{\partial G_i}{\partial n}
\,d\Gamma
\]

and

\[
G_{ij}
=
\int_{\Gamma_j}
G_i
\,d\Gamma
\]

Then the discretized BEM equation becomes

\[
\boxed{
H\Phi=GQ
}
\]

where

\[
\Phi=
\begin{bmatrix}
\phi_1\\
\phi_2\\
\vdots\\
\phi_{N_E}
\end{bmatrix}
\]

and

\[
Q=
\begin{bmatrix}
q_1\\
q_2\\
\vdots\\
q_{N_E}
\end{bmatrix}
\]

---

# 7. Gauss Quadrature

For regular boundary elements, the required integrals are evaluated numerically using Gauss-Legendre quadrature.

The standard Gauss integration formula is

\[
\boxed{
\int_{-1}^{1}f(z)\,dz
\approx
\sum_{k=1}^{N_G}
w_k f(z_k)
}
\]

where

- \(z_k\): Gauss points,
- \(w_k\): Gauss weights,
- \(N_G\): number of Gauss points.

The code generates them with

```python
gp, gw = np.polynomial.legendre.leggauss(N_GAUSS)
```

For example,

```python
N_GAUSS = 10
```

uses ten Gauss points for each regular element integral.

---

# 8. Mapping the Gauss Point to the Physical Element

Gauss quadrature is defined on

\[
-1\le z\le1
\]

but a physical boundary element extends from point \(\mathbf{a}\) to \(\mathbf{b}\).

The mapping is

\[
\boxed{
\mathbf{x}(z)
=
\frac12
\left[
(1-z)\mathbf{a}
+
(1+z)\mathbf{b}
\right]
}
\]

which appears in the code as

```python
x = 0.5 * (
    (1 - z) * a
    +
    (1 + z) * b
)
```

The Jacobian of this transformation is

\[
\boxed{
J=\frac{L_e}{2}
}
\]

and therefore

\[
d\Gamma=J\,dz
\]

The physical boundary integral becomes

\[
\int_{\Gamma_e}
f(\mathbf{x})\,d\Gamma
=
\int_{-1}^{1}
f(\mathbf{x}(z))J\,dz
\]

and Gauss quadrature gives

\[
\boxed{
\int_{\Gamma_e}
f\,d\Gamma
\approx
\sum_{k=1}^{N_G}
w_k f(z_k)J
}
\]

---

# 9. Element Integration Function

The function

```python
integrate_element(source, a, b, nvec)
```

calculates two quantities:

\[
I_G
=
\int_{\Gamma_e}G\,d\Gamma
\]

and

\[
I_H
=
\int_{\Gamma_e}
\frac{\partial G}{\partial n}
\,d\Gamma
\]

For each Gauss point,

```python
rvec = x - source
r = np.linalg.norm(rvec)
```

computes

\[
\mathbf{r}
=
\mathbf{x}-\boldsymbol{\xi}
\]

and

\[
r=|\mathbf{r}|
\]

The fundamental solution is then evaluated:

```python
G = -(1.0 / (2 * np.pi)) * np.log(r)
```

corresponding to

\[
G=-\frac{1}{2\pi}\ln r
\]

The normal derivative is

```python
dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)
```

corresponding to

\[
\frac{\partial G}{\partial n}
=
-\frac{1}{2\pi}
\frac{\mathbf{r}\cdot\mathbf{n}}{r^2}
\]

Finally,

```python
IG += w * G * J
IH += w * dGdn * J
```

implements

\[
I_G
\approx
\sum_k w_k G_k J
\]

and

\[
I_H
\approx
\sum_k
w_k
\left(
\frac{\partial G}{\partial n}
\right)_k
J
\]

---

# 10. Singular Self Element

A special case occurs when the source point lies on the same boundary element being integrated.

Then

\[
r\rightarrow0
\]

and the Green function becomes singular.

For the \(H\) matrix, the Cauchy Principal Value treatment produces the jump term

\[
\boxed{
H_{ii}=\frac12
}
\]

for a smooth boundary.

The code uses

```python
H[i, j] = 0.5
```

when

```python
i == j
```

For the logarithmic singularity in the \(G\) matrix, the self-element integral can be evaluated analytically for a constant straight element:

\[
\boxed{
G_{ii}
=
\frac{L_i}{2\pi}
\left[
1-\ln\left(\frac{L_i}{2}\right)
\right]
}
\]

which is implemented as

```python
Gmat[i, j] = (
    L[j]
    / (2 * np.pi)
    * (
        1.0
        - np.log(L[j] / 2.0)
    )
)
```

This avoids attempting ordinary Gauss quadrature directly at the singularity.

---

# 11. Assembly of the Global BEM Matrices

The program loops over all source elements \(i\) and integration elements \(j\).

```python
for i in range(NE):
    source = mid[i]

    for j in range(NE):
```

Each pair \((i,j)\) generates one matrix coefficient.

Thus,

\[
H=
\begin{bmatrix}
H_{11}&H_{12}&\cdots&H_{1N}\\
H_{21}&H_{22}&\cdots&H_{2N}\\
\vdots&\vdots&\ddots&\vdots\\
H_{N1}&H_{N2}&\cdots&H_{NN}
\end{bmatrix}
\]

and similarly for \(G\).

Unlike a typical FEM stiffness matrix, BEM matrices are generally **dense**, because every boundary element influences every collocation point through the Green function.

---

# 12. Boundary Condition

The prescribed boundary condition is

\[
\boxed{\phi=x}
\]

Since constant boundary elements are used, the prescribed value is evaluated at each element midpoint:

```python
phi_boundary = mid[:, 0]
```

Therefore,

\[
\phi_j=x_j
\]

for every boundary element.

The unknown is

\[
q_j=
\frac{\partial\phi}{\partial n}
\]

The BEM equation is

\[
H\Phi=GQ
\]

so

\[
\boxed{
Q=G^{-1}H\Phi
}
\]

The program solves this as

```python
q_boundary = np.linalg.solve(
    Gmat,
    H @ phi_boundary
)
```

Using `np.linalg.solve` is preferred over explicitly computing the matrix inverse.

---

# 13. Exact Boundary Flux

Since

\[
\phi=x
\]

the exact gradient is

\[
\nabla\phi=
\begin{bmatrix}
1\\
0
\end{bmatrix}
\]

The normal derivative is

\[
q=
\nabla\phi\cdot\mathbf{n}
\]

Therefore,

\[
\boxed{
q=n_x
}
\]

This is implemented as

```python
q_exact = normal[:, 0]
```

Physically:

- left wall: \(q=-1\),
- right wall: \(q=+1\),
- bottom wall: \(q=0\),
- top wall: \(q=0\).

The RMS boundary error is

\[
\boxed{
e_q
=
\sqrt{
\frac{1}{N_E}
\sum_{j=1}^{N_E}
(q_j-q_j^{exact})^2
}
}
\]

implemented by

```python
error_q = np.sqrt(
    np.mean(
        (q_boundary - q_exact)**2
    )
)
```

---

# 14. Interior Potential Reconstruction

An important feature of BEM is that the interior does not need to be meshed.

Once the boundary values \(\phi\) and \(q\) are known, the potential at any interior point \(\boldsymbol{\xi}\) can be reconstructed using

\[
\boxed{
\phi(\boldsymbol{\xi})
=
\int_\Gamma
\left[
Gq
-
\phi\frac{\partial G}{\partial n}
\right]
d\Gamma
}
\]

Because the point is inside the domain rather than on the boundary, the coefficient \(c\) is effectively 1 in this representation.

The code evaluates this as

```python
value += (
    IG * q_boundary[j]
    -
    IH * phi_boundary[j]
)
```

for every boundary element.

This is why BEM can display an internal potential field even though only the boundary was discretized.

The internal points are **evaluation points**, not additional unknown BEM nodes.

---

# 15. Interior Grid

For visualization, the program creates a set of interior field points:

```python
xs = np.linspace(0.05, 0.95, 25)
ys = np.linspace(0.05, 0.95, 25)
```

The endpoints are kept slightly away from the boundary to avoid near-singular integration difficulties during this simple demonstration.

A regular Cartesian grid is formed:

```python
X, Y = np.meshgrid(xs, ys)
```

For each point,

```python
Phi[iy, ix] = interior_phi(point)
```

evaluates the boundary integral representation.

The resulting matrix `Phi` contains

\[
\phi(x_i,y_j)
\]

throughout the interior.

---

# 16. Velocity Field

For potential flow,

\[
\boxed{
\mathbf{u}=\nabla\phi
}
\]

Therefore,

\[
u=
\frac{\partial\phi}{\partial x}
\]

and

\[
v=
\frac{\partial\phi}{\partial y}
\]

The code obtains these derivatives numerically:

```python
dphi_dy, dphi_dx = np.gradient(
    Phi,
    ys,
    xs
)
```

and then

```python
U = dphi_dx
V = dphi_dy
```

Hence,

\[
\boxed{
U\approx\frac{\partial\phi}{\partial x}
}
\]

and

\[
\boxed{
V\approx\frac{\partial\phi}{\partial y}
}
\]

For this example, the expected values are

\[
U\approx1
\]

and

\[
V\approx0
\]

throughout the domain.

---

# 17. Exact Interior Solution

The exact potential is

\[
\phi^{exact}=x
\]

which is represented by

```python
phi_exact = X
```

The maximum absolute error is

\[
\boxed{
e_\phi
=
\max
\left|
\phi^{BEM}-\phi^{exact}
\right|
}
\]

and is calculated with

```python
error_phi = np.max(
    np.abs(
        Phi - phi_exact
    )
)
```

The program also reports

```python
np.mean(U)
```

and

```python
np.mean(V)
```

which should approach

\[
\boxed{
\overline{U}\rightarrow1,
\qquad
\overline{V}\rightarrow0
}
\]

as the discretization becomes sufficiently accurate.

---

# 18. Visualization

The potential field is drawn using

```python
plt.contourf(
    X,
    Y,
    Phi,
    levels=20
)
```

This gives a contour map of

\[
\phi(x,y)
\]

The velocity field is drawn using

```python
plt.quiver(
    X[::2, ::2],
    Y[::2, ::2],
    U[::2, ::2],
    V[::2, ::2]
)
```

The expected result is a set of nearly horizontal arrows pointing in the positive \(x\)-direction.

Conceptually,

```text
→ → → → → → →

→ → → → → → →

→ → → → → → →

→ → → → → → →
```

---

# 19. Complete Computational Flow

The complete BEM calculation can be summarized as

\[
\boxed{
\nabla^2\phi=0
}
\]

\[
\Downarrow
\]

Green's identity

\[
\Downarrow
\]

\[
c\phi
+
\int_\Gamma
\phi\frac{\partial G}{\partial n}
\,d\Gamma
=
\int_\Gamma
Gq\,d\Gamma
\]

\[
\Downarrow
\]

Boundary discretization

\[
\Gamma
=
\sum_e \Gamma_e
\]

\[
\Downarrow
\]

Gauss quadrature

\[
\int_{\Gamma_e}(\cdots)d\Gamma
\approx
\sum_k
w_k(\cdots)_kJ
\]

\[
\Downarrow
\]

Matrix assembly

\[
\boxed{
H\Phi=GQ
}
\]

\[
\Downarrow
\]

Solve boundary unknowns

\[
\boxed{
Q=G^{-1}H\Phi
}
\]

\[
\Downarrow
\]

Interior reconstruction

\[
\phi(\boldsymbol{\xi})
=
\int_\Gamma
\left[
Gq
-
\phi\frac{\partial G}{\partial n}
\right]
d\Gamma
\]

\[
\Downarrow
\]

Velocity

\[
\boxed{
\mathbf{u}=\nabla\phi
}
\]

---

# 20. Why This Example Is Useful

This simple problem has an exact solution, so it is useful for checking whether the BEM implementation is correct.

It demonstrates the essential BEM ideas without additional complexity:

- only the boundary is discretized,
- the Green function transfers boundary information into the domain,
- Gauss quadrature evaluates regular boundary integrals,
- singular self-integrals require special treatment,
- the resulting BEM matrices are dense,
- interior values are reconstructed after the boundary solution is known.

---

# 21. Important Numerical Notes

## 21.1 Backslash syntax in Python

Avoid writing a continuation backslash with trailing spaces, such as

```python
dGdn = -(1.0 / (2 * np.pi)) * \ 
        np.dot(rvec, nvec) / (r * r)
```

because Python may generate a syntax error.

A safer form is

```python
dGdn = (
    -(1.0 / (2 * np.pi))
    * np.dot(rvec, nvec)
    / (r * r)
)
```

---

## 21.2 Near-singular integration

The present example treats exactly singular self-elements analytically.

However, if an interior evaluation point lies extremely close to a boundary element, ordinary Gauss quadrature can also lose accuracy even though the integral is not mathematically singular.

More advanced BEM programs may use:

- adaptive quadrature,
- element subdivision,
- special near-singular integration,
- singularity subtraction.

---

## 21.3 Constant elements

The current program uses constant elements because they are easy to understand.

More accurate BEM formulations may use linear or quadratic boundary elements:

\[
\phi(\xi)
=
\sum_jN_j(\xi)\phi_j
\]

and

\[
q(\xi)
=
\sum_jN_j(\xi)q_j
\]

with shape functions \(N_j\).

---

# 22. Requirements

The script requires

```text
numpy
matplotlib
```

Install them with

```bash
pip install numpy matplotlib
```

---

# 23. Running the Program

Save the Python code as, for example,

```text
bem_laplace_flow.py
```

and run

```bash
python bem_laplace_flow.py
```

The program prints:

- boundary-flux RMS error,
- maximum potential error,
- mean \(u\)-velocity,
- mean \(v\)-velocity.

It also displays the potential contour and velocity vectors.

---

# 24. Expected Physical Result

Because

\[
\phi=x
\]

the velocity should be

\[
\boxed{
\mathbf{u}=(1,0)
}
\]

Therefore:

- potential contours should be approximately vertical,
- velocity arrows should point horizontally to the right,
- \(U\) should be approximately 1,
- \(V\) should be approximately 0.

---

# 25. Main Concept

The most important idea of the example is

\[
\boxed{
\text{Solve only on the boundary}
\quad\Longrightarrow\quad
\text{reconstruct the solution anywhere inside}
}
\]

This is the defining computational idea behind the Boundary Element Method for Laplace-type problems.

