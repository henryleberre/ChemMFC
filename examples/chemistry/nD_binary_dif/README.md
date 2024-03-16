# nD Binary Diffusion

This example case showcases chemical diffusion within a binary mixture, replicating results from A. L. Sánchez et al. [1].

$$
\frac{\partial (\rho Y_k)}{\partial t} = D\nabla \cdot (\rho \nabla Y_k) - \nabla\cdot(\rho Y_k \vec{u})
$$

## Problem Descriptions

Equation of state:

$$ Y(\rho) = \frac{\rho_a(\rho_b-\rho)}{\rho(\rho_b-\rho_a)} $$

Half & Half Domain:

$$ \frac{\rho-\rho_b}{\rho_a-\rho_b} = \frac{1}{2} \left[1-\text{erf}\left(\frac{r}{2\sqrt{Dt}}\right)\right] $$

Disk & Sphere:

$$ \frac{\rho-\rho_b}{\rho_a-\rho_b} = \frac{1}{2} \left[\text{erf}\left(\frac{r+a}{2\sqrt{Dt}}\right)-\text{erf}\left(\frac{r-a}{2\sqrt{Dt}}\right)\right] $$

## Implementation

Using the above equations, we derive analytical expressions for the density and mass fraction scalar fields as functions of space and time.

#### Line, Rectangle, and Cuboid

$$
\begin{align*}
\rho\left(r,t\right) &= \rho_b + \frac{\rho_a-\rho_b}{2} \left[1-\text{erf}\left(\frac{r}{2\sqrt{Dt}}\right)\right] \\
Y_1(r,t) &= \frac{\rho_a\left[\rho_b-\rho(r,t)\right]}{\rho(r,t)\left[\rho_b-\rho_a\right]} \\
Y_2(r,t) &= 1 - Y_1(r,t)
\end{align*}
$$

#### Line, Disk, and Sphere

$$
\begin{align*}
\rho\left(r,t\right) &= \rho_b + \frac{\rho_a-\rho_b}{2} \left[\text{erf}\left(\frac{r+a}{2\sqrt{Dt}}\right)-\text{erf}\left(\frac{r-a}{2\sqrt{Dt}}\right)\right] \\
Y_1(r,t) &= \frac{\rho_a\left[\rho_b-\rho(r,t)\right]}{\rho(r,t)\left[\rho_b-\rho_a\right]} \\
Y_2(r,t) &= 1 - Y_1(r,t)
\end{align*}
$$

## Results

## References

[1] A. L. Sánchez, M. Vera, and A. Liñán, “Exact solutions for transient mixing of two gases of different densities,” Physics of Fluids, vol. 18, no. 7, p. 078102, Jul. 2006, doi: https://doi.org/10.1063/1.2221349.
‌