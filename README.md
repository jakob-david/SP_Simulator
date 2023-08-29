# SP_Simulator

## Numerical solution of the Schrödinger-Poisson equation 

### Schrödinger Poisson equation
The normal Schrödinger equation describes how a wave function of a particle behaves with time and is given in equation 1 where m is the mass and ℏ the reduced Planck constant. (Griffiths & Schroeter, 2018)

```math
i\hbar \frac {\partial u}{\partial t} =-\frac{\hbar^2}{2m}\frac{\partial^2u}{\partial x^2}+Vu
```

The equation only gets a physical meaning if it is squared and integrated as shown in equation 2. After that it can be interpreted as a probability to find the particle at a given point in space. (Griffiths & Schroeter, 2018)

```math
\int_{a}^{b}{\left|u(x,t)\right|^2dx}=p
```
Where $p$ is the probability of finding the particle between $a$ and $b$, at time $t$.

In general, there are three cases of this equation. The first is the free case with V=0, the second is the linear case is V=V(x) where the potential is only dependent on space and the non-linear case with V=V(x,u) (Mauser & Stimming, 2021). Since in the last case the potential is also dependent on the solution the whole equation becomes non-linear. The last non-linear case will not be discussed here. However, the written code can treat such potentials as well. 

The Schrödinger-Poisson equation modifies equation 1 by adding a Newtonian gravitational potential. Hereby, the probability density is interpreted as a mass density. The full PDE is given in equation 3 and 4. Such a system can also be used to describe things like exotic dark matter. (Mocz, 2023)

```math
i\hbar\frac{\partial u}{\partial t}=-\frac{\hbar^2}{2m}\frac{\partial^2u}{\partial x^2}+Vu+m\Phi u
```
```math
\Phi=\ \frac{1}{\nabla^2}4\pi Gm\left|u\right|^2
```

### Scaling 

Since the focus of this project lies only on solving the system numerically, the equations will be scaled by introducing a new dimension-less parameter $\varepsilon=\frac{\hbar}{m}$ and define a new potential with $V_{new}=\ \frac{V_{old}\ }{m}$. Furthermore, for the sake of simplicity, it will be assumed that  $\varepsilon\defeq1$ and $4\pi G\defeq1$. The new system is shown in equation 5 (Mauser & Stimming, 2021) and 6. 

```math
i\frac{\partial u}{\partial t}=-\frac{1}{2}\frac{\partial^2u}{\partial x^2}+Vu+\Phi u
```
```math
\Phi=\ \frac{1}{\nabla^2}\left|u\right|^2
```