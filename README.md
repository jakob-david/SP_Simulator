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

Since the focus of this project lies only on solving the system numerically, the equations will be scaled by introducing a new dimension-less parameter $\varepsilon=\frac{\hbar}{m}$ and define a new potential with $V_{new}=\ \frac{V_{old}\ }{m}$. Furthermore, for the sake of simplicity, it will be assumed that  $\varepsilon=1$ and $4\pi G=1$. The new system is shown in equation 5 (Mauser & Stimming, 2021) and 6. 

```math
i\frac{\partial u}{\partial t}=-\frac{1}{2}\frac{\partial^2u}{\partial x^2}+Vu+\Phi u
```
```math
\Phi=\ \frac{1}{\nabla^2}\left|u\right|^2
```

### Split-Operator Method 

The chosen method for solving the Schrödinger-Poisson equation numerically is a so-called Split-Operator Method. At first, it will only be introduced for the scaled Schrödinger equation which is shown in equation 7 (Mauser & Stimming, 2021). Later, it will be adapted to simulate the whole Schrödinger-Poisson system.

```math
i\frac{\partial u}{\partial t}=-\frac{1}{2}\frac{\partial^2u}{\partial x^2}+Vu
```
In a first step the Hamiltonian of the system will be split into a component in position space ${\hat{H}}_x=V$ and into a component in momentum space ${\hat{H}}_k=-\frac{1}{2}\frac{\partial^2}{\partial x^2}$. Next, a somewhat general solution to the system is assumed which is displayed in equation 8. It is also presumed that the system is simulated by concatenating a series of small timesteps $∆t$. (Schloss, 2022)

```math
u(x,t+∆t)=[e^{-i\hat{H}∆t}] u\left(x,t\right)=[e^{i(\hat{H}_x+\hat{H}_k)∆t}]u(x,t)
```

Using the Baker-Campbell-Housdorff theorem this solution can be split. The thus retrieved formula is shown in equation 9. (Schloss, 2022)

```math
u(x,t+∆t)=[e^{-i\hat{H}_x∆t}e^{-i\hat{H}_k∆t}e^{[i\hat{H}_x,i\hat{H}_k]∆t^2}]u(x,t)
```

Since the cross-term scales with ∆t2 in contrast to ∆t its contribution will be small and can therefore be ignored. Because of this, the overall error of this method is $Ο(∆t^{2})$ (Schloss, 2022).

Next, it can be observed that the part in position space and the part in momentum space can be dealt separately. The part in position space can be calculated relatively easily since it is just a simple multiplication. The part in momentum space is a little bit more difficult since it contains the operator $\frac{\partial^2}{\partial x^2}$. To circumvent this problem, the Fourier Transform is used, because using it will change the differential operator into a simple multiplication. The process is shown in equation 10, where k denotes the wave number according to the current position $x$ and $\mathcal{F}$ denotes the Fourier Transform. (Schloss, 2022)

```math
\frac{\partial^2u_x}{\partial x^2}=\mathcal{F}^{-1}\left[\frac{{k_x}^2}{2}\mathcal{F}\left[u_x\right]\right]
```

Because, the Fourier Transform is computationally inefficient, the Fast Fourier Transform (FFT) will be used. One complete timestep of the algorithm is displayed in equation 11, where $\mathcal{F}$ now denotes the FFT, ${\hat{U}}_x=e^{-i\hat{H}_x∆t}$ and ${\hat{U}}_k=e^{-i\hat{H}_k∆t}$. (Schloss, 2022)

```math
u(x,t+∆t)=[\hat{U}_k(∆t)\mathcal{F}[\hat{U}_x(∆t)u(x,t)]] +O(∆t^{2})
```

Now the algorithm for simulation the Schrödinger Equation is finished and can be easily adapted. Sine the Newtonian term also uses a differential operator, the FFT is used as well as is shown below (Johansson, 2010)

```math
V_{int}=\frac{1}{\nabla^2}\left|u_x\right|^2=\mathcal{F}^{-1}\left[\frac{1}{{k_x}^2}\mathcal{F}\left[\left|u_x\right|^2\right]\right]
```

Due to the fact that a division by $0$ is not allowed, one must choose a sufficiently small number greater than $0$ in the case that the wave number is $0$. The result of equation 12 can be thought of as a somewhat internal potential due to the Newtonian gravitational force. The routine in equation 11 will basically stay the same. The only thing that changes is the operator ${\hat{H}}_x$ which now will be $\hat{H}_x=V+V_int$ instead of just $V$ (Johansson, 2010).

