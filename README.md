# SP_Simulator

## Numerical solution of the Schrödinger-Poisson equation 

### Schrödinger Poisson equation
The normal Schrödinger equation describes how a wave function of a particle behaves with time and is given in equation 1 where m is the mass and ℏ the reduced Planck constant. (Griffiths & Schroeter, 2018)

$i\hbar \frac {\partial u}{\partial t} =-\frac{\hbar^2}{2m}\frac{\partial^2u}{\partial x^2}+Vu$

The equation only gets a physical meaning if it is squared and integrated as shown in equation 2. After that it can be interpreted as a probability to find the particle at a given point in space. (Griffiths & Schroeter, 2018)

\begin{align}
    g &= \int_a^b f(x)dx \label{eq1}\tag{1} \\
    a &= b + c \label{eq2}\tag{2}
\end{align}