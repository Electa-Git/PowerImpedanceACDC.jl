# Initialization

As the ABCD formulation is a linear representation of the power system,
nonlinear descriptions of the components such as power converters must
be linearized around an operating point. This operating point is
determined in the initialisation by solving the power flow equations
representing the combined AC/DC system. This initialization is done automatically within the package using the [`PowerImpedanceACDC.power_flow`](@ref) function.

```@docs
PowerImpedanceACDC.power_flow
```
To do so, the network is initialized using the optimal power flow tool
[^1] implemented as a Julia package, which can be found
in the [PowerModelsACDC](https://github.com/Electa-Git/PowerModelsACDC.jl) repository. The
package relies on the power flow models developed for the MatACDC
simulator [^2], which extends Matpower
[^3] AC power system models with the DC
representations and with power converters.

As a result, the constructed power system is divided into AC and DC
systems and the converters. It contains AC and DC branches and buses,
converters, generators, loads, shunts and storage elements. Components
implemented in this electromagnetic stability simulator are represented
using their equivalent models for the purpose of the power flow
analysis.

## AC and DC branches

AC and DC branches represent three-phase AC and DC connections between
buses respectively. Branches are grouped inside AC or DC grids (zones).
AC branches are defined with parameters described in
[^3], while DC branches parameters are given in
[^2].

For the purpose of modeling the system components, the model of the AC
branch as provided in [^3] is depicted in the following figure. Beside the shunt admittance
$j\frac{b_c}{2}$,the Julia package [PowerModelsACDC](https://github.com/Electa-Git/PowerModelsACDC.jl) supports the admittance as
$\frac{g_c}{2} + j\frac{b_c}{2}$. The full expression for the AC
admittance parameters is given by the equation:
$$ \textbf{Y}_{ac} = \begin{bmatrix}
    \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right) \, \frac{1}{\tau^2} & -\frac{y_s}{\tau \, \exp(-j\theta_{\rm shift})} \\
    -\frac{y_s}{\tau \, \exp(-j\theta_{\rm shift})} & \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right) 
    \end{bmatrix}.
$$ It should be noted that the AC network is
considered as balanced for the power flow solution, and thus, all the
components have diagonal matrix models, with equal elements on the
matrix diagonal. 

![Matpower AC branch model
[^3].](pictures/power_flow/branch_ac.png)

A DC branch is modeled with its equivalent series resistance
[^2]. For the power flow calculation, some components
are modeled as AC and DC branches and their models are described in
detail in this subsection.

### Impedance model

As described, an [Impedance](@ref) is modeled using ABCD
parameters: $$\begin{bmatrix}
    \textbf{A} & \textbf{B} \\
    \textbf{C} & \textbf{D}
    \end{bmatrix} = 
    \begin{bmatrix}
    \textbf{I} & \textbf{Z} \\
    \textbf{0} & \textbf{I}
    \end{bmatrix}.$$ In the case of the DC impedance, all matrices are
of size $1 \times 1$, while three-phase impedances are of size
$3 \times 3$.

The DC branch model is then given as a Thevenin equivalent series
impedance $r = \Re\{Z\}$. The AC branch is modeled as an ideal
transformer with $\tau = 1$ and $\theta_{\rm shift} = 0$, and with
$r_s = \Re\{\textbf{Z}(j\omega)\left\langle 1,1 \right\rangle\}$,
$x_s = \Im\{\textbf{Z}(j\omega)\left\langle 1,1 \right\rangle\}$,
$g_c = 0$ and $b_c = 0$.

It should be noted that this only refers to the branches interconnecting
several nodes in the system. AC shunt impedances should be treated
separately as a shunt component. DC shunt impedances can only be added
as DC loads.

### Transformer model

Since the [Transformer](@ref) model considered previously cannot be easily represented as the model for an AC branch, Y parameters are extracted from
the ABCD parameters.

In the case of DC branches, since ABCD parameters are each of size
$1 \times 1$ (i.e. scalars), the tap value can be determined as
$\tau = \sqrt{\frac{A}{D}}$, while the series impedance is obtained as
$r = \Re\{\frac{B}{\tau}\}$.

In the case of AC networks and three-phase transformers, using the
assumption of a balanced system, the submatrices
$\textbf{Y}\left\langle 1:3, 1:3 \right\rangle$,
$\textbf{Y}\left\langle 1:3, 4:6 \right\rangle$,
$\textbf{Y}\left\langle 4:6, 1:3 \right\rangle$ and
$\textbf{Y}\left\langle 4:6, 4:6 \right\rangle$ are diagonal. Thus, it
is sufficient to use a single diagonal value from each submatrix. Then,
$\textbf{Y}\left\langle 1,1 \right\rangle = \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right) \, \frac{1}{\tau^2}$,
$\textbf{Y}\left\langle 1,4 \right\rangle = \textbf{Y}\left\langle 4,1 \right\rangle = -\frac{y_s}{\tau \, \exp(-j\theta_{\rm shift})}$
and
$\textbf{Y}\left\langle 4,4 \right\rangle = \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right)$.
The following expressions are derived: $$\begin{aligned}
\nonumber &\tau = \sqrt{\frac{\textbf{Y}\left\langle 4,4 \right\rangle}{\textbf{Y}\left\langle 1,1 \right\rangle}}, \quad
\nonumber \theta_{\rm shift} = 0, \\
\nonumber &y_s = -\textbf{Y}\left\langle 1,4 \right\rangle \, \tau \, \exp(-j\theta_{\rm shift}) , \\
\nonumber &y_c = \textbf{Y}\left\langle 4,4 \right\rangle - y_s, \\
\nonumber &r_s = \Re\left\{\frac{1}{y_s}\right\}, \quad x_s = \Im\left\{\frac{1}{y_s}\right\}, \\
\nonumber &g_c = \Re\{y_c\}, \quad b_c = \Im\{y_c\}.
\end{aligned}$$

### Transmission line model

A transmission line (OHL, cable, cross-bonded cable or mixed OHL-cable)
is represented using its nominal $\pi$-model depicted in the next figure.
$$ 
\nonumber \textbf{Z}(j\omega) = \textbf{Y}_c^{-1} \sinh(\mathbf{\Gamma} l), \\
\textbf{Y}(j\omega) = \textbf{Y}_c \tanh(\mathbf{\Gamma} l).
$$

![Nominal $\pi$-model of the transmission
line.](pictures/transmission_line/nominal-PI.pdf)

For the DC case, the shunt admittance is not considered, while the
branch resistance is equal to $r = \Re\{Z(0)\}$.

For the balanced AC transmission line, the impedance and admittance
matrices are diagonal. It can be chosen as
$Z(j\omega) = \textbf{Z}(j\omega) \left\langle 1,1 \right\rangle$ and
$Y(j\omega) = \textbf{Y}\left\langle 1,1 \right\rangle$. Then, the AC
branch model is given by: $$\begin{aligned}
\nonumber & \tau = 0, \quad \theta_{\rm shift} = 0, \\
\nonumber & r_s = \Re\{Z(j\omega)\}, \quad x_s = \Im\{Z(j\omega)\}, \\
& g_c = \Re\{Y(j\omega)\}, \quad b_c = \Im\{Y(j\omega)\}.
\end{aligned}$$



## Power converter

A power converter is, in accordance with [^2], modeled
together with its phase reactor, filter and transformer. In order to
match the constructed MMC model in [Multi modular converter](@ref) only
the reactor is considered.

Losses of the converter are calculated in the form of
$P_{loss} = a + b I_c + c I_c^2$. Since the switches are modeled as
ideal, the model implementation supports only a constant value
$c = \frac{R_{arm}}{2}$.

![Power flow model of the power
converter.](pictures/power_flow/converter_pf.png)

Depending on the actual realisation of the converter's controls, the
parameters of the converter can be set as a DC voltage controlling or an
active power controlling converter.

## References

[^1]: H. Ergun, J. Dave, D. Van Hertem, and F. Geth, "Optimal power flow for AC/DC grids: Formulation, convex relaxation, linear approximation, and implementation," IEEE Transactions on Power Systems, vol. 34, no. 4, pp. 2980-2990, July 2019, doi: 10.1109/TPWRS.2019.2897835.

[^2]: J. Beerten, "MATACDC 1.0 User's manual," Department Electrical Engineering, University of Leuven, 2012. [Online]. Available: https://www.esat.kuleuven.be/electa/teaching/matacdc/MatACDCManual

[^3]: R. D. Zimmerman and C. E. Murillo-Sánchez, "MATPOWER 6.0 user's manual," PSERC, Tempe, AZ, USA, 2016.


