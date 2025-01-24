The main components of the hybrid AC and DC power system are DC and AC
grid equivalents (represented as sources), impedances, transformers,
transmission lines and cables, breakers, FACTS, shunt components and
converters (MMC converters, two level converters, etc.). In this
chapter, the representation of each of the components as a multiport
network using ABCD parameters is derived and presented, as implemented
in the tool.

## Impedance {#sec:impedance}

An impedance can be defined between $n$ input ports (nodes) and $n$
output ports (nodes). An impedance is represented as an $n \times n$
matrix $\mathbf{Z}$: $$\mathbf{Z} = \begin{bmatrix}
Z_{11} & Z_{12} & \cdots & Z_{1n} \\
\vdots & \vdots & \ddots & \vdots \\
Z_{n1} & Z_{n2} & \cdots & Z_{nn}
\end{bmatrix}$$ An example of an impedance with two input ports and two
output ports is given in Fig.
[1](#fig_impedance_model){reference-type="ref"
reference="fig_impedance_model"}.

![Model of the 2 input ports/2 output ports
impedance.](pictures/impedance/impedance_model.pdf){#fig_impedance_model}

To represent impedances as multiport components with ABCD parameters,
the following equation constructed using the Modified Nodal Analysis
(MNA) approach [@ho1975modified] needs to be solved. $$\small
\underbrace{\begin{bmatrix}
\underset{i}{\operatorname{diag}} \left\{\sum\limits_{j, Z_{ij} \neq 0} \frac{1}{Z_{ij}}\right\}_{n \times n} & \vline & \operatorname{diag}\{-1\}_{n \times n} \\
\hline 
\mathbf{N}_{1, n \times n} & \vline & \mathbf{0}_{n \times n} 
\end{bmatrix}}_{\mathbf{M}_1} \times
\begin{bmatrix}
\mathbf{V}_p \\
\mathbf{I}_p
\end{bmatrix} = 
\underbrace{\begin{bmatrix}
\mathbf{N}_{2, n \times n} & \vline & \mathbf{0}_{n \times n}   \\
\hline 
-\underset{i}{\operatorname{diag}} \left\{\sum\limits_{j, Z_{ji} \neq 0} \frac{1}{Z_{ji}}\right\}_{n \times n} & \vline & \operatorname{diag}\{-1\}_{n \times n}
\end{bmatrix}}_{\mathbf{M}_2} \times
\begin{bmatrix}
\mathbf{V}_s \\
\mathbf{I}_s
\end{bmatrix},$$ where matrices $\mathbf{N}_1$ and $\mathbf{N}_2$
consist of $n$ rows with $n$ columns with entries at the position
$(i,j)$ equal to $-\frac{1}{Z_{ji}}$ and $\frac{1}{Z_{ij}}$, for
$Z_{ij} \neq 0$ and $Z_{ji} \neq 0$ (where $i$ represents row and $j$
column in impedance matrix), respectively. The solution of the previous
system is given as $\mathbf{M}_1^{-1}\mathbf{M}_2$ if $\mathbf{M}_1$ is
invertible matrix, or by determining LU decomposition and reduced row
echelon form if it is not invertible.

For example, for the circuit depicted in Fig.
[1](#fig_impedance_model){reference-type="ref"
reference="fig_impedance_model"}, the equations would be:
$$\begin{bmatrix}
\frac{1}{Z_{11}} + \frac{1}{Z_{12}} & 0 & \vline & -1 & 0 \\
0 & \frac{1}{Z_{21}} + \frac{1}{Z_{22}} & \vline & 0 & -1 \\
\hline 
-\frac{1}{Z_{11}} & -\frac{1}{Z_{21}} & \vline & 0 & 0 \\
-\frac{1}{Z_{12}} & -\frac{1}{Z_{22}} & \vline & 0 & 0
\end{bmatrix} \times
\begin{bmatrix}
V_{11} \\
V_{12} \\
I_{11} \\
I_{12}
\end{bmatrix} = 
\begin{bmatrix}
\frac{1}{Z_{11}} & \frac{1}{Z_{12}} & \vline & 0 & 0 \\
\frac{1}{Z_{21}} & \frac{1}{Z_{22}}  & \vline & 0 & 0 \\
\hline 
-\left(\frac{1}{Z_{11}} + \frac{1}{Z_{21}}\right) & 0 & \vline & -1 & 0 \\
0 & -\left(\frac{1}{Z_{12}} + \frac{1}{Z_{22}}\right) & \vline & 0 & -1
\end{bmatrix} \times
\begin{bmatrix}
V_{21} \\
V_{22} \\
I_{21} \\
I_{22}
\end{bmatrix}.$$

## Transformer {#sec:transformer}

A transformer is modeled as in Fig.
[2](#fig_transformer_model){reference-type="ref"
reference="fig_transformer_model"}. The parameters of the transformer
can be defined explicitly or determined from on-site test data as
described in [@martinez2017power]. On-site test data are typically
presented in the form of open and short-circuit values of the primary
side voltage $V_1$ and current $I_1$ and secondary side voltage $V_2$
and current $I_2$, along with the core power losses $P_{1,core}$ and
winding power losses $P_{1,winding}$. The open and short-circuit test
should be performed on the secondary side.

Then the parameters from Fig.
[2](#fig_transformer_model){reference-type="ref"
reference="fig_transformer_model"} can be estimated as:
$$\begin{aligned}
\nonumber R_{ps} &= \dfrac{P_1^{short}}{(I_1^{short})^2}, \qquad L_{ps} &= \dfrac{Q_1^{short}}{\omega (I_1^{short})^2}, \\
\nonumber R_m &= \dfrac{(V_1^{open})^2}{P_1^{open}}, \qquad L_m &= \dfrac{(V_1^{open})^2}{\omega Q_1^{open}}, \\
\nonumber n &= \dfrac{V_1^{open}}{V_2^{open}}, & \\
\nonumber R_p &= \dfrac{R_{ps}}{2}, \qquad L_p &= \dfrac{L_{ps}}{2}, \\
R_s &= \dfrac{R_{ps}}{2n^2}, \qquad L_s &= \dfrac{L_{ps}}{2n^2},
\end{aligned}$$ knowing that
$Q_1^{o,s} = \sqrt{(V_1^{o,s}I_1^{o,s})^2 - {P_1^{o,s}}^2}$.

ABCD multiport parameters are then estimated as
[@bayo2018control; @wu2014impact]: $$\begin{bmatrix}
A & B \\
C & D
\end{bmatrix} = \mathbf{Y}_{turn} \times \left(\mathbf{Z}^p_{winding} \times \mathbf{Y}_{iron} \times \mathbf{N}_{tr} \times \mathbf{Z}^s_{winding} \parallel \mathbf{Z}_{stray}\right) \times \mathbf{Y}_{turn},
\label{eq_transformer_abcd}$$ where $\mathbf{Y}_{turn} = \begin{bmatrix}
1 & 0 \\
sC_t & 1
\end{bmatrix}$, $\mathbf{Z}^p_{winding} = \begin{bmatrix}
1 & sL_p + R_p \\
0 & 1
\end{bmatrix}$, $\mathbf{Y}_{iron} = \begin{bmatrix}
1 & 0 \\
\frac{1}{sL_m} + \frac{1}{R_m} & 1
\end{bmatrix}$, $\mathbf{Z}^s_{winding} = \begin{bmatrix}
1 & sL_s + R_s \\
0 & 1
\end{bmatrix}$, $Z_{stray} = \begin{bmatrix}
1 & \frac{1}{sC_{stray}} \\
0 & 1
\end{bmatrix}$ and $\mathbf{N}_{tr} = \begin{bmatrix}
n & 0 \\
0 & \frac{1}{n}
\end{bmatrix}$, with $n$ the turn ratio.

![Model of the
transformer.](pictures/transformer/transformer_model.pdf){#fig_transformer_model}

Three-phase transformers can be either in the YY and $\Delta$Y
configuration, where each of the three single-phase transformers is
represented by its equivalent from Fig.
[2](#fig_transformer_model){reference-type="ref"
reference="fig_transformer_model"}.

-   The YY configuration is derived from the equation
    [\[eq_transformer_abcd\]](#eq_transformer_abcd){reference-type="eqref"
    reference="eq_transformer_abcd"}, such that $$\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix} = 
    \begin{bmatrix}
    \operatorname{diag}\{A\}_{3 \times 3} & \operatorname{diag}\{B\}_{3 \times 3}  \\
    \operatorname{diag}\{C\}_{3 \times 3} & \operatorname{diag}\{D\}_{3 \times 3}   
    \end{bmatrix}.$$

-   The $\Delta$Y configuration is more complex and it is modeled using
    following equations. The inner primary and secondary stages of the
    transformer (i.e. all the components except the parasitic
    capacitances and the load impedance) are given by: $$\begin{bmatrix}
    A & B \\
    C & D
    \end{bmatrix}_{inner} = \mathbf{Z}^p_{winding} \times \mathbf{Y}_{iron} \times \mathbf{N}_{tr} =
    \begin{bmatrix}
    n + n \ (sL_p+R_p) \, \left(\frac{1}{sL_m} + \frac{1}{R_m}\right) & \frac{sL_p + R_p}{n} \\
    n \, \left(\frac{1}{L_m} + \frac{1}{R_m}\right) & \frac{1}{n}
    \end{bmatrix}.
    \label{eq_z_inner}$$ The $\Delta$Y configuration transforms voltages
    from the delta side $v^{a,b,c}_p$ to the wye side voltages
    $v^{a,b,c}_s$ as $v^{a,b,c}_p = \sqrt{3} \, v^{a,b,c}_s$, while the
    currents are related as: $$\mathbf{i}_p^{a,b,c} = \begin{bmatrix}
    \frac{1}{\sqrt{3}} & 0 & -\frac{1}{\sqrt{3}} \\
    -\frac{1}{\sqrt{3}}  & \frac{1}{\sqrt{3}} & 0 \\
    0 & -\frac{1}{\sqrt{3}} & \frac{1}{\sqrt{3}} 
    \end{bmatrix} \times \mathbf{i}_s^{a,b,c} .$$

    Using the previous voltage/current relations and ABCD representation
    of the inner transfer function in Eq.
    [\[eq_z_inner\]](#eq_z_inner){reference-type="eqref"
    reference="eq_z_inner"}, the inner impedance can be obtained as:
    $$\mathbf{Z}_{inner} =
    \begin{bmatrix}
    \operatorname{diag}\{A_{inner} \sqrt{3}\}_{3 \times 3} & \vline & \operatorname{diag}\{\frac{B_{inner}}{\sqrt{3}}\}_{3 \times 3} \\
    \hline
    \operatorname{diag}\{C_{inner} \sqrt{3}\}_{3 \times 3} & \vline &
    \begin{array}{ccc}
     \frac{D_{inner}}{\sqrt{3}} & 0 & -\frac{D_{inner}}{\sqrt{3}} \\
    -\frac{D_{inner}}{\sqrt{3}} & \frac{D_{inner}}{\sqrt{3}} & 0 \\
    0 & - \frac{D_{inner}}{\sqrt{3}} & \frac{D_{inner}}{\sqrt{3}} \\ 
    \end{array}
    \end{bmatrix}.$$ The transformer is eventually represented using
    ABCD parameters as: $$\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix} = (\mathbf{Y}_{turn} \times (\mathbf{Z}_{inner} \parallel \mathbf{Z}_{stray}) \times \mathbf{Y}_{turn}).$$

## Autotransformer

This type of model may be expanded into a multi-winding transformer,
e.g., a three-winding transformer. As an example, the positive and
negative, and zero sequence impedance of an autotransformer with YNa0(d)
configuration is shown in
Fig [3](#fig:threewindingtrf){reference-type="ref"
reference="fig:threewindingtrf"}. In Fig.
 [3](#fig:threewindingtrf){reference-type="ref"
reference="fig:threewindingtrf"}, H, X and Y refer to the high-voltage,
low-voltage and tertiary voltage side respectively. The per unit leakage
impedances may be obtained from the per unit leakage impedances
$Z_{\text{HX}}$, $Z_{\text{HY}}$ and $Z_{\text{XY}}$, as obtained using
the short-circuit test, and impedance to ground $Z_g$
as [@anderson1995analysis]: $$\begin{bmatrix}
        Z_{\text{X}}\\
        Z_{\text{Y}}\\
        Z_{\text{Z}}\\
    \end{bmatrix} = 
    \frac{1}{2}\begin{bmatrix}
        1 & 1 & -1\\
        1 & -1 & 1\\
        -1 & 1 & 1\\
    \end{bmatrix}
        \begin{bmatrix}
        Z_{\text{HX}}\\
        Z_{\text{HY}}\\
        Z_{\text{XY}}\\
    \end{bmatrix}\text{, and}$$ $$\begin{bmatrix}
        Z_{\text{X0}}\\
        Z_{\text{Y0}}\\
        Z_{\text{n0}}\\
    \end{bmatrix} = 
    \frac{1}{2}\begin{bmatrix}
        1 & 1 & -1 & \frac{n-1}{n}\\
        1 & -1 & 1 & -\frac{n-1}{n^2}\\
        -1 & 1 & 1 & \frac{1}{n}\\
    \end{bmatrix}
        \begin{bmatrix}
        Z_{\text{HX}}\\
        Z_{\text{HY}}\\
        Z_{\text{XY}}\\
        6Z_{\text{g}}\\
    \end{bmatrix},$$ where $n$ is the winding transformation ratio. The
phase (or physical) domain model may be derived from these sequence
impedances using the Fortescue transform.

<figure id="fig:threewindingtrf">
<p>    </p>
<figcaption>Three-winding transformer model for autotransformer with
YNa0(d) configuration in positive and negative sequence (a) and zero
sequence (b).</figcaption>
</figure>

## Transmission line

The ABCD model parameters of a transmission line can be defined as
[@castellanos1997full; @morched1999universal] $$\begin{bmatrix}
\mathbf{A} & \mathbf{B} \\
\mathbf{C} & \mathbf{D}
\end{bmatrix} =
\begin{bmatrix}
\cosh(\Gamma l) & \mathbf{Y}_c^{-1}\sinh(\Gamma l) \\
\mathbf{Y}_c \sinh(\Gamma l) &  \cosh(\Gamma l)
\end{bmatrix}
\label{eq_tl_abcd}$$ where
$\mathrm{\Gamma} = \sqrt{\mathbf{Z}\mathbf{Y}}$ and
$\mathbf{Y}_c = \mathbf{Z}^{-1} \, \Gamma$, and $l$ standing for the
line or cable length. The used formula is based on the
frequency-dependent phase domain model.

### Overhead line

Based on realizations of the transmission line as defined in PSCAD, see
Fig. [4](#fig_tl_pscad){reference-type="ref" reference="fig_tl_pscad"},
five possible realizations are defined as:

1.  flat (horizontal, which presents flat configuration without ground
    wires),

2.  vertical,

3.  delta (for at least three phases),

4.  offset (for at least three phases),

5.  concentric (for at least three phases).

Besides, the conductor positions can be added manually as absolute
$(x,y)$ positions.

Thus, the simulator enables the creation of overhead lines with the
**conductors** having the properties [@martinez2017power] presented in
the table [1](#tab:ohl:conductors){reference-type="ref"
reference="tab:ohl:conductors"}.

::: {#tab:ohl:conductors}
           Symbol           Meaning
  ------------------------- ---------------------------------------------------------------------
            $n_b$           number of bundles (or a number of phases)
          $n_{sb}$          number of subconductors per bundle
          $y_{bc}$          height of the lowest bundle above ground
       $\Delta y_{bc}$      vertical offset between bundles
       $\Delta x_{bc}$      horizontal offset between the lowest bundles
   $\Delta \Tilde{x}_{bc}$  horizontal offset in the case of concentric and offset organization
          $d_{sb}$          distance between closest subconductors with equidistant
                            concentric organization (symmetric)
          $d_{sag}$         maximal sag offset
            $r_c$           radius of the conductor
          $R_{dc}$          DC resistance of the conductor
            $g_c$           shunt conductance of the conductor
         $\mu_{r,c}$        relative permeability of the conductor
          positions         added manually
        organization        can be flat, vertical, concentric, delta and offset

  : OHL conductor parameters
:::

**Ground wires** have the properties presented in table
[2](#tab:ohl:groundwires){reference-type="ref"
reference="tab:ohl:groundwires"}.

::: {#tab:ohl:groundwires}
      Symbol     Meaning
  -------------- ------------------------------------------------------------------
      $n_g$      number of ground wires
   $\Delta x_g$  relative horizontal distance between ground wires
   $\Delta y_g$  relative vertical between ground wires and the lowest conductors
      $r_g$      radius of the ground wire
   $d_{g,sag}$   ground wire sag
    $R_{g,dc}$   DC resistance of the ground wires
   $\mu_{r,g}$   relative permeability of the ground wire

  : OHL groundwires' properties
:::

![PSCAD overhead line organization: 1) single conductor with groundwire;
2) single conductor; 3) 2 conductors flat; 4) 3 conductors flat; 5) 3
conductors delta; 6) 3 conductors horizontal; 7) 3 conductors vertical;
8) 3 conductors concentric; 9) 3 conductors
offset.](pictures/transmission_line/tl.png){#fig_tl_pscad
width="\\linewidth"}

The transmission line model is constructed using the procedure from
[@martinez2017power; @manitoba2003research]. The overhead transmission
line consists of $n_b$ including sub-conductors, stranding, etc. and
$n_g$ ground wires.

Each line/conductor positioned as $x_c$ relatively starting from the
central tower position and $y_c$ vertically, measured from ground, with
the sag at the midpoint between towers $d_{sag}$, see Fig.
[5](#fig_tl){reference-type="ref" reference="fig_tl"}a. Thus, the
modified vertical position is used in calculations as
$\hat{y}_c = y_c - \frac{2}{3} \, d_{sag}$. Conductor is formed using
$n_{sb}$ sub-conductors grouped in the bundle, where all sub-conductors
are grouped using symmetrical equidistant pattern with the distance
between the two nearest sub-conductors being $d_{bc}$, or a bundle
spacing. Using conductor position, the position of the each
sub-conductor can be estimated. Knowing the angle between two
sub-conductors on the circle and its radius $$\begin{aligned}
\nonumber \varphi = \frac{360^\circ}{n_{sb}}, \\
r = \frac{d_{sb}}{2 \sin(\varphi/2)}, 
\end{aligned}$$ the position can be estimated starting from the angle
$\varphi_s = \frac{\pi}{2}$ if the number of sub-conductors is odd, or
from $\varphi_s = \frac{\pi + \varphi}{2}$ for an even number of
sub-conductors, as follows: $$\begin{aligned}
\nonumber x_{bc} = x_c + r \cos(\varphi_s + k \, \varphi), \\
y_{bc} = y_c + r \sin(\varphi_s + k \, \varphi) - \frac{2}{3} \, d_{sag}, 
\end{aligned}$$ for $k \in \{1, 2, ..., n_{bc}\}$. If the number of
sub-conductors is equal to one, its position is given by
$(x_c, \hat{y}_c)$. Each conductor is characterized with the relative
permeability of the material $\mu_r$, the conductor dc resistance
$R_{dc}$ and the radius $r_i$.

![Overhead line modelling: a) tower and relative conductor positions; b)
sub-conductor
bundle.](pictures/transmission_line/transmission_line.pdf){#fig_tl}

Ground wires are modeled similarly, represented with their relative
position $(x_g,y_g)$, radius $r_g$, dc resistance $R_{gdc}$ and relative
permeability of the material $\mu_r$.

Earth parameters are given with permeability $\mu_e$, permittivity
$\epsilon_e$ and conductivity $\rho_e$.

In order to represent the transmission line using ABCD parameters, it is
necessary to calculate series impedance and shunt admittance matrices
[@martinez2017power]. Both matrices are of the size $n \times n$, where
$n = \sum\limits_{i=1}^{n_c} n^i_{bc} + n_g$. The impedance matrix has
the following form: $$\begin{aligned}
\mathbf{Z} = \operatorname{diag}(Z_i) + 
\begin{bmatrix}
Z_{0,11} & \cdots & Z_{0,1n} \\
\vdots & \ddots & \vdots \\
Z_{0,n1} & \cdots & Z_{0,nn}
\end{bmatrix}
\end{aligned}$$ where
$Z_i = \frac{m\rho_i}{2\pi r_i} \, \coth(0.733mr_i) + \frac{0.3179 \rho_i}{\pi r_i^2}$
for the $i$-th conductor/sub-conductor/ground wire and $r_i$ is its
radius, resistivity $\rho_i = R^i_{dc} \, \pi r_i^2$ and
$m = \sqrt{j\omega \,\frac{\mu_0 \mu_{r,i}}{\rho_i}}$; The components
$Z_{0, ij} = \frac{j\omega \, \mu_0}{2\pi} \, \log\left( \frac{\hat{D}_{ij}}{d_{ij}}\right)$
for $$\begin{aligned}
\nonumber &d_{ij} = \left\{\begin{array}{ll}
\sqrt{(x_i-x_j)^2 + (y_i-y_j)^2}, & \quad i \neq j, \\
r_i, & \quad i = j,
\end{array} \right. \\
\nonumber &D_{ij} = \left\{\begin{array}{ll}
\sqrt{(x_i-x_j)^2 + (y_i+y_j)^2}, & \quad i \neq j, \\
2y_i, & \quad i = j,
\end{array} \right. \\
\nonumber &\hat{D}_{ij} = \sqrt{(y_i + y_j + 2d_e)^2 + (x_i-x_j)^2}, \\
&d_e = \sqrt{\frac{1}{j\omega \, \mu_e (\sigma_e + j\omega \, \epsilon_e)}}.
\end{aligned}$$

The shunt admittance is a matrix formed as
$$\mathbf{Y} = s \, \mathbf{P}^{-1} + \mathbf{G}$$ from matrix
$\mathbf{P}$ with its components
$\mathbf{P}_{ij} = \frac{1}{2\pi\epsilon_0} \, \log\left( \frac{D_{ij}}{d_{ij}}\right)$
and $\mathbf{G} = \operatorname{diag}\{g_c\}$.

### Cable

The cable groups are implemented focusing on the available
configurations of the cables available in PSCAD, where a cable can be
included either placed inside the pipe, so called pipe-type cables, or
placed underground. Cables are usually coaxial with up to 4 layers of
both conductors and insulators.

Cables can be insulated or pipe-type coaxial cables. At the moment, only
a group of coaxial cables is implemented in the package. A cable group
consists of $n$ cables, each one have maximum three conducting layers
and three insulation layers, as can be seen from Fig.
[6](#fig_cable){reference-type="ref" reference="fig_cable"}. The
conducting layers of the cable are denoted as core, sheath and armor.
Between the conducting layers, there are insulators, except for the last
conductor where the insulator is not a strict necessity, but it is
common. For each conductor the following set of parameters is given:
$r^c_i$ and $r^c_o$ as conductor inner and outer radius in meters,
conductor relative permeability $\mu^c_r$ and conductor resistivity
$\rho_c$ (in \[$\Omega$m\]). The insulator is described using the
following parameters: $r^i_i$ and $r^i_o$ are the insulator inner and
outer radius in meters, $\epsilon^i$ is the insulator relative
permittivity and $\mu^i_r$ the insulator relative permeability.

![Coaxial cable.](pictures/transmission_line/cable.pdf){#fig_cable}

Additionally, the configuration parameters can be modified by adding two
semiconducting layers in the insulator 1, and implementing the sheath
consisting of the wire screen and outer sheath layer. In that case, the
procedure described in [@wu2014impact] is applied.

-   Conductor surface impedance

    A hollow conductor surface impedance is given by: $$\begin{aligned}
    \nonumber Z_{aa} = \frac{\rho_c m}{2\pi r^c_i} \, \coth(m(r^c_o - r^c_i)) + \frac{\rho^c}{2\pi r^c_i \, (r^c_i + r^c_o)} \; \left[\frac{\Omega}{\text{m}}\right]  & \mbox{for inner surface}, \\
    \nonumber Z_{bb} = \frac{\rho^c m}{2\pi r^c_o} \, \coth(m(r^c_o - r^c_i)) + \frac{\rho^c}{2\pi r^c_o \, (r^c_i + r^c_o)} \; \left[\frac{\Omega}{\text{m}}\right]  & \mbox{for outer surface}, \\
    Z_{ab} =  \frac{\rho^c m}{\pi (r^c_i + r^c_o)} \, \operatorname{csch}(m(r^c_o - r^c_i)) \; \left[\frac{\Omega}{\text{m}}\right],
    \end{aligned}$$ where $m = \sqrt{j\omega \mu^c_r}$ For a non-hollow
    conductor, the outer surface impedance is $$\begin{aligned}
    Z_{bb} = \frac{\rho^c m}{2\pi r^c_o} \, \coth(0.733 mr^c_o) + \frac{0.3179 \rho^c}{\pi {r^c_o}^2} \; \left[\frac{\Omega}{\text{m}}\right].
    \end{aligned}$$

-   The insulator layer between two conductors has an impedance
    $$Z_i = \frac{j\omega \mu_0 \mu^i_r}{2\pi} \, \log\left(\frac{r^i_o}{r^i_i}\right).$$

-   The earth return impedance of the cable and mutual between cables is
    $$\begin{aligned}
    Z_g = \frac{j\omega \, \mu_g}{2\pi} \, \left(-\log\left(\frac{\gamma m D}{2}\right) + \frac{1}{2} - \frac{2}{3} \, mH\right),
    \end{aligned}$$ for $$\begin{aligned}
    \nonumber &D = \left\{ \begin{array}{ll}
    \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2} & \quad \mbox{for cables } i \neq j, \\
    r_i & \quad \mbox{radius of the cable }i,
    \end{array} \right. \\
    \nonumber &H = \left\{ \begin{array}{ll}
    y_i + y_j & \quad \mbox{for cables } i \neq j, \\
    2y_i & \quad \mbox{for the cable }i,
    \end{array} \right. \\
    \end{aligned}$$ and $\gamma \approx 0.5772156649$ being Euler's
    constant.

According to [@ametani1980general; @de1987computation], one cable is
represented with its series impedance $\mathbf{Z}_{ii}$ matrix. Each
matrix $\mathbf{Z}_{ii}$ has the size $n_c \times n_c$ and its entries
for $j \in \{1, \ldots, n_c-1\}$ are given by $$\begin{aligned}
\nonumber \mathbf{Z}_{ii} \left\langle j,j \right\rangle = Z^j_{bb} + Z^j_i + Z^{j+1}_{aa}, \\ 
\nonumber \mathbf{Z}_{ii} \left\langle j,j+1 \right\rangle = Z_{ii}\left\langle j+1,j \right\rangle = - Z^{j+1}_{ab}, \\ 
\mathbf{Z}_{ii} \left\langle n_c, n_c \right\rangle = Z^{n_c}_{bb} + Z^{n_c}_i + Z^{ii}_g, 
\end{aligned}$$ and otherwise the matrix entries are 0.

Mutual surface impedances between the cables are given by matrix
$\mathbf{Z}_{ij}$ having all components equal to $Z^{ij}_g$.

The shunt admittance matrix can be estimated as
$\mathbf{Y} = s\mathbf{P}^{-1}$ and matrix $P$, which has the form
$$\mathbf{P} = \begin{bmatrix}
\mathbf{P}_{11} & \mathbf{P}_{12} & \cdots & \mathbf{P}_{1n} \\
\vdots & \ddots & & \vdots \\    
\mathbf{P}_{n1} & \mathbf{P}_{n2} & \cdots & \mathbf{P}_{nn}
\end{bmatrix}.$$ Matrices $\mathbf{P}_{ii}$ have components
$$\mathbf{P}_{ii} = \begin{bmatrix}
P_c + P_s + P_a & P_s + P_a & P_a \\
P_s + P_a & P_s + P_a & P_a \\
P_a & P_a & P_a
\end{bmatrix} + 
\begin{bmatrix}
P_{ii} & P_{ii} & P_{ii} \\
P_{ii} & P_{ii} & P_{ii} \\
P_{ii} & P_{ii} & P_{ii} 
\end{bmatrix},$$ where $P_{c,s,a}$ belong respectively to core, shield
and armor insulators and have the following values:
$P = \frac{\log(r_o/r_i)}{2\pi\epsilon}$ and
$P_{ii} = \frac{\log(2h_i/r)}{2\pi\epsilon_0}$ is a earth return.
Matrices $\mathbf{P}_{ij}$, for $i \neq j$,have all components equal to
$P_{ij} = \frac{\log(D_2/D_1)}{2\pi\epsilon_0}$, where
$D_1 = \sqrt{(x_i-x_j)^2 + (y_i-y_j)^2}$ and
$D_2 = \sqrt{(x_i-x_j)^2 + (y_i+y_j)^2}$ [@ametani1980general].

As it is valid to assume that sheath and armor are grounded, it is
allowed to use Kron reduction. Using Kron reduction, as proposed in
[@de2011effects; @rivas2002calculation], when applied to the matrices
$Y$ and $Z$ a compact shunt admittance and series impedance model is
obtained. For determining the corresponding ABCD parameters, the same
procedure is used as for the transmission line from equation
[\[eq_tl_abcd\]](#eq_tl_abcd){reference-type="eqref"
reference="eq_tl_abcd"}.

### Cross-bonded cables

Cables are cross-bonded in order to reduce sheath circulating currents.
The cross-bonding is made by transposing sheaths of the cable sections.
As in [@wu2014impact], this transposition can be made in ABCD domain.

Cross-bonded cables consists of minor sections as in Fig.
[7](#fig:cable:crossbonding){reference-type="ref"
reference="fig:cable:crossbonding"}, where the smaller sections
(referred to as minor sections) are marked with J, K and L. Minor
sections are then grouped into bigger (major) sections, for which all
the cable layers except the core are short connected to ground. Thus,
the ABCD parameters of the major section can be estimated using Kron
elimination.

![Cross-bonded
cable.](pictures/transmission_line/crossbonding.png){#fig:cable:crossbonding
width="\\linewidth"}

The procedure for determining the ABCD parameters of the whole
cross-bonded cable is as follows:

Let us assume that the ABCD parameters of each major section are marked
as $ABCD^r_\eta$. Thus, the equivalent cable ABCD parameters are given
by $ABCD = \prod\limits_{\eta = 1}^{n} ABCD^r_\eta$, where $n$ is the
number of major sections.

The equivalent ABCD parameters of one major section can be estimated as
follows:

-   Determine the ABCD parameters of the each minor section inside the
    major: $ABCD_{\eta, i}$, for $i \in \{1,m\}$ and $m$ is the number
    of the minor sections inside $\eta$ major section.

-   Reorganize the ABCD matrix for each minor section as:
    $M_{\eta, i} = \mathbf{R} \, ABCD_{\eta,i} \, \mathbf{R}^{-1}$. The
    matrix $\mathbf{R}$ is a transposition matrix that sorts voltages
    and currents from the form:
    $[V_{1,c}, \ V_{1,s}, \ V_{2,c}, \ \ldots, \ I_{1,c}, \ I_{1,s}, \ I_{2,c}, \ \ldots]^T$
    into
    $[V_{1,c}, \ V_{2,c}, \ V_{3,c}, \ V_{1,s} \ldots, \ I_{1,c}, \ I_{2,c}, \ I_{3,c}, \ \ldots]^T$.
    Basically, it groups first all core cable voltages, then all sheath
    voltages, \...

-   Apply transposition from A-B-C to C-A-B for all minor sections
    except for the first:
    $\mathbf{M}_{CB} \, \mathbf{T} \, \mathbf{M}_{\eta,i} \, \mathbf{T}^{-1}$,
    where $\mathbf{M}_{CB}$ introduces sheath cross-bonding losses.
    Matrix $\mathbf{M}_{CB}$ is the identity matrix except for the
    indices that belong to interconnections of sheath voltages and
    currents. For example, assuming that $n_c = 3$ (this is the number
    of the cables), the sheath is the second layer of the total $n_l$
    layers and thus\
    $\mathbf{M}_{CB} \left\langle n_c+1:2n_c, n_c*n_l + n_c + 1:n_c*n_l + 2n_c \right\rangle = \operatorname{diag}\{2Z_{CB}\}_{n_c \times n_c}$.
    The impedance $Z_{CB}$ presents the impedance from nonideal bonding.

-   Apply the ABCD reduction introduced in [@wu2014impact] and described
    in section
    [\[sec:abcd:reduction\]](#sec:abcd:reduction){reference-type="ref"
    reference="sec:abcd:reduction"}.

According to the previous description, the ABCD parameters of the major
section are given by:
$$ABCD_{\eta} = M_{\eta,1} \times \prod\limits_{i=2}^m \mathbf{M}_{CB} \, (\mathbf{T} \, \mathbf{M}_{\eta,i} \, \mathbf{T}^{-1})$$

### Mixed OHL-cables

Mixed OHL-cable components contain OHLs and cable sections.. Each OHL
and cable section is characterized individually and a complete 'mixed
OHL-cable' component is presented with the equivalent ABCD
representation.

This ABCD representation has the form:
$$ABCD = \prod\limits_{\eta=1}^n ABCD_\eta,$$ where $ABCD_\eta$ are the
ABCD parameters of an OHL or cable section, while $n$ is the total
number of sections.

## AC and DC grid equivalents

AC and DC grid equivalents can be modelled as either ideal AC and DC
sources respectively, or including an equivalent impedance (e.g.
short-circuit impedance to model the system strength) . These sources
are described using the following relations: $$\begin{aligned}
\nonumber \mathbf{V}_p &=& \mathbf{V}_s + \mathbf{Z} \, \mathbf{I}_s +  \mathbf{V}, \\
\mathbf{I}_p &=& \mathbf{I}_s,
\end{aligned}$$ for $\mathbf{V}$ being the vector of the voltage
source's voltages, $\mathbf{Z}$ the series equivalent impedance as
diagonal matrix with the values:
$\mathbf{Z} = \operatorname{diag}\{ Z_s\}$. As explained in
[\[sec_multiport\]](#sec_multiport){reference-type="ref"
reference="sec_multiport"} with $\mathbf{I}_p$ and $\mathbf{V}_p$ are
presented input voltage source currents and voltages, while with
$\mathbf{I}_s$ and $\mathbf{V}_s$ the output currents and voltages.

For the estimation of the equivalent impedance of the network,
independent voltage sources are short-circuited, which means that in
this case, the grid ABCD parameters can be represented as an identity
matrix. Additionally, the internal grid impedance can be added as an
serial connection of the impedance and voltage source. The ABCD
parameters of the equivalent network are now given by: $$\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix} = 
    \begin{bmatrix}
    \mathbf{I} & \mathbf{Z} \\
    \mathbf{0} & \mathbf{I}
    \end{bmatrix}.$$

## MMC {#sec:mmc}

Voltage source converters (VSC) are often implemented as modular
multiterminal converters (MMC). For example for HVDC applications, the
MMC has become the de-facto industry standard nowadays. They are used
for a fast and efficient conversion of energy. In the code, they are
represented as inverters or rectifiers, which do not only consider the
impedance as seen from either AC or DC side, but also consider the
coupling between the two side. Hence, they have two DC side pins and
three AC side pins. The converter is represented with its admittance
matrix, which is incorporated with the ABCD system of equations.

![MMC.](pictures/mmc/mmc.pdf){#fig_mmc}

### MMC model

An MMC is depicted in Fig. [8](#fig_mmc){reference-type="ref"
reference="fig_mmc"}. The variables from Fig.
[8](#fig_mmc){reference-type="ref" reference="fig_mmc"} are defined for
all three phases, $j \in \{ a,b,c \}$. The sets of submodules are
represented by their averaged equivalent, and thus, the following
equations for voltages and currents can be written: $$\begin{aligned}
\nonumber v_{Mj}^{U,L} = m^{U,L}_j v_{Cj}^{U,L}, \\
i_{Mj}^{U,L} = m^{U,L}_j i^{U,L}_j,
\end{aligned}$$ where $m^{U,L}_j$ are the corresponding insertion
indices.

Using $\Sigma-\Delta$ nomenclature, the variables can be represented as:
$$\begin{aligned}
\nonumber i^\Delta_j &=& i^U_j - i^L_j, \quad i^\Sigma_j = \dfrac{i^U_j + i^L_j}{2}, \\
\nonumber v_{Cj}^\Delta &=& \dfrac{v_{Cj}^U - v_{Cj}^L}{2}, \quad v_{Cj}^\Sigma = \dfrac{v_{Cj}^U + v_{Cj}^L}{2}, \\
\nonumber m^\Delta_j &=& m^U_j - m^L_j, \quad m^\Sigma_j = m^U_j + m^L_j, \\
\nonumber v_{Mj}^\Delta &=& \dfrac{-v_{Mj}^U + v_{Mj}^L}{2} = -\dfrac{m^\Delta_j v_{Cj}^\Sigma + m^\Sigma_j v_{Cj}^\Delta}{2}, \\
\nonumber v_{Mj}^\Sigma &=& \dfrac{v_{Mj}^U + v_{Mj}^L}{2} = \dfrac{m^\Sigma_j v_{Cj}^\Sigma + m^\Delta_j v_{Cj}^\Delta}{2},
\end{aligned}$$ In order to obtain the differential equations in the dqz
frame, Park's transformation is used to represent the system in a set of
several rotating frames, each of which is related to a different angular
frequency. The Park's transformation is as defined in
[\[eq:park\]](#eq:park){reference-type="eqref" reference="eq:park"}. It
should be noted that zero sequence of the $\Sigma$ components
corresponds to the DC currents and voltages. On the other hand Z (being
$Zd$ and $Zq$) incorporates third harmonic into $\Delta$ currents and
voltages modeling.

For the purpose of the modeling, the MMC converter is represented using
12 differential equations for the state variables
[@bergna2018generalized; @bergna2018pi]: $$\begin{aligned}
\nonumber \dfrac{d i^\Delta_d}{dt} &=& -\dfrac{v^G_d - v^\Delta_{Md} + R^{ac}_{eq} i^\Delta_d + \omega L^{ac}_{eq} i^\Delta_q}{L^{ac}_{eq}}, \\
\nonumber \dfrac{d i^\Delta_q}{dt} &=& -\dfrac{v^G_q - v^\Delta_{Mq} + R^{ac}_{eq} i^\Delta_q - \omega L^{ac}_{eq} i^\Delta_d}{L^{ac}_{eq}}, \\
 \nonumber \dfrac{d i^\Sigma_d}{dt} &=& -\dfrac{v^\Sigma_{Md} + R_{arm} i^\Sigma_d - 2 \omega L_{arm}i^\Sigma_q}{L_{arm}}, \\
 \nonumber \dfrac{d i^\Sigma_q}{dt} &=& -\dfrac{v^\Sigma_{Mq} + R_{arm} i^\Sigma_q + 2 \omega L_{arm}i^\Sigma_d}{L_{arm}}, \\
\nonumber \dfrac{d i^\Sigma_z}{dt} &=& -\dfrac{v^\Sigma_{Mz} - \frac{v_{dc}}{2} + R_{arm} i^\Sigma_z}{L_{arm}}, \\
\nonumber \dfrac{d v^\Delta_{Cd}}{dt} &=& \dfrac{N}{2C_{arm}} \, \left( i^\Sigma_z m^\Delta_d - \frac{i^\Delta_q m^\Sigma_q}{4} + i^\Sigma_d \ \left(\frac{m^\Delta_d}{2} + \frac{m^\Delta_{Zd}}{2}\right) - i^\Sigma_q \ \left(\frac{m^\Delta_q}{2} + \frac{m^\Delta_{Zq}}{2} \right) \right. \\
\nonumber &&\left. + i^\Delta_d \ \left( \frac{m^\Sigma_d}{4} + \frac{m^\Sigma_z}{2} \right) - 2\omega C_{arm} \ v_{Cq}^\Delta \right), \\
\nonumber \dfrac{d v^\Delta_{Cq}}{dt} &=& -\dfrac{N}{2C_{arm}} \, \left( \frac{i^\Delta_d m^\Sigma_q}{4} - i^\Sigma_z m^\Delta_q + i^\Sigma_q \ \left(\frac{m^\Delta_d}{2} - \frac{m^\Delta_{Zd}}{2}\right) + i^\Sigma_d \ \left(\frac{m^\Delta_q}{2} - \frac{m^\Delta_{Zq}}{2} \right) \right. \\
\nonumber &&\left. + i^\Delta_q \ \left( \frac{m^\Sigma_d}{4} - \frac{m^\Sigma_z}{2} \right) - 2\omega C_{arm} \ v_{Cd}^\Delta \right), \\
\nonumber \dfrac{d v^\Delta_{CZd}}{dt} &=& -\dfrac{N}{8C_{arm}} \, \left( i^\Delta_d m^\Sigma_d + 2i^\Sigma_d m^\Delta_d + i^\Delta_q m^\Sigma_q + 2 i^\Sigma_q m^\Delta_q  + 4i^\Sigma_z m^\Delta_{Zd} \right) - 3\omega v^\Delta_{CZq}, \\
\nonumber \dfrac{d v^\Delta_{CZq}}{dt} &=& -\dfrac{N}{8C_{arm}} \, \left( i^\Delta_q m^\Sigma_d + 2 i^\Sigma_d m^\Delta_q - i^\Delta_d m^\Sigma_q - 2 i^\Sigma_q m^\Delta_d  + 4i^\Sigma_z m^\Delta_{Zq} \right) + 3\omega v^\Delta_{CZd}, \\
\end{aligned}$$ $$\begin{aligned}
\nonumber \dfrac{d v_{Cd}^\Sigma}{dt} &=& \dfrac{N}{2C_{arm}} \, \left(i^\Sigma_d m^\Sigma_z + i^\Sigma_z m^\Sigma_d + i^\Delta_d \ \left(\frac{m^\Delta_d}{4} + \frac{m^\Delta_{Zd}}{4}\right) - i^\Delta_q \ \left(\frac{m^\Delta_q}{4} - \frac{m^\Delta_{Zq}}{4}\right) \right) + 2\omega C_{arm} v_{Cq}^\Sigma, \\
\nonumber \dfrac{d v_{Cq}^\Sigma}{dt} &=& -\dfrac{N}{2C_{arm}} \, \left(i^\Delta_q \ \left(\frac{m^\Delta_d}{4} - \frac{m^\Delta_{Zd}}{4}\right) - i^\Sigma_z m^\Sigma_q + i^\Delta_d \ \left(\frac{m^\Delta_q}{4} + \frac{m^\Delta_{Zq}}{4}\right) - i^\Sigma_q m^\Sigma_z \right) + 2\omega C_{arm} v_{Cd}^\Sigma, \\
\dfrac{d v_{Cz}^\Sigma}{dt} &=& -\dfrac{N}{8C_{arm}} \, \left(i^\Delta_d m^\Delta_d + i^\Delta_q m^\Delta_q +
2i^\Sigma_d m^\Sigma_d + 2i^\Sigma_q m^\Sigma_q + 4i^\Sigma_z m^\Sigma_z \right),
\end{aligned}$$ where $L^{ac}_{eq} = L_f + \frac{L_{arm}}{2}$ and
$R^{ac}_{eq} = R_f + \frac{R_{arm}}{2}$. The state variables are
$\mathbf{x} = [\mathbf{i}^\Delta_{dq}, \ \mathbf{i}^\Sigma_{dqz}, \ \mathbf{v}^\Delta_{CdqZ},\ \mathbf{v}_{Cdqz}^\Sigma]^T$.
The 12 algebraic relations used for determining 7 voltages
$[v^\Delta_{Md}, v^\Delta_{Mq}, v^\Delta_{MZd}, v^\Delta_{MZq}, v^\Sigma_{Md}, v^\Sigma_{Mq}, v^\Sigma_{Mz}]$
and insertion indeices
$\begin{bmatrix} m^\Delta_d, & m^\Delta_q, & m^\Delta_{Zd}, & m^\Delta_{Zq}, &
m^\Sigma_d, & m^\Sigma_q, & m^\Sigma_z \end{bmatrix}^T$ are given as:
$$\begin{aligned}
\nonumber v_{Md}^\Delta &=& \dfrac{m^\Delta_q v^\Sigma_{Cq}}{4} - \dfrac{m^\Delta_d v^\Sigma_{Cz}}{2} - \dfrac{m^\Delta_d v^\Sigma_{Cd}}{4} - \dfrac{m^\Delta_{Zd} v^\Sigma_{Cd}}{4} + \dfrac{m^\Delta_{Zq} v^\Sigma_{Cq}}{4} - \dfrac{m^\Sigma_d v^\Delta_{Cd}}{4} - \dfrac{m^\Sigma_z v^\Delta_{Cd}}{2} + \dfrac{m^\Sigma_q v^\Delta_{Cq}}{4} \\
\nonumber && -\dfrac{m^\Sigma_d v^\Delta_{CZd}}{4} + \dfrac{m^\Sigma_q v^\Delta_{CZq}}{4}, \\
\nonumber v^\Delta_{Mq} &=& \dfrac{m^\Delta_d v^\Sigma_{Cq}}{4} + \dfrac{m^\Delta_q v^\Sigma_{Cd}}{4} - \dfrac{m^\Delta_q v^\Sigma_{Cz}}{2} - \dfrac{m^\Delta_{Zd} v^\Sigma_{Cq}}{4} - \dfrac{m^\Delta_{Zq} v^\Sigma_{Cd}}{4} + \dfrac{m^\Sigma_d v^\Delta_{Cq}}{4} + \dfrac{m^\Sigma_q v^\Delta_{Cd}}{4} - \dfrac{m^\Sigma_z v^\Delta_{Cq}}{2} \\
\nonumber && -\dfrac{m^\Sigma_d v^\Delta_{CZq}}{4} - \dfrac{m^\Sigma_q v^\Delta_{CZd}}{4}, \\
\nonumber v^\Delta_{MZd} &=& -\dfrac{m^\Delta_d v^\Sigma_{Cd}}{4} - \dfrac{m^\Delta_q v^\Sigma_{Cq}}{4} - \dfrac{m^\Delta_{Zd} v^\Sigma_{Cz}}{2} - \dfrac{m^\Sigma_d v^\Delta_{Cd}}{4} - \dfrac{m^\Sigma_q v^\Delta_{Cq}}{4} - \dfrac{m^\Sigma_z v^\Delta_{Zd}}{2}, \\
\nonumber v^\Delta_{MZq} &=& -\dfrac{m^\Delta_d v^\Sigma_{Cq}}{4} - \dfrac{m^\Delta_q v^\Sigma_{Cd}}{4} - \dfrac{m^\Delta_{Zq} v^\Sigma_{Cz}}{2} - \dfrac{m^\Sigma_d v^\Delta_{Cq}}{4}+ \dfrac{m^\Sigma_q v^\Delta_{Cd}}{4} - \dfrac{m^\Sigma_z v^\Delta_{Zq}}{2}, \\
\nonumber v^\Sigma_{Md} &=& \dfrac{m^\Delta_d v^\Delta_{Cd}}{4} - \dfrac{m^\Delta_q v^\Delta_{Cq}}{4} + \dfrac{m^\Delta_d v^\Delta_{CZd}}{4} + \dfrac{m^\Delta_{Zd} v^\Delta_{Cd}}{4} + \dfrac{m^\Delta_q v^\Delta_{Zq}}{4} +
\dfrac{m^\Delta_{Zq} v^\Delta_{Cq}}{4} + \dfrac{m^\Sigma_d v^\Sigma_{Cz}}{2} + \dfrac{m^\Sigma_z v^\Sigma_{Cd}}{2}, \\
\nonumber v^\Sigma_{Mq} &=& \dfrac{m^\Delta_q v^\Delta_{Zd}}{4} - \dfrac{m^\Delta_q v^\Delta_{Cd}}{4} - \dfrac{m^\Delta_d v^\Delta_{Zq}}{4} - \dfrac{m^\Delta_d v^\Delta_{Cq}}{4} + \dfrac{m^\Delta_{CZd} v^\Delta_{Cq}}{4} -
\dfrac{m^\Delta_{Zq} v^\Delta_{Cd}}{4} + \dfrac{m^\Sigma_q v^\Sigma_{Cz}}{2} + \dfrac{m^\Sigma_z v^\Sigma_{Cq}}{2}, \\
v^\Sigma_{Mz} &=& \dfrac{m^\Delta_d v^\Delta_{Cd}}{4} + \dfrac{m^\Delta_q v^\Delta_{Cq}}{4} + \dfrac{m^\Delta_{Zd} v^\Delta_{CZd}}{4} + \dfrac{m^\Delta_{Zq} v^\Delta_{CZq}}{4} + \dfrac{m^\Sigma_d v^\Sigma_{Cd}}{4} + \dfrac{m^\Sigma_q v^\Sigma_{Cq}}{4} + \dfrac{m^\Sigma_z v^\Sigma_{Cz}}{2},
\end{aligned}$$ $$\begin{bmatrix}
m^\Delta_d \\ 
m^\Delta_q \\ 
m^\Delta_{Zd} \\ 
m^\Delta_{Zq} \\
m^\Sigma_d \\
m^\Sigma_q \\
m^\Sigma_z
\end{bmatrix} =
\dfrac{2}{v_{dc}} \ 
\begin{bmatrix}
-v^{\Delta}_{Md, ref} \\
-v^{\Delta}_{Mq, ref} \\
-v^{\Delta}_{MZd, ref} \\
-v^{\Delta}_{MZq, ref} \\
v^{\Sigma}_{Md, ref} \\
v^{\Sigma}_{Mq, ref} \\
v^{\Sigma}_{Mz, ref}
\end{bmatrix}.$$

The set of the previous 12 differential equations and the set of
algebraic equations are accompanied with the 7 equations for the
reference values of the voltages
$[\mathbf{v}^{\Delta }_{MdqZ, ref}, \mathbf{v}^{\Sigma}_{Mdqz, ref}]$.
The reference voltages are given as zero by default, except for the
value of $v^{\Sigma}_{Cz, ref} = \frac{v_{dc}}{2}$.

### Operating point

The converter's operating point can be defined manually or derived as a
result of solving the power flow equations of the interconnected power
system (including the converter's steady-state characteristics). In both
situations, the following fields should be present:

  ---------------------- ------------------------------------------------------
   $P_{min}$, $P_{max}$   minimum and maximum active AC power of the converter
           $P$             power flow estimated or predefined active AC power
   $Q_{min}$, $Q_{max}$            minimum and maximum reactive power
           $Q$             power flow estimated or predefined reactive power
         $P_{dc}$             power flow estimated or predefined DC power
         $V_{DC}$                              DC voltage
     $V_m$, $\theta$             amplitude and phase of the AC voltage
  ---------------------- ------------------------------------------------------

Using previous fields, the converter's operating point is estimated by
solving a set of linear differential equations to obtain converter's
steady-state. As a reference for the MMC's initial operating point using
values obtained using power flow, we define: $$\begin{aligned}
\nonumber i^{\Delta C}_{d, ref} = \frac{2}{3} \, \frac{(v^{GC}_d P + v^{GC}_q Q)}{{v^{GC}_d}^2 + {v^{GC}_q}^2}, \\
\nonumber i^{\Delta C}_{q, ref} = \frac{2}{3} \, \frac{(v^{GC}_q P - v^{GC}_d Q)}{{v^{GC}_d}^2 + {v^{GC}_q}^2}, \\
\nonumber i^\Sigma_{z, ref} = \frac{P_{dc}}{3 V_{DC}}, \\
\nonumber P_{ac, ref} = P, \\
\nonumber Q_{ar, ref} = Q, \\
\nonumber v_{dc, ref} = V_{DC}, \\
W^\Sigma_{z, ref} = \frac{3C_{arm} V_{DC}^2}{N}.
\end{aligned}$$

### Control implementations

For the PI controls in the dqz frame additional equations have been
developed
[@bergna2018generalized; @bergna2018pi; @sakinci2019generalized]. The
different controllers are considered to be tuned using a pole placement
method. Also the PLL is implemented using a PI controller structure as
in [@freytes2017analyse].

#### Phase locked loop (PLL)

The PLL is used to synchronize the converter's internal controller
frequency, used to control the currents in a rotating frame, to the grid
frequency. All converter variables are mapped to the dqz frame using the
same Park's transformation without a phase shift.

According to Fig. [9](#fig_mmc_pll){reference-type="ref"
reference="fig_mmc_pll"} the following equations for PLL can be written.
$$\begin{aligned}
    \nonumber \dfrac{d \xi_{pll}}{dt} &=& -v^{G,C}_q, \\
    \nonumber \dfrac{d \theta}{dt} &=& \Delta \omega, \\
    \nonumber \Delta \omega &=& -K_{p,pll} \ v^{G,C}_q + K_{i,pll} \ \xi_{pll}, \\
    \omega_C &=& \Delta \omega + \omega_0.
    
\end{aligned}$$

It should be noted that the output current control and circulating
current control are implemented in the converter's reference frame.
Thus, the rotation should be applied to the corresponding variables
before the control law is applied. However, since the converter's
internal dynamics is analyzed in the grid's reference frame, the output
of the mentioned controls should be restored to the grid's reference
frame using the inverse of the rotation matrix, given by:
$$T(\theta) = \begin{bmatrix}
        \cos(\theta) & -\sin(\theta) \\
        \sin(\theta) & \cos(\theta)
        \end{bmatrix},
        \label{eq_rot}$$ while its inverse is:
$$T^{-1}(\theta) = \begin{bmatrix}
        \cos(\theta) & \sin(\theta) \\
        -\sin(\theta) & \cos(\theta)
        \end{bmatrix}.
        \label{eq_irot}$$

![PLL implementation.](pictures/mmc/pll.pdf){#fig_mmc_pll}

The mapping of $i^\Delta_{dq, ref}$ from the grid's to the converter's
reference frame is done by: $$\begin{bmatrix}
        i^{\Delta C}_{d,ref} \\
        i^{\Delta C}_{q,ref}
        \end{bmatrix} = T(\theta) \ 
        \begin{bmatrix}
        i^{\Delta}_{d,ref} \\ 
        i^{\Delta}_{q,ref}
        \end{bmatrix},$$ Also currents $i^\Delta_{d}$ and $i^\Delta_{q}$
are mapped to: $$\begin{bmatrix}
        i^{\Delta C}_{d} \\ 
        i^{\Delta C}_{q} 
        \end{bmatrix}= T(\theta) \ 
        \begin{bmatrix}
        i^{\Delta}_{d} \\ 
        i^{\Delta}_{q}
        \end{bmatrix}.$$

Similarly: $$\begin{bmatrix}
        i^{\Sigma C}_{d,ref} \\
        i^{\Sigma C}_{q,ref}
        \end{bmatrix} =
        T(-2\theta) \ 
        \begin{bmatrix}
        i^\Sigma_{d,ref} \\
        i^\Sigma_{q,ref}
        \end{bmatrix}, \qquad 
        \begin{bmatrix}
        i^{\Sigma C}_{d} \\
        i^{\Sigma C}_{q}
        \end{bmatrix} =
        T(-2\theta) \ 
        \begin{bmatrix}
        i^\Sigma_{d} \\
        i^\Sigma_{q}
        \end{bmatrix}.$$

#### DC voltage control

DC voltage control (DCC) provides the reference value for
$i^{\Delta C}_{d,ref}$, depending of the variation of $v_{dc}$. The
control law provides the following equations $$\begin{aligned}
    \frac{d v_{dc}}{dt} &=& \frac{N}{6 C_{arm}} \left(i_{dc} - 3i^\Sigma_z \right) , \\
    \frac{d\xi_{v_{dc}}}{dt} &=& v_{dc, ref} - v_{dc}, \\
    \nonumber i^{\Delta C}_{d, ref} &=& - K_{p,dc} \, (v_{dc, ref} - v_{dc}) - K_{i,dc}  \, \xi_{v_{dc}}.
    
\end{aligned}$$

![DC voltage control.](pictures/mmc/dc_control.pdf){#fig_mmc_dc}

#### Output current control (OCC)

defines the reference values for the output currents $i^\Delta_{d,ref}$
and $i^\Delta_{q,ref}$ given in the grid reference frame. This control
method adds several equations: $$\begin{aligned}
    \nonumber \dfrac{d \xi^\Delta_{d}}{dt} &=& i^{\Delta C}_{d,ref} - i^{\Delta C}_d, \\
    \nonumber \dfrac{d \xi^\Delta_{q}}{dt} &=& i^{\Delta C}_{q,ref} - i^{\Delta C}_q, \\
    \nonumber v^{\Delta C}_{Md, ref} &=& K_{i,occ} \xi^\Delta_d + K_{p,occ} \ (i^{\Delta C}_{d, ref} - i^{\Delta C}_d) + \omega_C L^{ac}_{eq} i^\Delta_q + v^{G,C}_d, \\
     v^{\Delta C}_{Mq, ref} &=& K_{i,occ} \xi^\Delta_q + K_{p,occ} \ (i^{\Delta C}_{q, ref} - i^{\Delta C}_q) - \omega_C L^{ac}_{eq} i^{\Delta C}_d + v^{G,C}_q.
    
\end{aligned}$$ Voltages $v^{\Delta C}_{Md, ref}$ and
$v^{\Delta C}_{Mq, ref}$ are used in grid's reference frame for further
calculations: $$\begin{bmatrix}
        v^{\Delta}_{Md, ref} \\
        v^{\Delta}_{Mq, ref}
        \end{bmatrix} = 
        T^{-1}(\theta) \ 
        \begin{bmatrix}
        v^{\Delta C}_{Md, ref} \\
        v^{\Delta C}_{Mq, ref}
        \end{bmatrix}.$$

If the controller is defined only using bandwidth $\omega_n$ and $\zeta$
(instead of $K_p$ and $K_i$), the proportional and integral gains are
tuned as: $$\begin{aligned}
    K_{i,occ} &=& L^{ac}_{eq} \ \omega_n^2, \\
    K_{p,occ} &=& 2 \zeta \omega_n L^{ac}_{eq} - R^{ac}_{eq}.
    
\end{aligned}$$

![OCC implementation.](pictures/mmc/occ.pdf){#fig_mmc_occ}

#### Circulating current control (CCC)

The CCC is constructed to set the circulating current to it reference,
which is considered to be $i^\Sigma_{d,ref} = 0$,
$i^\Sigma_{q,ref} = 0$.

The equations added by the CCC are: $$\begin{aligned}
    \nonumber \dfrac{d \xi^\Sigma_d}{dt} &=& i^{\Sigma C}_{d,ref} - i^{\Sigma C}_{d}, \\
    \nonumber \dfrac{d i^\Sigma_{q}}{dt} &=& i^{\Sigma C}_{q,ref} - i^{\Sigma C}_q, \\
    \nonumber v_{Md, ref}^{\Sigma C} &=& -K_{i,ccc} \ \xi^\Sigma_d - K_{p,ccc} \ (i^{\Sigma C}_{d,ref} - i^{\Sigma C}_d) + 2\omega_C L_{arm} i^{\Sigma C}_q, \\
    v_{Mq, ref}^{\Sigma C} &=& -K_{i,ccc} \xi^\Sigma_q - K_{p,ccc} \ (i^{\Sigma C}_{q,ref} - i^{\Sigma C}_q) - 2\omega_C L_{arm} i^{\Sigma C}_d.
    
\end{aligned}$$ To return to the grid's reference frame, the following
transformation is applied: $$\begin{bmatrix}
        v^{\Sigma}_{Md, ref} \\
        v^{\Sigma}_{Mq, ref}
        \end{bmatrix} = 
        T^{-1}(-2\theta) \ 
        \begin{bmatrix}
        v^{\Sigma C}_{Md, ref} \\
        v^{\Sigma C}_{Mq, ref}
        \end{bmatrix}.$$

The proportional and integral gains are tuned as: $$\begin{aligned}
    K_i &=& L_{arm} \ \omega_n^2, \\
    K_p &=& 2 \zeta \omega_n L_{arm} - R_{arm}.
    
\end{aligned}$$

![CCC implementation.](pictures/mmc/ccc.pdf){#fig_mmc_ccc}

#### Energy control and zero current control

The energy control is built around the "zero" energy and as a result, it
provides a reference value for the 'zero' current $i^\Sigma_{z,ref}$.
The energy controller involves the following equations, as visible from
Fig. [13](#fig_mmc_ec_zcc){reference-type="ref"
reference="fig_mmc_ec_zcc"}a. $$\begin{aligned}
    \nonumber W^\Sigma_z &=& \dfrac{3C_{arm}}{2N} \ ({v_{Cd}^\Delta}^2 + {v_{Cq}^\Delta}^2 + {v_{CZd}^\Delta}^2 + {v_{CZq}^\Delta}^2 + {v_{Cd}^\Sigma}^2 + {v_{Cq}^\Sigma}^2 + 2{v_{Cz}^\Sigma}^2), \\
    \nonumber \dfrac{d \xi_{W^\Sigma_z}}{dt} &=& W^\Sigma_{z,ref} - W^\Sigma_z, \\
    \nonumber P_{ac} &=& \frac{3}{2} \ (v^{G,C}_d i^{\Delta C}_d + v^{G,C}_q i^{\Delta C}_q), \\
    i^\Sigma_{z,ref} &=& \frac{K_{p,ec} \ (W^\Sigma_{z,ref} - W^\Sigma_z) + K_{i,ec} \xi_{W^\Sigma_z} + P_{ac}}{3v_{dc}}.
    
\end{aligned}$$

Additionally, the zero current control (ZCC) sets the zero current to
the desired value. The implementation of this control is depicted in
Fig. [13](#fig_mmc_ec_zcc){reference-type="ref"
reference="fig_mmc_ec_zcc"}b. It can work without the energy controller.
$$\begin{aligned}
    \nonumber \dfrac{d \xi^\Sigma_z}{dt} &=& i^\Sigma_{z,ref} - i^\Sigma_z, \\
    v_{Mz, ref}^\Sigma &=& \frac{v_{dc}}{2} - K_{p,zcc} \ (i^\Sigma_{z,ref} - i^\Sigma_z) - K_{i,zcc} \ \xi^\Sigma_z.
    
\end{aligned}$$ The tuning of the ZCC employs the same principles as for
CCC.

![Energy control and ZCC
implementation.](pictures/mmc/energy_zcc.pdf){#fig_mmc_ec_zcc
width="\\linewidth"}

#### Active and reactive power control

An outer control loop for the control of the active and reactive power
can be added, see Fig. [14](#fig_mmc_power){reference-type="ref"
reference="fig_mmc_power"}. These control loops are used to successfully
estimate the AC currents $i^\Delta_{dq,ref}$. The control loops operate
according to the following equations: $$\begin{aligned}
    \nonumber P_{ac} &=& \frac{3}{2} \ (v^{G,C}_d i^{\Delta C}_d + v^{G,C}_q i^{\Delta C}_q), \\
    \nonumber Q_{ac} &=& \frac{3}{2} \ (-v^{G,C}_d i^{\Delta C}_q + v^{G,C}_q i^{\Delta C}_d), \\
    \nonumber \dfrac{d \xi_{P_{ac}}}{dt} &=& P_{ac,ref} - P_{ac}, \\
    \nonumber \dfrac{d \xi_{Q_{ac}}}{dt} &=& Q_{ac,ref} - Q_{ac}, \\
    \nonumber i^{\Delta C}_{d,ref} &=& K_p^{P_{ac}} \ (P_{ac,ref} - P_{ac}) + K_i^{P_{ac}} \xi_{P_{ac}}, \\
    i^{\Delta C}_{q,ref} &=& -K_p^{Q_{ac}} \ (Q_{ac,ref} - Q_{ac}) - K_i^{Q_{ac}} \xi_{Q_{ac}}.
    
\end{aligned}$$

![Active and reactive power control
implementation.](pictures/mmc/power.pdf){#fig_mmc_power}

### Steady-state solution and admittance model

The previous system of differential and algebraic equations is solved
for the equilibrium using Julia package NLsolve. After determining the
equilibrium, the system is represented as a multi-input multi-output
system (MIMO), where the variables\
${\mathbf{x} = \begin{bmatrix}
\mathbf{i}^\Delta_{dq} & \mathbf{i}^\Sigma_{dqz} & \mathbf{v}^\Delta_{CdqZ} & \mathbf{v}^\Sigma_{Cdqz} \end{bmatrix}}$
represent the state-variables, whereas the input vector is given as
$\mathbf{u} = \begin{bmatrix}
v_{dc} & v^G_d & v^G_q
\end{bmatrix}$. In order to obtain transfer functions from the input to
output, which is defined as $\mathbf{y} = \begin{bmatrix}
3i^\Sigma_z & i^\Delta_d & i^\Delta_q \end{bmatrix}$, the previous
equations are rewritten to satisfy the following form: $$\begin{aligned}
\nonumber \dot{\mathbf{x}}(t) &=& \mathbf{A}_{MIMO} \mathbf{x}(t) + \mathbf{B}_{MIMO} \mathbf{u}(t), \\
\mathbf{y}(t) &=& \mathbf{C}_{MIMO} \mathbf{x}(t) + \mathbf{D}_{MIMO} \mathbf{u}(t).
\label{eq_mimo_sys}
\end{aligned}$$

The corresponding matrices $\mathbf{A}_{MIMO}$, $\mathbf{B}_{MIMO}$,
$\mathbf{C}_{MIMO}$ and $\mathbf{D}_{MIMO}$ are determined as Jacobians
around the equilibrium for the state variables and inputs. The Jacobian
is determined using the Julia package ForwardDiff
[@RevelsLubinPapamarkou2016]. Applying the Laplace transform, the
previous system of equations
[\[eq_mimo_sys\]](#eq_mimo_sys){reference-type="eqref"
reference="eq_mimo_sys"} transforms to: $$\begin{aligned}
\nonumber s\mathbf{X}(s) &=& \mathbf{A}_{MIMO} \mathbf{X}(s) + \mathbf{B}_{MIMO} \mathbf{U}(s), \\
\mathbf{Y}(s) &=& \mathbf{C}_{MIMO} \mathbf{X}(s) + \mathbf{D}_{MIMO} \mathbf{U}(s).
\label{eq_mimo_sys_laplace}
\end{aligned}$$ The MIMO transfer function is thus given by:
$$\mathbf{Y}_{MMC}(s) = \mathbf{Y}(s)\mathbf{U}(s)^{-1} = \mathbf{C}_{MIMO} \ (s\mathbf{I}-\mathbf{A}_{MIMO})^{-1} \ \mathbf{B}_{MIMO} + \mathbf{D}_{MIMO}.
    \label{eq:mmc:mimo}$$

The thus obtained matrix transfer function form the following matrix of
admittances $$\mathbf{Y}_{MMC}(s) = 
    \begin{bmatrix}
    Y_{zz} & Y_{zd} & Y_{zq} \\
    Y_{dz} & Y_{dd} & Y_{dq} \\
    Y_{qz} & Y_{qd} & Y_{qq} \\
    \end{bmatrix}$$ that connect vector of currents
$[i_{dc}(s), \ i^\Delta_d(s), \ i^\Delta_q(s) ]^T$ with the voltages
vector $[v_{dc}(s), \ v^G_d(s), \ v^G_q(s)]^T$.

If we model the MMC as in Fig. [15](#fig_mmc_model){reference-type="ref"
reference="fig_mmc_model"}, the converter can be represented as a one
input, two output component.

![MMC block model.](pictures/mmc/mmc_simplified.pdf){#fig_mmc_model}

As ABCD parameters can be defined properly only for the multiport
networks with the same number of the input and the output pins, the
equations which are solved for the MMC are given with the matrix
$\mathbf{Y}_{MMC}$.

For the case when the converter controls the DC voltage, the vector of
inputs and outputs are: $\mathbf{u} = \begin{bmatrix}
i_{dc} & v^G_d & v^G_q
\end{bmatrix}$ and $\mathbf{y} = \begin{bmatrix}
v_{dc} & i^\Delta_d & i^\Delta_q 
\end{bmatrix}$. Then, the system is also represented as MIMO
[\[eq:mmc:mimo\]](#eq:mmc:mimo){reference-type="eqref"
reference="eq:mmc:mimo"}, but in order to determine
$\mathbf{Y}_{MMC}(s)$, a transformation is applied to change the
positions of $v_{dc}(s)$ and $i_{dc}(s)$.

## Shunt reactor {#sec:shunt_reactor}

The simulator allows for the implementation of single-phase and
three-phase shunt reactors. It is possible to account for the layering
of the component, and the winding of each phase is constructed by
connecting all layers in series, as represented in
Figure [16](#fig_shunt_reactor_model){reference-type="ref"
reference="fig_shunt_reactor_model"}(a). The equivalent circuit of the
single-phase model is given in
Figure [16](#fig_shunt_reactor_model){reference-type="ref"
reference="fig_shunt_reactor_model"}(b), where inductances, resistances,
and parasitic capacitances are modeled as lumped components.

![(a) Layers configuration and (b) equivalent circuit diagram for the
shunt reactor
[@wu2014impact]](pictures/shunt_reactor/shunt_reactor.PNG){#fig_shunt_reactor_model}

A shunt reactor is characterised by the following set of parameters
[@wu2014impact Section 2.5]. For a three-phase shunt reactor, it is
assumed that each phase is characterised by the same parameters. The
three-phase shunt reactor can be connected in Wye or in Delta.

  -------------- ------------------------------------------------------------------------------
      $pins$                                  the number of phases
       $N$                                    the number of layers
      $L_k$                  the series inductance of layer $k$, with $k=1, ..., N$
      $R_k$                  the series resistance of layer $k$, with $k=1, ..., N$
      $C_k$                the cross-over resistance of layer $k$, with $k=1, ..., N$
   $C_{k\_k-1}$   the inter-layer capacitance between layers $k$ and $k-1$, with $k=2, ..., N$
    $C_{1\_E}$               the capacitance between layer 1 and the earthed screen
  -------------- ------------------------------------------------------------------------------

If inductances, resistances, cross-over capacitances and inter-layer
capacitances are known for each layer, they can be provided as vectors
of values. Otherwise, total or average values can be specified, in which
case the layer values are obtained $\forall k$ as: $$\begin{aligned}
    L_k &= \frac{L_\tx{tot}}{N} \\
    R_k &= \frac{R_\tx{tot}}{N} \\
    C_k &= C_\tx{CO, avg} \\
    C_{k\_k-1} &= C_\tx{IL, avg}
\end{aligned}$$ where CO stands for *cross-over* and IL stands for
*inter-layer*.

### Calculation of the ABCD matrix

The procedure for determining the ABCD matrix is as follows. First, the
matrix is obtained for an N-port component, assuming that the layers are
disconnected (N pins on the $p$ side and N pins on the $q$ side). The
connection of the layers in series is done at a later stage by means of
boundary conditions. The corresponding $2N$-by-$2N$ ABCD matrix is
obtained as:
$$\mathbf{ABCD} = \mathbf{K}_{CIL,q-side} \cdot \mathbf{K}_{RLC} \cdot \mathbf{K}_{CIL,p-side}$$
where the matrices are defined as (see
Fig. [17](#fig_shunt_reactor_calc){reference-type="ref"
reference="fig_shunt_reactor_calc"}):

-   for the overall N-port component with all layers disconnected;

-   for the N-port component comprising all p-side inter-layer
    capacitors;

-   the N-port component comprising the RL series components in parallel
    with the cross-over capacitors;

-   for the N-port component comprising all q-side inter-layer
    capacitors.

![Calculation example for
$N=3$](pictures/shunt_reactor/shunt_reactor_ABCD_calculations.jpg){#fig_shunt_reactor_calc
width="\\linewidth"}

#### Matrix $\mathbf{K}_{RLC}$ 

Matrix $\mathbf{K}_{RLC}$ can be obtained by writing the following set
of equations: $$\begin{cases}
                I_{b1} = I_{a1} \\
                I_{b2} = I_{a2} \\
                \vdots \\
                I_{bN} = I_{aN} \\
                U_{b1} = U_{a1} - ((sL_1+R_1)^{-1} + sC_1)^{-1}I_{a1} \\
                U_{b2} = U_{a2} - ((sL_2+R_2)^{-1} + sC_2)^{-1}I_{a2} \\
                \vdots \\
                U_{bN} = U_{aN} - ((sL_N+R_N)^{-1} + sC_N)^{-1}I_{aN}
            \end{cases}$$ which results in an ABCD matrix of the form:
$$\mathbf{K}_{RLC} = 
            \begin{bmatrix}
            \mathbf{I} & \mathbf{Z} \\
            \mathbf{0} & \mathbf{I}
            \end{bmatrix},$$ with in particular the block element
$\mathbf{Z}$: $$\mathbf{Z} = 
            \begin{bmatrix}
            -((sL_1+R_1)^{-1} + sC_1)^{-1} & 0 & \cdots & 0 \\
            0 & -((sL_2+R_2)^{-1} + sC_2)^{-1} & \cdots & 0 \\
            \vdots & \vdots & \ddots & \vdots \\
            0 & 0 & \cdots & -((sL_N+R_N)^{-1} + sC_N)^{-1}
            \end{bmatrix}.$$

#### Matrices $\mathbf{C}_{IL-LE,p}$ and $\mathbf{C}_{IL-LE,q}$ 

Matrices $\mathbf{C}_{IL-LE,p}$ and $\mathbf{C}_{IL-LE,q}$ can be
obtained by writing the following set of equations, for example for
matrix $\mathbf{C}_{IL-LE,p}$: $$\begin{cases}
                U_{a1} = U_{p1} \\
                U_{a2} = U_{p2} \\
                \vdots \\
                U_{aN} = U_{pN} \\
                I_{a1} = I_{p1} - s\frac{1}{2}C_{1\_E}(U_{p1}-0     ) + s\frac{1}{2}C_{2\_1}(U_{p2}-U_{p1}) \\
                I_{a2} = I_{p2} - s\frac{1}{2}C_{2\_1}(U_{p2}-U_{p1}) + s\frac{1}{2}C_{3\_2}(U_{p3}-U_{p2}) \\
                \vdots \\
                I_{aN} = I_{pN} - s\frac{1}{2}C_{N\_N-1}(U_{pN}-U_{p(N-1)})
            \end{cases}$$ which results in matrices of the form:
$$\mathbf{K}_{CIL,q-side} = 
                \begin{bmatrix}
                \mathbf{I} & \mathbf{0} \\
                \mathbf{Y}_q & \mathbf{I}
                \end{bmatrix},\qquad 
                \mathbf{K}_{CIL,p-side} = 
                \begin{bmatrix}
                \mathbf{I} & \mathbf{0} \\
                \mathbf{Y}_p & \mathbf{I}
                \end{bmatrix}.$$ with in particular the block elements
$\mathbf{Y}_p$ and $\mathbf{Y}_q$:
$$\scalemath{0.8}{\mathbf{Y}_p = \mathbf{Y}_q = 
            s\frac{1}{2}
            \begin{bmatrix}
            -C_{1\_E} - C_{2\_1} & C_{2\_1} & 0 & 0 & \cdots & 0 & 0 & 0  \\
            C_{2\_1} & -C_{2\_1} - C_{3\_2} & C_{3\_2} & 0 & \cdots & 0 & 0 & 0  \\
            0 & C_{3\_2} & -C_{3\_2} - C_{4\_3} & C_{4\_3} & \cdots & 0 & 0 & 0  \\
            \vdots & & & & \ddots & & & \vdots \\
            0 & 0 & 0 & 0 & \cdots & C_{N-1\_N-2} & -C_{N-1\_N-2}-C_{N\_N-1} & C_{N\_N-1} \\
            0 & 0 & 0 & 0 & \cdots & 0 & C_{N\_N-1} & -C_{N\_N-1}
            \end{bmatrix}.}$$

### Application of the boundary conditions

The application of the boundary conditions allows to define the
connections of the layers in series, such as: $$\begin{gathered}
        p_1 \longleftrightarrow p_2 \\
        q_2 \longleftrightarrow q_3 \\
        p_3 \longleftrightarrow p_4 \\
        q_4 \longleftrightarrow q_5 \\
        \vdots
    
\end{gathered}$$ which results in the following set of conditions:
$$\begin{cases}
            U_{p1} = U_{p2} \\
            U_{q2} = U_{q3} \\
            U_{p3} = U_{p4} \\
            U_{q4} = U_{q5} \\
            \vdots \\
            I_{p1} + I_{p2} = 0 \\
            I_{q2} + I_{q3} = 0 \\
            I_{p3} + I_{p4} = 0 \\
            I_{q4} + I_{q5} = 0 \\
            \vdots
        \end{cases}$$ These conditions can be imposed to the system by a
series of rows and columns operations. A transformation from **ABCD** to
**Y** makes the application of the boundary conditions easier, as it
gathers all voltages in the input vector and all currents in the output
vector. $$\begin{bmatrix}
            \mathbf{U}_q \\
            \mathbf{I}_q
        \end{bmatrix} = 
        \begin{bmatrix}
        \mathbf{A} & \mathbf{B} \\
        \mathbf{C} & \mathbf{D}
        \end{bmatrix}
        \begin{bmatrix}
            \mathbf{U}_p \\
            \mathbf{I}_p
        \end{bmatrix}
        \longrightarrow
        \begin{bmatrix}
            \mathbf{I}_p \\
            \mathbf{I}_q
        \end{bmatrix} = 
        \begin{bmatrix}
        -\mathbf{B}^{-1}\mathbf{A} & \mathbf{B}^{-1} \\
        \mathbf{C} -\mathbf{D}\mathbf{B}^{-1}\mathbf{A} & \mathbf{D}\mathbf{B}^{-1}
        \end{bmatrix}
        \begin{bmatrix}
            \mathbf{U}_p \\
            \mathbf{U}_q
        \end{bmatrix}$$ The equality of voltages $U_x$ and $U_y$ is
expressed by adding column $y$ to column $x$, removing column $y$ and
replacing $U_x$ by $U_z$. This voltage is an intermediate voltage of the
series connection which is not relevant and will be eliminated later.\
The fact that the sum of two currents $I_x$ and $I_y$ is equal to zero
is expressed by adding row $y$ to row $x$, removing row $y$ and
replacing $I_x$ by 0.

### Reduction of the system

The system is reduced by eliminating all intermediate voltages from the
voltage vector and all zero entries from the current vector. This is
done by first exchanging the intermediate voltages with the zero
entries. The procedure to do so is explained in [@wu2014impact]. After
the exchange, the columns corresponding to zero currents and the rows
corresponding to intermediate voltages can simple be removed.

Eventually, we obtain an admittance description of the form:
$$\begin{bmatrix}
            \mathbf{I}_i \\
            \mathbf{I}_o
        \end{bmatrix} = 
        \begin{bmatrix}
        Y_{11} & Y_{12} \\
        Y_{21} & Y_{22}
        \end{bmatrix}
        \begin{bmatrix}
            \mathbf{U}_i \\
            \mathbf{U}_o
        \end{bmatrix} \label{eqn:shunt_reactor_admittance_representation}$$

### Single-phase and three-phase connections

![Single and three-phase
connections](pictures/shunt_reactor/three_phase_shunt_reactor.jpg){#fig:three_phase_shunt_reactor}

For all connections presented in
Fig. [18](#fig:three_phase_shunt_reactor){reference-type="ref"
reference="fig:three_phase_shunt_reactor"}, we obtain an ABCD
representation in the form of: $$\begin{bmatrix}
            \mathbf{U}_o \\
            \mathbf{I}_o
        \end{bmatrix} = 
        \begin{bmatrix}
        \mathbf{I} & \mathbf{0} \\
        \mathbf{Y} & \mathbf{I}
        \end{bmatrix}
        \begin{bmatrix}
            \mathbf{U}_i \\
            \mathbf{I}_i
        \end{bmatrix}$$ where $\mathbf{I}$ and $\mathbf{0}$ are the
identity and zero matrices of adequate size; $\textbf{U}$ and
$\textbf{I}$ are scalar for the single phase connection and vectors for
the three-phase configurations, e.g.: $$\mathbf{U}_o = 
        \begin{bmatrix}
        U_{Ao} \\
        U_{Bo} \\
        U_{Co}
        \end{bmatrix}$$ Finally, we have the following, where $Y_{11}$
and $Y_{12}$ are as defined in
Eq. [\[eqn:shunt_reactor_admittance_representation\]](#eqn:shunt_reactor_admittance_representation){reference-type="ref"
reference="eqn:shunt_reactor_admittance_representation"}:

-   Single phase: $$\mathbf{Y} = -Y_{11}$$

-   Three-phase delta: $$\mathbf{Y} =
                    \begin{bmatrix}
                        Y_{12}-Y_{11} & -Y_{12} & Y_{11} \\
                        Y_{11} & Y_{12}-Y_{11} & -Y_{12} \\
                        -Y_{12} & Y_{11} & Y_{12}-Y_{11}
                    \end{bmatrix}$$

-   Three-phase wye grounded: $$\mathbf{Y} =
                    \begin{bmatrix}
                        -Y_{11} & 0 & 0 \\
                        0 & -Y_{11} & 0 \\
                        0 & 0 & -Y_{11}
                    \end{bmatrix}$$

-   Three-phase wye ungrounded: $$\mathbf{Y} = \frac{1}{3}Y_{11}
                    \begin{bmatrix}
                        -2 & 1 & 1 \\
                        1 & -2 & 1 \\
                        1 & 1 & -2
                    \end{bmatrix}$$

\@inproceedingsreveyrand2018multiport, title=Multiport conversions
between S, Z, Y, h, ABCD, and T parameters, author=Reveyrand, Tibault,
booktitle=2018 International Workshop on Integrated Nonlinear Microwave
and Millimetre-wave Circuits (INMMIC), pages=1--3, year=2018,
organization=IEEE

\@articlefrickey1994conversions, title=Conversions between S, Z, Y, H,
ABCD, and T parameters which are valid for complex source and load
impedances, author=Frickey, Dean A, journal=IEEE Transactions on
microwave theory and techniques, volume=42, number=2, pages=205--211,
year=1994, publisher=IEEE

\@articledorfler2012kron, title=Kron reduction of graphs with
applications to electrical networks, author=Dorfler, Florian and Bullo,
Francesco, journal=IEEE Transactions on Circuits and Systems I: Regular
Papers, volume=60, number=1, pages=150--163, year=2012, publisher=IEEE

\@articlexu2005harmonic, title=Harmonic resonance mode analysis,
author=Xu, Wilsun and Huang, Zhenyu and Cui, Yu and Wang, Haizhen,
journal=IEEE Transactions on Power Delivery, volume=20, number=2,
pages=1182--1190, year=2005, publisher=IEEE

\@articleho1975modified, title=The modified nodal approach to network
analysis, author=Ho, Chung-Wen and Ruehli, Albert and Brennan, Pierce,
journal=IEEE Transactions on circuits and systems, volume=22, number=6,
pages=504--509, year=1975, publisher=IEEE

\@articlewang2014modeling, title=Modeling and analysis of harmonic
stability in an AC power-electronics-based power system, author=Wang,
Xiongfei and Blaabjerg, Frede and Wu, Weimin, journal=IEEE Transactions
on Power Electronics, volume=29, number=12, pages=6421--6432, year=2014,
publisher=IEEE

\@booksood2006hvdc, title=HVDC and FACTS controllers: applications of
static converters in power systems, author=Sood, Vijay K, year=2006,
publisher=Springer Science & Business Media

\@bookpadiyar1990hvdc, title=HVDC power transmission systems: technology
and system interactions, author=Padiyar, KR, year=1990, publisher=New
Age International

\@inproceedingsmiddlebrook1978design, title=Design techniques for
preventing input-filter oscillations in switched-mode regulators,
author=Middlebrook, RD and others, booktitle=Proc. Powercon, volume=5,
pages=A3--1, year=1978

\@articlewang2018harmonic, title=Harmonic Stability in Power
Electronic-Based Power Systems: Concept, Modeling, and Analysis,
author=Wang, Xiongfei and Blaabjerg, Frede, journal=IEEE Transactions on
Smart Grid, volume=10, number=3, pages=2858--2870, year=2018,
publisher=IEEE

\@bookmilano2010power, title=Power system modelling and scripting,
author=Milano, Federico, year=2010, publisher=Springer Science &
Business Media

\@articlestamatiou2017analytical, title=Analytical derivation of the
DC-side input admittance of the direct-voltage controlled modular
multilevel converter, author=Stamatiou, Georgios and Beza, Mebtu and
Bongiorno, Massimo and Harnefors, Lennart, journal=IET Generation,
Transmission & Distribution, volume=11, number=16, pages=4018--4030,
year=2017, publisher=IET

\@ARTICLEharnefors2007input, author=L. Harnefors and M. Bongiorno and S.
Lundberg, journal=IEEE Transactions on Industrial Electronics,
title=Input-Admittance Calculation and Shaping for Controlled
Voltage-Source Converters, year=2007, volume=54, number=6,
pages=3323-3334, keywords=electric admittance;power
convertors;input-admittance calculation;controlled voltage-source
converters;admittance shaping;Shape control;Voltage
control;Admittance;Matrix converters;Power system dynamics;Control
systems;Power system stability;Pulse width modulation converters;Power
conversion;Frequency conversion;Converter control;dissipativeness;grid
interaction;oscillations;passivity;stability,
doi=10.1109/TIE.2007.904022, ISSN=1557-9948, month=Dec,

\@bookkundur1994power, title=Power system stability and control,
author=Kundur, Prabha and Balu, Neal J and Lauby, Mark G, volume=7,
year=1994, publisher=McGraw-hill New York

\@phdthesisbayo2018control, title=Control Interactions in Power Systems
with Multiple VSC HVDC Converters, school=KU Leuven, author=Bayo Salas,
A, year=2018

\@articleji2019harmonic, title=Harmonic Stability Analysis of MMC-Based
DC System using DC Impedance Model, author=Ji, Ke and Tang, Guangfu and
Yang, Jie and Li, Yunfeng and Liu, Dong, journal=IEEE Journal of
Emerging and Selected Topics in Power Electronics, year=2019,
publisher=IEEE

\@articleagbemuko2019dynamic, title=Dynamic modelling and interaction
analysis of multi-terminal VSC-HVDC grids through an impedance-based
approach, author=Agbemuko, Adedotun J and Domı́nguez-Garcı́a, José Luis
and Prieto-Araujo, Eduardo and Gomis-Bellmunt, Oriol,
journal=International Journal of Electrical Power & Energy Systems,
volume=113, pages=874--887, year=2019, publisher=Elsevier

\@articlelyu2018harmonic, title=Harmonic State-Space Based Small-Signal
Impedance Modeling of a Modular Multilevel Converter With Consideration
of Internal Harmonic Dynamics, author=Lyu, Jing and Zhang, Xin and Cai,
Xu and Molinas, Marta, journal=IEEE Transactions on Power Electronics,
volume=34, number=3, pages=2134--2148, year=2018, publisher=IEEE

\@inproceedingsroose2018impedance, title=Impedance-Based DC Side
Stability Assessment of VSC-HVDC Systems with Control Time Delay,
author=Roose, T and Bayo-Salas, A and Beerten, J, booktitle=2018 20th
European Conference on Power Electronics and Applications (EPE'18 ECCE
Europe), pages=1--10, year=2018

\@articlebessegato2018method, title=A method for the calculation of the
ac-side admittance of a modular multilevel converter, author=Bessegato,
Luca and Harnefors, Lennart and Ilves, Kalle and Norrga, Staffan,
journal=IEEE transactions on power electronics, volume=34, number=5,
pages=4161--4172, year=2018, publisher=IEEE

\@inproceedingssun2018adequacy, title=Adequacy Analysis of Overhead Line
Model for Harmonic Stability Analysis of Grid-Connected Voltage-Source
Converters, author=Sun, Yin and Wu, L and Wang, XF and de Jong, Erik CW
and Blåbjerg, Frede and Cuk, Vladimir and Cobben, JFG, booktitle=2018
IEEE 19th Workshop on Control and Modeling for Power Electronics
(COMPEL), pages=1--8, year=2018, organization=IEEE

\@articleliu2018oscillatory, title=An Oscillatory Stability Criterion
Based on the Unified $dq$-Frame Impedance Network Model for Power
Systems With High-Penetration Renewables, author=Liu, Huakun and Xie,
Xiaorong and Liu, Wei, journal=IEEE Transactions on Power Systems,
volume=33, number=3, pages=3472--3485, year=2018, publisher=IEEE

\@bookerickson2001fundamentals, title=Fundamentals of power electronics,
author=Erickson, Robert W and Maksimović, Dragan, year=2001,
publisher=Springer Science & Business Media

\@bookmartinez2017power, title=Power system transients: parameter
determination, author=Martinez-Velasco, Juan A, year=2017, publisher=CRC
press

\@articlemanitoba2003research, title=Research Centre Inc,
author=Manitoba, HVDC, journal=PSCAD/EMTDC user's guide, year=2003

\@articlecastellanos1997full, title=Full frequency-dependent
phase-domain transmission line model, author=Castellanos, F and Marti,
JR, journal=IEEE transactions on Power Systems, volume=12, number=3,
pages=1331--1339, year=1997, publisher=IEEE

\@articlemorched1999universal, title=A universal model for accurate
calculation of electromagnetic transients on overhead lines and
underground cables, author=Morched, Atef and Gustavsen, Bjorn and
Tartibi, Manoocher, journal=IEEE Transactions on Power Delivery,
volume=14, number=3, pages=1032--1038, year=1999, publisher=IEEE

\@articlede2011effects, title=Effects of conductor counter-transposition
on the positive-sequence impedance and losses of cross-bonded cables,
author=De Léon, Francisco and Márquez-Asensio, Manuel L and
Álvarez-Cordero, Gabriel, journal=IEEE Transactions on Power Delivery,
volume=26, number=3, pages=2060--2063, year=2011, publisher=IEEE

\@articlerivas2002calculation, title=Calculation of frequency-dependent
parameters of power cables: Matrix partitioning techniques,
author=Rivas, Richard A and Marti, Jose R, journal=IEEE transactions on
power delivery, volume=17, number=4, pages=1085--1092, year=2002,
publisher=IEEE

\@articleametani1980general, title=A general formulation of impedance
and admittance of cables, author=Ametani, A, journal=IEEE transactions
on power apparatus and systems, number=3, pages=902--910, year=1980,
publisher=IEEE

\@articlewu2014impact, title=Impact of EHV/HV underground power cables
on resonant grid behavior, author=Wu, Lei and others, journal=Eindhoven:
Eindhoven University of Technology, year=2014, publisher=Citeseer

\@articlevan1986efficient, title=An efficient algorithm for cable
theory, applied to blowfly photoreceptor cells and LMC's, author=Van
Hateren, JH, journal=Biological cybernetics, volume=54, number=4-5,
pages=301--311, year=1986, publisher=Springer

\@articlede1987computation, title=Computation of cable impedances based
on subdivision of conductors, author=de Arizon, Paloma and Dommel,
Hermann W, journal=IEEE Transactions on Power Delivery, volume=2,
number=1, pages=21--27, year=1987, publisher=IEEE

\@articlebergna2018generalized, title=Generalized voltage-based
state-space modeling of modular multilevel converters with constant
equilibrium in steady state, author=Bergna-Diaz, Gilbert and Freytes,
Julian and Guillaud, Xavier and D'Arco, Salvatore and Suul, Jon Are,
journal=IEEE Journal of Emerging and Selected Topics in Power
Electronics, volume=6, number=2, pages=707--725, year=2018,
publisher=IEEE

\@phdthesisfreytes2017analyse, title=Analyse de stabilité en petit
signaux des Convertisseurs Modulaires Multiniveaux et application à
l'étude d'interopérabilité des MMC dans les Réseaux HVDC,
author=Freytes, Julian, year=2017, school=Ecole centrale de Lille

\@articlebergna2018pi, title=PI Passivity-based Control and Performance
Analysis of MMC Multi-Terminal HVDC Systems, author=Bergna-Diaz, Gilbert
and Zonetti, Daniele and Sanchez, Santiago and Ortega, Romeo and
Tedeschi, Elisabetta, journal=IEEE Journal of Emerging and Selected
Topics in Power Electronics, year=2018, publisher=IEEE

\@articlesakinci2019generalized, title=Generalized Dynamic Phasor
Modeling of the MMC for Small-Signal Stability Analysis, author=Sakinci,
Özgür Can and Beerten, Jef, journal=IEEE Transactions on Power Delivery,
volume=34, number=3, pages=991--1000, year=2019, publisher=IEEE

\@articleRevelsLubinPapamarkou2016, title = Forward-Mode Automatic
Differentiation in Julia, author = Revels, J. and Lubin, M. and
Papamarkou, T., journal = arXiv:1607.07892 \[cs.MS\], year = 2016, url =
https://arxiv.org/abs/1607.07892

\@ARTICLEergun2019optimal, author=H. Ergun and J. Dave and D. Van Hertem
and F. Geth, journal=IEEE Transactions on Power Systems, title=Optimal
Power Flow for AC/DC Grids: Formulation, Convex Relaxation, Linear
Approximation, and Implementation, year=2019, volume=34, number=4,
pages=2980-2990, keywords=AC-DC power convertors;approximation
theory;HVDC power convertors;HVDC power transmission;power grids;power
transmission control;reactive power control;AC-DC grids;linear
approximation;active power control capabilities;reactive power control
capabilities;HVDC converter stations;power systems;ancillary
services;optimal power flow model;convex relaxation
formulation;parameterized ac-dc converter model;common ac optimal power
flow formulations;dc nodes;converter station technologies;ac
nodes;ancillary security;open-source tool;Mathematical model;HVDC
transmission;AC-DC power converters;Numerical
models;Inductors;Impedance;Linear approximation;HVDC
transmission;flexible ac transmission systems;power system analysis
computing, doi=10.1109/TPWRS.2019.2897835, ISSN=0885-8950, month=July,

\@articlebeerten2012matacdc, title=MATACDC 1.0 User's manual,
author=Beerten, Jef, journal=Department Electrical Engineering,
University of Leuven, url =
https://www.esat.kuleuven.be/electa/teaching/matacdc/MatACDCManual,
year=2012

\@articleergun2018powermodelsacdc, title=PowerModelsACDC repository,
author=Ergun, H. and Geth, F. and Van Hertem, D., journal=Department
Electrical Engineering, University of Leuven, url =
https://github.com/hakanergun/PowerModelsACDC.jl, year=2018

\@articlezimmerman2016matpower, title=MATPOWER 6.0 user's manual,
author=Zimmerman, Ray D and Murillo-Sánchez, Carlos E, journal=PSERC:
Tempe, AZ, USA, year=2016

\@articlelevron2017modeling, title=Modeling power networks using dynamic
phasors in the dq0 reference frame, author=Levron, Yoash and Belikov,
Juri, journal=Electric Power Systems Research, volume=144,
pages=233--242, year=2017, publisher=Elsevier

\@articlebelikov2018integration, title=Integration of long transmission
lines in large-scale dq0 dynamic models, author=Belikov, Juri and
Levron, Yoash, journal=Electrical Engineering, volume=100, number=2,
pages=1219--1228, year=2018, publisher=Springer

\@articledarco2019time, title=Time-Invariant State-Space Model of an AC
Cable by DQ-Representation of Frequency-Dependent PI-Sections,
author=D'Arco, Salvatore and Suul, Jon Are and Beerten, Jef,
journal=Proceedings of the 2019 IEEE PES Innovative Smart Grid
Technologies Europe, year=2019, publisher=IEEE

\@articlerim1990transformers, title=Transformers as equivalent circuits
for switches: General proofs and DQ transformation-based analyses,
author=Rim, Chun T and Hu, Dong Y and Cho, Gyu H, journal=IEEE
Transactions on Industry Applications, volume=26, number=4,
pages=777--785, year=1990, publisher=IEEE

\@bookanderson1995analysis, title=Analysis of faulted power systems,
author=Anderson, Paul M and Anderson, Paul M, volume=445, year=1995,
publisher=IEEE press New York
