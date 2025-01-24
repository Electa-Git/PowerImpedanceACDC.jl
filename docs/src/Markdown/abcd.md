Classical circuit theory applied to power systems relies on the
description of the system using admittance matrices or hybrid matrices
[^1] ; @kundur1994power]. Commonly, the circuit equations
are solved using admittance/hybrid matrices applying Kirchhoff's laws
and Ohm's law relying on the Modified Nodal Analysis (MNA) approach for
the components described using linear models. The admittance based
representation is also gaining popularity for the assessment of harmonic
stability in systems with power electronic components
[@stamatiou2017analytical; @bessegato2018method].

Although the admittance representation of the system and its components
has a simple definition and a physical dimension, system and component
configurations exist without impedance or admittance parameters defined.
In Fig. [1](#fig_examples){reference-type="ref"
reference="fig_examples"}a, where the port voltage and current have
subscript "$p$" for an input port and "$s$" for the output port, a case
is depicted where the series impedance cannot be represented using
impedance parameters because of the open connection at the input and at
the output port. Similarly, the example from
Fig. [1](#fig_examples){reference-type="ref" reference="fig_examples"}b
shows the case of a shunt admittance. It cannot be described using
admittance parameters, since the short connection between the ports
would give an infinite value for the interconnection admittance. A
hybrid port representation could be applied in these cases, but its
usage for determining input and output impedance of the network is
everything but intuitive.

![Examples of the circuit in which is not defined: a) impedance matrix;
b) admittance matrix.](pictures/abcd/examples.pdf){#fig_examples
width="50%"}

To overcome the challenge of nonexistent $Z$ or $Y$ parameter
representation, we propose a generalised algorithm for representing the
power system, and its constituting components using multiport ABCD
parameters instead (Fig.
[2](#fig_multiport_network){reference-type="ref"
reference="fig_multiport_network"}). Multiport networks can include
polyphase AC networks, multi-pole DC networks, etc. The motivation for
the choice of the ABCD parameters' system representation stems from the
fact that ABCD parameters provide a direct connection between the
voltages and currents at the input ports, and voltages and currents at
the output ports.

![Polyphase power system using multiport ABCD
parameters.](pictures/abcd/multiport_network.pdf){#fig_multiport_network
width="50%"}

As explained, admittance-based representations of networks are becoming
a promising approach gaining popularity to be used in impedance-based
stability assessments to investigate harmonic stability with power
electronic converters [@liu2018oscillatory], e.g. for Voltage Source
Converter High Voltage Direct Current (VSC HVDC) systems
[@ji2019harmonic; @bayo2018control]. The use of ABCD parameters for such
a stability assessment, however, is only recently starting to get
attention in literature. Besides [@bayo2018control], where the method
was proposed to build small network models for analyzing high-frequency
interactions in the kHz-range, recently, also [@sun2018adequacy]
assessed interactions, but within the bandwidth of the converter
controllers (up to several 100 Hz), building on the work from
[@wu2014impact] to create a network equivalent containing a
frequency-dependent model of a single overhead line.

What has been missing so far, however, is a generalised modeling
framework that allows automatically constructing an equivalent
impedance, including both active and passive components, at any node in
the network. Therefore, a systematically derived modeling framework
using an ABCD representation for the multiport power system and its
components has been developed. Each component is represented using ABCD
parameters, and its corresponding model is presented in detail in this
report. The model of the DC side impedance of a state-of-the-art VSC
HVDC-based MMC is given and the complete system modeling is presented
for a two-terminal MMC-based VSC HVDC system. The report also summarises
how to use the ABCD parameters for a harmonic system stability analysis.

## ABCD parameters basics

The system is represented as an interconnection of components. To
simplify the calculation of the transfer functions and/or input and
output impedances, each component is modeled as a multiport network as
depicted in [Figure 2](#fig:closing_impedance) : 

<div style="text-align: center;">
    <img src="C:/Users/asaad/.julia/dev/hvdcstability.jl\docs\src\pictures\abcd\closing_impedance.jpg" alt="Closing Impedance" width="500">
    <p><em>Fig. 2: Closing Impedance diagram.</em></p>
</div>
The input voltages and currents are
vectors denoted as $\mathbf{V}_p$ and $\mathbf{I}_p$, while the output
voltages and currents are $\mathbf{V}_s$ and $\mathbf{I}_s$. Generally,
a multiport network has the same number of input and output ports, and
thus, the dimensions of the vectors are the same, denoted as $n$.

A multiport network can be represented with ABCD parameters, where each
of parameters $\mathbf{A}$, $\mathbf{B}$, $\mathbf{C}$ and $\mathbf{D}$
represent $n \times n$ matrices and 
```math
\begin{bmatrix}
\mathbf{V}_p \\
\mathbf{I}_p
\end{bmatrix} =
\begin{bmatrix}
\mathbf{A} & \mathbf{B} \\
\mathbf{C} & \mathbf{D}
\end{bmatrix} \times
\begin{bmatrix}
\mathbf{V}_s \\
\mathbf{I}_s
\end{bmatrix}.
```
As with components in an electrical power systems, their multiport ABCD
parameter representations can be interconnected. Two possible
connections are series and parallel connections.

-   The series connection of two multiport networks is depicted in Fig.
    [3](#fig_series_connection){reference-type="ref"
    reference="fig_series_connection"}. The ABCD multiport
    representation is especially desirable for this type of connection
    because the new parameters are determined in a simple matter as
    follows. $$\begin{bmatrix}
    \mathbf{V}_p \\
    \mathbf{I}_p
    \end{bmatrix} = \underbrace{\begin{bmatrix}
    \mathbf{A}_1 & \mathbf{B}_1 \\
    \mathbf{C}_1 & \mathbf{D}_1
    \end{bmatrix} \times
    \begin{bmatrix}
    \mathbf{A}_2 & \mathbf{B}_2 \\
    \mathbf{C}_2 & \mathbf{D}_2
    \end{bmatrix} }_{\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix}}
    \times
    \begin{bmatrix}
    \mathbf{V}_s \\
    \mathbf{I}_s
    \end{bmatrix}$$

    ![Series connected multiport
    networks.png](pictures/abcd/serial_connection.pdf){#fig_series_connection}

-   The parallel connection, depicted in Fig.
    [4](#fig_parallel_connection){reference-type="ref"
    reference="fig_parallel_connection"}, is more complex to calculate.
    In the case of nonzero matrices $\mathbf{B}_1$ and $\mathbf{B}_2$,
    the parallel connection is represented as: $$\small
    \begin{bmatrix}
    \mathbf{V}_p \\
    \mathbf{I}_p
    \end{bmatrix} = \underbrace{\begin{bmatrix}
    (\mathbf{B}_1^{-1} + \mathbf{B}_2^{-1})^{-1} \, (\mathbf{B}_1^{-1}\mathbf{A}_1 + \mathbf{B}_2^{-1}\mathbf{A}_2) & (\mathbf{B}_1^{-1} + \mathbf{B}_2^{-1})^{-1} \\
    \mathbf{C}_1 + \mathbf{C}_2 + (\mathbf{D}_2 - \mathbf{D}_1) \, (\mathbf{B}_1 + \mathbf{B}_2)^{-1} \, (\mathbf{A}_1 - \mathbf{A}_2) & \mathbf{D}_1 + (\mathbf{D}_2 - \mathbf{D}_1) \, (\mathbf{B}_1 + \mathbf{B}_2)^{-1} \, \mathbf{B}_1
    \end{bmatrix}}_{\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix}}
    \times
    \begin{bmatrix}
    \mathbf{V}_s \\
    \mathbf{I}_s
    \end{bmatrix}.$$ The formula is also valid for networks whose
    matrices are of dimension $1 \times 1$, i.e. two-port networks. If
    some of the matrices cannot be inverted, the previous equation
    becomes: $$\small
    \begin{bmatrix}
    \mathbf{V}_p \\
    \mathbf{I}_p
    \end{bmatrix} = \underbrace{\begin{bmatrix}
    \mathbf{A}_i & \mathbf{0} \\
    \mathbf{C}_1 + \mathbf{C}_2 + (\mathbf{D}_2 - \mathbf{D}_1) \,\mathbf{B}_j^{-1} \, (\mathbf{A}_1 - \mathbf{A}_2) & \mathbf{D}_i 
    \end{bmatrix}}_{\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix}}
    \times
    \begin{bmatrix}
    \mathbf{V}_s \\
    \mathbf{I}_s
    \end{bmatrix},$$ where $i,j \in \{1,2\}$ and $i$ denotes the
    invertible matrix $\mathbf{B}_i$ with $j \neq i$.

    ![Parallel connected multiport
    networks.](pictures/abcd/parallel_connection.pdf){#fig_parallel_connection}

It should be noted that the term 'multiport' network encompasses both
single component with one or more input and output ports, as well as a
subnetwork with defined input and output 'ports' (or nodes).

## Determining the input/output impedance of the network

Let us assume that every output port, represented with the voltage
$V_{si}$ and the current $I_{si}$, is closed with an impedance $Z_{ti}$.
Then we can write $V_{si} = Z_{ti} \, I_{si}$, or in matrix form:
$$\mathbf{V}_s = \mathbf{Z}_t \odot \mathbf{I}_s = \widetilde{\mathbf{Z}}_t \times \mathbf{I}_s,$$
with $\odot$ denoting the Hadamard product, $\mathbf{Z}_t$ the
corresponding closing impedance column vector, and
$\widetilde{\mathbf{Z}}_t = \operatorname{diag}\{\mathbf{Z}_t\}$ (see
Fig. [5]  (#fig:abcd:closing_impedance){reference-type="ref"
reference="fig:abcd:closing_impedance"}).

The input impedance can be then rewritten from: $$\begin{bmatrix}
\mathbf{V}_p \\
\mathbf{I}_p
\end{bmatrix} =
\begin{bmatrix}
\mathbf{A} & \mathbf{B} \\
\mathbf{C} & \mathbf{D}
\end{bmatrix} \times 
\begin{bmatrix}
\widetilde{\mathbf{Z}}_t \times \mathbf{I}_s \\
\mathbf{I}_s
\end{bmatrix},$$ as
$$\mathbf{Z}_{p} = (\mathbf{A}\times \widetilde{\mathbf{Z}}_t + \mathbf{B}) \times (\mathbf{C}\times \widetilde{\mathbf{Z}}_t + \mathbf{D})^{-1}.$$

Similarly, by closing the input ports with a diagonal impedance
$\widetilde{\mathbf{Z}}_t$, the impedance as seen from the output ports
can be estimated as:
$$\mathbf{Z}_s = (\widetilde{\mathbf{Z}}_t \times \mathbf{C} - \mathbf{A})^{-1} \times (\widetilde{\mathbf{Z}}_t \times \mathbf{D} - \mathbf{B}).$$

![Closing impedance on the output
side.](pictures/abcd/closing_impedance.png){#fig:abcd:closing_impedance
width="45%"}

## Determination of the combined system ABCD parameters

The ABCD parameters can be used to determine the impedance "visible"
from the desired node or a component port. The power system can take
forms with different number of component 'pins', e.g. (but not limited
to) a three-phase AC system using $abc$ (or phase) parameters with
components having 3 input and output pins, a positive sequence
equivalent (components with 2 input and output pins); or a monopolar
(represented by 1 input and output pin) or bipolar DC system (2 input
and 2 output pins). To determine the impedance "visible" from the
desired node/nodes in the system, the partition of the system is formed
recursively, containing only the nodes and the components included in
the path between the desired nodes. The example of the obtained
subsystem is depicted in Fig. [6](#fig_system_abcd){reference-type="ref"
reference="fig_system_abcd"}, where the nodes denoted as $V_{p1,s}$,
$V_{p2,s}$ and $V_{s1,s}$, $V_{s2,s}$ represent $2 \times 2$ port
system, respectively.

![Model of the polyphase
subsystem.](pictures/abcd/system_abcd2.pdf){#fig_system_abcd}

As can be seen from Fig. [6](#fig_system_abcd){reference-type="ref"
reference="fig_system_abcd"}, the subsystem between input and output
nodes contains $m$ components, where each component
$j \in \{1, \ldots, m\}$ is represented with $p^j_i$ inputs and $p^j_o$
outputs. The subsystem also contains a total number of $n_n$ nodes, of
which $n_o$ nodes are output nodes (denoted as $V_{s1,s}$ and $V_{s2,s}$
in Fig. [6](#fig_system_abcd){reference-type="ref"
reference="fig_system_abcd"}) or ground nodes. Subsystem has $n_c$ input
currents/voltages (denoted as $I_{p1,s}$ and $I_{p2,s}$ in Fig.
[6](#fig_system_abcd){reference-type="ref"
reference="fig_system_abcd"}).

Let us assume the following naming convention. The $i$th component has
input voltages and currents denoted as $\mathbf{V}_{pi}$ and
$\mathbf{I}_{pi}$ (positive currents enter the component and exit the
node), and output voltages and currents $\mathbf{V}_{si}$ and
$\mathbf{I}_{si}$ (positive currents exit the component and enter the
node). $\mathbf{I}_i$ are the input subsystem currents entering the
subsystem and exiting the input nodes, and $\mathbf{I}_0$ are the
currents through ground(s) and the output pins, and enter the output
nodes.

Now, the set of $n_n + \sum\limits_{i=1}^m 2p^i_p$ equations can be
written:

-   $n_n$ equations for each of the nodes inside the network denoting
    the sum of the currents entering and leaving the node. By
    convention, a current is considered positive when it leaves the
    node.

-   $\sum\limits_{i = 1}^m p^i_p$ equations giving ABCD component
    relationships between input voltages $\mathbf{V}_{pi}$ and output
    voltages $\mathbf{V}_{si}$, and currents $\mathbf{I}_{si}$.

-   $\sum\limits_{i = 1}^m p^i_p$ equations giving ABCD component
    relationships between input currents $\mathbf{I}_{pi}$ and output
    voltages $\mathbf{V}_{si}$, and currents $\mathbf{I}_{si}$.

The unknown variables are $n_v = n_n - n_o$ node voltages, $n_c$ input
currents denoted as $\mathbf{I}_i$ and
$\sum\limits_{i = 1}^m p^i_p + p^i_s$ component currents
$\mathbf{I}_{pi}$ and $\mathbf{I}_{si}$.

The complete set of equations is written in a matrix form and consists
of $n_n + \sum\limits_{i = 1}^m 2p^i_p$ equations with
${n_v + n_c + \sum\limits_{i = 1}^m (p^i_p + p^i_s)}$ variables and
matrix of outputs with the size
$\left(n_n + \sum\limits_{i = 1}^m 2p^i_p\right) \times 2n_o$. It is:
$$\mathbf{M} \times \mathbf{X} = \mathbf{N} \times \mathbf{Y}\,\,\,
(8)$$ where the matrices $\mathbf{M}$ and $\mathbf{N}$
consist of numerical and symbolic coefficients, vector
$\mathbf{X} = \begin{bmatrix}
V_1 & \cdots & V_{n_v} & I_{i1} & \cdots & I_{i n_c} & \mathbf{I}_{p1} & \mathbf{I}_{s1} & \cdots & \mathbf{I}_{pm} & \mathbf{I}_{sm}
\end{bmatrix}$ consists of the unknown variables and vector
${\mathbf{Y} = \left.{\begin{bmatrix} \mathbf{V}_{0j}, \mathbf{I}_{0j}\end{bmatrix}^T}\right|_{j=1}^{n_o}}$
of the output and ground voltages and currents. The solution of the
previous system [\[eq_system\]](#eq_system){reference-type="eqref"
reference="eq_system"} is given as reduced row echelon (or gaussian
elimination) form of concatenated matrices $[\mathbf{M}, \mathbf{N}]$.

## Transformations between abc to dqz frames

The complete AC power system is modeled in $dqz$-frame to fit the
developed power converter model. For that purpose, it was necessary to
apply $abc$ to $dqz$ transformation.

In order to transform three-phase voltages and currents from the
stationary $abc$-frame to the rotating $dqz$-frame, Park's
transformation defined as
$$\textbf{P}_{\omega_0}(t) = \frac{2}{3} \ \begin{bmatrix}
    \cos(\omega_0 t) & \cos\left(\omega_0 t-\frac{2\pi}{3}\right) & \cos\left(\omega_0 t-\frac{4\pi}{3}\right) \\
    \sin(\omega_0 t) & \sin\left(\omega_0 t-\frac{2\pi}{3}\right) & \sin\left(\omega_0 t-\frac{4\pi}{3}\right) \\
    \frac{1}{2} & \frac{1}{2} & \frac{1}{2}
    \end{bmatrix},
    \label{eq:park}$$ is employed. The inverse Park's transformation is
given as
$$\textbf{P}^{-1}_{\omega_0}(t) = \frac{3}{2} \, \textbf{P}^T_{\omega_0}(t) + \frac{1}{2} \,
    \begin{bmatrix}
    0 & 0 & 1 \\
    0 & 0 & 1 \\
    0 & 0 & 1
    \end{bmatrix}.
    \label{eq:ipark}$$

For the transformation of the admittance from the $abc$- to the
$dqz$-frame, the following theorem can be formulated.

There has been reported work [@levron2017modeling] that applies the
transformation formula only for symmetrical systems. The formula is also
successfully applied for modeling overhead lines in the $dq$-frame in
[@belikov2018integration; @darco2019time].

The following theorem presents derivation of the transformation from
$abc$-frame to $dq$-frame for the only one admittance $3 \times 3$. In
the case of the multiport parameters representation the same formula can
be used for the each of four $3 \times 3$ submatrices giving
interconnections between inputs and outputs.

::: theorem
Every $3 \times 3$ admittance in the abc-domain $\mathbf{Y}(j\omega)$
can be transformed to the $dq$-domain without loss of generality as
$\mathbf{Y}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{Y}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{Y}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}$,
where $\mathbf{a}$ is a transformation matrix defined as
$$\textbf{a} = \begin{bmatrix}
    1 & \exp(j\varphi) & \exp(2j\varphi) \\
    j & j \exp(j\varphi) & j\exp(2j\varphi) \\
    0 & 0 & 0
    \end{bmatrix},$$ for $\varphi = \frac{2\pi}{3}$.
:::

::: proof
*Proof.* Currents and voltages in the $dqz$-frame are related to the
currents and voltages in the $abc$-frame as: $$\begin{aligned}
\mathbf{i}_{dqz}(t) = \mathbf{P}_{\omega_0}(t) \, \mathbf{i}_{abc}(t), \\
\mathbf{v}_{dqz}(t) = \mathbf{P}_{\omega_0}(t) \, \mathbf{v}_{abc}(t),
\end{aligned}$$ and vice versa, $abc$ quantities can be transformed to
their $dqz$ equivalents as: $$\begin{aligned}
\mathbf{i}_{abc}(t) = \mathbf{P}^{-1}_{\omega_0}(t) \, \mathbf{i}_{dqz}(t), \\
\mathbf{v}_{abc}(t) = \mathbf{P}^{-1}_{\omega_0}(t) \, \mathbf{v}_{dqz}(t),
\end{aligned}$$ with $\omega_0$ being the angular frequency of the
rotation frame. In the previous equations, multiplications become
convolutions in the spectral domain after applying the Fourier
transform, see Appendix

[\[sec:appendix:fourier\]](#sec:appendix:fourier){reference-type="ref"
reference="sec:appendix:fourier"}.

In the spectral domain, the relation between the currents and voltages
in the $abc$-frame can be written as:
$$\mathbf{I}_{abc}(j\omega) = \mathbf{Y}(j\omega)\mathbf{V}_{abc}(j\omega).$$
This equation can be further used as:
$$\mathbf{I}_{dqz}(j\omega) = \frac{1}{2\pi} \, \mathbf{P}_{\omega_0}(j\omega) \ast (\mathbf{Y}(j\omega) \, \mathbf{V}_{abc}(j\omega)) = \frac{1}{2\pi} \, \mathbf{P}_{\omega_0}(j\omega) \ast \left(\mathbf{Y}(j\omega) \,
    \frac{1}{2\pi} \,\mathbf{P}^{-1}_{\omega_0}(j\omega) \ast \mathbf{V}_{dqz}(j\omega)\right),$$
where the Fourier transform applied to Park's transform gives:
$${P}_{\omega_0}(j\omega) = \frac{2\pi}{3} \, \left( 
        \textbf{a} \, \delta(\omega + \omega_0) + \overline{\textbf{a}} \, \delta(\omega - \omega_0) + \textbf{c} \, \delta(\omega)
        \right) \,
        (12)$$ where $\overline{\textbf{a}}$ denotes
conjugate of the matrix $\textbf{a}$ and 
$$\begin{aligned}
    \textbf{a} = \begin{bmatrix}
    1 & \exp(j\varphi) & \exp(2j\varphi) \\
    j & j \exp(j\varphi) & j\exp(2j\varphi) \\
    0 & 0 & 0
    \end{bmatrix}, \\
    \textbf{c} = \begin{bmatrix}
    0 & 0 & 0 \\
    0 & 0 & 0 \\
    1 & 1 & 1
    \end{bmatrix},   
\end{aligned}$$
 for $\varphi = \frac{2\pi}{3}$. Similarly, the Fourier
transform applied to the inverse Park's transform is given by:
$$\textbf{P}^{-1}_{\omega_0}(j\omega) = \pi \, (\textbf{a}^T \, \delta(\omega + \omega_0) + \overline{\textbf{a}}^T \, \delta(\omega - \omega_0) + 2 \textbf{c}^T \, \delta(\omega))
        $$ Denoting
$\mathbf{G}(j\omega) = \frac{1}{2\pi} \,\mathbf{P}^{-1}_{\omega_0}(j\omega) \ast \mathbf{V}_{dqz}(j\omega)$,
one can obtain: $$\begin{aligned}
\nonumber \mathbf{G}(j\omega) &=& \frac{1}{2} \, (\mathbf{a}^T \, \delta(\omega + \omega_0) + \overline{\mathbf{a}}^T \, \delta(\omega - \omega_0) + 2 \mathbf{c}^T \, \delta(\omega)) \ast \mathbf{V}_{dqz}(j\omega) = \\
& =& \frac{1}{2} \, (\mathbf{a}^T \, \mathbf{V}_{dqz}(j(\omega + \omega_0)) + \overline{\mathbf{a}}^T \, \mathbf{V}_{dqz}(j(\omega - \omega_0)) + 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(\omega)).
\end{aligned}$$ Similarly, considering that
$\mathbf{H}(j\omega) = \mathbf{Y}(j\omega) \, \mathbf{G}(j\omega)$, we
can write: $$\begin{aligned}
\nonumber \mathbf{I}_{dqz}(j\omega) &=& \frac{1}{3} \, (\mathbf{a} \, \delta(\omega + \omega_0) + \overline{\mathbf{a}} \, \delta(\omega - \omega_0) + \mathbf{c} \, \delta(\omega)) \ast
    \mathbf{H}(j\omega)  = \\
\nonumber & = &\frac{1}{3} \, (\mathbf{a} \, \mathbf{H}(j(\omega + \omega_0)) + \overline{\mathbf{a}} \, \mathbf{H}(j(\omega - \omega_0)) + \mathbf{c} \, \mathbf{H}(j\omega)) = \\
\nonumber & = & \frac{1}{3} \, (\mathbf{a} \, \mathbf{Y}(j(\omega + \omega_0)) \, \mathbf{G}(j(\omega + \omega_0)) + \overline{\mathbf{a}} \, \mathbf{Y}(j(\omega - \omega_0)) \, \mathbf{G}(j(\omega - \omega_0)) + \mathbf{c} \, \mathbf{Y}(j\omega) \, \mathbf{G}(j\omega)) = \\
\nonumber & = & \frac{1}{6} \, \mathbf{a}\, \mathbf{Y}(j(\omega+\omega_0) \, (\mathbf{a}^T \, \mathbf{V}_{dqz}(j(\omega+2\omega_0)) + \overline{\mathbf{a}}^T \, \mathbf{V}_{dqz}(j\omega) + 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j(\omega+\omega_0))) + \\
\nonumber  && +   \overline{\mathbf{a}}\, \mathbf{Y}(j(\omega-\omega_0) \, (\mathbf{a}^T \, \mathbf{V}_{dqz}(j\omega) + \overline{\mathbf{a}}^T \, \mathbf{V}_{dqz}(j(\omega-2\omega_0)) + 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j(\omega-\omega_0))) + \\
\nonumber & & +   \mathbf{c} \, \mathbf{Y}(j\omega) \, (\mathbf{a}^T \, \mathbf{V}_{dqz}(j(\omega + \omega_0)) + \overline{\mathbf{a}}^T \, \mathbf{V}_{dqz}(j(\omega-\omega_0)) + 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j\omega))).
\end{aligned}$$

Since, we need to represent everything using $d$- and $q$-components,
the corresponding expressions are zero in $dq$ except for the zero
value. These expressions are: $$\begin{aligned}
{(\mathbf{a}\, \mathbf{Y}(j(\omega+\omega_0) \, 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j(\omega+\omega_0)))_{dq} = \mathbf{0}_{2 \times 2}}, \\
(\overline{\mathbf{a}}\, \mathbf{Y}(j(\omega-\omega_0) \, 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j(\omega-\omega_0)))_{dq} = \mathbf{0}_{2 \times 2}, \\
(\mathbf{c} \, \mathbf{Y}(j\omega) \, (\mathbf{a}^T \, \mathbf{V}_{dqz}(j(\omega + \omega_0)) + \overline{\mathbf{a}}^T \, \mathbf{V}_{dqz}(j(\omega-\omega_0)) + 2 \mathbf{c}^T \, \mathbf{V}_{dqz}(j\omega)))_{dq} = \mathbf{0}_{2 \times 2}.
\end{aligned}$$ Then, $$\begin{aligned}
\nonumber \mathbf{I}_{dq}(j\omega) &=& \frac{1}{6} \, (\mathbf{a} \, \mathbf{Y}(j(\omega+\omega_0)) \, \overline{\mathbf{a}}^T + \overline{\mathbf{a}} \, \mathbf{Y}(j(\omega-\omega_0)) \, \mathbf{a}^T)_{dq} \, \mathbf{V}_{dq}(j\omega) +  \\
& & +  \frac{1}{6} \, (\mathbf{a} \, \mathbf{Y}(j(\omega+\omega_0)) \, \mathbf{a}^T)_{dq} \, \mathbf{V}_{dq}(j(\omega+2\omega_0) +  \frac{1}{6} \, (\overline{\mathbf{a}} \, \mathbf{Y}(j(\omega-\omega_0)) \, \overline{\mathbf{a}}^T)_{dq} \, \mathbf{V}_{dq}(j(\omega-2\omega_0)).
\end{aligned}$$

Since the rotation to the $dqz$-frame applies spectral symmetry (by
multiplying with sine and cosine functions) around $\omega = 0$,
$\omega = \omega_0$ and $\omega = -\omega_0$,
$\mathbf{V}_{dq}(j(\omega-2\omega_0)) = \mathbf{V}_{dq}(j(\omega+2\omega_0)) = \mathbf{V}_{dq}(j\omega)$.
Finally,
$$\mathbf{I}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{Y}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{Y}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq} \, \mathbf{V}_{dq}(j\omega)$$
and
$$\mathbf{Y}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{Y}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{Y}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}.$$ ◻
:::
It should be noted that in the case of $dqz$-representation, the
zero-sequence can be determined directly from the formula
[\[eq:theorem:idqz\]](#eq:theorem:idqz){reference-type="eqref"
reference="eq:theorem:idqz"}.

One example of the short circuit impedance of the three-phase overhead
line is depicted in Fig.
[7](#fig:ohl_transformation){reference-type="ref"
reference="fig:ohl_transformation"}, where the impedance can be seen
before and after the application of the transformation. The obtained
diagrams correspond to the waveforms obtained in [@darco2019time].

<figure id="fig:ohl_transformation">
<embed src="pictures/abcd/ohl_no_transformation.pdf" />
<p>(a)</p>
<embed src="pictures/abcd/ohl_transformation.pdf" />
<p>(b)</p>
<figcaption>Short circuit impedance of the three phase overhead line:
(a) without applied transformation; (b) with applied
transformation.</figcaption>
</figure>

The obtained formula for the admittance can be checked on a few
well-known examples:

-   A three-phase inductor set described by the formula:
    $L \, \dot{\mathbf{i}}_{abc} = \mathbf{v}_{abc}$, by using formula
    [\[eq:abcd:transformation_formula\]](#eq:abcd:transformation_formula){reference-type="eqref"
    reference="eq:abcd:transformation_formula"} transforms to
    $$\mathbf{Y}_{dq}(j\omega) = \begin{bmatrix}
            \frac{-j\omega}{L(\omega^2 -\omega_0^2)} & \frac{\omega_0}{L(\omega^2 -\omega_0^2)}\\
            -\frac{\omega_0}{L(\omega^2 -\omega_0^2)} & \frac{-j\omega}{L(\omega^2 -\omega_0^2)}
            \end{bmatrix},$$ which corresponds to the results given in
    [@rim1990transformers].

-   A three-phase capacitor set described as
    $C \dot{\mathbf{v}}_{abc} = \mathbf{i}_{abc}$ gives
    $$\mathbf{Y}_{dq}(j\omega) = \begin{bmatrix}
            j\omega C & \omega_0 C \\
            -\omega_0 C & j\omega C
            \end{bmatrix},$$ which corresponds to the results given in
    [@rim1990transformers].

In the case of ABCD parameters (given in the $abc$ domain), the same
transformation can be applied to every matrix $\mathbf{A}$,
$\mathbf{B}$, $\mathbf{C}$ and $\mathbf{D}$. The new matrices ABCD
parameters in the $dq$ domain are: $$\begin{aligned}
\mathbf{A}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{A}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{A}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}, \\
\mathbf{B}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{B}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{B}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}, \\
\mathbf{C}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{C}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{C}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}, \\
\mathbf{D}_{dq}(j\omega) = \frac{1}{3} \, ((\mathbf{a} \, \mathbf{D}(j(\omega+\omega_0))  + \overline{\mathbf{a}} \, \mathbf{D}(j(\omega-\omega_0))) \, \Re\{\mathbf{a}\}^T)_{dq}.
\end{aligned}$$

## Transformation from bipolar to equivalent monopolar representation

As power converters are modeled in this package as a three pins
components, where one pin corresponds to the DC side connection and two
pins are used for the AC-side connection represented in the $dq$-frame,
it is necessary to represent the DC network with a single-pin
components. Bipolar DC components are then reduced to their $1 \times 1$
equivalent representation.

Bipolar DC components are represented by means of ABCD parameters as:
$$\begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix} =
    \begin{bmatrix}
    a_{11} & a_{12} & b_{11} & b_{12} \\
    a_{21} & a_{22} & b_{21} & b_{22} \\
    c_{11} & c_{12} & d_{11} & d_{12} \\
    c_{21} & c_{22} & d_{21} & d_{22}
    \end{bmatrix},$$ and $$\begin{bmatrix}
    v_{p1} \\
    v_{p2} \\
    i_{p1} \\
    i_{p2}
    \end{bmatrix} =
    \begin{bmatrix}
    \mathbf{A} & \mathbf{B} \\
    \mathbf{C} & \mathbf{D}
    \end{bmatrix} \times
    \begin{bmatrix}
    v_{s1} \\
    v_{s2} \\
    i_{s1} \\
    i_{s2}
    \end{bmatrix}.$$

For balanced bipolar DC networks connected to power converters or DC
sources, the following relation is valid:
$v_{s1} = -v_{s2} = \frac{v_s}{2}$ and $i_{s1} = -i_{s2} = i_s$. At the
input of the DC network component, it is known that
$v_p = v_{p1} - v_{p2}$ and $i_{p1} = -i_{p2} = i_p$. Substituting these
relationships, the equation becomes: $$\begin{bmatrix}
    v_p \\
    i_p
    \end{bmatrix} =
    \begin{bmatrix}
    \frac{a_{11}+a_{22}-a_{12}-a_{21}}{2} & \frac{b_{11}+b_{22}-b_{12}-b_{21}}{2} \\
    \frac{c_{11}+c_{22}-c_{12}-c_{21}}{2} & \frac{d_{11}+d_{22}-d_{12}-d_{21}}{2}
    \end{bmatrix} \times
    \begin{bmatrix}
    v_s \\
    i_s
    \end{bmatrix},$$ which represents an equivalent single-pin model.

## Reduction of the ABCD matrix {#sec:abcd:reduction}

In the case of cross-bonded cable when the outer conducting layers are
grounded, it is necessary to reduce the ABCD parameters matrix. For that
purpose it can be applied ABCD matrix reduction formula.

The overall ABCD matrix is divided into parts: matrix part with the
superscript 11 should be kept (e.g. core layer of the cable), 22 should
be removed (e.g. belongs to the sheath and armor) and 12 and 21 are
their interconnections.

$$\begin{bmatrix}
        \textbf{A} & \textbf{B} \\
        \textbf{C} & \textbf{D}
        \end{bmatrix} =
        \begin{bmatrix}
        \textbf{A}_{11} & \textbf{A}_{12} & \textbf{B}_{11} & \textbf{B}_{12} \\
        \textbf{A}_{21} & \textbf{A}_{22} & \textbf{B}_{21} & \textbf{B}_{22} \\
        \textbf{C}_{11} & \textbf{C}_{12} & \textbf{D}_{11} & \textbf{D}_{12} \\
        \textbf{C}_{21} & \textbf{C}_{22} & \textbf{D}_{21} & \textbf{D}_{22} \\
        \end{bmatrix}$$

The new, reduced matrix is obtained by applying the formula:
```math
\begin{aligned}
    \widetilde{\textbf{A}} = \textbf{A}_{11} - (\textbf{A}_{12} \textbf{Z}_s + \textbf{B}_{12}) \, \textbf{E} \, (\textbf{Z}_p \, \textbf{C}_{21} + \textbf{A}_{21}), \\
    \widetilde{\textbf{B}} = \textbf{B}_{11} - (\textbf{A}_{12} \textbf{Z}_s + \textbf{B}_{12}) \, \textbf{E} \, (\textbf{Z}_p \, \textbf{D}_{21} + \textbf{B}_{21}), \\
    \widetilde{\textbf{C}} = \textbf{C}_{11} - (\textbf{C}_{12} \textbf{Z}_s + \textbf{D}_{12}) \, \textbf{E} \, (\textbf{Z}_p \, \textbf{C}_{21} + \textbf{A}_{21}), \\
    \widetilde{\textbf{D}} = \textbf{D}_{11} - (\textbf{C}_{12} \textbf{Z}_s + \textbf{D}_{12}) \, \textbf{E} \, (\textbf{Z}_p \, \textbf{D}_{21} + \textbf{B}_{21}), \\
    \textbf{E} = ((\textbf{A}_{22} \, \textbf{Z}_s + \textbf{B}_{22}) + \textbf{Z}_p \, (\textbf{C}_{22} \, \textbf{Z}_s + \textbf{D}_{22}))^{-1}, 
\end{aligned}
 ```
 for $\textbf{Z}_p$ and $\textbf{Z}_s$ being diagonal
quadratic matrices representing closing loads (impedances) of the pins
that should be reduced from the input and output side, respectively.

Newly defined ABCD parameters are: $$\begin{bmatrix}
        \widetilde{\textbf{A}} &  \widetilde{\textbf{B}} \\
        \widetilde{\textbf{C}} & \widetilde{\textbf{D}}
        \end{bmatrix} .$$

## Conclusion

Although ABCD parameters can only be properly defined when the number of
input and output nodes (voltages and currents) are the same, this
multiport representation has multiple advantages:

-   The input and output multiport impedance can be found directly.

-   There is the unique representation of the each multiport network
    using ABCD parameters. For instance, ABCD parameters are defined
    even in cases where the admittance matrix does not exist, e.g., in
    case of a infinite shunt admittance.

-   ABCD parameters operate with voltages and currents and thus, the
    values inside the ABCD matrix have a clear physical dimension and
    "meaning". This cannot be said for H (hybrid) parameters, which is
    usually used for RF and microelectronics simulations.

-   There is a unique relationship between multiport $Z$, $Y$, $H$, $S$
    and ABCD multiport parameters
    [@reveyrand2018multiport; @frickey1994conversions].

    By the definition Z parameters, or impedance parameters, provide
    relation between voltages and currents of the multiport network.
    $$\begin{bmatrix}
            \textbf{V}_p \\
            \textbf{V}_s 
            \end{bmatrix} = \textbf{Z} \times
            \begin{bmatrix}
            \textbf{I}_p \\
            \textbf{I}_s
            \end{bmatrix}$$ Similarly, Y parameters, or admittance
    parameters, give relation between currents and voltages.
    $$\begin{bmatrix}
            \textbf{I}_p \\
            \textbf{I}_s 
            \end{bmatrix} = \textbf{Y} \times
            \begin{bmatrix}
            \textbf{V}_p \\
            \textbf{V}_s
            \end{bmatrix}$$ Hybrid parameters are defined as
    $$\begin{bmatrix}
            \textbf{V}_p \\
            \textbf{I}_s 
            \end{bmatrix} = \textbf{H} \times
            \begin{bmatrix}
            \textbf{I}_p \\
            \textbf{V}_s
            \end{bmatrix}$$ and S parameters like scattering parameters,
    are defined in terms of incident $\textbf{a}$ and reflected
    $\textbf{b}$ waves: $$\begin{bmatrix}
            \textbf{a}_p \\
            \textbf{b}_p 
            \end{bmatrix} = \textbf{S} \times
            \begin{bmatrix}
            \textbf{a}_s \\
            \textbf{b}_s
            \end{bmatrix}.$$

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

[^1]: F. Milano, Power System Modelling and Scripting. Springer Science & Business Media, 2010.

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
