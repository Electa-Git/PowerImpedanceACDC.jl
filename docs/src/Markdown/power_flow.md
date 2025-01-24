As the ABCD formulation is a linear representation of the power system,
nonlinear descriptions of the components such as power converters must
be linearized around an operating point. This operating point is
determined in the initialisation by solving the power flow equations
representing the combined AC/DC system.

To do so, the network is initialized using the optimal power flow tool
[@ergun2019optimal] implemented as a Julia package, which can be found
in the PowerModelsACDC repository [@ergun2018powermodelsacdc]. The
package relies on the power flow models developed for the MatACDC
simulator [@beerten2012matacdc], which extends Matpower
[@zimmerman2016matpower] AC power system models with the DC
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
[@zimmerman2016matpower], while DC branches parameters are given in
[@beerten2012matacdc].

For the purpose of modeling the system components, the model of the AC
branch as provided in [@zimmerman2016matpower] is depicted in Fig.
[1](#fig:power_flow:ac_branch){reference-type="ref"
reference="fig:power_flow:ac_branch"}. Beside the shunt admittance
$j\frac{b_c}{2}$,the Julia package PowerModelsACDC
[@ergun2018powermodelsacdc] supports the admittance as
$\frac{g_c}{2} + j\frac{b_c}{2}$. The full expression for the AC
admittance parameters is given by the equation:
$$\textbf{Y}_{ac} = \begin{bmatrix}
    \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right) \, \frac{1}{\tau^2} & -\frac{y_s}{\tau \, \exp(-j\theta_{\rm shift})} \\
    -\frac{y_s}{\tau \, \exp(-j\theta_{\rm shift})} & \left(y_s + \frac{g_c}{2} + j\frac{b_c}{2} \right) 
    \end{bmatrix}.
    \label{eq:ac_branch}$$ It should be noted that the AC network is
considered as balanced for the power flow solution, and thus, all the
components have diagonal matrix models, with equal elements on the
matrix diagonal.

![Matpower AC branch model
[@zimmerman2016matpower].](pictures/power_flow/branch_ac.png){#fig:power_flow:ac_branch
width="75%"}

A DC branch is modeled with its equivalent series resistance
[@beerten2012matacdc]. For the power flow calculation, some components
are modeled as AC and DC branches and their models are described in
detail in this subsection.

### Impedance

As described in Subsection
[\[sec:impedance\]](#sec:impedance){reference-type="ref"
reference="sec:impedance"}, an impedance is modeled using ABCD
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

### Transformer

Since the transformer model considered in
[\[sec:transformer\]](#sec:transformer){reference-type="ref"
reference="sec:transformer"} cannot be easily represented as the model
in Fig. [1](#fig:power_flow:ac_branch){reference-type="ref"
reference="fig:power_flow:ac_branch"}, Y parameters are extracted from
the ABCD parameters, as described in the Appendix in
Eq. [\[eq:abcd_to_y\]](#eq:abcd_to_y){reference-type="eqref"
reference="eq:abcd_to_y"}.

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

### Transmission line

A transmission line (OHL, cable, cross-bonded cable or mixed OHL-cable)
is represented using its nominal $\pi$-model depicted in Fig.
[2](#fig:transmission_line:equivalent_pi){reference-type="ref"
reference="fig:transmission_line:equivalent_pi"}, where
$$\begin{aligned}
\nonumber \textbf{Z}(j\omega) = \textbf{Y}_c^{-1} \sinh(\mathbf{\Gamma} l), \\
\textbf{Y}(j\omega) = \textbf{Y}_c \tanh(\mathbf{\Gamma} l).
\label{eq:tl_distributed}
\end{aligned}$$

![Nominal $\pi$-model of the transmission
line.](pictures/transmission_line/nominal-PI.pdf){#fig:transmission_line:equivalent_pi
width="50%"}

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

## Shunt components

Shunt reactors and capacitors are defined with their admittance value as
$y = g_s + jb_s$ [@zimmerman2016matpower].

## Generators

In the current version of the simulator, generators are simplified to
ideal three-phase AC sources, and are defined as reference buses.

## Power converter

A power converter is, in accordance with [@beerten2012matacdc], modeled
together with its phase reactor, filter and transformer. In order to
match the constructed MMC model in Section
[\[sec:mmc\]](#sec:mmc){reference-type="ref" reference="sec:mmc"}, only
the reactor is considered.

Losses of the converter are calculated in the form of
$P_{loss} = a + b I_c + c I_c^2$. Since the switches are modeled as
ideal, the model implementation supports only a constant value
$c = \frac{R_{arm}}{2}$.

![Power flow model of the power
converter.](pictures/power_flow/converter_pf.png){#fig:power_flow:converter
width="\\linewidth"}

Depending on the actual realisation of the converter's controls, the
parameters of the converter can be set as a DC voltage controlling or an
active power controlling converter.

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
