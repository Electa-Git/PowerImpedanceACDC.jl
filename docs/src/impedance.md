# Impedance

An impedance can be defined between $n$ input ports (nodes) and $n$ output ports (nodes). An impedance is represented as an $n \times n$matrix $\mathbf{Z}$: 
```math
\mathbf{Z} = \begin{bmatrix}Z_{11} & Z_{12} & \cdots & Z_{1n} \\\vdots & \vdots & \ddots & \vdots \\
Z_{n1} & Z_{n2} & \cdots & Z_{nn}\end{bmatrix}
``` 
An example of an impedance with two input ports and twooutput ports is given in Fig.[1].

![Model of the 2 input ports/2 output ports
impedance.](C:/Users/asaad/.julia/dev/hvdcstability.jl/docs/src/pictures/impedance/impedance_model.png)
*Fig.1: Model of the 2 input ports/2 output ports impedance.*

To represent impedances as multiport components with ABCD parameters,the following equation constructed using the Modified Nodal Analysis
(MNA) approach [^1] needs to be solved. 
```math
\small
\underbrace{\begin{bmatrix}
\underset{i}{\operatorname{diag}} \left\{\sum\limits_{j, Z_{ij} \neq 0} \frac{1}{Z_{ij}}\right\}_{n \times n} & \big| & \operatorname{diag}\{-1\}_{n \times n} \\
\hline 
\mathbf{N}_{1, n \times n} & \big| & \mathbf{0}_{n \times n} 
\end{bmatrix}}_{\mathbf{M}_1} \times
\begin{bmatrix}
\mathbf{V}_p \\
\mathbf{I}_p
\end{bmatrix} = 
\underbrace{\begin{bmatrix}
\mathbf{N}_{2, n \times n} & \big| & \mathbf{0}_{n \times n}   \\
\hline 
-\underset{i}{\operatorname{diag}} \left\{\sum\limits_{j, Z_{ji} \neq 0} \frac{1}{Z_{ji}}\right\}_{n \times n} & \big| & \operatorname{diag}\{-1\}_{n \times n}
\end{bmatrix}}_{\mathbf{M}_2} \times
\begin{bmatrix}
\mathbf{V}_s \\
\mathbf{I}_s
\end{bmatrix},
```
where matrices $\mathbf{N}_1$ and $\mathbf{N}_2$consist of $n$ rows with $n$ columns with entries at the position$(i,j)$ equal to $-\frac{1}{Z_{ji}}$ and $\frac{1}{Z_{ij}}$, for$Z_{ij} \neq 0$ and $Z_{ji} \neq 0$ (where $i$ represents row and $j$ column in impedance matrix), respectively. The solution of the previoussystem is given as $\mathbf{M}_1^{-1}\mathbf{M}_2$ if $\mathbf{M}_1$ is invertible matrix, or by determining LU decomposition and reduced row echelon form if it is not invertible.

For example, for the circuit depicted in Fig.[1], the equations would be:
```math
\begin{bmatrix}
\frac{1}{Z_{11}} + \frac{1}{Z_{12}} & 0 & \big| & -1 & 0 \\
0 & \frac{1}{Z_{21}} + \frac{1}{Z_{22}} &\big|e & 0 & -1 \\
\hline 
-\frac{1}{Z_{11}} & -\frac{1}{Z_{21}} & \big| & 0 & 0 \\
-\frac{1}{Z_{12}} & -\frac{1}{Z_{22}} & \big| & 0 & 0
\end{bmatrix} \times
\begin{bmatrix}
V_{11} \\
V_{12} \\
I_{11} \\
I_{12}
\end{bmatrix} = 
\begin{bmatrix}
\frac{1}{Z_{11}} & \frac{1}{Z_{12}} & \big| & 0 & 0 \\
\frac{1}{Z_{21}} & \frac{1}{Z_{22}}  & \big| & 0 & 0 \\
\hline 
-\left(\frac{1}{Z_{11}} + \frac{1}{Z_{21}}\right) & 0 & \big| & -1 & 0 \\
0 & -\left(\frac{1}{Z_{12}} + \frac{1}{Z_{22}}\right) & \big| & 0 & -1
\end{bmatrix} \times
\begin{bmatrix}
V_{21} \\
V_{22} \\
I_{21} \\
I_{22}
\end{bmatrix}.
```
## Impedance definition
```@meta
CurrentModule = PowerImpedanceACDC
```

```@docs
PowerImpedanceACDC.impedance
```
## References
[^1]: C.-W. Ho, A. Ruehli, and P. Brennan, "The modified nodal approach to network analysis," IEEE Transactions on Circuits and Systems, vol. 22, no. 6, pp. 504-509, 1975.