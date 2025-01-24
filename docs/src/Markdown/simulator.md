The simulator is implemented in the Julia programming language and
consists of multiple structures organised within a specific hierarchy.
The top class consists of a structure named `Network`, which is
"responsible" for the creation of the complete system. The structure
`Network` consists of the sub-classes `Nets` and `Elements`. Each `Net`
is a set of connections between elements, and the set `Elements`
presents a collection of all the elements in the `Network` (system).

An `Element` can be made as a special power system component, or as a
composite element, which can contain a number of components and their
connections. The composite element is by default defined as a `Network`.

The simulator is implemented using the following Julia packages:
SymEngine, LinearAlgebra, NLsolve and ForwardDiff
[@RevelsLubinPapamarkou2016] for symbolic and numerical calculations;
DelimitedFiles and FileIO for reading and writing in files;
DataStructures and Parameters for data types support; PGFplotsX for the
visualisation. For the steady-state initialisation, the following two
power flow packages are used, namely PowerModels and PowerModelsACDC as
well as Ipopt and JuMP.

## Files organization

All package-related files are organized in the folder `/src/Network`.

-   `HVDCstability.jl` - generates the complete simulator package. It
    calls all used packages and links all generated Julia files.

-   `compat.jl` - Julia file with the definitions that implement the
    compatibility between different Julia versions.

-   `globals.jl` - defines all global constants and functions used
    inside the package (empty at the moment).

-   `Network.jl` - generates the top structure `Network` and all its
    functions.

-   `Components` - folder containing files with the component
    definitions:

    -   `AbstractElement.jl` - generates the structure `Element`, which
        presents the type of each component in the power system.

    -   `element_types.jl` - contains paths to all defined components
        and it links them to the structure `Element.`

    -   `converter` - folder containing converter definitions:

        -   `controller.jl` - abstract structure `Controller`, which can
            be used for the desired controller definition. So far, only
            the derived `PI_controller` is implemented.

        -   `converter.jl` - defines the abstract type `Converter` and
            the set of functions that can be generally applied to the
            wide range of power converters' definitions.

        -   `MMC.jl` - creates a specific type of the power converter
            named MMC and the functions that can be only applied to
            MMCs.

    -   `impedance` - folder containing an impedance definition and the
        functions related to impedances:

        -   `impedance.jl` - implements `Impedance` structure and its
            functions.

    -   `shunt_reactor` - folder with the shunt reactor definitions:

        -   `shunt_reactor.jl` - implementation of the component
            `Shunt_reactor` and its functions.

    -   `source` - folder with source definitions:

        -   `source.jl` - defines the abstract type `Source` and the set
            of functions that can be generally applied to the wide range
            of voltage source definitions.

        -   `dc_source.jl` - defines the functions for the DC voltage
            source.

        -   `ac_source.jl` - defines the functions for the AC voltage
            source.

        -   TO DO - Can be added: current source definition.

    -   `transformer` - folder containing transformer and
        autotransformer definitions:

        -   `transformer` - defines the structure `Transformer`, which
            can be single-phase and three-phase, at this moment both in
            a YY and $\Delta$Y configuration. NOTE: Can be extended with
            other transformer implementations and/or realisations.

        -   `autotransformer.jl` - implements the `Autotransformer`
            structure that supports only one autotransformer
            implementation. For a more general definition than the one
            provided in these notes, the model should be extended.

    -   `transmission_line` - folder with transmission line
        descriptions:

        -   `transmission_line.jl` - general `Transmission_line`
            structure and functions.

        -   `overhead_line.jl` - the `Overhead_line` structure, derived
            from the `Transmission_line` structure, and the functions
            specific to OHL.

        -   `cable.jl` - the `Cable` structure derived from the
            `Transmission_line` structure, and the function specific to
            a cable.

        -   `crossbonded_cable.jl` - the `Crossbonded_cable` structure,
            derived from the\
            `Transmission_line` structure, and the function specific to
            a cross-bonded cable.

        -   `mixed_OHL_cable.jl` - the `Mixed_OHL_cable` structure,
            derived from the `Transmission_line` structure, and the
            function specific for this structure.

    -   `tools` - folder with the tools for the `Element` structure:

        -   `plot.jl` - implements bode plotting.

        -   `abcd_parameters.jl` - various functions used for ABCD
            parameter calculations and manipulations.

        -   `kron.jl` - Kron elimination for ABCD and Y parameters.

        -   `tools.jl` - merges all files inside the `tools` folder in
            order to create a proper file linking. Every time the new
            tool is generated, its path has to be added to this file.

    <!-- -->

    -   `Solvers` - folder containing power system solving
        functionalities:

        -   `solvers.jl` - file containing the path of each solving
            possibility. It represents the linking file for all the
            other files in the folder `Solvers`.

        -   `determine_impedance.jl` - implementation of the impedance
            determination function.

        -   `stability.jl` - implementation of the function used to
            determine the feedback transfer function of the
            interconnection between a power converter and the rest of
            the network.

        -   `make_abcd.jl` - makes an ABCD representation of the
            subnetwork.

        -   `make_y.jl` - makes a Y representation of the subnetwork.

-   `GUI` - Place for the GUI, to be implemented.

Package test files in folder `/test/`:

-   `runtests.jl` - initialization of the packages used for tests and
    inclusion of the test files.

-   `tests.jl` - collection of the tests to be run together with the
    values for comparison.

-   `tests` - folder containing various tests.

Documentation which is generated automatically in folder `/docs`:

-   `make.jl` - file that defines how the documentation should be
    translated from `.md` files to HTML. Also defines how HTML files
    should be connected and in which `.git` account is deployed.

-   `README.md` - read me file for the github.

-   `src` - folder containing `.md` files with the descriptions of the
    functions which are to be added in the HTML documentation.

Additional files:

-   `gen_pr.jl` - generates the `Project.toml` file from the written
    `REQUIRE` file, which contains the names of packages used inside
    this package.

-   `Project.toml` - a file that is automatically generated by calling
    `gen_pr.jl`. `REQUIRE` - all packages called inside this package
    have to be named inside this file, each in the new line (for example
    package: `PowerModels`, etc.). See the file for more information.

-   `.travis.yml` - Travis file used for checking implemented package
    functionalities by [travis-ci.com](travis-ci.com){.uri}.

-   `LICENSE`

-   `HVDCstability.pdf` - user and development manual.

-   `README.md` - github read me file.

-   `.gitignore` - files and folders that should not be uploaded to
    github.

## Network structure

As introduced, a `Network` consists of a structure defining the
collection of all components and their interconnections. The components
are defined as type `Element` and the nodes of the system are of type
dictionary (`Dict`). A `Symbol` is used to refer to each node, which
consists of the set of element pins to which it is connected.

It should be noted that the designator `gnd` is chosen as a universal
symbol for the ground. It is not necessary to initialise it inside the
network, since it is generated directly from the provided designator. In
order to provide decoupling between different circuit parts, symbols
`gnd` with a desirable suffix can be chosen (e.g. `gnd1`, `gnd2`, etc.).
Using the same designator `gnd` for multiple components means that those
components are short connected through the ground.

A `Network` is initialized with the macro `@network` as follows:

        @network begin #= ... =# end

where `#= ... =#` denotes the set of expressions describing the network.
It provides a simple domain-specific language to describe networks. The
begin/end block can hold element definitions of the form
`refdes = elementfunc(params)`, where `refdes` is the user-defined
symbol that refers to the constructed element, `elementfunc` is the
function called with the list of parameters `params` for the creation of
the element. Connection specifications are in the following form. (Note
that $==$ can also be used in place of $\longleftrightarrow$.)

``` {mathescape="" commandchars="\\\\\\{\\}"}
refdes[pin1] $\longleftrightarrow$ refdes2[pin2] $\longleftrightarrow$ MyNode
```

The part $\longleftrightarrow \text{MyNode}$ is optional, but it gives
an opportunity to give a symbolic name to the specified connection.

The following code can be used to construct a network with an impedance
(a pure resistance in this case) and a voltage source.

``` {mathescape="" commandchars="\\\\\\{\\}"}
net = @network begin
    src = dc\_source(V = 5)
    r = impedance(z = 1000, pins = 1)
    src[1.1] $\longleftrightarrow$ r[1.1]
    src[2.1] $\longleftrightarrow$ r[2.1]
end
```

Alternatively, connection specifications can be given after an element
specification, separated by commas. In that case, the mention of
`refdes` may be omitted, defaulting to the current element, as follows:

``` {mathescape="" commandchars="\\\\\\{\\}"}
@network begin
    src = dc\_source(V = 5)
    r = impedance(z = 1000, pins = 1), src[1.1] $\longleftrightarrow$ [1.1], src[2.1] $\longleftrightarrow$ [2.1]
end
```

Finally, a connection endpoint may simply be in the form of the symbolic
`netname`, to connect to a named net. (Such named nets are created as
needed.)

``` {mathescape="" commandchars="\\\\\\{\\}"}
@network begin
    src = dc\_source(V = 5), [2.1] $\longleftrightarrow$ gnd
    r = impedance(z = 1000, pins = 1), [1.1] $\longleftrightarrow$ src[1.1], [2.1] $\longleftrightarrow$ gnd
end
```

### Add and delete components

Different power system components are added as type `Element` in the
`Network`. Functions that add and delete components are as follows.

-   `add!(n::Network, elem::Element)`

    Adds the element `elem` to the network `n`, creating and returning a
    new, unique reference designator [for the
    element]{style="color: blue"} (not for the network, whose designator
    remains the same). The pins of the element are left disconnected.

-   `add!(n::Network, designator::Symbol, elem::Element)`

    Adds the element `elem` to the network `n` with the reference
    designator `designator`, leaving its pins unconnected. If the
    network already contains an element with the same designator, it is
    removed first.

-   `delete!(n::Network, designator::Symbol)`

    Deletes the element that corresponds to this designator from the
    network `n` (disconnecting all its pins). The element is completely
    removed from the memory.

-   `connect!(n::Network, pins::Union``Symbol,Tuple``Symbol,Any``...)`

    Connects the given pins (or named nets - set of nodes) to each other
    in the network `n`. Named nets are given as `Symbol`s, pins are
    given as `Tuple{Symbols,Any}`, where the first entry is the
    reference `designator` of an element in `n`, and the second entry is
    the pin name. For convenience, the latter is automatically converted
    to a `Symbol` as needed.

-   `disconnect!(n::Network, p::Tuple``Symbol,Symbol``)`

    Disconnects the given pin `p` from anything else in the network `n`.
    The pin is given as a `Tuple{Symbols,Any}`, where the first entry is
    the reference `designator` of an element in `n`, and the second
    entry is the pin name. For convenience, the latter is automatically
    converted to a `Symbol` as needed. Note that if e.g. three pins
    `p1`, `p2`, and `p3` are connected then `disconnect!(n, p1)` will
    disconnect `p1` from `p2` and `p3`, but leave `p2` and `p3`
    connected to each other.

-   `composite_element(subnet::Network, input_pins::Array{Any},`\
    `output_pins::Array{Any})`

    Creates a network element from the (sub-)network `net`. The
    `input_pins` and `output_pins` define input and output nodes of the
    element.

Using functions `add!` and `connect!`.

``` {mathescape="" commandchars="\\\\\\{\\}"}
network = Network()
add!(network, :r, impedance(z = 1e3, pins = 1))
add!(network, :src, dc\_source(V = 5))
connect!(network, (:src, 2.1), (:r, 2.1), :gnd) # connect to gnd node
```

Creating element from the network.

``` {mathescape="" commandchars="\\\\\\{\\}"}
my\_network = @network begin
   r1 = impedance(z = 10e3, pins = 1)
   r2 = impedance(z = 10e3, pins = 1), [1.1] == r1[2.1]
   c = impedance(z = 10e3, pins = 1), [1.1] == r2[1.1], [2.1] == r2[2.1]
   src = dc\_source(V = 5), [1.1] == r1[1.1], [2.1] == r2[2.1]
end
composite\_element(my\_network, Any[(:r2, Symbol(1.1))], Any[(:r2, Symbol(2.1))])
```

### Additional checks

After the construction of a network `n`, at the end of the macro
`n = ``network begin #= ... =# end`, the program calls the following two
functions to check if the network is well connected and to construct the
ABCD parameter equivalents for all the components. The following
functions are only called internally.

-   `check_lumped_elements(network :: Network)`

    This function checks if the network is well connected, e.g. that
    there is no lumped (i.e. disconnected) pins in the network. If the
    network is not well connected, the user receives an error message.

    **Note:** If in the future use and development this function is not
    needed anymore, it can be removed from the calling at the end of
    `@network` function. The directive that should be deleted then is
    `push!(ccode.args, :(check_lumped_elements(network)))`

-   `power_flow(network :: Network)`

    This function is called in order to set the operating points of the
    converters. It generates a dictionary with the syntax used in
    [@ergun2018powermodelsacdc] and it solves the power flow problem
    using the Julia package PowerModelsACDC [@ergun2018powermodelsacdc].
    The results are used to update the operating point of the
    converters.

    If there is not a converter in the network, the function is not
    called.

### Impedance determination and stability assessment

-   ` function determine_impedance(network::Network; input_pins :: Array{Any},`\
    `output_pins :: Array{Any}, elim_elements :: Array{Symbol},`\
    `omega_range = (-3, 5, 100), , parameters_type = :ABCD)`

    For the selected multiport, it is possible to determine the
    impedance using the procedure described in Section
    [\[sec_multiport\]](#sec_multiport){reference-type="ref"
    reference="sec_multiport"}. A multiport is depicted in Fig.
    [\[fig_multiport_network\]](#fig_multiport_network){reference-type="ref"
    reference="fig_multiport_network"}.

    Input or outputs pins can be connected to some elements, which
    should not be considered for the impedance estimation. Those
    elements are listed as symbols in `elim_elements`.

    The function generates the impedance as seen from the port defined
    by input and output pins. The impedance is calculated numerically at
    each frequency point along a user-defined frequency range. The
    number of frequency points is user-defined or else set to 1000 by
    default.

    Additionally, it can be chosen how the network will be solved: using
    ABCD or Y parameters. Defining `parameters_type` as `:ABCD` sets the
    solving using ABCD parameters, while setting it to `:Y` defines the
    usage of Y parameters. In these cases, the function
    `determine_impedance` internally calls either function `make_abcd`
    or `make_y`.

    The specification for the determination of an impedance is given in
    the example, where the network consists of a DC voltage source and a
    cable.

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    net = @network begin
        vs = dc\_source(V = 500e3)
        c = cable(length = 100e3, positions = [(0,1)], earth\_parameters = (1,1,1),
        C1 = Conductor($r_o$ = 24.25e-3, $\rho$ = 1.72e-8), 
        C2 = Conductor($r_i$ = 41.75e-3, $r_o$ = 46.25e-3, $\rho$ = 22e-8),
        C3 = Conductor($r_i$ = 49.75e-3, $r_o$ = 60.55e-3, $\rho$ = 18e-8, $\mu_r$ = 10),
        I1 = Insulator($r_i$ = 24.25e-3, $r_o$ = 41.75e-3, $\epsilon_r$ = 2.3),
        I2 = Insulator($r_i$ = 46.25e-3, $r_o$ = 49.75e-3, $\epsilon_r$ = 2.3),
        I3 = Insulator($r_i$ = 60.55e-3, $r_o$ = 65.75e-3, $\epsilon_r$ = 2.3))
        vs[1.1] $\longleftrightarrow$ c[1.1] $\longleftrightarrow$ Node1
        vs[2.1] $\longleftrightarrow$ c[2.1] $\longleftrightarrow$ gnd
    end
    ```

    To determine the impedance visible from the voltage source `vs`, the
    following command should be called:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    imp, omega = determine\_impedance(net, elim\_elements = [:vs], 
            input\_pins = Any[:Node1], output\_pins = Any[:gnd],
            omega\_range = (-1,6,10000))
    ```

    The impedance is determined inside the network `net`, from the
    element `vs` (without including the element `vs`) and the port
    defined with `input_pins` as array consisting of `Node1` and the
    `output_pins` containing array with `gnd`. The impedance is
    estimated over the frequency range $10^{-1}$ \[rad/s\] to $10^{6}$
    \[rad/s\] in 10000 points.

    The function returns two complex arrays: the first array being the
    impedance array and the second array containing the angular
    frequencies at which the impedance is calculated.

-   `function check_stability(net :: Network, mmc :: Element,`\
    `direction :: Symbol = :dc)`

    This function determines two impedances inside the network, from
    which it forms the feedback transfer function. It allows "cutting"
    the power network next to the converter on its dc or ac side
    (determined by `direction` parameter). Afterwards, it checks the
    impedance $Z_{conv}$, obtained by 'looking' into the converter as
    well as the other impedance $Z_h$ as seen from the converter when
    'looking' into the rest of the system.

    Using the previous two impedances, the feedback transfer function is
    estimated as $Z_h  Y_{conv}$.

    The impedances are calculated for the angular frequencies whose
    range is defined by `omega_range`.

### Saving and plotting data

Any component in the network has its own data representation. Component
data can be saved and plotted using the following functions.

-   `save_data(element :: Element, file_name :: String;`\
    `omega_range = (-3, 5, 1000), scale = :log)`

    This function saves component specific data in a csv textual file.
    The data is saved as frequency dependent. Thus, an additional
    parameter `omega_range` provides the possibility to manually add the
    frequency range for saving the data. The scale can be given as
    logarithmic with `:log` and linear with `:lin`.

-   `plot_data(element :: Element; omega_range = (-3, 5, 1000), scale = :log)`

    This function is a component-defined function for plotting data, it
    differs from component to component. It plots the component defined
    data with the desired frequency range. An additional parameter
    `omega_range` provides the possibility to manually add the frequency
    scale for saving data. The scale can be given as logarithmic with
    `:log` and linear with `:lin`.

-   `bode(transfer_function :: Array{Any}; omega_range = (-3, 5, 100),`\
    `titles :: Array{String} = [""], omega = [], axis_type = :loglog,`\
    `save_data = false)`

    This function is used for plotting the transfer function or
    frequency dependent data as a Bode plot. The function takes
    frequency points in the form of `omega_range` or as mapped values
    `omega`. For a nice display, the labels can be given as strings by
    calling the parameter `titles`.

    It can be specified how to plot data:

    -   `:loglog` - logarithmic frequency scale and logarithmic
        impedance in dB;

    -   `:linlog` - linear frequency scale and logarithmic impedance in
        dB;

    -   `:loglin` -logarithmic frequency scale and linear impedance
        (magnitude);

    -   `:logrealimag` - logarithmic frequency scale and real/imaginary
        part;

    -   `:linrealimag` - linear frequency scale and real/imaginary part.

    Data can be saved by setting `save_data = true`.

## Component definition

### Impedance

An impedance is constructed as an element by calling the function:\
`impedance(;z :: Union{Int, Float64, Basic, Array{Basic}} = 0, pins :: Int = 0)`.

The function creates an impedance with the specified number of
input/output pins. The number of input pins is equal to the number of
output pins. The impedance expression `exp` has to be given in Ohms and
can have both numerical and symbolic values (example: `z = s-2`). Pins
are named: 1.1, 1.2, \..., 1.n and 2.1, 2.2, \..., 2.n, where n is a
number of input/output pins.

In the case of a $1 \times 1$ impedance, the parameter `z` has only one
value. Example: `impedance(z = 1000, pins = 1)`.

If the impedance is multiport (i.e. if it has a number of input pins and
output pins greater than 1), then its value is given as an array `[val]`
with one, $n$ or $n \times n$ number of values. When the arrway as
length 1, then the impedance is defined with only diagonal nonzero
values equal to val. When the array has a length equal to the number of
pins, the impedance has only diagonal nonzero values equal to the values
in the array. When the array is of size $n^2 \times 1$, the impedance
matrix is of size $n \times n$ and all its values are defined
accordingly to the array.\

``` {mathescape="" commandchars="\\\\\\{\\}"}
# 3$\times$3 impedance with diagonal values equal to s
impedance(z = [s], pins = 3)     
# 3$\times$3 impedance with diagonal values equal 2, s, 0.5s, respectively
impedance(z = [2,s,s/2], pins = 3) 
# 2$\times$2 impedance with all values defined
impedance(z = [1,s,3,4], pins = 2) 
```

In this example, `s` is the Laplace operator. The impedance can be added
in the form of a transfer function, which can be either a polynomial or
a nonlinear function of `s` and even a noninteger powers of `s`, e.g. a
pure time delay may be represented by an exponential function of `s`.

The code below can be used to define the network with its impedances, as
depicted in Fig.
[\[fig:impedance_example\]](#fig:impedance_example){reference-type="ref"
reference="fig:impedance_example"}, and to generate its bode plot.

``` {mathescape="" commandchars="\\\\\\{\\}"}
using SymEngine

s = symbols("s")
net = @network begin
    vs = dc\_source(V = 5)
    z1 = impedance(z = s+2, pins = 1)
    z2 = impedance(z = s, pins = 1)
    z3 = impedance(z = s, pins = 1)

    vs[1.1] $\longleftrightarrow$ z1[1.1] $\longleftrightarrow$ Node1
    z1[2.1] $\longleftrightarrow$ z2[1.1] $\longleftrightarrow$ z3[1.1]
    vs[2.1] $\longleftrightarrow$ z2[2.1] $\longleftrightarrow$ z3[2.1] $\longleftrightarrow$ gnd
end
imp, omega = determine\_impedance(net, elim\_elements = [:vs],
                            input\_pins = Any[:Node1], output\_pins= Any[:gnd])
bode(imp, omega = omega)
```

<figure id="fig:impedance_example_bode">
<embed src="pictures/impedance/test_example_1.pdf"
style="width:65.0%" />
<embed src="pictures/impedance/impedance_test.pdf"
style="width:80.0%" />
<figcaption>Magnitude and a phase of the equivalent
impedance.</figcaption>
</figure>

The equivalent impedance magnitude and phase are depicted in Fig.
[1](#fig:impedance_example_bode){reference-type="ref"
reference="fig:impedance_example_bode"}, which corresponds to the fact
that the impedance seen from the DC source is equal to $2+1.5s$.

### Transformer

A transformer element is created using the function
`transformer(;args...)`, which creates a single phase transformer or a
three-phase transformer in a `YY` or $\Delta$`Y` configuration.

The component is defined using the following structure.

``` {mathescape="" commandchars="\\\\\\{\\}"}
@with\_kw mutable struct Transformer
    value :: Array\{Basic\} = []         # ABCD value

    pins :: Int = 1                    # marks single or three phase
    organization :: Symbol = :YY       # three phase organization (:YY or :$\Delta$Y)

    $\omega$ :: Union\{Int, Float64\} = 2*$\pi$*50 # rated frequency in [Hz]
    $V_1^o$ :: Union\{Int, Float64\} = 0      # open circuit primary voltage [V]
    $V_1^s$ :: Union\{Int, Float64\} = 0      # short circuit primary voltage [V]
    $I_1^o$ :: Union\{Int, Float64\} = 0      # open circuit primary current [V]
    $I_1^s$ :: Union\{Int, Float64\} = 0      # short circuit primary current [V]
    $P_1^o$ :: Union\{Int, Float64\} = 0      # open circuit losses on primary side [W]
    $P_1^s$ :: Union\{Int, Float64\} = 0      # short circuit losses on primary side [W]
    $V_2^o$ :: Union\{Int, Float64\} = 0      # open circuit secondary voltage [V]
    $V_2^s$ :: Union\{Int, Float64\} = 0      # short circuit secondary voltage [V]

    $n$  :: Union\{Int, Float64\} = 0       # turn ratio
    $L_p$ :: Union\{Int, Float64\} = 0      # primary side inductance [H]
    $R_p$ :: Union\{Int, Float64\} = 0      # primary side resistance [$\Omega$]
    $R_s$ :: Union\{Int, Float64\} = 0      # secondary side resistance [$\Omega$]
    $L_s$ :: Union\{Int, Float64\} = 0      # secondary side inductance [H]
    $L_m$ :: Union\{Int, Float64\} = 0      # magnetising inductance [H]
    $R_m$ :: Union\{Int, Float64\} = 0      # magnetising resistance [$\Omega$]
    $C_t$ :: Union\{Int, Float64\} = 0      # turn-to-turn capacitance [F]
    $C_s$ :: Union\{Int, Float64\} = 0      # stray capacitance [F]
end
```

The pins are defined as : 1.1, 2.1 for single phase transformers and as:
1.1, 1.2, 1.3, 2.1, 2.2, 2.3 for three-phase transformers.\
For a circuit containing a single-phase transformer and an inductive
load, the circuit definition is given in the following listing.

``` {mathescape="" commandchars="\\\\\\{\\}"}
s = symbols("s")
net = @network begin
    vs = dc\_source(V = 5e3)
    t = transformer($V_1^o$ = 2.4e3, $V_1^s$ = 51.87, $V_2^o$ = 240, $P_1^o$ = 171.1, $P_1^s$ = 642.1,
                    $I_1^o$ = 0.48, $I_1^s$ = 20.83, $C_s$ = 12e-6, $C_t$ = 7e-6)
    z = impedance(pins = 1, z = 25e-3*s)

    vs[1.1] $\longleftrightarrow$ t[1.1] $\longleftrightarrow$ Node1
    t[2.1] $\longleftrightarrow$ z[1.1]
    vs[2.1] $\longleftrightarrow$ z[2.1] $\longleftrightarrow$ gnd
end
```

The impedance seen from the voltage source is depicted in Fig.
[2](#fig:transformer_example_bode){reference-type="ref"
reference="fig:transformer_example_bode"}.

![Magnitude and a phase of the equivalent impedance of the network
containing a
transformer.](pictures/transformer/transformer_example.pdf){#fig:transformer_example_bode
width="80%"}

### Overhead line

An overhead line is defined with a structure and added as a field
`element_value` inside an Element. It is called with the function:
`overhead_line(;args...)`.

The function `overhead_line` generates the element `elem` with the
`element_value` of the type `Transmission_line`. Arguments should be
given according to the `Overhead_line` fields:

1.  length - line length \[m\]

2.  conductors - defined in the following structure:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    struct Conductors
        # number of bundles (phases)
        $n^b$ :: Int = 1            
        # number of subconductors per bundle
        $n^{sb}$ :: Int = 1                
        # height above the ground of the lowest bundle  [m]
        $y^{bc}$ :: Union\{Int, Float64\}  = 0     
        # vertical offset between the bundles [m]
        $\Delta y^{bc}$ :: Union\{Int, Float64\} = 0    
        # horizontal offset between the lowest bundles  [m]
        $\Delta x^{bc}$ :: Union\{Int, Float64\} = 0    
        # horizontal offset in group of bundles    [m]
        $\Delta \widetilde{x}^{bc}$ :: Union\{Int, Float64\} = 0  
        # sag offset    [m]
        $d^{sag}$ :: Union\{Int, Float64\} = 0   
        # subconductor spacing (symmetric)  [m]
        $d^{sb}$ :: Union\{Int, Float64\} = 0 
        # conductor radius  [m]
        $r^c$ :: Union\{Int, Float64\} = 0    
        # DC resistance for the entire conductor [$\Omega$/m]
        $R^{dc}$ :: Union\{Int, Float64\} = 0
        # shunt conductance
        $g^c$ :: Union\{Int, Float64\} = 1e-11   
        # relative conductor permeability
        $\mu_r^c$ :: Union\{Int, Float64\} = 1     
        # add absolute positions manually
        positions :: Tuple\{Vector\{Union\{Int, Float64\}\}, 
                        Vector\{Union\{Int, Float64\}\}\} = ([],[])   
        # organization can be :flat, :vertical, :delta, :concentric, :offset
        organization :: Symbol = Symbol()
    end
    ```

3.  groundwires - defined in the following structure:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    struct Groundwires
        # number of groundwires (typically 0 or 2)
        $n^g$ :: Int = 0                 
        # horizontal offset between groundwires [m]
        $\Delta x^g$ :: Union\{Int, Float64\} = 0  
        # vertical offset between the lowest conductor and groundwires  [m]
        $\Delta y^g$ :: Union\{Int, Float64\} = 0 
        # ground wire radius  [m]
        $r^g$ :: Union\{Int, Float64\} = 0    
        # sag offset [m]
        $d^{gsag}$ ::  Union\{Int, Float64\} = 0    
        # groundwire DC resistance [$\Omega$/m]
        $R^{gdc}$ :: Union\{Int, Float64\} = 0    
        # relative groundwire permeability
        $\mu_r^g$ :: Union\{Int, Float64\} = 1       
        # add absolute positions manually
        positions :: Tuple\{Vector\{Union\{Int, Float64\}\}, 
                    Vector\{Union\{Int, Float64\}\}\} = ([],[])    
        eliminate :: Bool = true
    end
    ```

4.  `earth_parameters` - with default value (1,1,1) and defining
    $(\mu_r, \epsilon_r, \rho)$ in units of (\[\], \[\], \[$\Omega$m\])

``` {mathescape="" commandchars="\\\\\\{\\}"}
overhead\_line(length = 227e3, earth\_parameters = (1,1,100),
    conductors = Conductors($n^b$ = 2, $n^{sb}$ = 2, organization = :flat,
    $R^{dc}$ = 0.06266, $r^c$ = 0.01436, $y^{bc}$ = 27.5, $\Delta x^{bc}$ = 11.8, $d^{sb}$ = 0.4572, $d^{sag}$ = 10),
    groundwires = Groundwires($n^g$ = 2, $\Delta x^g$ = 6.5, $\Delta y^g$ = 7.5, $R^{gdc}$ = 0.9196, 
    $r^g$ = 0.0062, $d^{gsag}$ = 10))
```

The short-circuit impedance for the transmission line defined in the
previous listing is plotted in Fig.
[3](#fig_example_tl){reference-type="ref" reference="fig_example_tl"}.

![Overhead line at the DC side
example.](pictures/transmission_line/transmission_line_example.pdf){#fig_example_tl
width="\\linewidth"}

### Cable

A cable is defined with a structure and added as a field `element_value`
inside an Element. It is called with the function: `cable(;args...)`.

The function `cable()` generates the element `elem` with the
`element_value` of the type `Cable`. Arguments should be given according
to the `Cable` fields:

1.  length - length of the cable in \[m\]

2.  `earth_parameters` - with default value (1,1,1) and defining
    $(\mu_r\_earth, \epsilon_r\_earth, \rho\_earth)$ in units of (\[\],
    \[\], \[$\Omega$m\]), i.e. ground (earth) relative premeability,
    relative permittivity and ground resistivity

3.  conductors - dictionary with the key symbol being: C1, C2, C3 or C4,
    and the value given with the structure `Conductor`. If the sheath
    consists of a metallic screen and a sheath, then a screen must be
    added with a key symbol SC and a sheath with key symbol C2.

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    struct Conductor
        $r_i$ :: Union\{Int, Float64\} = 0           inner radius
        $r_o$ :: Union\{Int, Float64\} = 0            # outer radius
        $\rho$  :: Union\{Int, Float64\} = 0            # conductor resistivity [$\Omega$m]
        $\mu_r$  :: Union\{Int, Float64\} = 1           # relative permeability
        
        $A$ :: Union\{Int, Float64\} = 0             # nominal area
    end
    ```

4.  insulators - a dictionary with as key symbols: I1, I2, I3 and I4,
    and the value given with the structure `Insulator`. For the 2nd
    insulator, the semiconducting layers can be added by specifying the
    outer radius of the inner semiconducting layer and the inner radius
    of the outer semiconducting layer.

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    struct Insulator
        $r_i$ :: Union\{Int, Float64\} = 0               # inner radius
        $r_o$ :: Union\{Int, Float64\} = 0               # outer radius
        $\epsilon_r$ :: Union\{Int, Float64\} = 1               # relative permittivity
        $\mu_r$ :: Union\{Int, Float64\} = 1               # relative permeability
        
        # inner semiconductor outer radius
        $a$ :: Union\{Int, Float64\} = 0    
        # outer semiconductor inner radius
        $b$ :: Union\{Int, Float64\} = 0                
    end
    ```

5.  positions - given as an array in (x,y) format

6.  type - symbol representing aerial or underground cable

7.  configuration - symbol with two possible values: coaxial (default)
    and pipe-type

8.  eliminate - indicator weather the shaeth and armor layers should be
    grounded

``` {mathescape="" commandchars="\\\\\\{\\}"}
cable(length = 100e3, positions = [(0,1)], 
    C1 = Conductor($r_o$ = 24.25e-3, $\rho$ = 1.72e-8), 
    C2 = Conductor($r_i$ = 41.75e-3, $r_o$ = 46.25e-3, $\rho$ = 22e-8),
    C3 = Conductor($r_i$ = 49.75e-3, $r_o$ = 60.55e-3, $\rho$ = 18e-8, $\mu_r$ = 10),
    I1 = Insulator($r_i$ = 24.25e-3, $r_o$ = 41.75e-3, $\epsilon_r$ = 2.3),
    I2 = Insulator($r_i$ = 46.25e-3, $r_o$ = 49.75e-3, $\epsilon_r$ = 2.3),
    I3 = Insulator($r_i$ = 60.55e-3, $r_o$ = 65.75e-3, $\epsilon_r$ = 2.3))
```

A Bode plot is provided in Fig.
[4](#fig_example_cable){reference-type="ref"
reference="fig_example_cable"} for the short-connected cable defined in
the previous listing with a length of 100 km.

![Cable example for line length 100
km.](pictures/transmission_line/cable_example.pdf){#fig_example_cable
width="\\linewidth"}

### MMC {#sec:mmc_simulator}

An MMC is implemented as an `Element` using the function
`mmc(;args...)`. The field `element_value` of the element is defined
according to the following structure:

``` {mathescape="" commandchars="\\\\\\{\\}"}
$\omega_0$ :: Union\{Int, Float64\} = 100*$\pi$ # nominal angular frequency

P :: Union\{Int, Float64\} = -10              # active power [MW]
Q :: Union\{Int, Float64\} = 3                # reactive power [MVA]
P\_dc :: Union\{Int, Float64\} = 100           # DC power [kW]
P\_min :: Union\{Float64, Int\} = -100         # min active power output [MW]
P\_max :: Union\{Float64, Int\} = 100          # max active power output [MW]
Q\_min :: Union\{Float64, Int\} = -50          # min reactive power output [MVA]
Q\_max :: Union\{Float64, Int\} = 50           # max reactive power output [MVA]

$\theta$ :: Union\{Int, Float64\} = 0
$V_m$ :: Union\{Int, Float64\} = 333             # AC voltage, amplitude [kV]
$V^{dc}$ :: Union\{Int, Float64\} = 640            # DC-bus voltage [kV]

$L_{arm}$ :: Union\{Int, Float64\}  = 50e-3        # arm inductance [H]
$R_{arm}$ :: Union\{Int, Float64\}  = 1.07         # equivalent arm resistance [$\Omega$]
$C_{arm}$ :: Union\{Int, Float64\}  = 10e-3        # capacitance per submodule [F]
$N$ :: Int = 400                             # number of submodules per arm

$L_r$ :: Union\{Int, Float64\}  = 60e-3          # inductance of the phase reactor [H]
$R_r$ :: Union\{Int, Float64\}  = 0.535          # resistance of the phase reactor [$\Omega$]

# used inside the functions
controls :: OrderedDict\{Symbol, Controller\} = OrderedDict\{Symbol, Controller\}()
equilibrium :: Array\{Union\{Int, Float64\}\} = [0]
A :: Array\{Union\{Int, Float64\}\} = [0]
B :: Array\{Union\{Int, Float64\}\} = [0]
C :: Array\{Union\{Int, Float64\}\} = [0]
D :: Array\{Union\{Int, Float64\}\} = [0]
```

The constructed MMC has two pins on the AC side: 2.1 and 2.2, and one
pin on its DC-side: 1.1 and 2.2. The component is described using two
ABCD parameters with a matrix of size $4 \times 4$.

The controls are defined as `PI_control` and the keyword `occ` is for
output current control, `ccc` for circulating current control, `zcc` for
zero current control, `power` for active and reactive current control,
`energy` for zero energy control and `dc` for the DC voltage control.\

The following example demonstrates the implementation of an MMC with
output current control, circulating current control and power and energy
controls. It is assumed that $P = 1000 \text{ MW}$, $Q = 0$,
$V_m = 320 \text{ kV}$ and $\theta = 0$.

``` {mathescape="" commandchars="\\\\\\{\\}"}
mmc(energy = PI\_control($K_p$ = 120, $K_i$ = 400),
    occ = PI\_control($\zeta$ = 0.7, bandwidth = 1000),
    ccc = PI\_control($\zeta$ = 0.7, bandwidth = 300),
    zcc = PI\_control($\zeta$ = 0.7, bandwidth = 300),
    power = PI\_control($K_p$ = 2.0020e-07, $K_i$ = 1.0010e-04))
```

The admittances of the MMC can be displayed using the function
`plot_data` and are presented in Figs.
[5](#fig_mmc_ac){reference-type="ref" reference="fig_mmc_ac"},
[6](#fig_mmc_dc_side){reference-type="ref" reference="fig_mmc_dc_side"}
and [7](#fig_mmc_acdc){reference-type="ref" reference="fig_mmc_acdc"}.

![AC side admittance of the MMC.](pictures/mmc/ac_side.pdf){#fig_mmc_ac
width="65%"}

![DC side impedance of the
MMC.](pictures/mmc/dc_side.pdf){#fig_mmc_dc_side width="65%"}

![Admittance interconnection between AC and DC
side.](pictures/mmc/acdc_side.pdf){#fig_mmc_acdc width="65%"}

In order to demonstrate the functionality of the simulator, multiple
diagrams are generated for the same test example, but with base values
$P = 1000 \text{ MVA}$, $Q = 0$, $V_{m} = 333 \text{ kV}$ and
$\theta = 0$.

In the case when the impedance between two out of three pins denoted as
$x$ and $y$ ($x,y \in \{d,q,z\}$) is estimated, while the third pin is
short connected to the ground, the expected impedance is
$Z_{eq} = \frac{1}{Y_{xx}}$, where
$Y_{xx} \in \{ Y_{dd}, Y_{qq}, 3Y_{zz}\}$. The same result is obtained
in the simulation and depicted in Fig.
[8](#fig:mmc_example_sc_impedance){reference-type="ref"
reference="fig:mmc_example_sc_impedance"}.

For the case when one pin is considered as "open" and the impedance is
determined between other pins, the simulated results are depicted in
Fig. [9](#fig:mmc_example_oc_impedance){reference-type="ref"
reference="fig:mmc_example_oc_impedance"}. The impedance between pins
$x$ ($x = d,q$) and $z$ for the short connected third pin, the expected
impedance is:
$$Z_{eq} = \dfrac{1}{3\left(Y_{zz} - \frac{Y_{zx} Y_{xz}}{Y_{xx}} \right)}.$$
Impedance between pins $x$ and $y$ for $x,y = d,q$ with open DC pin $z$
is $$Z_{eq} = \dfrac{1}{Y_{xx} - \frac{Y_{zx} Y_{xz}}{Y_{zz}}}.$$

<figure id="fig:mmc_example_sc_impedance">
<embed src="pictures/mmc/impedance_dq.pdf" />
<p>(a)</p>
<embed src="pictures/mmc/impedance_zz.pdf" />
<p>(c)</p>
<embed src="pictures/mmc/impedance_qd.pdf" />
<p>(b)</p>
<figcaption>Impedance between two pins when the third is short
connected: (a) between d and q pins; (b) between q and d pins; and (c)
between dc pin and short connected d and q pins.</figcaption>
</figure>

<figure id="fig:mmc_example_oc_impedance">
<embed src="pictures/mmc/impedance_dq_z_open.pdf" />
<p>(a)</p>
<embed src="pictures/mmc/impedance_zd_q_open.pdf" />
<p>(c)</p>
<embed src="pictures/mmc/impedance_qd_z_open.pdf" />
<p>(b)</p>
<embed src="pictures/mmc/impedance_zq_d_open.pdf" />
<p>(d)</p>
<figcaption>Impedance between two pins when the third is open
connection: (a) between d and q pins; (b) between q and d pins; (c)
between z and d; and (d) between z and q.</figcaption>
</figure>

An example for the simulation of the MMC with an incorporated PLL is
given with the following listing.

``` {mathescape="" commandchars="\\\\\\{\\}"}
mmc(energy = PI\_control($K_p$ = 120, $K_i$ = 400),
    occ = PI\_control($\zeta$ = 0.7, bandwidth = 1000),
    ccc = PI\_control($\zeta$ = 0.7, bandwidth = 300),
    zcc = PI\_control($\zeta$ = 0.7, bandwidth = 300),
    power = PI\_control($K_p$ = 2.0020e-07, $K_i$ = 1.0010e-04),
    pll = PI\_control($K_p$ = 2e-3, $K_i$ = 2))
```

The admittances of the MMC can be displayed using the function
`plot_data` and are presented in Figs.
[10](#fig_mmc_ac_pll){reference-type="ref" reference="fig_mmc_ac_pll"},
[11](#fig_mmc_dc_side_pll){reference-type="ref"
reference="fig_mmc_dc_side_pll"} and
[12](#fig_mmc_acdc_pll){reference-type="ref"
reference="fig_mmc_acdc_pll"}.

![AC side admittance of the MMC with
PLL.](pictures/mmc/ac_side_pll.pdf){#fig_mmc_ac_pll width="65%"}

![DC side impedance of the MMC with PLL
incorporated.](pictures/mmc/dc_side_pll.pdf){#fig_mmc_dc_side_pll
width="65%"}

![Admittance interconnection between AC and DC side of the MMC with PLL
applied.](pictures/mmc/acdc_side_pll.pdf){#fig_mmc_acdc_pll width="65%"}

### Testing of the point-to-point hybrid power system

The operation of the complete system is tested on the hybrid power
system example depicted in Fig. [13](#fig:ptp_mmc){reference-type="ref"
reference="fig:ptp_mmc"} with the parameter values given in the table
[1](#table:mmc_example){reference-type="ref"
reference="table:mmc_example"}. Gen 1 and gen 2 are modelled as infinite
buses and it is assumed that the power flows from gen 1 to gen 2. It
should be noted that the external grids are modeled as infinite buses
(if their internal impedance is $Z_{int} = 0$), and otherwise as
Thevenin equivalent (that is, if $Z_{inf} \neq 0$).

![Point-to-point HVDC-based power
system.](pictures/mmc/point_to_point-marked.png){#fig:ptp_mmc
width="\\linewidth"}

By specifying and solving the power flow problem for the hybrid power
system [13](#fig:ptp_mmc){reference-type="ref" reference="fig:ptp_mmc"},
the operating point of the converter is obtained. The power flow problem
is formed under the assumption that gen 1 provides the active power
$P = 100 \text{ MW}$ and reactive power $Q = 0$, and that the AC voltage
is balanced, only contains a fundamental frequency component, with the
amplitude $V_m = 320 \text{ kV}$. Gen 2 takes the same amount of active
and reactive power and provides the AC voltage with the same amplitude.

The solution of the power flow specifies that MMC 1 preserves DC voltage
at $V_{DC} = 320 \text{ kV}$, while the AC voltage has amplitude
$V_m = 336.09 \text{ kV}$ with the phase shift
$\theta = -0.054 \text{ rad}$. It takes $P = -100.26 \text{ MW}$ active
power and $Q = 0$ reactive power, while it gives to the DC power network
$P_{dc} = 100.23 \text{ MW}$. Similarly, MMC 2 preserves AC power as
$P = 100 \text{ MW}$ and $Q = -219.76 \text{ MVA}$. Voltages on the DC
and AC side are defined as: $V_{DC} = 319.95 \text{ kV}$,
$V_m = 312.32 \text{ kV}$ and $\theta = 0.0485 \text{ rad}$,
respectively.

::: {#table:mmc_example}
  **OHL1 and OHL2 for 320 kV**                    
  ----------------------------------------------- ---------------------------------------------------------------------
  length                                          200 km
  organization                                    flat
  conductor relative horizontal distance          $\Delta x_{bc} = 10 \text{ m}$
  relative vertical position                      $y_{bc} = 30 \text{ m}$
  conductor sag                                   $d_{sag} = 10 \text{ m}$
  conductors and subconductors                    $n_c = 3$, $n_{sb} = 1$
  conductor radius                                $r_c = 0.015 \text{ m}$
  conductor resistivity                           $R_{dc} = 0.063 \, \Omega\text{m}$
  groundwire horizontal distance                  $\Delta x_g = 6.5 \text{ m}$
  groundwire vertical distance                    $\Delta y_g = 7.5 \text{ m}$
  groundwire sag                                  $d_{gsag} = 10 \text{ m}$
  groundwire radius                               $r_g = 0.0062 \text{ m}$
  groundwire resistivity                          $R_{gdc} = 0.92 \text{ m}$
  conductance                                     $g_c = 1 \cdot 10^{-11} \, \text{S}/\text{m}$
  **Cable 320 kV**                                
  length                                          100 km
  number of cables                                $n = 2$
  positions                                       $x = -0.5 \text{ m}, 0.5 \text{ m}$, $y = 1 \text{ m}$
  core (Cu) radius                                $r^{c1}_i = 0, \ r^{c1}_o = 24.25 \text{ mm}$
  insulator 1 (XLPE) radius                       $r^{i1}_i = r^{c1}_o, r^{i1}_o = 41.75 \text{ mm}$
  sheath (Pb) radius                              $r^{c2}_i = r^{i1}_o, r^{c2}_o = 46.25 \text{ mm}$
  insulator 2 (XLPE) radius                       $r^{i2}_i = r^{c2}_o, r^{i2}_o = 49.75 \text{ mm}$
  armor (steel) radius                            $r^{c3}_i = r^{i2}_o, r^{c3}_o = 60.55 \text{ mm}$
  insulator 3 (XLPE)                              $r^{i3}_i = r^{c3}_o, r^{i3}_o = 65.75 \text{ mm}$
  **MMC 1 and MMC 2**                             
  rated power of the system                       $S_r = 100 \text{ MVA}$
  AC-grid voltage                                 $V_{m} = 320 \text{ kV}$
  DC-bus voltage                                  $V_{DC} = 320 \text{ kV}$
  line frequency                                  $f_0 = 50 \text{ Hz}$
  number of submodules per arm                    $N = 401$
  submodule capacitance                           $C_{arm} = 10 \text{ mF}$
  arm inductance                                  $L_{arm} = 50 \text{ mH}$
  arm resistance                                  $R_{arm} = 1.07 \, \Omega$
  phase reactor inductance                        $L_r = 60 \text{ mH}$
  phase reactor resistance                        $R_r = 0.535 \, \Omega$
  circulating current control                     $\zeta = 0.7$, $\omega_n = 1000 \, \frac{\text{rad}}{\text{s}}$
  output current control                          $\zeta = 0.7$, $\omega_n = 300 \, \frac{\text{rad}}{\text{s}}$
  DC voltage control (only MMC1)                  $K_{p,dc} = 0.01$, $K_{i,dc} = 2$
  active and reactive power control (only MMC2)   $K_{p, PQ} = 2.002 \cdot 10^{-7}$, $K_{i,PQ} = 1.001 \cdot 10^{-4}$
  zero current control                            $\zeta = 0.7$, $\omega_n = 300 \, \frac{\text{rad}}{\text{s}}$
  zero energy control                             $K_{p,ec} = 120$, $K_{i,ec} = 400$

  : System example: point-to-point HVDC-based power system
:::

As it is marked in Fig. [13](#fig:ptp_mmc){reference-type="ref"
reference="fig:ptp_mmc"}, we will check the stability of the system by
"cutting" it at the DC side of the converter MMC 1 and finding the
feedback transfer function.By running the `check_stability` routine, the
values of the $Z_{MMC1}$, $Z_{eq}$ and $Y_{MMC1} Z_{eq}$ over the
frequency range are obtained, see Fig.
[14](#fig:bode_ptp_mmc){reference-type="ref"
reference="fig:bode_ptp_mmc"}. As it can be seen from the Fig.
[14](#fig:bode_ptp_mmc){reference-type="ref"
reference="fig:bode_ptp_mmc"}, the feedback transfer function reaches 0
dB after 100 $\frac{\text{rad}}{\text{s}}$, but the phase margin is
positive, proving that the system is stable.

This simulation example runs for 61.88 s in the Julia programming
language under a Windows 10 operating system, on the CPU Intel Core
i7-8565U 1.99 GHz and with 16 GB RAM memory.

![Bode plot of the $Z_{MMC1}$, $Z_{eq}$ and
$Y_{MMC1} Z_{eq}$.](pictures/mmc/bode_ptp_mmc.pdf){#fig:bode_ptp_mmc
width="\\linewidth"}

## Quick start - Implementation of a new component

1.  Make a new folder in `/Network/Components/new_folder` with the name
    of the new component;

2.  Make a file in the folder with the name `new_component_name.jl`;

3.  Implement the component as the struct type. The struct can contain
    any number of required fields with default initial values. The
    following lines provide universal struct format, whose contructor is
    generated automatically by the compiler. The constructor is then
    called by `New_component_name()` and it gives as output one instance
    of the component with the initialized all parameters initialized to
    the default values.

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    export new\_component\_name
        
        @with\_kw mutable struct New\_component\_name
            field :: Type\_field = initial\_value
        end
    ```

4.  Make a new function called `new_component_name`:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    """
        Comments and description of the function
        """
        function new\_component(;args... )
        
            # Insert here the code that builds the ABCD representation 
            # or the equivalent component parameters 
            # used later for ABCD representation of the component
                
            elem = New\_component\_name()
            element = Element(element\_value = elem, input\_pins = nb\_input\_pins,
                output\_pins = nb\_output\_pins)
            return element
        end
    ```

5.  Make the functions for the determination of ABCD and Y parameters:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    function eval\_abcd(elem :: New\_component\_name, s :: Complex)
        
            # Number of instructions for evaluating ABCD parameters
        
            return ABCD
        end
    ```

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    function eval\_y(elem :: New\_component\_name, s :: Complex)
        
            # Number of instructions for evaluating Y parameters
        
            return Y
        end
    ```

    Naturally, one function may call the other, depending on the way the
    ABCD and Y matrix are built.\
    Also, the pin configuration (1.1, 1.2, \... 2.1, 2.2, \...) that
    will be used to interconnect the element within a network will
    automatically follow the order in which the variables (voltages or
    currents) appear in the input and output vectors of the ABCD
    representation. This is straightforward for single-phase (SISO)
    components. For three-phase (MIMO) components that might be
    unbalanced, the recommended order is that of the phases (a, b, c).
    There might be cases with different numbers of input and output
    pins, such as three-winding transformers, in which case the output
    variables could be considered in the order (a1, b1, c1, a2, b2, c2).

6.  If the component is intended to be added to the power flow, then the
    two functions for setting power flow parameters are added:

    ``` {mathescape="" commandchars="\\\\\\{\\}"}
    function make\_power\_flow\_ac!(elem :: New\_component\_name, 
            dict :: Dict\{String, Any\}, global\_dict :: Dict\{String, Any\})
                
            # Insert here the code to fill-in the powerflow dictionary.
            # Example:
            key = length(dict["new\_powerflow\_component"])
            dict["new\_powerflow\_component"][string(key)]["parameter"] = value
            # where new\_powerflow\_component, parameter and value must be adapted.
        end
            
        function make\_power\_flow\_dc!(elem :: New\_component\_name, 
                dict :: Dict\{String, Any\}, global\_dict :: Dict\{String, Any\})
                
            # Insert here the code to fill-in the powerflow dictionary.
        end
    ```

7.  Add functions for plotting and saving component related data:

    `function save_data(comp :: New_component_type, file_name :: String, omegas)`,

    `function plot_data(comp :: New_component_type, omegas)`;

8.  In the file `element_types`, add a new line
    `include("new_folder/new_component_name.jl")` to link it with the
    power system structures (see for example definition of the component
    `Transformer`).

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
