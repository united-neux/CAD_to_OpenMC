---
title: "CAD_to_OpenMC: from CAD design to particle transport"
tags:
    - Python
    - Monte Carlo
    - Model preparation
    - CAD
    - Mesh

authors:
    - name: Erik Bergbäck Knudsen
      orcid: 0000-0003-1827-9188
      corresponding: true
      affiliation: 1
      equal contribution: true
    - name: Lorenzo Chierici
      affiliation: 2

affiliations:
    - name: United Neux, Denmark
      index: 1
    - name: Copenhagen Atomics A/S, Denmark
      index: 2

date: XX September 2024
bibliography: bibliography.bib
---

# Summary
<!-- Here's a citation @knudsen2013mcxtrace , or in a different format [@knudsen2013mcxtrace].-->
We present CAD_to_OpenMC - a practical tool for allowing geometries designed in
CAD-systems to be used in Monte Carlo particle transport simulation systems.
Written in Python, and easily installable through pip and/or conda channels,
CAD_to_openMC does so by importing a geometry in the form of a step-file, a
common export format among many CAD systems, identifying objects within that
file, and running a triangularization algorithm on all surfaces contained within
the geometry.
The geometry is written to an HDF5 file in a format defined by the
DAGMC-library, which is usable by the main particle-transport packages.


# Statement of need
Most Monte Carlo neutronics, and general particle transport software packages (e.g., OpenMC, GEANT4, MCNP, Fluka) [@openmc_2013;
@geant4_2003; @mcnp_2022; @fluka_2005; @fluka_2014] are ubiquitous in nuclear science.
Almost all of these use Constructive Solid Geometry (CSG) to describe the geometry in which computations are to be performed. CSG is computationally efficient but also tricky to work with, in that all object boundaries must be constructed from analytical descriptions of geometrical
primaries (cylinders, planes, 2nd-order surfaces, etc.).
Most engineering design, however, is performed using CAD tools, which are
optimized for ease of design.

These facts form a compelling argument for bridging the gap between these
worlds. So far, a number of solutions have been proposed which can mostly be put
into two categories.

1. tools that convert CAD models into CSG [@du_inversecsg_2018; @catalan_geouned_2024], and
2. tools that allow ray-tracing on CAD geometries, or close derivatives thereof.

In the second category, one of the most widespread open source solutions is DAGMC [@wilson_acceleration_2010],
a library with couplings to OpenMC, MCNP, fluka, and Geant4.
DAGMC is proven to be fast and reliable, yet one hurdle remains: DAGMC
requires a discretized (meshed) approximation of the CAD geometry. In general, it is not a trivial requirement.

CAD software packages are, for the most
part, proprietary packages with little focus on being interchangeable. One
upshot of this is a notable lack of overlapping standard formatting between
packages, the exception being the STEP-format, which is backed by an
international standard [@stepStd_2024]. Thus, a tool should support this format if
it is not to be tied to any particular CAD engine.

The CAD_to_OpenMC tool eases the process of generating a meshed description
of a CAD-generated geometry (in the form of a step file) ready for inclusion in
transport codes through DAGMC.

While other active projects exist that target similar problems (e.g., cad_to_dagmc [@cadtodagmc], and stellarmesh [@stellarmesh]),
CAD_to_OpenMC is designed with generality in mind. This is defined as:
First, it is aimed at
working for all step-geometries, with no assumptions on geometry. Second, it must be relatively easy
to add code for a new backend should one wish to do so. Hence the code structure with one frontend and several backend classes.
Last, it must be easy to extract, generate, and manipulate material tags from the underlying step model. As models become large and complex, maintaining a separate list of materials, in sequence, is seen as cumbersome and error-prone.

# Method
CAD_to_OpenMC uses OCCT [@occt3d] to interact with CAD geometries, i.e. step-files. Once the geometry has been imported, one of
several meshing backends may be called to create a discretized version of the
geometry. In the end, the generated geometry is exported into the DAGMC-expected format.

Similar to most CAD systems, OCCT uses a hierarchical boundary representation model (BREP).
In a simplified picture, objects consist of a set of bounded surfaces (faces), formed by a set of boundaries (edges),
which themselves consist of line segments connecting points or simply mathematical descriptions of curves.
More specifically, the BREP-description used in OCCT also includes notions of shells (the boundary of a volume), loops (a circular set of edges), etc.
The BREP concept fundamentally different from CSG, where objects are formed by
boolean operations on halfspaces. For example, the OCCT BREP-description allows several more types of operations for forming a geometry,
such as extrusions and rotations.

Thus, a need for conversion always exists, regardless
of whether one discretizes the geometry or not.



## Triangularization / Surface Meshing
<!-- Describe the imprinting and why it is very important especially for this case -->
When a curved geometry is to be discretized, the result is an approximation.
This leads to potential problems with overlaps between regions if

a. objects are close to each other, such as for a cylindrical can with thin walls
b. objects share a surface such as the case for a liquid inside a can.
c. when objects have surfaces that touch each other.

The first case puts a constraint on the absolute tolerance of the
discretization, i.e., triangles have to be small enough not to cause crossing
surfaces.

In terms of case b, the discretization process must only process a shared surface once.
his means that objects cannot be independently processed. CAD_to_OpenMC handles this by
generating a hash code for each surface upon
processing. Each time a surface is encountered, the code is checked, and if the
surface has already been encountered, it will be reused.

In the last case (c), imprinting has to be performed. This is the process
where the boundary curves of one surface are projected onto another, splitting it into two or more sub surfaces.

<!-- meshing backends -->
### Meshing backends
CAD_to_OpenMC is made to be flexible in terms of the
algorithms, backends, used to generate triangularized surfaces.
The present release supports a range of backends: `{'gmsh','stl','stl2','db'}`.
In most cases `stl2` performs well. It uses a very basic
algorithm, with no restriction on triangle aspect ratio, but also produces
discretized geometry files of moderate size. In practice, the aspect ratio does
not pose problem for particle transport but other applications are more sensitive.

In some cases, the simple algorithm fails to create a watertight model. In such
cases, we recommend using the `db` backend if available. The `gmsh` backend consistently
produces a practical model, but often has the problem that it is memory hungry. For example,
neither the ARE nor the MSRE (see below) models could be run this way on our
available hardware (64 GB workstation).


### Material tags
<!-- The way material tags are extracted here -->
If needed, CAD_to_OpenMC can use the CAD-generated part names as material tags.
The default behavior is to use the first part of the part name as a material
tag. This may be changed by supplying a dictionary, as
the tags-argument. Here, the keys are regular
expressions, and the values are the material tags to use. Unmatched parts are tagged as vacuum (the default) by may be tagged by extracte part names. 
Thu,s a subset of the parts may be retagged.
The below example shows how to tag all parts with the name `wall` in them with `concrete`
and all parts ending with bellows with `steel`.
```python
tag_dict = {".*bellows" : "steel", ".*wall.*" : "concrete"}
```

<!--
here's a test equation:
$$x= \int_a^{b+23} \frac{1}{1+gx} \mathrm{d}x \label{eq:integral}$$
which is nice

I am now going to refer to the equation: \autoref{eq:integral}.
-->



## Implicit complement
CAD_to_OpenMC can add an implicit complement to the output file. That is, the material tag
applied to any part of the geometry *not* claimed by any
CAD part. This is done by simply assigning the name of the material as a string (with a suffix of `\_comp\`)
to the implicit complement attribute of the base Assembly object.
That string gets picked up by DAGMC system and is used for the unclaimed volume.

# Results
We have chosen three reactor models as test systems. A tabletop reactor and two full-scale molten salt reactors. The former
(GODIVA IV) model is included in the ICSBEP-benchmark project [@icsbep_2022] as case HEU-MET-FAST-086; the latter two were part of the molten salt reactor program at ORNL.

## GODIVA IV
This reactor consists of a cylindrical core (\autoref{fig:GIV_CAD}) held by three sets of clamps set. Additionally, the core has three vertical holes into which control and burst rods may be inserted from below. The rods themselves are similar in composition to the
fuel elements.

The benchmark includes five experimental cases, which differ in terms of control- and burst-rod positions.
Three geometries are described:

1. a detailed model which includes as-built things (curved clamp etc.),
2. simplified model, where all surfaces are along the principal axis, all corners 90 deg. etc., and
3. cylindrically symmetric model.

The benchmark reports only experimental results for the two first, but contains MCNP geometries for
all three in [@icsbep_2022]. The 3rd model is created in order to use some legacy analysis tools which are purely 1D/2D.

Corresponding to \autoref{fig:GIV_CAD}, \autoref{fig:GIV_meshed} shows the discretized version of the two reactor models used for our further analysis.

![CAD-drawings of the Godiva IV reactor, detailed (left) and simplified (right) versions. Note the rectangular clamps and supports in the simplified version. Visible are also the set of control and burst rods.\label{fig:GIV_CAD}](figs/GIV_both.pdf){#GIV_CAD width="50%"}

![Discretized Godiva IV-models, detailed (left) and simplified (right) versions.\label{fig:GIV_meshed}](figs/both_meshed.pdf){#GIV_meshed width=50%}

![Part-by-part comparison between volume calculations using stochastic volume estimation in OpenMC compared with direct volume calculations reported by CAD-software for the detailed model (green) and simplified benchmark model (magenta). The indicated intervals are the computational error margins, almost exclusively stemming from the estimated error in the stochastic volume computation. The additional black circles denote the relative difference between the extracted benchmark CSG-model (as computed by stochastic volume runs in OpenMC) and the CAD-model.\label{fig:voldiff}](figs/both_rel_voldiff.png)

Figure \ref{fig:voldiff} shows
differences in volumes between the discretized models and the exact CAD model
for the various objects in the model. We used the stochastic volume estimator of OpenMC.
[@openmc_2013] for the discretized models.

Generally, differences in volume have a much bigger influence on the neutronics
of a reactor than do small boundary changes (with constant volume). Hence, this
is a useful measure for performance. The errors found (fig. \ref{fig:voldiff}) are dominated by
the error in the stochastic volume estimator, not the volume error itself. This
is evidenced by the very small error in burst- and control-rod volume for the
detailed model, which were run with smaller tolerances. Additionally, we have used the benchmark CSG-model
to verify a few volumes in the detailed model. A CSG model generally yields shorter runtimes, which allows tighter
tolerance while also remaining practical. We find that the CSG-benchmark deviates slightly from the CAD-model constructed
from drawings, suggesting there is an underlying discrepancy internally in the benchmark, which may have to be addressed.

The considered cases have different settings for the control- and
burst-rods each (see tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos}).

|case | CR 1 top | CR 2 top |BR top | $k_{eff}$ CAD | $k_{eff}$ CSG | $k_{eff}$ Lit.|
|-----|------|-----|-----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.98026       | 0.98187       | 0.9865        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98101       | 0.98185       | 0.9867        |
| 3   | -0.493   | -3.794   | 0.0   | 0.98124       | 0.98297       | 0.9878        |
| 4   | -0.469   | -0.447   | -2.970| 0.98745       | 0.98359       | 0.9883        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98706       | 0.98844       | 0.9933        |
Table: Control rod (CR) and burst rod (BR) positions for the 5 cases of the Godiva-IV benchmark/simplified model from HEU-MET-FAST-086
[@icsbep_2022; @goda2021]\label{tab_bm_giv_rod_pos}. Measures in inches
withdrawn from the fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC denotes Monte Carlo, whereas Lit. denotes numbers from the benchmark.

|case | CR 1 top | CR 2 top |BR top | $k_{eff}$ CAD | $k_{eff}$ CSG | $k_{eff}$ Lit.|
|-----|------|-----|-----|----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.97905      | 0.98303       | 0.9880        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98390      | 0.98275       | 0.9880        |
| 3   | -0.493   | -3.794   | 0.0   | 0.97885      | 0.98330       | 0.9887        |
| 4   | -0.469   | -0.447   | -2.970| 0.98352      | 0.98426       | 0.9897        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98390      | 0.98969       | 0.9945        |
Table: Control rod (CR) and burst rod (BR) positions for the 5 cases of the detailed Godiva-IV model from HEU-MET-FAST-08
[@icsbep_2022; @goda2021]\label{tab_det_giv_rod_pos}. Measures in inches
withdrawn from the fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC denotes Monte Carlo, whereas Lit. denotes numbers from the benchmark report.

To get a model in which the rods can be moved, the core geometry, burst, and control-rods were discretized individually
The full reactor model is then assembled as an OpenMC geometry.

Two minor adjustments ($\approx 1$ mm) had to be made to the stated measurements; locking bolts and intermediate sub-assembly plate, for the model to fit. Both edits were ~= 1 mm and did not affect the results. We assume the errors are misprints in the drawings.

It is clear from tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos} that the agreement between model and
experiment is not perfect, yet the difference between CAD-based and CSG-based
models is small (generally on the order of 5e-3). We take this as proof that the geometry generation
works as intended.


## Molten Salt Reactors
![Rendering of triangularized models of the ARE- and MSRE-reactors as generated using CAD_to_OpenMC. The MSRE model (left) includes both reactor core, liquid fuel contained therein, graphite moderator stringers, as well as thermal shielding and reactor enclosure. The reactor pit _is_ included in the model, but we have excluded it from the image to make the core more visible. The ARE model (right) includes the set of three safety/shim rods and the regulating rod in the center.\label{fig:msreAre}](figs/msre_are.png)


\autoref{fig:msreAre} (left) shows pictures of the meshed MSRE and ARE geometries, including reactor enclosure control rods, etc.
The Aircraft Reactor Experiment (ARE) and Molten Salt Reactor Experiment (MSRE) were two molten salt reactor experiments carried out at Oak Ridge National Laboratory; ARE in November '59 and MSRE was run between '65 and 70.
In the confines of this article, this pair serves as examples of complex reactor geometries that can be handled by CAD_to_OpenMC. Very detailed CAD models and example scripts used to compute these numbers are available for the set, which we have used as inputs [@msreData; @areData].
\autoref{tab_bmark} tabulates the $k_{eff}$-values computed using a materials composition set to mimic the reported values as closely as possible.
It is clear that, similar to the case for the GODIVA IV device, the modelled
values are not in complete agreement with the reported ones, yet this may be
explained by possible discrepancies between the drawings in accessible reports and the actual experiment. Any engineering realities not written down in these reports are likely lost.
This is particularly true for ARE, where the discrepancy is largest, and details are scarce.
It should be noticed that the calculated keff for the MSRE case has better agreement than any other MSRE criticality benchmarks known to the authors.

| model\label{tab_bmark} | $k_{eff}$ | $k_{eff,lit}$ |
|-------|-----------|-------|
| ARE   | 1.06582   | 1.0   |
|MSRE   | 1.00549  | 1.0|

Table: Criticality numbers, $k_{eff}$, for a two molten salt reactors
as computed by the combination of CAD-model, CAD_to_OpenMC, and OpenMC, using
ENDF v8.0 cross-sections. Also tabled are corresponding (albeit trivial) experimental values found in
literature,$k_{eff,lit}$ [@cottrell_are_operation_1955; @robertson_msre_1965].

<!-- meshing timings?-->


# Discussion and Conclusion
We submit that the tool presented is a convenient tool for making CAD geometries available for Monte Carlo particle transport. By utilizing the DAGMC-layer, the resulting geometries are not restricted to OpenMC, but in fact may be used also in MCNP, fluka, etc. Experience has shown that a particularly useful feature is to extract tags from CAD-defined parts and interpret them as material tags for transport. This enables a consistent material naming scheme throughout the entire modeling procedure.

Finally, as noted, other active projects exist targeting the problem. In the interest of efficiency and resource management, there are active efforts aiming unification, which may bear fruit in coming releases.


# References
