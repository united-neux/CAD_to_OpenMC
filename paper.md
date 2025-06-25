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
We present `CAD_to_OpenMC`: a practical tool for allowing geometries designed in
CAD systems to be used in Monte Carlo particle transport simulation systems.
Written in Python, and easily installable through pip and/or conda channels,
`CAD_to_openMC` does so by importing a geometry in the form of a step-file, a
common export format among many CAD systems, identifying objects within that
file, and running a triangularization algorithm on all surfaces contained within
the geometry.
The geometry is written to an HDF5 file in a format defined by the
DAGMC library, which is usable by the main particle-transport packages.


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
a library with couplings to OpenMC, MCNP, Fluka, and Geant4.
DAGMC is proven to be fast and reliable, yet one hurdle remains: DAGMC
requires a discretized (meshed) approximation of the CAD geometry. In general, it is not a trivial requirement.

Most CAD software packages are, though often proprietary packages, able to write to the STEP-format, which is backed by an
international standard [@stepStd_2024]. For this reason, the STEP format was chosen as the input format.

The `CAD_to_OpenMC` tool eases the process of generating a meshed description
of a CAD-generated geometry (in the form of a step file) ready for inclusion in
transport codes through DAGMC.

# Method
`CAD_to_OpenMC` uses OCCT [@occt3d] to interact with CAD geometries, i.e., step-files. Once the geometry has been imported, one of
several meshing backends may be called to create a discretized version of the
geometry. In the end, the generated geometry is exported into the DAGMC-expected format.

OCCT uses a hierarchical boundary representation model (BREP).
In a simplified picture, objects consist of a set of bounded surfaces (faces), formed by a set of boundaries (edges),
which themselves consist of line segments connecting points or simply mathematical descriptions of curves.
More specifically, the BREP-description used in OCCT also includes notions of shells (the boundary of a volume), loops (a circular set of edges), etc.
The BREP concept is fundamentally different from CSG, where objects are formed by
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

In terms of case b, each surface must only be discretized once. `CAD_to_OpenMC` handles this by
filling a hash table for surfaces upon
processing, enabling reuse without reprocessing.

In the last case (c), imprinting has to be performed. This is the process
where the boundary curves of one surface are projected onto another, splitting it into two or more sub surfaces.

<!-- meshing backends -->
### Meshing backends
`CAD_to_OpenMC` is made to be flexible in terms of the
algorithms, backends, used to generate triangularized surfaces.
The present release supports a range of backends: `{'gmsh','stl','stl2','db'}`.
In most cases `stl2` performs well, but sometimes fails to produce a watertight mesh.
In such cases, we recommend using the `db` backend if available. Further, the `gmsh` backend consistently
produces a practical model, but is memory hungry.

### Material tags
<!-- The way material tags are extracted here -->
`CAD_to_OpenMC` can use the CAD part names as material tags.
The default is to use the first part of the part name as a material
tag. This may be changed by supplying a dictionary as
the tags-argument. Here, the keys are regular
expressions, and the values are the material tags to use. Unmatched parts are tagged as vacuum (the default) but may be tagged by extracted part names. 
The below example shows how to tag all parts with the name `wall` in them with `concrete`
and all parts ending with bellows with `steel`.
```python
tag_dict = {".*bellows" : "steel", ".*wall.*" : "concrete"}
```


## Implicit complement
`CAD_to_OpenMC` can add an implicit complement to the output file. That is, the material tag
applied to any part of the geometry *not* claimed by any
CAD part. This is done by simply assigning the name of the material as a string (with a suffix of `\_comp\`)
to the implicit complement attribute of the base Assembly object.
That string gets picked up by DAGMC system and is used for the unclaimed volume.

# Results
We have chosen a reactor model as a test system: A tabletop reactor (GODIVA IV) model from the ICSBEP-benchmark project [@icsbep_2022] as case HEU-MET-FAST-086.

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

![Part-by-part comparison between volume calculations using stochastic volume estimation in OpenMC and volumes reported by CAD-software for the detailed (green) and simplified (magenta) benchmark models. The intervals are the computational error margins, almost exclusively stemming from the estimated error in the stochastic volume computation. The additional black circles denote the relative difference between the extracted benchmark CSG-model (as computed by stochastic volume runs in OpenMC) and the CAD-model.\label{fig:voldiff}](figs/both_rel_voldiff.png)

Figure \ref{fig:voldiff} shows
volume differences between discretized and exact CAD models
for the various objects in the model.

Differences in volume often influence the neutronics
of a reactor more than boundary shifts. Hence, this
is a useful measure for performance. The errors found (Fig. \ref{fig:voldiff}) are dominated by
the error in the stochastic volume estimator, not the volume error itself. This
is evidenced by the very small error in burst- and control-rod volume for the
detailed model, which were run with smaller tolerances. Additionally, we have used the benchmark CSG-model
to verify a few volumes in the detailed model. A CSG model generally yields shorter runtimes, which allows tighter
tolerances while remaining practical. The CSG-benchmark deviates slightly from the CAD-model constructed
from drawings, suggesting an underlying discrepancy internally in the benchmark.

The considered cases have different settings for the control- and
burst-rods (see Tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos}).

|case | CR 1 top | CR 2 top |BR top | $k_{\textrm{eff}}$ CAD | $k_{\textrm{eff}}$ CSG | $k_{\textrm{eff}}$ Lit.|
|-----|------|-----|-----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.98026       | 0.98187       | 0.9865        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98101       | 0.98185       | 0.9867        |
| 3   | -0.493   | -3.794   | 0.0   | 0.98124       | 0.98297       | 0.9878        |
| 4   | -0.469   | -0.447   | -2.970| 0.98745       | 0.98359       | 0.9883        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98706       | 0.98844       | 0.9933        |
Table: Control-rod (CR) and burst-rod (BR) positions for the five cases of the Godiva-IV benchmark/simplified model from HEU-MET-FAST-086
[@icsbep_2022; @goda2021]\label{tab_bm_giv_rod_pos}. Measures in inches
withdrawn from the fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC denotes Monte Carlo, whereas Lit. denotes numbers from the benchmark.

|case | CR 1 top | CR 2 top |BR top | $k_{\textrm{eff}}$ CAD | $k_{\textrm{eff}}$ CSG | $k_{\textrm{eff}}$ Lit.|
|-----|------|-----|-----|----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.97905      | 0.98303       | 0.9880        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98390      | 0.98275       | 0.9880        |
| 3   | -0.493   | -3.794   | 0.0   | 0.97885      | 0.98330       | 0.9887        |
| 4   | -0.469   | -0.447   | -2.970| 0.98352      | 0.98426       | 0.9897        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98390      | 0.98969       | 0.9945        |
Table: Control-rod (CR) and burst-rod (BR) positions for the five cases of the detailed Godiva-IV model from HEU-MET-FAST-08
[@icsbep_2022; @goda2021]\label{tab_det_giv_rod_pos}. Measures in inches
withdrawn from the fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC denotes Monte Carlo, whereas Lit. denotes numbers from the benchmark report.

Two minor adjustments ($\approx 1$ mm) had to be made to the stated measurements; locking bolts and intermediate sub-assembly plate, for the model to fit. Both edits were ~1 mm and did not affect the results. We assume the errors are misprints in the drawings.

It is clear from Tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos} that the agreement between model and
experiment is not perfect, yet the difference between CAD-based and CSG-based
models is small (generally on the order of $5 \times 10^{-3}$). We take this as proof that the geometry generation
works as intended.


# Discussion and Conclusion
We submit that the tool presented is a convenient tool for making CAD geometries available for Monte Carlo particle transport. By utilizing the DAGMC-layer, the resulting geometries are not restricted to OpenMC, but in fact may be used also in MCNP, Fluka, etc. Experience has shown that a particularly useful feature is to extract tags from CAD-defined parts and interpret them as material tags for transport. This enables a consistent material naming scheme throughout the entire modeling procedure.

Finally, as noted, other active projects exist targeting the problem. In the interest of efficiency and resource management, there are active efforts aiming unification, which may bear fruit in coming releases.


# References
