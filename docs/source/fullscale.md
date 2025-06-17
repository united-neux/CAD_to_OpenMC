---
bibliography: bibliography.bib
---
```{bibliography} bibliography.bib
```

# Full scale runs
Below we show two more examples where this CAD_to_OpenMC has been used in
conjunction with OpenMC to arrive at full-scale simulations of molten salt
reactors, comparing the results with experiment. This provides an examples of
not only the workings of OpenMC, but also a view into the deatils of modelling
(or lack thereof). It has often been remarked that triangularized geometry is
not good for criticality. We find however that we get just as close to
experiment with this kind of geometry as with a CSG-geometry, suggesting that
the geometrical differences stemming from discretiztion errors are smaller than differences between (reported) model and reality. 

## Molten Salt Reactors
:::{figure} ../images/msre_are.png
:name: fig:msreAre

Rendering of triangularized models of the ARE- and MSRE-reactors as generated
using CAD_to_OpenMC. The MSRE model (left) includes both reactor core, liquid
fuel contained therein, graphite moderator stringers, as well as thermal
shielding and reactor enclosure. The reactor pit _is_ included in the model,
but we have excluded it from the image to make the core more visible. The ARE
model (right) includes the set of three safety/shim rods and the regulating rod
in the center.
:::

{numref}`fig:msreAre` shows pictures of the meshed MSRE and ARE geometries,
including reactor enclosure, control rods, etc.
The Aircraft Reactor Experiment (ARE) {footcite:p}`cottrell_are_operation_1955`
and Molten Salt Reactor Experiment (MSRE){footcite:p}`robertson_msre_1965` were
two molten salt reactor experiments carried out at Oak Ridge National
Laboratory; ARE in November '59 and MSRE was run between '65 and '70.
In the confines of this project, this pair serves as examples of complex
reactor geometries that can be handled by CAD_to_OpenMC. Very detailed CAD
models and example scripts used to compute these numbers are available for the
set, which we have used as inputs {footcite:p}`msreData,areData`.
{numref}`tabbmark` tabulates the {math}`k_{eff}`-values computed using a
materials composition set to mimic the reported values as closely as possible.
It is clear that, similar to the case for the GODIVA IV device, the modelled
values are not in complete agreement with the reported ones, yet this may be
explained by possible discrepancies between the drawings in accessible reports
and the actual experiment. Any engineering realities not written down in these
reports are likely and unfortunately lost.
This is particularly true for ARE, where the discrepancy is largest, and
details are scarce.
It should be noticed that the calculated keff for the MSRE case has better
agreement than any other MSRE criticality benchmarks known to the authors.

:::{table} Criticality numbers, {math}`k_{eff}`, for two molten salt reactors as computed by the combination of CAD-model, CAD_to_OpenMC, and OpenMC, using ENDF v8.0 cross-sections. Also tabled are corresponding (albeit trivial) experimental values found in literature,$k_{eff,lit}$
:name: tabbmark
 
| model | {math}`k_{eff}` | {math}`k_{eff,lit}` |
|-------|-----------|-------|
| ARE   | 1.06582   | 1.0   |
|MSRE   | 1.00549  | 1.0|

:::

:::{footbibliography}
:::
<!-- meshing timings?-->

