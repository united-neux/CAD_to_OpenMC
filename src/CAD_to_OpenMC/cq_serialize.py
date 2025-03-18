# SPDX-FileCopyrightText: 2025 E. B. Knudsen <erik@united-neux.eu>
# SPDX-License-Identifier: MIT

import copyreg
from io import BytesIO

import cadquery as cq
import OCP


def _inflate_shape(data: bytes):
    with BytesIO(data) as bio:
        return cq.Shape.importBrep(bio)


def _reduce_shape(shape: cq.Shape):
    with BytesIO() as stream:
        shape.exportBrep(stream)
        return _inflate_shape, (stream.getvalue(),)


def _inflate_transform(*values: float):
    trsf = OCP.gp.gp_Trsf()
    trsf.SetValues(*values)
    return trsf


def _reduce_transform(transform: OCP.gp.gp_Trsf):
    return _inflate_transform, tuple(
        transform.Value(i, j) for i in range(1, 4) for j in range(1, 5)
    )


def register():
    """
    Registers pickle support functions for common CadQuery and OCCT objects.
    """

    for cls in (
        cq.Edge,
        cq.Compound,
        cq.Shell,
        cq.Face,
        cq.Solid,
        cq.Vertex,
        cq.Wire,
    ):
        copyreg.pickle(cls, _reduce_shape)

    copyreg.pickle(cq.Vector, lambda vec: (cq.Vector, vec.toTuple()))
    copyreg.pickle(OCP.gp.gp_Trsf, _reduce_transform)
    copyreg.pickle(
        cq.Location, lambda loc: (cq.Location, (loc.wrapped.Transformation(),))
    )
