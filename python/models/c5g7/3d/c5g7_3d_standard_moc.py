# THIS MODEL HAS NOT BEEN TESTED YET. MAY NOT BE COMPLETE OR CORRECT.

import gc
import sys

from um2 import (
    initialize,
    finalize,
    get_c5g7_materials,
    MPACTModel,
    MESH_QUAD,
    MESH_QUADRATIC_QUAD,
)

from um2.tools.utils import string_to_lattice


def build_c5g7_3d():
    model = MPACTModel()

    materials = get_c5g7_materials()
    (
        uo2,
        mox43,
        mox70,
        mox87,
        fiss_chamber,
        guide_tube,
        moderator,
    ) = materials

    for mat in materials:
        model.add_material(mat)

    pin_pitch = 1.26
    radius = 0.54

    rings = [3, 2]
    radii = [radius, 0.62]
    num_azi = 8

    cyl_pin_mesh_type = MESH_QUADRATIC_QUAD
    cyl_pin_id = model.add_cylindrical_pin_mesh(
        pin_pitch,
        radii,
        rings,
        num_azi,
        2,
    )

    xy_extents = (pin_pitch, pin_pitch)

    rect_pin_mesh_type = MESH_QUAD
    rect_pin_id = model.add_rectangular_pin_mesh(xy_extents, 5, 5)

    # 6 cylindrical coarse cells:
    # first 24 faces are fuel/material i, remaining 24 are moderator material 6.
    for i in range(6):
        mat_ids = [i] * 24 + [6] * 24
        model.add_coarse_cell(
            xy_extents,
            cyl_pin_mesh_type,
            cyl_pin_id,
            mat_ids,
        )

    # 1 rectangular moderator coarse cell:
    mat_ids = [6] * 25
    model.add_coarse_cell(
        xy_extents,
        rect_pin_mesh_type,
        rect_pin_id,
        mat_ids,
    )

    # RTMs
    for i in range(7):
        model.add_rtm([[i]])

    uo2_lattice = string_to_lattice("""
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 5 0 0 5 0 0 5 0 0 0 0 0
        0 0 0 5 0 0 0 0 0 0 0 0 0 5 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 5 0 0 5 0 0 5 0 0 5 0 0 5 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 5 0 0 5 0 0 4 0 0 5 0 0 5 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 5 0 0 5 0 0 5 0 0 5 0 0 5 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 5 0 0 0 0 0 0 0 0 0 5 0 0 0
        0 0 0 0 0 5 0 0 5 0 0 5 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    """)

    mox_lattice = string_to_lattice("""
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
        1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1
        1 2 2 2 2 5 2 2 5 2 2 5 2 2 2 2 1
        1 2 2 5 2 3 3 3 3 3 3 3 2 5 2 2 1
        1 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2 1
        1 2 5 3 3 5 3 3 5 3 3 5 3 3 5 2 1
        1 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1
        1 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1
        1 2 5 3 3 5 3 3 4 3 3 5 3 3 5 2 1
        1 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1
        1 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1
        1 2 5 3 3 5 3 3 5 3 3 5 3 3 5 2 1
        1 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2 1
        1 2 2 5 2 3 3 3 3 3 3 3 2 5 2 2 1
        1 2 2 2 2 5 2 2 5 2 2 5 2 2 2 2 1
        1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    """)

    h2o_lattice = [[6 for _ in range(17)] for _ in range(17)]

    model.add_lattice(uo2_lattice)
    model.add_lattice(mox_lattice)
    model.add_lattice(h2o_lattice)

    model_height = 214.2
    num_slices = 10
    num_fuel_slices = 9 * num_slices // 10
    z_slices = [i * model_height / num_slices for i in range(num_slices + 1)]

    lattice_ids = [2] * num_slices
    lattice_ids[:num_fuel_slices] = [0] * num_fuel_slices
    model.add_assembly(lattice_ids, z_slices)

    lattice_ids = [2] * num_slices
    lattice_ids[:num_fuel_slices] = [1] * num_fuel_slices
    model.add_assembly(lattice_ids, z_slices)

    lattice_ids = [2] * num_slices
    model.add_assembly(lattice_ids, z_slices)

    core = string_to_lattice("""
        0 1 2
        1 0 2
        2 2 2
    """)

    model.add_core(core)

    model.write("c5g7.xdmf")


if __name__ == "__main__":
    initialize()
    try:
        build_c5g7_3d()
        gc.collect()
    finally:
        finalize()