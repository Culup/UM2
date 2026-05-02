import gc
import sys

from um2 import (
    initialize,
    finalize,
    get_c5g7_materials,
    MPACTModel,
    add_cylindrical_pin_lattice_2d,
    overlay_coarse_grid,
    set_global_mesh_size,
    generate_mesh,
    MESH_QUADRATIC_TRI,
)

from um2.tools.utils import string_to_lattice

def build_c5g7_2d(num_coarse_cells):
    radius = 0.54
    pin_pitch = 1.26
    assembly_pitch = 21.42

    materials = get_c5g7_materials()
    uo2, mox43, mox70, mox87, fiss_chamber, guide_tube, moderator = materials

    xy_extents = [(pin_pitch, pin_pitch)] * 6
    pin_radii = [[radius]] * 6
    pin_materials = [
        [uo2],
        [mox43],
        [mox70],
        [mox87],
        [fiss_chamber],
        [guide_tube],
    ]

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

    add_cylindrical_pin_lattice_2d(
        uo2_lattice, xy_extents, pin_radii, pin_materials, 0.0, 2.0 * assembly_pitch
    )
    add_cylindrical_pin_lattice_2d(
        uo2_lattice, xy_extents, pin_radii, pin_materials, 1.0 * assembly_pitch, 1.0 * assembly_pitch
    )
    add_cylindrical_pin_lattice_2d(
        mox_lattice, xy_extents, pin_radii, pin_materials, 0.0, 1.0 * assembly_pitch
    )
    add_cylindrical_pin_lattice_2d(
        mox_lattice, xy_extents, pin_radii, pin_materials, 1.0 * assembly_pitch, 2.0 * assembly_pitch
    )
    
    model = MPACTModel()
    for mat in materials:
        model.add_material(mat)
    
    model.add_coarse_grid(3.0 * assembly_pitch, 3.0 * assembly_pitch, num_coarse_cells, num_coarse_cells)
    overlay_coarse_grid(model, moderator)

    set_global_mesh_size(pin_pitch / 4)
    generate_mesh(MESH_QUADRATIC_TRI)

    model.export_model("c5g7_2d.inp", "c5g7_2d.xdmf", True)


if __name__ == "__main__":
    initialize()
    try:
        if len(sys.argv) != 2:
            raise SystemExit("Usage: python c5g7_2d_variable_grid.py num_coarse_cells")
        
        num_coarse_cells = int(sys.argv[1])

        build_c5g7_2d(num_coarse_cells)
        gc.collect()
    finally:
        finalize()