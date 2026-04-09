from .config cimport Int, Float, MatID
from libc.stdint cimport uint8_t
from libc.stdlib cimport malloc, free
from ._um2 cimport (
    UM2_Color_t,
    UM2_MeshType_t,
    um2SizeOfInt,
    um2SizeOfFloat,
    um2Initialize,
    um2Finalize,
    um2NewMPACTModel,
    um2DeleteMPACTModel,
    um2ReadMPACTModel,
    um2MPACTModelAddMaterial,
    um2MPACTModelAddCoarseGrid,
    um2OverlayCoarseGrid,
    um2AddCylindricalPinLattice2D,
    um2SetMeshFieldFromKnudsenNumber,
    um2GenerateMesh,
    um2MPACTModelImportCoarseCellMeshesAndWrite,
    um2MPACTExportModel,
    um2MPACTNumCoarseCells,
    um2MPACTNumRTMs,
    um2MPACTNumLattices,
    um2MPACTNumAssemblies,
    um2MPACTCoreNumCells,
    um2MPACTAssemblyNumCells,
    um2MPACTLatticeNumCells,
    um2MPACTRTMNumCells,
    um2NewMaterial,
    um2DeleteMaterial,
    um2MaterialNumNuclides,
    um2MaterialSetName,
    um2MaterialSetColor,
    um2MaterialSetTemperature,
    um2MaterialSetDensity,
    um2MaterialGetName,
    um2MaterialGetColor,
    um2MaterialGetTemperature,
    um2MaterialGetDensity,
    um2MaterialNumDensities,
    um2MaterialNumDensity,
    um2MaterialZaids,
    um2MaterialZaid,
    um2MaterialAddNuclide,
    um2MaterialPopulateXSec,
    um2MaterialSetUO2,
    um2MaterialSetH2O,
    um2MaterialSetZirc4,
    um2MaterialSetSS304,
    um2NumC5G7Materials,
    um2GetC5G7Materials,
)

MESH_INVALID = UM2_MESH_INVALID
MESH_TRI = UM2_MESH_TRI
MESH_QUAD = UM2_MESH_QUAD
MESH_TRIQUAD = UM2_MESH_TRIQUAD
MESH_QUADRATIC_TRI = UM2_MESH_QUADRATIC_TRI
MESH_QUADRATIC_QUAD = UM2_MESH_QUADRATIC_QUAD
MESH_QUADRATIC_TRIQUAD = UM2_MESH_QUADRATIC_TRIQUAD

def size_of_int():
    cdef int n
    um2SizeOfInt(&n)
    return n


def size_of_float():
    cdef int n
    um2SizeOfFloat(&n)
    return n

def initialize():
    um2Initialize()

def finalize():
    um2Finalize()

cdef class MPACTModel:
    cdef void * model_ptr

    def __cinit__(self):
        self.model_ptr = NULL
        um2NewMPACTModel(&self.model_ptr)

    def close(self):
        if self.model_ptr != NULL:
            um2DeleteMPACTModel(self.model_ptr)
            self.model_ptr = NULL

    def __dealloc__(self):
        if self.model_ptr != NULL:
            um2DeleteMPACTModel(self.model_ptr)
            self.model_ptr = NULL

    def read_from_file(self, path):
        um2ReadMPACTModel(self.model_ptr, path.encode('utf-8'))

    def add_material(self, Material material):
        um2MPACTModelAddMaterial(self.model_ptr, material.material_ptr)

    def add_coarse_grid(self, Float width, Float height, Int nx, Int ny):
        um2MPACTModelAddCoarseGrid(self.model_ptr, width, height, nx, ny)

    def import_coarse_cell_meshes_and_write(self, mesh_path, model_path, bint write_knudsen_data=True):
        um2MPACTModelImportCoarseCellMeshesAndWrite(self.model_ptr,
                                                    mesh_path.encode("utf-8"),
                                                    model_path.encode("utf-8"),
                                                    1 if write_knudsen_data else 0)

    def export_model(self, mesh_path, model_path, bint write_knudsen_data=True):
        um2MPACTExportModel(self.model_ptr,
                            mesh_path.encode("utf-8"),
                            model_path.encode("utf-8"),
                            1 if write_knudsen_data else 0)

    def num_coarse_cells(self):
        cdef Int n
        um2MPACTNumCoarseCells(self.model_ptr, &n)
        return n

    def num_rtms(self):
        cdef Int n
        um2MPACTNumRTMs(self.model_ptr, &n)
        return n

    def num_lattices(self):
        cdef Int n
        um2MPACTNumLattices(self.model_ptr, &n)
        return n

    def num_assemblies(self):
        cdef Int n
        um2MPACTNumAssemblies(self.model_ptr, &n)
        return n

    def core_num_cells(self):
        cdef Int nx
        cdef Int ny
        um2MPACTCoreNumCells(self.model_ptr, &nx, &ny)
        return (nx, ny)

    def assembly_num_cells(self, Int asy_id):
        cdef Int nx
        um2MPACTAssemblyNumCells(self.model_ptr, asy_id, &nx)
        return nx

    def lattice_num_cells(self, Int lat_id):
        cdef Int nx
        cdef Int ny
        um2MPACTLatticeNumCells(self.model_ptr, lat_id, &nx, &ny)
        return (nx, ny)

    def rtm_num_cells(self, Int rtm_id):
        cdef Int nx
        cdef Int ny
        um2MPACTRTMNumCells(self.model_ptr, rtm_id, &nx, &ny)
        return (nx, ny)

def overlay_coarse_grid(MPACTModel model, Material material):
    um2OverlayCoarseGrid(model.model_ptr, material.material_ptr)

def set_mesh_field_from_knudsen_number(MPACTModel model,
                                       int dim,
                                       double target_kn,
                                       double fuel_mfp_threshold=-1.0,
                                       double fuel_mfp_scale=-1.0,
                                       double abs_mfp_threshold=-1.0,
                                       double abs_mfp_scale=-1.0):
    um2SetMeshFieldFromKnudsenNumber(
        dim,
        model.model_ptr,
        target_kn,
        fuel_mfp_threshold,
        fuel_mfp_scale,
        abs_mfp_threshold,
        abs_mfp_scale,
    )

def generate_mesh(UM2_MeshType_t mesh_type):
    um2GenerateMesh(mesh_type)

def add_cylindrical_pin_lattice_2d(lattice_ids,
                                   xy_extents,
                                   pin_radii,
                                   pin_materials,
                                   Float origin_x,
                                   Float origin_y):
    cdef Int ny, nx, num_pin_types
    cdef Int i, j, k, idx
    cdef Int total_radii = 0
    cdef Int total_materials = 0

    cdef Int * lattice_buf = NULL
    cdef Float * extents_buf = NULL
    cdef Float * radii_flat_buf = NULL
    cdef Int * radii_offsets_buf = NULL
    cdef void ** mats_flat_buf = NULL
    cdef Int * material_offsets_buf = NULL

    if not lattice_ids:
        raise ValueError("lattice_ids cannot be empty")

    ny = len(lattice_ids)
    nx = len(lattice_ids[0])
    for row in lattice_ids:
        if len(row) != nx:
            raise ValueError("all lattice rows must have the same length")

    num_pin_types = len(xy_extents)
    if len(pin_radii) != num_pin_types:
        raise ValueError("pin_radii must have one list per pin type")
    if len(pin_materials) != num_pin_types:
        raise ValueError("pin_materials must have one list per pin type")

    for i in range(num_pin_types):
        if len(xy_extents[i]) != 2:
            raise ValueError("each xy_extent must have length 2")
        if len(pin_radii[i]) != len(pin_materials[i]):
            raise ValueError("each pin type must have matching radii/material counts")
        total_radii += len(pin_radii[i])
        total_materials += len(pin_materials[i])

    lattice_buf = <Int *>malloc(nx * ny * sizeof(Int))
    extents_buf = <Float *>malloc(2 * num_pin_types * sizeof(Float))
    radii_flat_buf = <Float *>malloc(total_radii * sizeof(Float))
    radii_offsets_buf = <Int *>malloc((num_pin_types + 1) * sizeof(Int))
    mats_flat_buf = <void **>malloc(total_materials * sizeof(void *))
    material_offsets_buf = <Int *>malloc((num_pin_types + 1) * sizeof(Int))

    if (lattice_buf == NULL or extents_buf == NULL or
        radii_flat_buf == NULL or radii_offsets_buf == NULL or
        mats_flat_buf == NULL or material_offsets_buf == NULL):
        free(lattice_buf)
        free(extents_buf)
        free(radii_flat_buf)
        free(radii_offsets_buf)
        free(mats_flat_buf)
        free(material_offsets_buf)
        raise MemoryError()

    try:
        idx = 0
        for j in range(ny):
            for i in range(nx):
                lattice_buf[idx] = <Int>lattice_ids[j][i]
                idx += 1

        for i in range(num_pin_types):
            extents_buf[2 * i] = <Float>xy_extents[i][0]
            extents_buf[2 * i + 1] = <Float>xy_extents[i][1]

        idx = 0
        radii_offsets_buf[0] = 0
        for i in range(num_pin_types):
            for j in range(len(pin_radii[i])):
                radii_flat_buf[idx] = <Float>pin_radii[i][j]
                idx += 1
            radii_offsets_buf[i + 1] = idx

        idx = 0
        material_offsets_buf[0] = 0
        for i in range(num_pin_types):
            for j in range(len(pin_materials[i])):
                if not isinstance(pin_materials[i][j], Material):
                    raise TypeError("pin_materials entries must be Material objects")
                mats_flat_buf[idx] = (<Material>pin_materials[i][j]).material_ptr
                idx += 1
            material_offsets_buf[i + 1] = idx

        um2AddCylindricalPinLattice2D(
            lattice_buf,
            nx,
            ny,
            extents_buf,
            radii_flat_buf,
            radii_offsets_buf,
            <const void * const *>mats_flat_buf,
            material_offsets_buf,
            num_pin_types,
            origin_x,
            origin_y)
    finally:
        free(lattice_buf)
        free(extents_buf)
        free(radii_flat_buf)
        free(radii_offsets_buf)
        free(mats_flat_buf)
        free(material_offsets_buf)

cdef uint8_t _checked_byte(int value):
    if value < 0 or value > 255:
        raise ValueError("color channel must be between 0 and 255")
    return <uint8_t>value

cdef class Color:
    cdef UM2_Color_t _color

    def __cinit__(self, int r=0, int g=0, int b=0, int a=255):
        self._color.r = _checked_byte(r)
        self._color.g = _checked_byte(g)
        self._color.b = _checked_byte(b)
        self._color.a = _checked_byte(a)

    @property
    def r(self):
        return self._color.r

    @r.setter
    def r(self, int value):
        self._color.r = _checked_byte(value)

    @property
    def g(self):
        return self._color.g

    @g.setter
    def g(self, int value):
        self._color.g = _checked_byte(value)

    @property
    def b(self):
        return self._color.b

    @b.setter
    def b(self, int value):
        self._color.b = _checked_byte(value)

    @property
    def a(self):
        return self._color.a

    @a.setter
    def a(self, int value):
        self._color.a = _checked_byte(value)

    def __repr__(self):
        return f"Color({self._color.r}, {self._color.g}, {self._color.b}, {self._color.a})"

cdef UM2_Color_t _as_um2_color_t(Color color):
    return color._color

cdef Color _from_um2_color_t(UM2_Color_t color):
    cdef Color py_color = Color.__new__(Color)
    py_color._color = color
    return py_color

cdef class Material:
    cdef void * material_ptr
    cdef bint owner

    @staticmethod
    cdef Material from_ptr(void * ptr, bint owner=True):
        cdef Material m = Material.__new__(Material)
        m.material_ptr = ptr
        m.owner = owner
        return m

    def __cinit__(self):
        self.material_ptr = NULL
        self.owner = False

    def __init__(self):
        if self.material_ptr == NULL:
            um2NewMaterial(&self.material_ptr)
            self.owner = True

    def __dealloc__(self):
        if self.owner and self.material_ptr != NULL:
            um2DeleteMaterial(self.material_ptr)
            self.material_ptr = NULL
    
    def num_nuclides(self):
        cdef Int n = um2MaterialNumNuclides(self.material_ptr)
        return n

    def set_name(self, name):
        um2MaterialSetName(self.material_ptr, name.encode('utf-8'))

    def set_color(self, Color color):
        um2MaterialSetColor(self.material_ptr, _as_um2_color_t(color))

    def set_temperature(self, Float temp):
        um2MaterialSetTemperature(self.material_ptr, temp)

    def set_density(self, Float density):
        um2MaterialSetDensity(self.material_ptr, density)

    def get_name(self):
        return um2MaterialGetName(self.material_ptr).decode('utf-8')

    def get_color(self):
        cdef UM2_Color_t color = um2MaterialGetColor(self.material_ptr)
        return _from_um2_color_t(color)

    def get_temperature(self):
        return um2MaterialGetTemperature(self.material_ptr)

    def get_density(self):
        return um2MaterialGetDensity(self.material_ptr)

    def num_densities(self):
        cdef Int n = um2MaterialNumNuclides(self.material_ptr)
        cdef Int i
        return [um2MaterialNumDensity(self.material_ptr, i) for i in range(n)]

    def num_density(self, Int i):
        return um2MaterialNumDensity(self.material_ptr, i)

    def zaids(self):
        cdef Int n = um2MaterialNumNuclides(self.material_ptr)
        cdef Int i
        return [um2MaterialZaid(self.material_ptr, i) for i in range(n)]

    def zaid(self, Int i):
        return um2MaterialZaid(self.material_ptr, i)

    def add_nuclide(self, Int zaid, Float num_density):
        um2MaterialAddNuclide(self.material_ptr, zaid, num_density)
    
    def populate_xsec(self):
        um2MaterialPopulateXSec(self.material_ptr)

    def set_uo2(self, Float wt_u235, Float wt_gad):
        um2MaterialSetUO2(self.material_ptr, wt_u235, wt_gad)

    def set_h2o(self):
        um2MaterialSetH2O(self.material_ptr)

    def set_zirc4(self):
        um2MaterialSetZirc4(self.material_ptr)
    
    def set_ss304(self):
        um2MaterialSetSS304(self.material_ptr)

def get_c5g7_materials():
    cdef Int n = um2NumC5G7Materials()
    cdef void ** mats = <void **>malloc(n * sizeof(void *))
    cdef Int i
    cdef list result = []

    if mats == NULL:
        raise MemoryError()

    try:
        um2GetC5G7Materials(mats)
        for i in range(n):
            result.append(Material.from_ptr(mats[i], owner=False))
    finally:
        free(mats)

    return result