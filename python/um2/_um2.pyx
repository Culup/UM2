from ._um2 cimport Int
from ._um2 cimport (
    um2SizeOfInt,
    um2SizeOfFloat,
    um2Initialize,
    um2Finalize,
    um2NewMPACTModel,
    um2DeleteMPACTModel,
    um2ReadMPACTModel,
    um2MPACTNumCoarseCells,
    um2MPACTNumRTMs,
    um2MPACTNumLattices,
    um2MPACTNumAssemblies,
    um2MPACTCoreNumCells,
    um2MPACTAssemblyNumCells,
    um2MPACTLatticeNumCells,
    um2MPACTRTMNumCells
)

cdef class MPACTModel:
    cdef void * model_ptr

    def __cinit__(self):
        um2NewMPACTModel(&self.model_ptr)

    def __dealloc__(self):
        um2DeleteMPACTModel(self.model_ptr)

    def read_from_file(self, path):
        um2ReadMPACTModel(path.encode('utf-8'), &self.model_ptr)

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

def size_of_int(int size):
        cdef int n
        um2SizeOfInt(&n)
        return n


def size_of_float(int size):
        cdef int n
        um2SizeOfFloat(&n)
        return n

def initialize():
    um2Initialize()

def finalize():
    um2Finalize()