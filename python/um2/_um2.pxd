from libc.stdint cimport int32_t

ctypedef int32_t Int

cdef extern from "um2c.h":
    void um2SizeOfInt(int * size)
    void um2SizeOfFloat(int * size)

    void um2Initialize()
    void um2Finalize()

    void um2NewMPACTModel(void ** model)
    void um2DeleteMPACTModel(void * model)
    void um2ReadMPACTModel(const char * path, void ** model)

    void um2MPACTNumCoarseCells(void * model, Int * n)
    void um2MPACTNumRTMs(void * model, Int * n)
    void um2MPACTNumLattices(void * model, Int * n)
    void um2MPACTNumAssemblies(void * model, Int * n)

    void um2MPACTCoreNumCells(void * model, Int * nx, Int * ny)
    void um2MPACTAssemblyNumCells(void * model, Int asy_id, Int * nx)
    void um2MPACTLatticeNumCells(void * model, Int lat_id, Int * nx, Int * ny)
    void um2MPACTRTMNumCells(void * model, Int rtm_id, Int * nx, Int * ny)

    
