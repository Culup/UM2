#pragma once

#include <um2/config.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  uint8_t r;
  uint8_t g;
  uint8_t b;
  uint8_t a;
} UM2_Color_t;

typedef enum {
  UM2_MESH_INVALID = 0,
  UM2_MESH_TRI = 3,
  UM2_MESH_QUAD = 4,
  UM2_MESH_TRIQUAD = 7,
  UM2_MESH_QUADRATIC_TRI = 6,
  UM2_MESH_QUADRATIC_QUAD = 8,
  UM2_MESH_QUADRATIC_TRIQUAD = 14
} UM2_MeshType_t;

//==============================================================================
// Data sizes 
//==============================================================================

void
um2SizeOfInt(int * size);

void
um2SizeOfFloat(int * size);

//==============================================================================
// Memory management
//==============================================================================

void
um2Malloc(void ** p, Int size);

void
um2Free(void * p);

//==============================================================================
// Initialization and finalization
//==============================================================================

void
um2Initialize();

void
um2Finalize();


//==============================================================================
// MPACT model
//==============================================================================

void
um2NewMPACTModel(void ** model);

void
um2DeleteMPACTModel(void * model);

void
um2ReadMPACTModel(void * model, char const * path);

// Model construction
// ------------------------------------------------------------------------------
void
um2MPACTModelAddMaterial(void * model, void const * material);

void
um2MPACTModelAddCoarseGrid(void * model, Float width, Float height, Int nx, Int ny);

void
um2MPACTModelImportCoarseCellMeshesAndWrite(void * model,
                                            char const * mesh_path,
                                            char const * model_path,
                                            int write_knudsen_data);

void
um2MPACTExportModel(void * model,
                      char const * mesh_path,
                      char const * model_path,
                      int write_knudsen_data);

Int
um2MPACTModelAddCylindricalPinMesh(void * model,
                                   Float pitch,
                                   Float const * radii,
                                   Int const * rings,
                                   Int n_radii,
                                   Int num_azi,
                                   Int mesh_order);

Int
um2MPACTModelAddRectangularPinMesh(void * model,
                                   Float width,
                                   Float height,
                                   Int nx,
                                   Int ny);

void
um2MPACTModelAddCoarseCell(void * model,
                           Float width,
                           Float height,
                           UM2_MeshType_t mesh_type,
                           Int mesh_id,
                           MatID const * mat_ids,
                           Int n_mat_ids);

void
um2MPACTModelAddRTM(void * model,
                    Int const * ids,
                    Int nx,
                    Int ny);

void
um2MPACTModelAddLattice(void * model,
                        Int const * ids,
                        Int nx,
                        Int ny);

void
um2MPACTModelAddAssembly(void * model,
                         Int const * lattice_ids,
                         Float const * z_slices,
                         Int n_lattice_ids);

void
um2MPACTModelAddCore(void * model,
                     Int const * ids,
                     Int nx,
                     Int ny);

void
um2MPACTModelWrite(void * model,
                   char const * path);

// Num
//------------------------------------------------------------------------------
void
um2MPACTNumCoarseCells(void * model, Int * n);

void
um2MPACTNumRTMs(void * model, Int * n);

void
um2MPACTNumLattices(void * model, Int * n);

void
um2MPACTNumAssemblies(void * model, Int * n);

// NumCells
//------------------------------------------------------------------------------
void
um2MPACTCoreNumCells(void * model, Int * nx, Int * ny);

void
um2MPACTAssemblyNumCells(void * model, Int asy_id, Int * nx);

void
um2MPACTLatticeNumCells(void * model, Int lat_id, Int * nx, Int * ny);

void
um2MPACTRTMNumCells(void * model, Int rtm_id, Int * nx, Int * ny);

// GetChild
//------------------------------------------------------------------------------
void
um2MPACTCoreGetChild(void * model, Int ix, Int iy, Int * child);

void
um2MPACTAssemblyGetChild(void * model, Int asy_id, Int ix, Int * child);

void
um2MPACTLatticeGetChild(void * model, Int lat_id, Int ix, Int iy, Int * child);

void
um2MPACTRTMGetChild(void * model, Int rtm_id, Int ix, Int iy, Int * child);

// CoarseCell
//------------------------------------------------------------------------------
void
um2MPACTCoarseCellNumFaces(void * model, Int cc_id, Int * num_faces);

void
um2MPACTCoarseCellWidth(void * model, Int cc_id, Float * width);

void
um2MPACTCoarseCellHeight(void * model, Int cc_id, Float * height);

void
um2MPACTCoarseCellFaceAreas(void * model, Int cc_id, Float * areas);

void
um2MPACTCoarseCellFaceContaining(void * model, Int cc_id, Float x, Float y, Int * face_id);

void
um2MPACTCoarseCellFaceCentroid(void * model, Int cc_id, Int face_id, Float * x, Float * y);

void
um2MPACTCoarseCellMaterialIDs(void * model, Int cc_id, MatID * mat_ids); 

void
um2MPACTIntersectCoarseCell(void * model, Int cc_id, Float origin_x, Float origin_y,
                            Float direction_x, Float direction_y, Float * intersections,
                            Int * n);

// RTM
//------------------------------------------------------------------------------
void
um2MPACTRTMWidth(void * model, Int rtm_id, Float * width);

void
um2MPACTRTMHeight(void * model, Int rtm_id, Float * height);

// Heights
//-----------------------------------------------------------------------------
void
um2MPACTCoarseCellHeights(void * model, Int * n, Int ** cc_ids, Float ** heights);

void
um2MPACTRTMHeights(void * model, Int * n, Int ** rtm_ids, Float ** heights);

void
um2MPACTLatticeHeights(void * model, Int * n, Int ** lat_ids, Float ** heights);

void
um2MPACTAssemblyHeights(void * model, Int asy_id, Float * heights);

void
um2MPACTCoarseCellFaceData(void * model, Int cc_id, Int * mesh_type, Int * num_vertices,
                           Int * num_faces, Float ** vertices, Int ** fv);

//==============================================================================
// GMSH Model Manipulation
//==============================================================================

void
um2OverlayCoarseGrid(void * model, void const * material);

void
um2AddCylindricalPinLattice2D(Int const * lattice_ids, Int nx, Int ny,
                              Float const * xy_extents,               // length = 2 * num_pin_types
                              Float const * pin_radii_flat,
                              Int const * pin_radii_offsets,          // length = num_pin_types + 1
                              void const * const * pin_materials_flat,
                              Int const * pin_material_offsets,       // length = num_pin_types + 1
                              Int num_pin_types,
                              Float origin_x,
                              Float origin_y);

void
um2SetMeshFieldFromKnudsenNumber(Int dim, void const * model,
                                 double kn_target, double mfp_threshold = -1.0,
                                 double mfp_scale = -1.0,
                                 double abs_mfp_threshold = -1.0, 
                                 double abs_mfp_scale = -1.0);

void
um2SetGlobalMeshSize(double const mesh_size);

void
um2GenerateMesh(UM2_MeshType_t mesh_type);

//==============================================================================
// Materials
//==============================================================================
void 
um2NewMaterial(void ** material);

void
um2DeleteMaterial(void * material);

Int
um2MaterialNumNuclides(void const * material);

void
um2MaterialSetName(void * material, char const * name);

void
um2MaterialSetColor(void * material, UM2_Color_t color);

void
um2MaterialSetTemperature(void * material, Float temperature);

void
um2MaterialSetDensity(void * material, Float density);

char const * 
um2MaterialGetName(void const * material);

UM2_Color_t
um2MaterialGetColor(void const * material);

Float
um2MaterialGetTemperature(void const * material);

Float
um2MaterialGetDensity(void const * material);

void
um2MaterialNumDensities(void const * material, Float * num_densities);

Float
um2MaterialNumDensity(void const * material, Int i);

void
um2MaterialZaids(void const * material, Int * zaids);

Int
um2MaterialZaid(void const * material, Int i);

void
um2MaterialAddNuclide(void * material, Int zaid, Float num_density);

void
um2MaterialPopulateXSec(void * material);

void
um2MaterialSetUO2(void * material, Float wt_u235, Float wt_gad);

void
um2MaterialSetH2O(void * material);

void 
um2MaterialSetSS304(void * material);

void
um2MaterialSetZirc4(void * material);

Int
um2NumC5G7Materials(void);

void
um2GetC5G7Materials(void ** materials);



#ifdef __cplusplus
}
#endif
