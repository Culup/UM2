#include <um2/config.hpp>
#include <um2/physics/material.hpp>
#include <um2/stdlib/string.hpp>
#include <um2/stdlib/vector.hpp>

#if UM2_USE_MPACT_XSLIBS
#  include <um2/common/cast_if_not.hpp>
#  include <um2/common/color.hpp>
#  include <um2/common/settings.hpp>
#  include <um2/physics/cross_section_library.hpp>
#endif

#include "../test_macros.hpp"

Float constexpr n_avo_barn = 0.602214076;

TEST_CASE(addNuclide)
{
  um2::Material m;
  m.addNuclide("H1", 1.0);
  m.addNuclide("He4", 1.0);
  m.addNuclide("U235", 1.0);
  m.addNuclide("U-238", 1.0);
  m.addNuclide("Cm-244", 1.0);

  Float constexpr h_wt = 1.00783;
  Float constexpr o_wt = 15.9949;
  Float constexpr h2o_wt = 2 * h_wt + o_wt;
  um2::Material h2o;
  h2o.setDensity(0.75);
  h2o.addNuclideWt("H1", 2 * h_wt / h2o_wt);
  h2o.addNuclideWt("O16", o_wt / h2o_wt);
  Float const h_num_density = 0.75 * 2 * 0.602214076 / h2o_wt;
  Float const o_num_density = h_num_density / 2;
  ASSERT_NEAR(h2o.numDensity(0), h_num_density, 1e-6);
  ASSERT_NEAR(h2o.numDensity(1), o_num_density, 1e-6);

  um2::Material h2o_atom;
  h2o_atom.setDensity(0.75);
  um2::Vector<um2::String> const symbols = {"H1", "O16"};
  um2::Vector<Float> const percents = {2.0 / 3.0, 1.0 / 3.0};
  h2o_atom.addNuclidesAtomPercent(symbols, percents);
  ASSERT_NEAR(h2o_atom.numDensity(0), h_num_density, 1e-6);
  ASSERT_NEAR(h2o_atom.numDensity(1), o_num_density, 1e-6);
}

#if UM2_USE_MPACT_XSLIBS

TEST_CASE(getXS)
{
  um2::XSLibrary const lib8(um2::settings::xs::library_path + "/" + um2::mpact::XSLIB_8G);
  um2::Material fuel;
  fuel.setName("Fuel");
  fuel.setDensity(castIfNot<Float>(10.42));
  fuel.setTemperature(castIfNot<Float>(565.0));
  fuel.setColor(um2::forestgreen);
  fuel.addNuclide("U235", castIfNot<Float>(1.0));
  fuel.addNuclide("O16", castIfNot<Float>(1.0));

  fuel.populateXSec(lib8);
  ASSERT_NEAR(fuel.xsec().t(0), 9, 1)
  ASSERT_NEAR(fuel.xsec().t(1), 14, 1)
}

#endif

TEST_CASE(setH2O)
{
  um2::Material water;
  water.setDensity(0.743); // g/cm^3
  water.setH2O();
  ASSERT(water.numNuclides() == 2);
  Float constexpr m_h1 = 1.00784;
  Float constexpr m_o16 = 15.9949146;
  Float const m_h2o = 2 * m_h1 + m_o16;
  Float const n_h2o = 0.743 * n_avo_barn / m_h2o;
  Float const n_h = 2 * n_h2o;
  Float const n_o = n_h2o;
  int iH = -1, iO = -1;
  for (int i = 0; i < water.numNuclides(); ++i) {
    if (water.zaid(i) == 1001) iH = i;
    if (water.zaid(i) == 8016) iO = i;
  }
  ASSERT(iH != -1);
  ASSERT(iO != -1);
  ASSERT_NEAR(water.numDensity(iH), n_h, 1e-6);
  ASSERT_NEAR(water.numDensity(iO), n_o, 1e-6);
}

TEST_CASE(setZirc4)
{
  um2::Material zirc4;
  zirc4.setDensity(6.55);
  zirc4.setZirc4();

  ASSERT(zirc4.numNuclides() == 23);

  int iZr90=-1, iSn120=-1, iFe56=-1, iCr52=-1;
  for (int i = 0; i < zirc4.numNuclides(); ++i) {
    int const z = zirc4.zaid(i);
    if (z == 40090) iZr90 = i;
    if (z == 50120) iSn120 = i;
    if (z == 26056) iFe56 = i;
    if (z == 24052) iCr52 = i;
    ASSERT(zirc4.numDensity(i) > 0);
  }

  ASSERT(iZr90 != -1);
  ASSERT(iSn120 != -1);
  ASSERT(iFe56 != -1);
  ASSERT(iCr52 != -1);

  Float constexpr m_zr = 91.222;
  Float constexpr wt_zr = 1.0 - (0.015 + 0.002 + 0.001); // 0.982
  Float const n_zr = wt_zr * 6.55 * n_avo_barn / m_zr;
  Float constexpr ab_zr90 = 51.47 / 100.0;

  ASSERT_NEAR(zirc4.numDensity(iZr90), ab_zr90 * n_zr, 1e-6);
}

TEST_CASE(setSS304)
{
  um2::Material ss304;
  ss304.setDensity(8.0);
  ss304.setSS304();

  ASSERT(ss304.numNuclides() == 24);

  int iFe56=-1, iCr52=-1, iNi58=-1, iC12 = -1, iMn55=-1, iSi28=-1, iP31=-1, iS32=-1;
  for (int i = 0; i < ss304.numNuclides(); ++i) {
    int const z = ss304.zaid(i);
    if (z == 26056) iFe56 = i;
    if (z == 24052) iCr52 = i;
    if (z == 28058) iNi58 = i;
    if (z == 6012) iC12 = i;
    if (z == 25055) iMn55 = i;
    if (z == 14028) iSi28 = i;
    if (z == 15031) iP31 = i;
    if (z == 16032) iS32 = i;
    ASSERT(ss304.numDensity(i) > 0);
  }

  ASSERT(iFe56 != -1);
  ASSERT(iCr52 != -1);
  ASSERT(iNi58 != -1);
  ASSERT(iC12 != -1);
  ASSERT(iMn55 != -1);
  ASSERT(iSi28 != -1);
  ASSERT(iP31 != -1);
  ASSERT(iS32 != -1);

  Float constexpr m_fe = 55.845;
  Float constexpr wt_fe = 1.0 - (0.19 + 0.0925 + 0.00080 + 0.02 + 0.010 + 0.00045 + 0.00030); // 0.68595
  Float const n_fe = wt_fe * 8.0 * n_avo_barn / m_fe;
  Float constexpr ab_fe56 = 91.754 / 100.0;

  ASSERT_NEAR(ss304.numDensity(iFe56), ab_fe56 * n_fe, 1e-6);
}

TEST_SUITE(Material)
{
  TEST(addNuclide);
  TEST(setH2O);
  TEST(setZirc4);
  TEST(setSS304);
#if UM2_USE_MPACT_XSLIBS
  TEST(getXS);
#endif
}

auto
main() -> int
{
  RUN_SUITE(Material);
  return 0;
}
