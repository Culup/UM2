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
  water.setDensity(0.726); // g/cm^3
  water.setH2O();
  Float const relativeDensity = 0.726 / 0.743;
  ASSERT(water.numNuclides() == 2);
  int iH = -1, iO = -1;
  for (int i = 0; i < water.numNuclides(); ++i) {
    if (water.zaid(i) == 1001)
      iH = i;
    if (water.zaid(i) == 8016)
      iO = i;
  }
  ASSERT(iH != -1);
  ASSERT(iO != -1);
  ASSERT_NEAR(water.numDensity(iH), 1.11915E-01 * relativeDensity, 1e-6);
  ASSERT_NEAR(water.numDensity(iO), 8.88085E-01 * relativeDensity, 1e-6);
}

TEST_CASE(setZirc4)
{
  um2::Material zirc4;
  zirc4.setDensity(6.56);
  zirc4.setZirc4();
  Float const relativeDensity = 6.56 / 6.56;
  ASSERT(zirc4.numNuclides() == 29);

  int iZr90 = -1, iSn120 = -1, iFe56 = -1, iCr52 = -1, iHf178 = -1;
  for (int i = 0; i < zirc4.numNuclides(); ++i) {
    int const z = zirc4.zaid(i);
    if (z == 40090)
      iZr90 = i;
    if (z == 50120)
      iSn120 = i;
    if (z == 26056)
      iFe56 = i;
    if (z == 24052)
      iCr52 = i;
    if (z == 72178)
      iHf178 = i;
    ASSERT(zirc4.numDensity(i) > 0);
  }

  ASSERT(iZr90 != -1);
  ASSERT(iSn120 != -1);
  ASSERT(iFe56 != -1);
  ASSERT(iCr52 != -1);
  ASSERT(iHf178 != -1);

  ASSERT_NEAR(zirc4.numDensity(iZr90), 4.98086E-01 * relativeDensity, 1e-6);
  ASSERT_NEAR(zirc4.numDensity(iSn120), 4.77153E-03 * relativeDensity, 1e-6);
  ASSERT_NEAR(zirc4.numDensity(iFe56), 1.92992E-03 * relativeDensity, 1e-6);
  ASSERT_NEAR(zirc4.numDensity(iCr52), 8.36988E-04 * relativeDensity, 1e-6);
  ASSERT_NEAR(zirc4.numDensity(iHf178), 2.71973E-05 * relativeDensity, 1e-6);
}

TEST_CASE(setSS304)
{
  um2::Material ss304;
  ss304.setDensity(8.0);
  ss304.setSS304();
  Float const relativeDensity = 8.0 / 8.0;
  ASSERT(ss304.numNuclides() == 17);

  int iFe56 = -1, iCr52 = -1, iNi58 = -1, iCnat = -1, iMn55 = -1, iSinat = -1, iP31 = -1;
  for (int i = 0; i < ss304.numNuclides(); ++i) {
    int const z = ss304.zaid(i);
    if (z == 26056)
      iFe56 = i;
    if (z == 24052)
      iCr52 = i;
    if (z == 28058)
      iNi58 = i;
    if (z == 6000)
      iCnat = i;
    if (z == 25055)
      iMn55 = i;
    if (z == 14000)
      iSinat = i;
    if (z == 15031)
      iP31 = i;
    ASSERT(ss304.numDensity(i) > 0);
  }

  ASSERT(iFe56 != -1);
  ASSERT(iCr52 != -1);
  ASSERT(iNi58 != -1);
  ASSERT(iCnat != -1);
  ASSERT(iMn55 != -1);
  ASSERT(iSinat != -1);
  ASSERT(iP31 != -1);

  ASSERT_NEAR(ss304.numDensity(iFe56), 6.28376E-01 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iCr52), 1.59029E-01 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iNi58), 6.38386E-02 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iCnat), 7.99365E-04 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iMn55), 2.00001E-02 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iSinat), 1.00002E-02 * relativeDensity, 1e-6);
  ASSERT_NEAR(ss304.numDensity(iP31), 4.50008E-04 * relativeDensity, 1e-6);
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
