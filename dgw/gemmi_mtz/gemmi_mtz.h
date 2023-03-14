
#ifndef DIALS_GEMMI_MTZ_H
#define DIALS_GEMMI_MTZ_H
#include <cctbx/miller.h>
#include <dials/array_family/reflection_table.h>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/mtz.hpp>
#include <gemmi/unitcell.hpp>
#include <iostream>
#include <set>

namespace dials {
namespace gemmi_mtz {

typedef cctbx::miller::index<> miller_index;

void create_mtz(const char *title, af::reflection_table reflections) {
  gemmi::Mtz mtz;
  mtz.title = title;

  mtz.history.emplace_back("From DIALS etc. ");
  mtz.cell = gemmi::UnitCell(10, 20, 30, 90, 90, 90);
  mtz.spacegroup = gemmi::find_spacegroup_by_number(19);
  mtz.add_base();
  const char *pxd[3] = {"DIALS", "XTAL", "FROMDIALS"};
  double wavelength = 1.0;
  mtz.datasets.push_back({1, pxd[0], pxd[1], pxd[2], mtz.cell, wavelength});
  mtz.add_column("M/ISYM", 'Y', 0, -1, false);
  mtz.add_column("BATCH", 'B', 0, -1, false);
  mtz.add_column("I", 'J', 0, -1, false);
  mtz.add_column("SIGI", 'Q', 0, -1, false);
  mtz.add_column("IPR", 'J', 0, -1, false);
  mtz.add_column("SIGIPR", 'Q', 0, -1, false);
  // mtz.add_column("BG", 'R', 0, -1, false);
  // mtz.add_column("SIGBG", 'R', 0, -1, false);
  // mtz.add_column("FRACTIONCALC", 'R', 0, -1, false);
  // mtz.add_column("XDET", 'R', 0, -1, false);
  // mtz.add_column("YDET", 'R', 0, -1, false);
  // mtz.add_column("ROT", 'R', 0, -1, false);
  // mtz.add_column("LP", 'R', 0, -1, false);
  // mtz.add_column("QE", 'R', 0, -1, false);

  // Get columns out of the table
  af::const_ref<miller_index> h = reflections["miller_index"];
  af::const_ref<int> batch = reflections["batch"];
  af::const_ref<double> isum = reflections["intensity.sum.value"];
  af::const_ref<double> sigisum = reflections["intensity.sum.variance"];
  af::const_ref<double> iprf = reflections["intensity.prf.value"];
  af::const_ref<double> sigiprf = reflections["intensity.prf.variance"];

  // Set MTZ size
  mtz.nreflections = (int)h.size();
  mtz.data.resize(mtz.columns.size() * h.size());

  // Prepare to loop over reflections
  std::set<int> frames;
  gemmi::UnmergedHklMover hkl_mover(mtz.spacegroup);
  size_t k = 0;
  for (std::size_t i = 0; i < h.size(); i++) {
    // std::cout << h[i] << std::endl;
    std::array<int, 3> hkl = {h[i][0], h[i][1], h[i][2]};

    int isym = hkl_mover.move_to_asu(hkl);

    for (size_t j = 0; j != 3; ++j)
      mtz.data[k++] = (float)hkl[j]; // HKL
    mtz.data[k++] = (float)isym;     // M/ISYM
    int frame = (float)batch[i];
    frames.insert(frame);
    mtz.data[k++] = (float)frame;      // BATCH
    mtz.data[k++] = (float)isum[i];    // I
    mtz.data[k++] = (float)sigisum[i]; // SIGI
    mtz.data[k++] = (float)iprf[i];    // IPR
    mtz.data[k++] = (float)sigiprf[i]; // SIGIPR
  }

  mtz.sort(5);
  // std::cout << "writing file" << std::endl;
  const char *output_path = "foo.mtz";
  mtz.write_to_file(output_path);
  return;
}

} // namespace gemmi_mtz
} // namespace dials

#endif // DIALS_GEMMI_MTZ_H
