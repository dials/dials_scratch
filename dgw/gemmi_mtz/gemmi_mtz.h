
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
  mtz.cell = gemmi::UnitCell(10, 20, 30, 90, 90, 90);    // dummy cell
  mtz.spacegroup = gemmi::find_spacegroup_by_number(19); // dummy sg
  mtz.add_base();
  const char *pxd[3] = {"DIALS", "XTAL", "FROMDIALS"};
  double wavelength = 1.0; // dummy wavelength
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
  af::const_ref<int> batch_col = reflections["batch"];
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
    std::array<int, 3> hkl = {h[i][0], h[i][1], h[i][2]};
    int isym = hkl_mover.move_to_asu(hkl);

    for (size_t j = 0; j != 3; ++j)
      mtz.data[k++] = (float)hkl[j]; // HKL
    mtz.data[k++] = (float)isym;     // M/ISYM
    int frame = (float)batch_col[i];
    frames.insert(frame);
    mtz.data[k++] = (float)frame;      // BATCH
    mtz.data[k++] = (float)isum[i];    // I
    mtz.data[k++] = (float)sigisum[i]; // SIGI
    mtz.data[k++] = (float)iprf[i];    // IPR
    mtz.data[k++] = (float)sigiprf[i]; // SIGIPR
  }

  // Prepare a batch header
  gemmi::Mtz::Batch batch;
  batch.set_dataset_id(1);

  // We don't set lbcell (refinement flags for unit cell),
  // because it's probably not used by any program anyway.
  // batch.ints[4] to [9] = left unset

  // We also skip jumpax, which is defined as:
  // reciprocal axis closest to principle goniostat axis E1
  // batch.ints[11] = left unset

  // We assume one crystal, one goniostat axis, one detector.
  batch.ints[12] = 1; // ncryst
  batch.ints[14] = 2; // ldtype 3D
  batch.ints[15] = 1; // jsaxs - goniostat scan axis number
  batch.ints[17] = 1; // ngonax - number of goniostat axes
  batch.ints[19] = 1; // ndet

  batch.set_cell(mtz.cell);                     // batch.floats[0] to [5]
  gemmi::Vec3 s0(-1, 0, 0);                     // dummy s0 vector
  gemmi::Mat33 U = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // dummy U mat
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      batch.floats[6 + 3 * i + j] = (float)U[j][i];

  batch.floats[21] =
      float(0.1); // dummy reflection width (full width) (degrees)
  // In the so-called "Cambridge" frame (as used by Mosflm),
  // the principal rotation axis is along z
  // and the incident beam S0 is along x.
  // Therefore, we set the rotation axis phi (scanax) to:
  batch.floats[38 + 2] = 1.f; // scanax = [0, 0, 1]

  double oscillation_range = 0.1; // dummy oscillation range
  batch.floats[47] = float(oscillation_range);
  // E1,E2,E3 vectors are the goniostat axes in Cambridge laboratory frame.
  // E1 is set to rotation axis. E2 and E3 are not set for ngonax==1.
  batch.floats[59 + 2] = 1.f; // e1 = scanax

  // Idealised source vector is -x ([-1 0 0]), antiparallel to beam.
  batch.floats[80 + 0] = -1.f; // source[0]

  // s0 source vector (including tilts), CMtz::MTZBAT::so in libccp4
  batch.floats[83] = (float)s0.x;
  batch.floats[84] = (float)s0.y;
  batch.floats[85] = (float)s0.z;

  batch.set_wavelength((float)wavelength); // batch.floats[86]
  // or batch.set_wavelength(iset.wavelength);
  // Detector geometry.
  double detector_distance = 100; // dummy dx[0]
  batch.floats[111] = (float)detector_distance;
  batch.floats[113] = 1.f;         // detlm[0][0][0]
  batch.floats[114] = (float)1024; // dummy NX
  batch.floats[115] = 1.f;
  batch.floats[116] = (float)1024; // dummy NY
  batch.axes.push_back("PHI");     // gonlab[0]

  double starting_angle = 0;
  int starting_frame = 1;
  for (int frame : frames) {
    batch.number = frame;
    double phistt =
        starting_angle + oscillation_range * (frame - starting_frame);
    batch.floats[36] = float(phistt);
    batch.floats[37] = float(phistt + oscillation_range); // phiend
    mtz.batches.push_back(batch);
  }

  mtz.sort(5);
  const char *output_path = "foo.mtz";
  mtz.write_to_file(output_path);
  return;
}

} // namespace gemmi_mtz
} // namespace dials

#endif // DIALS_GEMMI_MTZ_H
