#ifndef DIALS_SCRATCH_JMP_FITTING_HELPERS_H
#define DIALS_SCRATCH_JMP_FITTING_HELPERS_H

#include <map>
#include <vector>

class PixelList {
public:
  PixelList() {}

protected:
  std::map<vec3<int>, std::vector<std::size_t>> lookup_;
  std::map<std::size_t, std::vector<vec3<int>>> rlookup_;
};

#endif // DIALS_SCRATCH_JMP_FITTING_HELPERS_H
