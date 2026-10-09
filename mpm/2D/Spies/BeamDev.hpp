#ifndef BEAM_DEV_HPP
#define BEAM_DEV_HPP

#include "Spy.hpp"

#include <vector>

#include "vec2.hpp"

struct BeamDev : public Spy {
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  std::string filename;
  std::ofstream file;
  double KineticEnergy;
  size_t id_MP_end;
  vec2r pos0;
};

#endif /* end of include guard: BEAM_DEV_HPP */
