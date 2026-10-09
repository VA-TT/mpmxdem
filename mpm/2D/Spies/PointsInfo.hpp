#ifndef POINTS_INFO_HPP
#define POINTS_INFO_HPP

#include "Spy.hpp"

#include <vector>

#include "vec2.hpp"

struct PointsInfo : public Spy {
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  std::string filename;
  std::ofstream file;
  double KineticEnergy;
  size_t id_MP_end;
  std::vector<vec2r> pos0;
};

#endif /* end of include guard: POINTS_INFO_HPP */
