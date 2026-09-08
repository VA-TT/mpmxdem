#ifndef CHCL_DEM_MULTI_HPP
#define CHCL_DEM_MULTI_HPP

#include "CHCL_DEM.hpp"

#include <string>
#include <vector>

// Computationally homogenised law using several initial DEM RVEs.
// Each material point still owns an independent PBC3Dbox simulation. The
// initial RVE is selected pseudo-randomly from fileNames using MP.nb.
struct CHCL_DEM_Multi : public CHCL_DEM {
  std::vector<std::string> fileNames;

  std::string getRegistrationName() override;
  void read(std::istream& is) override;
  void write(std::ostream& os) override;
  void init(MaterialPoint& MP) override;
  std::string getFileNameFor(const MaterialPoint& MP) const;
};

#endif
