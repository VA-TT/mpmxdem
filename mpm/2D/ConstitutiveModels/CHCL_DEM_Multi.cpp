#include "CHCL_DEM_Multi.hpp"

#include "Core/MPMbox.hpp"
#include "Core/MaterialPoint.hpp"

#include <cstdint>
#include <stdexcept>

namespace {

// SplitMix64 finalizer: gives a reproducible pseudo-random RVE index for a
// material-point number without depending on construction/read order.
std::uint64_t mixMaterialPointNumber(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

}  // namespace

std::string CHCL_DEM_Multi::getRegistrationName() { return std::string("CHCL_DEM_Multi"); }

void CHCL_DEM_Multi::read(std::istream& is) {
  size_t count = 0;
  is >> count;
  if (count == 0) { throw std::runtime_error("CHCL_DEM_Multi requires at least one DEM RVE file"); }

  fileNames.resize(count);
  for (std::string& name : fileNames) {
    if (!(is >> name)) { throw std::runtime_error("CHCL_DEM_Multi: missing DEM RVE file name"); }
  }
}

void CHCL_DEM_Multi::write(std::ostream& os) {
  os << fileNames.size();
  for (const std::string& name : fileNames) { os << ' ' << name; }
  os << '\n';
}

std::string CHCL_DEM_Multi::getFileNameFor(const MaterialPoint& MP) const {
  if (fileNames.empty()) { throw std::runtime_error("CHCL_DEM_Multi has no DEM RVE file"); }
  const size_t index = static_cast<size_t>(mixMaterialPointNumber(MP.nb) % fileNames.size());
  return fileNames[index];
}

void CHCL_DEM_Multi::init(MaterialPoint& MP) {
  MP.isDoubleScale = true;
  MP.PBC = new PBC3Dbox;

  if (box->computationMode == true) {
    const std::string selectedFile = getFileNameFor(MP);
    MP.PBC->loadConf(selectedFile.c_str());

    // PBC3D and MPMbox use opposite stress sign conventions.
    MP.stress.xx = -MP.PBC->Sig.xx;
    MP.stress.xy = -MP.PBC->Sig.xy;
    MP.stress.yx = -MP.PBC->Sig.yx;
    MP.stress.yy = -MP.PBC->Sig.yy;
    MP.outOfPlaneStress = -MP.PBC->Sig.zz;
  }
}
