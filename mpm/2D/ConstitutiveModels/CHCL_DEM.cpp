#include "CHCL_DEM.hpp"

#include "Core/MPMbox.hpp"
#include "Core/MaterialPoint.hpp"

#include "PBC3D.hpp"
#include "MohrCoulomb.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <stdexcept>

double CHCL_DEM::shortestThickness(const mat9r& h) {
  // h stores lattice vectors as columns. V/|b x c| is the distance
  // between the two faces spanned by b,c; lengths of columns are not
  // thicknesses when the cell is sheared.
  const vec3r a(h.xx, h.yx, h.zx), b(h.xy, h.yy, h.zy), c(h.xz, h.yz, h.zz);
  const double volume = h.det();
  const double area = std::max({norm(cross(b,c)), norm(cross(c,a)), norm(cross(a,b))});
  if (!(volume > 0.0) || !(area > 0.0) || !std::isfinite(volume) || !std::isfinite(area)) {
    throw std::runtime_error("RVE replacement: invalid periodic cell geometry");
  }
  return volume / area;
}

bool CHCL_DEM::replaceIfNeeded(MPMbox& MPM, size_t p) {
  if (MPM.CHCL.replaceThicknessFactor <= 0.0) return false;
  MaterialPoint& mp = MPM.MP[p];
  if (!mp.PBC) throw std::runtime_error("RVE replacement: missing DEM cell");
  const auto it = MPM.models.find(MPM.CHCL.replacementModel);
  if (it == MPM.models.end() || dynamic_cast<MohrCoulomb*>(it->second) == nullptr) {
    throw std::runtime_error("CHCL.replaceRVE must refer to a declared MohrCoulomb model");
  }
  double rmax = mp.PBC->Rmax;
  if (!(rmax > 0.0)) {
    for (const auto& grain : mp.PBC->Particles) rmax = std::max(rmax, grain.radius);
  }
  if (!(rmax > 0.0)) throw std::runtime_error("RVE replacement: invalid maximum grain radius");
  const double dmax = 2.0 * rmax;
  mat4r previousInverse = MPM.prevF(p);
  previousInverse.inverse();
  const mat4r increment = mp.F * previousInverse;
  mat9r increment3D;
  increment3D.xx=increment.xx; increment3D.xy=increment.xy;
  increment3D.yx=increment.yx; increment3D.yy=increment.yy; increment3D.zz=1.0;
  const double current = shortestThickness(mp.PBC->Cell.h);
  const double predicted = shortestThickness(increment3D * mp.PBC->Cell.h);
  if (std::min(current,predicted) > MPM.CHCL.replaceThicknessFactor*dmax) return false;

  // Do not change density/mass/volume, total strain, stress or its zz/xz/yz
  // components. MC starts from this exact state, applying this step once.
  // Its in-plane return mapping still assumes zz is intermediate; no 3D
  // yield surface or retrospective DEM plastic-strain reconstruction is added.
  ConstitutiveModel* replacement = it->second;
  replacement->init(mp);
  mp.constitutiveModel = replacement;
  mp.isTracked = false; // no DEM configuration exists after replacement
  delete mp.PBC;
  mp.PBC = nullptr;
  // This physical model transition must remain visible even with a quiet logger.
  // Serialize the whole line because different MPs can convert in parallel.
  {
    static std::mutex notificationMutex;
    const std::lock_guard<std::mutex> lock(notificationMutex);
    fmt::print(stderr, "[RVE replacement] DEM -> Mohr-Coulomb: MP={} t={} s model={} "
                       "hmin/dmax={} predicted={} threshold={}\n",
               mp.nb, MPM.t, replacement->key, current/dmax, predicted/dmax,
               MPM.CHCL.replaceThicknessFactor);
    std::fflush(stderr);
  }
  replacement->updateStrainAndStress(MPM,p);
  return true;
}

std::string CHCL_DEM::getRegistrationName() { return std::string("CHCL_DEM"); }

CHCL_DEM::CHCL_DEM() {}

void CHCL_DEM::read(std::istream& is) { is >> fileName; }

void CHCL_DEM::write(std::ostream& os) { os << fileName << '\n'; }

// The elastic properties cannot be get that way, so, as a convention, -1 is returned
double CHCL_DEM::getYoung() { return -1.0; }
double CHCL_DEM::getPoisson() { return -1.0; }

void CHCL_DEM::init(MaterialPoint& MP) {
  MP.isDoubleScale = true;

  MP.PBC = new PBC3Dbox;
  if (box->computationMode == true) {
    MP.PBC->loadConf(fileName.c_str());

    // transfert the current stress
    // !!! (Sign convention is opposed between PBC3D and MPMbox) !!!
    MP.stress.xx = -MP.PBC->Sig.xx;
    MP.stress.xy = -MP.PBC->Sig.xy;
    MP.stress.yx = -MP.PBC->Sig.yx;
    MP.stress.yy = -MP.PBC->Sig.yy;
    MP.outOfPlaneStress = -MP.PBC->Sig.zz;
    // The DEM cell is 3D and its stress has no reason to be plane-strain.
    // Without these two the deviatoric invariant of the cell is truncated, and
    // with it the effective friction seen by the post-processing.
    MP.outOfPlaneShearXZ = -0.5 * (MP.PBC->Sig.xz + MP.PBC->Sig.zx);
    MP.outOfPlaneShearYZ = -0.5 * (MP.PBC->Sig.yz + MP.PBC->Sig.zy);
  }
}

void CHCL_DEM::updateStrainAndStress(MPMbox& MPM, size_t p) {
  if (replaceIfNeeded(MPM,p)) return;
  size_t* I = &(MPM.Elem[MPM.MP[p].e].I[0]);

  // Get the total strain increment from node velocities
  vec2r vn;
  mat4r dstrain;
  const vec2r *gNp = MPM.gradN(p);
  for (size_t r = 0; r < element::nbNodes; r++) {
    dstrain.xx += (MPM.nodes[I[r]].vel.x * gNp[r].x) * MPM.dt;
    dstrain.xy +=
        0.5 * (MPM.nodes[I[r]].vel.x * gNp[r].y + MPM.nodes[I[r]].vel.y * gNp[r].x) * MPM.dt;
    dstrain.yy += (MPM.nodes[I[r]].vel.y * gNp[r].y) * MPM.dt;
  }
  dstrain.yx = dstrain.xy;  // symmetic tensor

  MPM.MP[p].strain += dstrain;
  MPM.MP[p].deltaStrain = dstrain;

  mat4r prev_F_inv = MPM.prevF(p);
  prev_F_inv.inverse();
  mat4r Finc2D = MPM.MP[p].F * prev_F_inv;

  // remember here that MPM is 2D and DEM is 3D
  mat9r Finc3D;
  Finc3D.xx = Finc2D.xx;
  Finc3D.xy = Finc2D.xy;
  Finc3D.yx = Finc2D.yx;
  Finc3D.yy = Finc2D.yy;
  Finc3D.zz = 1.0;  // assuming plane strain

  mat9r SigAvg;

  // clang-format off
  MPM.MP[p].PBC->transform(
    Finc3D, 
    MPM.dt, 
	  MPM.CHCL.minDEMstep, 
    MPM.CHCL.rateAverage, 
    MPM.CHCL.criticalDEMTimeStepFactor, 
		SigAvg
  );
  // clang-format on

  // We force the symmetry of the stress tensor
  // because it's ok for the DEM side, but probably not for MPM
  SigAvg.symmetrize();

  // Stress
  // !!! (Sign convention is opposed between PBC3D and MPMbox) !!!
  MPM.MP[p].stress.xx = -SigAvg.xx;
  MPM.MP[p].stress.xy = -SigAvg.xy;
  MPM.MP[p].stress.yx = -SigAvg.yx;
  MPM.MP[p].stress.yy = -SigAvg.yy;
  MPM.MP[p].outOfPlaneStress = -SigAvg.zz;
  // SigAvg has just been symmetrized, so xz == zx and yz == zy here.
  MPM.MP[p].outOfPlaneShearXZ = -SigAvg.xz;
  MPM.MP[p].outOfPlaneShearYZ = -SigAvg.yz;
}
