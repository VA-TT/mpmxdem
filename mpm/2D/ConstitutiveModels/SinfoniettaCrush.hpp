#ifndef SINFONIETTA_CRUSH_HPP
#define SINFONIETTA_CRUSH_HPP

/**
 * @file SinfoniettaCrush.hpp
 * @brief Header file for the SinfoniettaCrush class.
 *
 * This file contains the declaration of the SinfoniettaCrush class, which is
 * a constitutive model for a material point.
 *
 * The SinfoniettaCrush model is an extension of the SinfoniettaClassica model,
 * which is used for the crushing of granular materials.
 *
 * The SinfoniettaCrush model is a simple model that can be used to model
 * the behavior of granular materials under a wide range of conditions.
 *
 * The class relies on various components such as the Young's modulus, Poisson's
 * ratio, and the rigidity of the material point.
 *
 * The file also defines constants and includes necessary headers for the
 * implementation of SinfoniettaCrush.
 */

#include "ConstitutiveModel.hpp"

#include "Rigidity.hpp"
#include "transitFunc.hpp"

struct SinfoniettaCrush : public ConstitutiveModel {
  double Young;
  double Poisson;
  Rigidity C;

  // See page 49 of Quentin's PhD
  double beta{0.0};    // non-associativity characterisation (beta = 3 => associated)
  double beta_p{0.0};  // plastic compliance (1/H)
  double kappa{0.0};   // shear hardening parameter (often 0)
  double varphi{0.0};  // friction angle (characteristic state)
  double pc0{0.0};     // pre-consolidation pressure

  double phiStar0{0.0};   // poro libérable initiale
  linearToPlateau bfunc;  //
  double ginf{0.0};
  double epv0{0.0};
  double l0{0.0};

  SinfoniettaCrush(double young = 200.0e6, double poisson = 0.2);
  std::string getRegistrationName();
  void read(std::istream& is);
  void write(std::ostream& os);
  void init(MaterialPoint& MP);
  void updateStrainAndStress(MPMbox& MPM, size_t p);
  double getYoung();
  double getPoisson();

 private:
  double func_f(mat9r Sigma, double q);
  double func_h(mat9r Sigma, double q, mat9r Ep);
  mat9r func_df_ds(mat9r Sigma, double q);
  mat9r func_dg_ds(mat9r Sigma, double q);
  double func_df_dq(mat9r Sigma, double q);

  const double epsilon_pressure{1e-13};
  const double yield_tol{1e-8};
  double z{0.0};  // depends on varphi
};

#endif /* end of include guard: SINFONIETTA_CRUSH_HPP */
