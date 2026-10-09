// Headless checks of the actual viewer calculations and color paths.
// Build: cmake --build BUILD --target check_viewer_formulas
// Run: BUILD/check_viewer_formulas
#define main viewer_main
#include "../See/see.cpp"
#undef main

#include <stdexcept>

static void near(double value, double expected) {
  if (std::abs(value - expected) > 1e-6) throw std::runtime_error("viewer formula mismatch");
}
static void colorMatches(size_t i, double value) {
  colorRGBA expected;
  colorTable.getRGB(value, &expected);
  const auto &actual = precompColors[i];
  if (actual.r != expected.r || actual.g != expected.g || actual.b != expected.b)
    throw std::runtime_error("viewer color uses a different formula");
}
int main() {
  INIT_TIMERS();
  Conf.MP.resize(3);
  SmoothedData.resize(3);
  for (size_t i=0; i<3; ++i) {
    auto &d = SmoothedData[i];
    d.stress.xx = d.stress.yy = d.outOfPlaneStress = -100.0 * (i+1);
    d.outOfPlaneShearXZ = 30.0 * (i+1);
    d.outOfPlaneShearYZ = 40.0 * (i+1);
    // Isotropic IN-PLANE expansion is deviatoric in 3D because Lzz=0.
    d.velGrad.xx = d.velGrad.yy = i+1;
    Conf.MP[i].F = mat4r::unit();
    Conf.MP[i].F.xx = std::exp(0.2);
    Conf.MP[i].F.yy = std::exp(-0.2);
    near(computePressure(i), 100.0*(i+1));
    near(computeShearRate(i), 2.0*(i+1)/std::sqrt(3.0));
    near(computeRheologyFriction(i), 0.5);
    near(computeMohrFriction(i), 0.5);
    near(computeEquivalentShearStrain(i), 0.4);
  }
  precomputeColors(2);
  for (size_t i=0;i<3;++i) colorMatches(i, 100.0*(i+1));
  precomputeColors(6);
  for (size_t i=0;i<3;++i) colorMatches(i, computeShearRate(i));
  Conf.MP.clear();
  SmoothedData.clear();
  precomputeColors(10); // An empty configuration must not dereference minmax.
  std::cout << "Viewer formula and color checks: PASS\n";
}
