#include "see.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <typeinfo>

#include "ConstitutiveModels/CHCL_DEM.hpp"
#include "ConstitutiveModels/CHCL_DEM_Multi.hpp"
#include "ConstitutiveModels/MohrCoulomb.hpp"
#include "Obstacles/Circle.hpp"
#include "Obstacles/Line.hpp"

float computeShearRate(size_t i) {
  const float Lxx = (float)SmoothedData[i].velGrad.xx;
  const float Lyy = (float)SmoothedData[i].velGrad.yy;
  const float Lxy = (float)SmoothedData[i].velGrad.xy;
  const float Lyx = (float)SmoothedData[i].velGrad.yx;
  return sqrtf((Lxx - Lyy) * (Lxx - Lyy) + (Lxy + Lyx) * (Lxy + Lyx));
}

float computeEquivalentShearStrain(size_t i) {
  const float exx = (float)Conf.MP[i].strain.xx;
  const float eyy = (float)Conf.MP[i].strain.yy;
  const float exy = 0.5f * (float)(Conf.MP[i].strain.xy + Conf.MP[i].strain.yx);
  return sqrtf((exx - eyy) * (exx - eyy) + 4.0f * exy * exy);
}

float percentile(std::vector<float> values, float fraction) {
  if (values.empty()) return 0.0f;
  const size_t index = std::min(values.size() - 1, (size_t)(fraction * (float)(values.size() - 1)));
  std::nth_element(values.begin(), values.begin() + index, values.end());
  return values[index];
}

float computePressure(size_t i) {
  const bool is3D = Conf.MP[i].isDoubleScale;
  const float sx  = (float)SmoothedData[i].stress.xx;
  const float sy  = (float)SmoothedData[i].stress.yy;
  const float sz  = is3D ? (float)SmoothedData[i].outOfPlaneStress : 0.0f;
  return is3D ? -(sx + sy + sz) / 3.0f : -(sx + sy) / 2.0f;
}

float loadRepresentativeDiameter(const std::string& fileName) {
  static std::map<std::string, float> cachedDiameters;
  const auto cached = cachedDiameters.find(fileName);
  if (cached != cachedDiameters.end()) return cached->second;

  PBC3Dbox sample;
  sample.loadConf(fileName.c_str());
  const float diameter = 2.0f * (float)sample.Rmean;
  cachedDiameters[fileName] = diameter;
  return diameter;
}

float computeRepresentativeDiameter(size_t i) {
  if (i < ADs.size() && ADs[i].Rmean > 0.0) return 2.0f * (float)ADs[i].Rmean;

  const auto* multi = dynamic_cast<const CHCL_DEM_Multi*>(Conf.MP[i].constitutiveModel);
  if (multi != nullptr) return loadRepresentativeDiameter(multi->getFileNameFor(Conf.MP[i]));

  const auto* chcl = dynamic_cast<const CHCL_DEM*>(Conf.MP[i].constitutiveModel);
  if (chcl != nullptr) return loadRepresentativeDiameter(chcl->fileName);

  return 0.0f;
}

float computeRheologyFriction(size_t i) {
  const bool is3D = Conf.MP[i].isDoubleScale;
  const float sx  = (float)SmoothedData[i].stress.xx;
  const float sy  = (float)SmoothedData[i].stress.yy;
  const float sxy = 0.5f * (float)(SmoothedData[i].stress.xy + SmoothedData[i].stress.yx);
  const float sz  = is3D ? (float)SmoothedData[i].outOfPlaneStress : 0.0f;

  const float p = computePressure(i);
  if (p <= 1e-6f) return 0.0f;

  const float sxx = sx + p;
  const float syy = sy + p;
  const float szz = is3D ? sz + p : 0.0f;
  const float sxz = is3D ? (float)SmoothedData[i].outOfPlaneShearXZ : 0.0f;
  const float syz = is3D ? (float)SmoothedData[i].outOfPlaneShearYZ : 0.0f;
  const float tau = sqrtf(0.5f * (sxx * sxx + syy * syy + szz * szz + 2.0f * (sxy * sxy + sxz * sxz + syz * syz)));

  return tau / p;
}

float computeInertialNumber(size_t i) {
  const float pressure  = computePressure(i);
  const float shearRate = computeShearRate(i);
  const float mass = (float)SmoothedData[i].rho * (float)Conf.MP[i].vol;
  const bool uses3DStress = Conf.MP[i].isDoubleScale;
  if (uses3DStress) {
    const float diameter = computeRepresentativeDiameter(i);
    if (diameter <= 0.0f) return 0.0f;
    return shearRate * sqrtf(mass / (diameter * (fabsf(pressure) + 1e-6f)));
  }
  return shearRate * sqrtf(mass / (fabsf(pressure) + 1e-6f));
}

float computePackingFraction(size_t i) {
  if (!Conf.MP[i].isDoubleScale) return 0.0f;

  const auto *model = dynamic_cast<const CHCL_DEM *>(Conf.MP[i].constitutiveModel);
  if (model == nullptr) return 0.0f;

  const auto *multi = dynamic_cast<const CHCL_DEM_Multi *>(model);
  const std::string sampleFile = multi ? multi->getFileNameFor(Conf.MP[i]) : model->fileName;
  static std::string loadedFile;
  static double grainDensity = 0.0;
  if (loadedFile != sampleFile) {
    PBC3Dbox sample;
    sample.loadConf(sampleFile.c_str());
    grainDensity = sample.density;
    loadedFile   = sampleFile;
  }

  if (grainDensity <= 0.0) return 0.0f;
  return (float)(SmoothedData[i].rho / grainDensity);
}

float computeMohrFriction(size_t i) {
  const bool is3D = Conf.MP[i].isDoubleScale;
  const float sx  = (float)SmoothedData[i].stress.xx;
  const float sy  = (float)SmoothedData[i].stress.yy;
  const float sxy = 0.5f * (float)(SmoothedData[i].stress.xy + SmoothedData[i].stress.yx);
  const float sz  = is3D ? (float)SmoothedData[i].outOfPlaneStress : 0.0f;

  // Largest Mohr circle, including sigma_zz for a 3D REV.
  const float radius2D   = sqrtf(0.25f * (sx - sy) * (sx - sy) + sxy * sxy);
  const float sigmaA     = -0.5f * (sx + sy) + radius2D;
  const float sigmaB     = -0.5f * (sx + sy) - radius2D;
  const float sigma1     = is3D ? std::max({sigmaA, sigmaB, -sz}) : std::max(sigmaA, sigmaB);
  const float sigma3     = is3D ? std::min({sigmaA, sigmaB, -sz}) : std::min(sigmaA, sigmaB);
  const float mohrCenter = 0.5f * (sigma1 + sigma3);
  const float mohrRadius = 0.5f * (sigma1 - sigma3);

  const auto *model     = dynamic_cast<const MohrCoulomb *>(Conf.MP[i].constitutiveModel);
  const float cohesion  = model == nullptr ? 0.0f : std::max(0.0f, (float)model->Cohesion);
  const float amplitude = hypotf(mohrCenter, cohesion);

  if (amplitude > 1e-10f && mohrRadius < amplitude) {
    const float phi = asinf(mohrRadius / amplitude) - atan2f(cohesion, mohrCenter);
    return sinf(std::max(0.0f, phi));
  }

  return 0.0f;
}

void printHelp() {
  std::cout << " B/S     | shear band contrast / equivalent shear strain" << std::endl;

  std::cout << std::endl;
  std::cout << "-KEY-----|--DESCRIPTION-----------------------------" << std::endl;
  std::cout << " +       | load next configuration file" << std::endl;
  std::cout << " -       | load previous configuration file" << std::endl;
  std::cout << " =       | fit the view" << std::endl;
  std::cout << " a       | MP colors = sig__/sig__ (key to be changed)" << std::endl;
  std::cout << " b       | show/hide colored gradient background" << std::endl;
  std::cout << " c       | show/hide MP contours" << std::endl;
  std::cout << " d       | show/hide node DOFs" << std::endl;
  // std::cout << " e       | ___" << std::endl;
  // std::cout << " f       | ___" << std::endl;
  std::cout << " g       | show/hide background grid" << std::endl;
  std::cout << " h       | print this help" << std::endl;
  std::cout << " i       | print info" << std::endl;
  // std::cout << " j       | ___" << std::endl;
  // std::cout << " k       | ___" << std::endl;
  // std::cout << " l       | ___" << std::endl;
  std::cout << " m       | show/hide Material Points" << std::endl;
  std::cout << " n       | go to conf number (input in console)" << std::endl;
  // std::cout << " o       | ___" << std::endl;
  // std::cout << " p       | ___" << std::endl;
  std::cout << " q       | quit" << std::endl;
  // std::cout << " r       | ___" << std::endl;
  std::cout << " s       | enable/disable MP representation with deformed box" << std::endl;
  // std::cout << " t       | ___" << std::endl;
  // std::cout << " u       | ___" << std::endl;
  // std::cout << " v       | ___" << std::endl;
  // std::cout << " w       | ___" << std::endl;
  std::cout << " x       | show/hide MP stress directions" << std::endl;
  // std::cout << " y       | ___" << std::endl;
  std::cout << " z       | make a single screenshot" << std::endl;
  std::cout << " Z       | make a series of screenshots" << std::endl;
  std::cout << "---------|--MP COLORS-------------------------------" << std::endl;
  std::cout << " 0       | none" << std::endl;
  std::cout << " 1       | velocity magnitude" << std::endl;
  std::cout << " 2       | pressure" << std::endl;
  std::cout << " 3       | density" << std::endl;
  std::cout << " 4       | sig_yy" << std::endl;
  std::cout << " 5       | DEM-cell damage (MPMxDEM only)" << std::endl;
  std::cout << " 6       | deps_q" << std::endl;
  std::cout << " 7       | volume variation" << std::endl;
  std::cout << " 8       | fx (with obstacles, for debug)" << std::endl;
  std::cout << " 9       | fy (with obstacles, for debug)" << std::endl;
  std::cout << "---------|--USE MENU FOR OTHER COLORS---------------" << std::endl;
  std::cout << std::endl;
  std::cout << " I/M/R   | inertial number / principal friction / rheology friction" << std::endl;
  std::cout << " p       | export rheology data" << std::endl;
}

void printInfo() {
  std::cout << "\nCurrent Conf = " << confNum << "\n\n";
  // to be continued
}

void printInertialNumbers() {
  char fileName[256];
  snprintf(fileName, 256, "rheology_conf%d.txt", confNum);
  std::ofstream outFile(fileName);

  if (!outFile.is_open()) {
    std::cout << "Error: Cannot open file " << fileName << " for writing.\n";
    return;
  }

  const bool useDEM =
      std::any_of(Conf.MP.begin(), Conf.MP.end(), [](const MaterialPoint &mp) { return mp.isDoubleScale; });
  outFile << "# MP_index  MP_x  MP_y  I  mu_mohr_sin_phi  mu_rheology_tau_over_p  packing_fraction_Phi ("
          << (useDEM ? "DEM" : "MPM") << " mode, smoothed data)\n";
  std::cout << "Exporting I, mu_mohr = sin(phi), mu_rheology = tau/p, packing fraction Phi\n";

  for (size_t i = 0; i < Conf.MP.size(); i++) {
    outFile << i << "  " << Conf.MP[i].pos.x << "  " << Conf.MP[i].pos.y << "  " << computeInertialNumber(i) << "  "
            << computeMohrFriction(i) << "  " << computeRheologyFriction(i) << "  " << computePackingFraction(i)
            << "\n";
  }

  outFile.close();
  std::cout << "Data saved to " << fileName << " (" << Conf.MP.size() << " MPs)\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case '0': {
    precomputeColors(0);
  } break;

  case '1': {
    precomputeColors(1);
  } break;

  case '2': {
    precomputeColors(2);
  } break;

  case '3': {
    precomputeColors(3);
  } break;

  case '4': {
    precomputeColors(4);
  } break;

  case '5': {
    if (!ADs.empty()) precomputeColors(5);
  } break;

  case '6': {
    precomputeColors(6);
  } break;

  case '7': {
    precomputeColors(7);
  } break;

  case '8': {
    precomputeColors(8);
  } break;

  case '9': {
    if (!ADs.empty()) precomputeColors(9);
  } break;

  case 'a': { // tmp: change the key
    precomputeColors(11);
  } break;

  case 'b': {
    show_background = 1 - show_background;
  } break;

  case 'c': {
    MP_contour = 1 - MP_contour;
  } break;

  case 'd': {
    show_node_dofs = 1 - show_node_dofs;
  } break;

  case 'g': {
    show_grid = 1 - show_grid;
  } break;

  case 'h': {
    printHelp();
  } break;

  case 'i': {
    printInfo();
  } break;

  case 'I': {
    precomputeColors(10);
  } break;

  case 'M': {
    precomputeColors(11);
  } break;

  case 'R': {
    precomputeColors(12);
  } break;

  case 'S': {
    precomputeColors(13);
  } break;

  case 'B': {
    precomputeColors(14);
  } break;

  case 'p': {
    printInertialNumbers();
  } break;

  case 'm': {
    show_MPs = 1 - show_MPs;
  } break;

  case 'n': {
    std::cout << "number: ";
    int num;
    std::cin >> num;
    try_to_readConf(num, Conf, confNum);
  } break;

  case 'q': {
    exit(0);
  } break;

  case 's': {
    MP_deformed_shape = 1 - MP_deformed_shape;
  } break;

  case 'x': {
    show_stress_directions = 1 - show_stress_directions;
  } break;

  case 'z': {
    std::cout << "image saved in 'oneshot.tga'\n";
    screenshot("oneshot.tga");
  } break;
  case 'Z': {
    // be carreful there's no way to stop this loop
    // if the process if too long
    while (try_to_readConf(confNum + 1, Conf, confNum)) {
      char name[256];
      snprintf(name, 256, "shot%d.tga", confNum);
      display();
      screenshot(name);
    }
    std::cout << "series of images saved in 'shot<n>.tga'\n";
  } break;

  case '-': {
    if (confNum > 0) try_to_readConf(confNum - 1, Conf, confNum);
  } break;

  case '+': {
    try_to_readConf(confNum + 1, Conf, confNum);
  } break;

  case '=': {
    fit_view();
    reshape(width, height);
  } break;

  } // end switch(Key)

  glutPostRedisplay();
}

void updateTextLine() {
  textZone.addLine("conf%d,  t = %0.4g s", confNum, Conf.t);
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    // display();
    glutPostRedisplay();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON: {
      int modifiers = glutGetModifiers();
      if (modifiers == GLUT_ACTIVE_SHIFT) mouse_mode = PAN;
      else if (modifiers == GLUT_ACTIVE_CTRL) mouse_mode = ZOOM;
      else mouse_mode = ROTATION;
    } break;
    case GLUT_MIDDLE_BUTTON: // will be removed (?)
      mouse_mode = ZOOM;
      break;
    }
  }
}

void motion(int x, int y) {
  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case ZOOM: {
    double ddy = (worldBox.ymax - worldBox.ymin) * dy;
    double ddx = (worldBox.xmax - worldBox.xmin) * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax += ddx;
    worldBox.ymin -= ddy;
    worldBox.ymax += ddy;
  } break;

  case PAN: {
    double ddx = (worldBox.xmax - worldBox.xmin) * dx;
    double ddy = (worldBox.ymax - worldBox.ymin) * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax -= ddx;
    worldBox.ymin += ddy;
    worldBox.ymax += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  // display();
  glutPostRedisplay();
}

void display() {
  glTools::clearBackground((bool)show_background);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_grid == 1) drawGrid();
  if (show_MPs == 1) drawMPs();
  if (show_obstacles == 1) drawObstacles();
  if (show_stress_directions == 1) drawStressDirections();

  if (color_option > 0) colorBar.show(width, height, colorTable);

  textZone.draw();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  worldBox.xmin = 0.0;
  worldBox.ymin = 0.0;
  worldBox.xmax = 0.0;
  worldBox.ymax = 0.0;
  for (size_t i = 0; i < Conf.nodes.size(); i++) {
    if (Conf.nodes[i].pos.x > worldBox.xmax) worldBox.xmax = Conf.nodes[i].pos.x;
    if (Conf.nodes[i].pos.y > worldBox.ymax) worldBox.ymax = Conf.nodes[i].pos.y;
  }
}

void reshape(int w, int h) {
  width          = w;
  height         = h;
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  double left    = worldBox.xmin;
  double right   = worldBox.xmax;
  double bottom  = worldBox.ymin;
  double top     = worldBox.ymax;
  double worldW  = right - left;
  double worldH  = top - bottom;
  double dW      = 0.2 * worldW;
  double dH      = 0.2 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;
  if (worldW >= worldH) {
    worldH = worldW / aspect;
    top    = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * aspect;
    right  = 0.5 * (left + right + worldW);
    left   = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  //glutPostRedisplay();
}

void drawGrid() {
  glColor4f(0.6f, 0.6f, 0.6f, 1.0f);
  glLineWidth(1.0f);

  std::vector<node> &nodes = Conf.nodes;
  size_t *I;
  for (size_t e = 0; e < Conf.Elem.size(); e++) {
    I = &(Conf.Elem[e].I[0]);
    glBegin(GL_LINE_LOOP);
    for (size_t r = 0; r < 4; r++) glVertex2d(nodes[I[r]].pos.x, nodes[I[r]].pos.y);
    glEnd();
  }

  if (show_node_dofs) {

    double s = 0.25 * 0.5 * (Conf.Grid.lx + Conf.Grid.ly);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glLineWidth(1.5f);
    for (size_t n = 0; n < nodes.size(); n++) {
      if (nodes[n].xfixed) {
        glBegin(GL_LINE_LOOP);
        glVertex2d(nodes[n].pos.x, nodes[n].pos.y);
        glVertex2d(nodes[n].pos.x - s, nodes[n].pos.y + s / 3);
        glVertex2d(nodes[n].pos.x - s, nodes[n].pos.y - s / 3);
        glEnd();
      }
      if (nodes[n].yfixed) {
        glBegin(GL_LINE_LOOP);
        glVertex2d(nodes[n].pos.x, nodes[n].pos.y);
        glVertex2d(nodes[n].pos.x - s / 3, nodes[n].pos.y - s);
        glVertex2d(nodes[n].pos.x + s / 3, nodes[n].pos.y - s);
        glEnd();
      }
    }
  }

  if (show_node_velocity_directions) {
    double s = 0.25 * sqrt(2) * (Conf.Grid.lx + Conf.Grid.ly);
    glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    for (size_t n = 0; n < nodes.size(); n++) {
      if (fabs(nodes[n].q.x) < 1e-10 && fabs(nodes[n].q.y) < 1e-10) { continue; }
      vec2r dir = nodes[n].q;
      dir.normalize();
      dir *= s;
      glVertex2d(nodes[n].pos.x, nodes[n].pos.y);
      glVertex2d(nodes[n].pos.x + dir.x, nodes[n].pos.y + dir.y);
    }
    glEnd();
  }
}

void precomputeColors(int n) {
  precompColors.clear();
  precompColors.resize(SmoothedData.size());

  if (n >= 0 && color_option != n) { color_option = n; }

  switch (color_option) {

  case 0: {
    for (size_t i = 0; i < SmoothedData.size(); i++) { precompColors[i].set(204, 204, 230, 255); }
  } break;

  case 1: {

    colorBar.setTitle("velocity magnitude");
    float vmax = 0.0f;
    float vmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float vel = (float)norm(SmoothedData[i].vel);
      if (vel > vmax) vmax = vel;
      if (vel < vmax) vmin = vel;
    }
    colorTable.setMinMax(0.0f, vmax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by velocity magnitude (vmin = " << vmin << ", vmax = " << vmax << ")\n";

    for (size_t i = 0; i < SmoothedData.size(); i++) {
      float vel = (float)norm(SmoothedData[i].vel);
      colorTable.getRGB(vel, &precompColors[i]);
    }
  } break;

  case 2: {

    colorBar.setTitle("pressure");
    float pmax = -std::numeric_limits<float>::max();
    float pmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float p = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      if (p > pmax) pmax = p;
      if (p < pmin) pmin = p;
    }
    colorTable.setMinMax(pmin, pmax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by pressure (pmin = " << pmin << ", pmax = " << pmax << ")\n";

    for (size_t i = 0; i < SmoothedData.size(); i++) {
      float p = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      colorTable.getRGB(p, &precompColors[i]);
    }
  } break;

  case 3: {

    colorBar.setTitle("density");
    float rhomax = 0.0f;
    float rhomin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float rho = (float)SmoothedData[i].rho;
      if (rho > rhomax) rhomax = rho;
      if (rho < rhomin) rhomin = rho;
    }
    colorTable.setMinMax(rhomin, rhomax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by density (rhomin = " << rhomin << ", rhomax = " << rhomax << ")\n";

    for (size_t i = 0; i < SmoothedData.size(); i++) {
      float rho = (float)SmoothedData[i].rho;
      colorTable.getRGB(rho, &precompColors[i]);
    }
  } break;

  case 4: {
    colorBar.setTitle("sig_xx");
    float pmax = -std::numeric_limits<float>::max();
    float pmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float p = (float)SmoothedData[i].stress.xx;
      if (p > pmax) pmax = p;
      if (p < pmin) pmin = p;
    }
    colorTable.setMinMax(pmin, pmax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by sig_xx (s_xx_min = " << pmin << ", s_xx_max = " << pmax << ")\n";

    for (size_t i = 0; i < SmoothedData.size(); i++) {
      float p = (float)SmoothedData[i].stress.xx;
      colorTable.getRGB(p, &precompColors[i]);
    }
  } break;

  case 5: {
    if (ADs.empty()) break;

    colorBar.setTitle("DEM-cell damage (DEM)");
    float Dmax = 0.0f;
    for (size_t i = 0; i < ADs.size(); i++) {
      if (ADsREF[i].NB == 0.0) continue;
      float D = 1.0f - (float)ADs[i].NB / (float)ADsREF[i].NB;
      if (D > Dmax) Dmax = D;
    }
    colorTable.setMinMax(0.0f, Dmax);
    colorTable.setTableID(6);
    colorTable.Rebuild();
    std::cout << "MP colored by DEM-cell damage [0, " << Dmax << "] \n";

    for (size_t i = 0; i < ADs.size(); i++) {
      if (ADsREF[i].NB == 0.0) continue;
      float D = 1.0f - (float)ADs[i].NB / (float)ADsREF[i].NB;
      colorTable.getRGB(D, &precompColors[i]);
    }
  } break;

  case 6: {
    colorBar.setTitle("equivalent shear strain rate");
    float depsqmax = -std::numeric_limits<float>::max();
    float depsqmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float d1    = (float)(SmoothedData[i].velGrad.xx - SmoothedData[i].velGrad.yy);
      float d2    = (float)(SmoothedData[i].velGrad.xy + SmoothedData[i].velGrad.yx);
      float depsq = (float)sqrt(d1 * d1 + d2 * d2);
      if (depsq > depsqmax) depsqmax = depsq;
      if (depsq < depsqmin) depsqmin = depsq;
    }
    colorTable.setMinMax(depsqmin, depsqmax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by equivalent shear strain rate (min = " << depsqmin << ", max = " << depsqmax
              << ")\n";

    for (size_t i = 0; i < SmoothedData.size(); i++) {
      float d1    = (float)(SmoothedData[i].velGrad.xx - SmoothedData[i].velGrad.yy);
      float d2    = (float)(SmoothedData[i].velGrad.xy + SmoothedData[i].velGrad.yx);
      float depsq = (float)sqrt(d1 * d1 + d2 * d2);
      colorTable.getRGB(depsq, &precompColors[i]);
    }
  } break;

  case 7: {
    colorBar.setTitle("volume variation");
    float volvarmax = -std::numeric_limits<float>::max();
    float volvarmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      if (fabs(MPREF[i].F.xx * MPREF[i].F.yy - MPREF[i].F.xy * MPREF[i].F.yx) == 0) continue;
      float volvar = (float)(fabs(Conf.MP[i].F.xx * Conf.MP[i].F.yy - Conf.MP[i].F.xy * Conf.MP[i].F.yx) /
                             fabs(MPREF[i].F.xx * MPREF[i].F.yy - MPREF[i].F.xy * MPREF[i].F.yx)) -
                     1.0f;
      if (volvar > volvarmax) volvarmax = volvar;
      if (volvar < volvarmin) volvarmin = volvar;
    }
    colorTable.setMinMax(volvarmin, volvarmax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by volvar (volvar_min = " << volvarmin << ", volvar_max = " << volvarmax << ")\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      if (fabs(MPREF[i].F.xx * MPREF[i].F.yy - MPREF[i].F.xy * MPREF[i].F.yx) == 0) continue;
      float volvar = (float)(fabs(Conf.MP[i].F.xx * Conf.MP[i].F.yy - Conf.MP[i].F.xy * Conf.MP[i].F.yx) /
                             fabs(MPREF[i].F.xx * MPREF[i].F.yy - MPREF[i].F.xy * MPREF[i].F.yx)) -
                     1.0f;
      colorTable.getRGB(volvar, &precompColors[i]);
    }
  } break;

  case 8: {
    colorBar.setTitle("fx");
    float fxmax = -std::numeric_limits<float>::max();
    float fxmin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float fx = (float)Conf.MP[i].contactf.x;
      if (fx > fxmax) fxmax = fx;
      if (fx < fxmin) fxmin = fx;
    }
    colorTable.setMinMax(fxmin, fxmax);
    colorTable.setTableID(2);
    colorTable.Rebuild();
    std::cout << "MP colored by fx (fx_min = " << fxmin << ", fx_max = " << fxmax << ")\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB((float)Conf.MP[i].contactf.x, &precompColors[i]); }
  } break;

  case 9: {
    colorBar.setTitle("fy");
    float fymax = -std::numeric_limits<float>::max();
    float fymin = std::numeric_limits<float>::max();
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float fy = (float)Conf.MP[i].contactf.y;
      if (fy > fymax) fymax = fy;
      if (fy < fymin) fymin = fy;
    }
    colorTable.setMinMax(fymin, fymax);
    colorTable.setTableID(2);
    colorTable.Rebuild();
    std::cout << "MP colored by fy (fy_min = " << fymin << ", fy_max = " << fymax << ")\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB((float)Conf.MP[i].contactf.y, &precompColors[i]); }
  } break;

  case 10: {
    const bool useDEM =
        std::any_of(Conf.MP.begin(), Conf.MP.end(), [](const MaterialPoint &mp) { return mp.isDoubleScale; });
    std::vector<float> inertialNumber(Conf.MP.size(), 0.0f);

    colorBar.setTitle(useDEM ? "Inertial number (MPMxDEM)" : "Inertial number (MPM)");
    for (size_t i = 0; i < Conf.MP.size(); i++) { inertialNumber[i] = computeInertialNumber(i); }

    const auto minmax = std::minmax_element(inertialNumber.begin(), inertialNumber.end());
    const float Imin  = *minmax.first;
    const float Imax  = *minmax.second;
    colorTable.setMinMax(Imin, Imax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by Inertial Number (I_min = " << Imin << ", I_max = " << Imax << ")\n";
    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB(inertialNumber[i], &precompColors[i]); }
  } break;

  case 11: {
    colorBar.setTitle("Mohr friction coefficient sin(phi)");
    std::vector<float> frictionCoefficient(Conf.MP.size(), 0.0f);

    for (size_t i = 0; i < Conf.MP.size(); i++) { frictionCoefficient[i] = computeMohrFriction(i); }

    const auto minmax = std::minmax_element(frictionCoefficient.begin(), frictionCoefficient.end());
    const float mumin = *minmax.first;
    const float mumax = *minmax.second;
    colorTable.setMinMax(mumin, mumax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by Mohr friction coefficient (mu_min = " << mumin << ", mu_max = " << mumax << ")\n";
    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB(frictionCoefficient[i], &precompColors[i]); }
  } break;

  case 12: {
    colorBar.setTitle("rheology mu = tau/p");
    std::vector<float> frictionCoefficient(Conf.MP.size(), 0.0f);

    for (size_t i = 0; i < Conf.MP.size(); i++) { frictionCoefficient[i] = computeRheologyFriction(i); }

    const auto minmax = std::minmax_element(frictionCoefficient.begin(), frictionCoefficient.end());
    const float mumin = *minmax.first;
    const float mumax = *minmax.second;
    colorTable.setMinMax(mumin, mumax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by mu from deviatoric stress (mu_min = " << mumin << ", mu_max = " << mumax << ")\n";
    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB(frictionCoefficient[i], &precompColors[i]); }
  } break;

  case 13: {
    colorBar.setTitle("equivalent total shear strain");
    std::vector<float> shearStrain(Conf.MP.size(), 0.0f);

    for (size_t i = 0; i < Conf.MP.size(); i++) { shearStrain[i] = computeEquivalentShearStrain(i); }

    const auto minmax = std::minmax_element(shearStrain.begin(), shearStrain.end());
    const float strainMin = *minmax.first;
    const float strainMax = *minmax.second;
    colorTable.setMinMax(strainMin, strainMax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by equivalent total shear strain (min = " << strainMin
              << ", max = " << strainMax << ")\n";
    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB(shearStrain[i], &precompColors[i]); }
  } break;

  case 14: {
    colorBar.setTitle("shear-band contrast: log10(shear rate)");
    std::vector<float> logShearRate(Conf.MP.size(), 0.0f);
    float maximumRate = 0.0f;
    for (size_t i = 0; i < Conf.MP.size(); i++) { maximumRate = std::max(maximumRate, computeShearRate(i)); }

    // A relative floor avoids log(0) without introducing a unit-dependent threshold.
    const float floorRate = std::max(1e-30f, maximumRate * 1e-8f);
    for (size_t i = 0; i < Conf.MP.size(); i++) {
      logShearRate[i] = log10f(std::max(computeShearRate(i), floorRate));
    }

    float displayMin = percentile(logShearRate, 0.05f);
    float displayMax = percentile(logShearRate, 0.95f);
    if (displayMax <= displayMin) displayMax = displayMin + 1.0f;
    colorTable.setMinMax(displayMin, displayMax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by log10(equivalent shear strain rate), clipped to percentiles 5-95 (min = "
              << displayMin << ", max = " << displayMax << ")\n";
    for (size_t i = 0; i < Conf.MP.size(); i++) { colorTable.getRGB(logShearRate[i], &precompColors[i]); }
  } break;

  default: {
    for (size_t i = 0; i < Conf.MP.size() /*SmoothedData.size()*/; i++) { precompColors[i].set(204, 204, 230, 255); }
  } break;
  } // end switch(color_option)
}

void setColor(int i, float alpha) {
  glColor4f((GLfloat)precompColors[i].r / 255.0f, (GLfloat)precompColors[i].g / 255.0f,
            (GLfloat)precompColors[i].b / 255.0f, alpha);
}

void drawStressDirections() {

  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

  double fac = 2.5e-6;

  mat4r S, V, D;
  for (size_t i = 0; i < Conf.MP.size(); ++i) {
    double xc = Conf.MP[i].pos.x;
    double yc = Conf.MP[i].pos.y;

    S    = SmoothedData[i].stress;
    S.xy = 0.5 * (S.xy + S.yx);
    S.yx = S.xy;
    S.sym_eigen(V, D);

    glBegin(GL_LINES);
    if (fabs(D.xx) > fabs(D.yy)) {
      glLineWidth(2.0f);
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glVertex2d(xc - fac * V.xx * fabs(D.xx), yc - fac * V.yx * fabs(D.xx));
      glVertex2d(xc + fac * V.xx * fabs(D.xx), yc + fac * V.yx * fabs(D.xx));
      glLineWidth(1.0f);
      glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
      glVertex2d(xc - fac * V.xy * fabs(D.yy), yc - fac * V.yy * fabs(D.yy));
      glVertex2d(xc + fac * V.xy * fabs(D.yy), yc + fac * V.yy * fabs(D.yy));
    } else {
      glLineWidth(1.0f);
      glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
      glVertex2d(xc - fac * V.xx * fabs(D.xx), yc - fac * V.yx * fabs(D.xx));
      glVertex2d(xc + fac * V.xx * fabs(D.xx), yc + fac * V.yx * fabs(D.xx));
      glLineWidth(2.0f);
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glVertex2d(xc - fac * V.xy * fabs(D.yy), yc - fac * V.yy * fabs(D.yy));
      glVertex2d(xc + fac * V.xy * fabs(D.yy), yc + fac * V.yy * fabs(D.yy));
    }

    glEnd();
  }
}

void drawMPs() {
  if (mouse_mode != NOTHING) {
    if (show_grid == 0) { drawGrid(); }
    return;
  }

  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.MP.size(); ++i) {

    double xc = Conf.MP[i].pos.x;
    double yc = Conf.MP[i].pos.y;
    double R  = 0.5 * Conf.MP[i].size;

    if (MP_deformed_shape == 1) {
      setColor((int)i);

      glBegin(GL_POLYGON);
      for (size_t r = 0; r < 4; r++) { glVertex2d(SmoothedData[i].corner[r].x, SmoothedData[i].corner[r].y); }
      glEnd();

      if (MP_contour == 1) {
        glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
        glBegin(GL_LINE_LOOP);
        for (size_t r = 0; r < 4; r++) { glVertex2d(SmoothedData[i].corner[r].x, SmoothedData[i].corner[r].y); }
        glEnd();
      }

    } else {
      setColor((int)i);

      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2d(xc + R * cos(angle), yc + R * sin(angle));
      }
      glEnd();

      if (MP_contour == 1) {
        glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
        glBegin(GL_LINE_LOOP);
        for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
          glVertex2d(xc + R * cos(angle), yc + R * sin(angle));
        }
        glEnd();
      }
    }
  }
}

void drawObstacles() {
  for (size_t iObst = 0; iObst < Conf.Obstacles.size(); iObst++) {
    if (Conf.Obstacles[iObst]->getRegistrationName() == "Circle") {
      Circle *C = static_cast<Circle *>(Conf.Obstacles[iObst]);

      glColor4f(0.5f, 0.0f, 0.0f, 0.1f);
      glLineWidth(1.0f);
      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2d(C->pos.x + C->R * cos(angle), C->pos.y + C->R * sin(angle));
      }
      glEnd();

      glColor4f(0.5f, 0.0f, 0.0f, 1.0f);
      glLineWidth(2.0f);
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2d(C->pos.x + C->R * cos(angle), C->pos.y + C->R * sin(angle));
      }
      glEnd();

    } else if (Conf.Obstacles[iObst]->getRegistrationName() == "Line") {
      Line *L = static_cast<Line *>(Conf.Obstacles[iObst]);

      glColor4f(0.5f, 0.0f, 0.0f, 0.1f);
      glLineWidth(1.0f);
      double w = (Conf.Grid.lx + Conf.Grid.ly) * 0.5;
      glBegin(GL_POLYGON);
      glVertex2d(L->pos.x, L->pos.y);
      glVertex2d(L->pos.x + L->len * L->udir.x, L->pos.y + L->len * L->udir.y);
      glVertex2d(L->pos.x + L->len * L->udir.x - w * L->n.x, L->pos.y + L->len * L->udir.y - w * L->n.y);
      glVertex2d(L->pos.x - w * L->n.x, L->pos.y - w * L->n.y);
      glEnd();

      glColor4f(0.5f, 0.0f, 0.0f, 1.0f);
      glLineWidth(2.0f);
      glBegin(GL_LINES);
      glVertex2d(L->pos.x, L->pos.y);
      glVertex2d(L->pos.x + L->len * L->udir.x, L->pos.y + L->len * L->udir.y);
      glEnd();
    }
  }
}

bool try_to_readConf(int num, MPMbox &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d.txt", num);
  if (fileTool::fileExists(file_name)) {
    OKNum = num;
    char co_file_name[256];
    snprintf(co_file_name, 256, "conf%d.txt_micro", num);
    readConf(file_name, co_file_name, CF);
    updateTextLine();
  } else {
    std::cout << file_name << " does not exist" << std::endl;
    return false;
  }
  return true;
}

void readConf(const char *file_name, const char *co_file_name, MPMbox &CF) {
  std::cout << "Read " << file_name << std::endl;
  CF.computationMode = false;
  CF.clean();
  CF.read(file_name);
  CF.postProcess(SmoothedData);
  precomputeColors(color_option);
  textZone.addLine("conf%d,  t = %0.4g s", CF.iconf, CF.t);
  if (fileTool::fileExists(co_file_name)) {
    std::cout << "  with additional data in file " << co_file_name << std::endl;
    readAdditionalData(co_file_name);
  }
  if (MPREF.empty()) { MPREF = CF.MP; }
}

void readAdditionalData(const char *fileName) {
  std::ifstream is(fileName);
  std::string lineTop;
  getline(is, lineTop);
  ADs.clear();
  additionalData D;
  while (is.good()) {
    is >> D.MP_x >> D.MP_y >> D.NInt >> D.NB >> D.TF >> D.FF >> D.Rmean >> D.Vmean >> D.VelMean >> D.VelMin >>
        D.VelMax >> D.VelVar >> D.Vsolid >> D.Vcell >> D.h_xx >> D.h_xy >> D.h_yx >> D.h_yy >> D.ReducedPartDistMean;
    ADs.push_back(D);
  }

  if (ADsREF.empty()) {
    std::cout << "REF ADDITIONAL DATA\n";
    for (size_t i = 0; i < ADs.size(); i++) { ADsREF.push_back(ADs[i]); }
  }
}

int screenshot(const char *filename) {
  // http://forum.devmaster.net/t/rendering-a-single-frame-to-a-file-with-opengl/12469/2

  // we will store the image data here
  unsigned char *pixels;
  // the thingy we use to write files
  FILE *shot;
  // we get the width/height of the screen into this array
  int screenStats[4];

  // get the width/height of the window
  glGetIntegerv(GL_VIEWPORT, screenStats);

  // generate an array large enough to hold the pixel data
  // (width*height*bytesPerPixel)
  pixels = new unsigned char[screenStats[2] * screenStats[3] * 3];
  // read in the pixel data, TGA's pixels are BGR aligned
  glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR, GL_UNSIGNED_BYTE, pixels);

  // open the file for writing. If unsucessful, return 1
  if ((shot = fopen(filename, "wb")) == NULL) return 1;

  // this is the tga header it must be in the beginning of
  // every (uncompressed) .tga
  unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // the header that is used to get the dimensions of the .tga
  // header[1]*256+header[0] - width
  // header[3]*256+header[2] - height
  // header[4] - bits per pixel
  // header[5] - ?
  unsigned char header[6] = {((unsigned char)(screenStats[2] % 256)),
                             ((unsigned char)(screenStats[2] / 256)),
                             ((unsigned char)(screenStats[3] % 256)),
                             ((unsigned char)(screenStats[3] / 256)),
                             24,
                             0};

  // write out the TGA header
  fwrite(TGAheader, sizeof(unsigned char), 12, shot);
  // write out the header
  fwrite(header, sizeof(unsigned char), 6, shot);
  // write the pixels
  fwrite(pixels, sizeof(unsigned char), screenStats[2] * screenStats[3] * 3, shot);

  // close the file
  fclose(shot);
  // free the memory
  delete[] pixels;

  // return success
  return 0;
}

void menu(int num) {
  switch (num) {

  case 0: {
    exit(0);
  } break;

  // Display options
  case 100: {
    show_MPs = 1 - show_MPs;
  } break;
  case 101: {
    show_background = 1 - show_background;
  } break;
  case 102: {
    MP_contour = 1 - MP_contour;
  } break;
  case 103: {
    MP_deformed_shape = 1 - MP_deformed_shape;
  } break;

  case 104: {
    show_stress_directions = 1 - show_stress_directions;
  } break;

  // Color Material Points
  case 200: {
    precomputeColors(0);
  } break;
  case 201: {
    precomputeColors(1);
  } break;
  case 202: {
    precomputeColors(2);
  } break;
  case 203: {
    precomputeColors(3);
  } break;
  case 204: {
    precomputeColors(4);
  } break;
  case 205: {
    if (!ADs.empty()) precomputeColors(5);
  } break;
  case 206: {
    precomputeColors(6);
  } break;
  case 207: {
    precomputeColors(7);
  } break;
  case 208: {
    precomputeColors(8);
  } break;
  case 209: {
    precomputeColors(9);
  } break;
  case 210: {
    precomputeColors(10);
  } break;
  case 211: {
    precomputeColors(11);
  } break;
  case 212: {
    precomputeColors(12);
  } break;
  case 213: {
    precomputeColors(13);
  } break;
  case 214: {
    precomputeColors(14);
  } break;

  // Grid informations
  case 300: {
    show_grid = 1 - show_grid;
  } break;
  case 301: {
    show_node_dofs = 1 - show_node_dofs;
  } break;
  case 302: {
    show_node_velocity_directions = 1 - show_node_velocity_directions;
  } break;

  default: {
    std::cout << "Not yet plugged!" << std::endl;
  }

  } // end switch

  glutPostRedisplay();  
}

void buildMenu() {
  int submenu100 = glutCreateMenu(menu); // Display options
  glutAddMenuEntry("Show/Hide Material Points", 100);
  glutAddMenuEntry("Show/Hide Background", 101);
  glutAddMenuEntry("Show/Hide MP Contours", 102);
  glutAddMenuEntry("Show/Hide MP Deformed Shape", 103);
  glutAddMenuEntry("Show/Hide MP Stress Directions", 104);

  int submenu200 = glutCreateMenu(menu); // Color Material Points
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Velocity Magnitude", 201);
  glutAddMenuEntry("Pressure, tr(Sig)/3", 202);
  glutAddMenuEntry("Density", 203);
  glutAddMenuEntry("Sig xx", 204);
  glutAddMenuEntry("DEM-cell damage (DEM)", 205);
  glutAddMenuEntry("Equivalent Shear Strain Rate (active band)", 206);
  glutAddMenuEntry("Volume Variation", 207);
  glutAddMenuEntry("Fx (with Obstacle)", 208);
  glutAddMenuEntry("Fy (with Obstacle)", 209);
  glutAddMenuEntry("Inertial Number (DEM/MPM)", 210);
  glutAddMenuEntry("mu from Principal Stresses", 211);
  glutAddMenuEntry("mu from Deviatoric Stress (rheology)", 212);
  glutAddMenuEntry("Equivalent Total Shear Strain", 213);
  glutAddMenuEntry("Shear-band Contrast (log, percentile)", 214);

  int submenu300 = glutCreateMenu(menu); // Grid informations
  glutAddMenuEntry("Show/Hide Grid", 300);
  glutAddMenuEntry("Imposed DoF", 301);
  glutAddMenuEntry("Velocity Directions", 302);

  int submenu400 = glutCreateMenu(menu); // Material Points informations
  glutAddMenuEntry("Bidon", 401);
  // ...

  // Main menu
  glutCreateMenu(menu);
  glutAddSubMenu("Display Options", submenu100);
  glutAddSubMenu("Color Material Points", submenu200);
  glutAddSubMenu("Grid Informations", submenu300);
  glutAddSubMenu("Material Points Informations", submenu400);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char *argv[]) {
  // The profiler has to be initialised here too, and not only in Runners/run.cpp:
  // its START_TIMER dereferences the 'current' timer, which INIT_TIMERS is what
  // allocates. Any function shared with mpmbox may carry a timer -- postProcess
  // reaches updateLiveNodeList, for instance -- and would otherwise crash the
  // viewer on a null pointer.
  INIT_TIMERS();

  Conf.computationMode = false;

  if (argc == 1) {
    confNum = 0;
    try_to_readConf(confNum, Conf, confNum);
  } else if (argc == 2) {
    confNum = 0;
    readConf(argv[1], "###", Conf);
  }

  mouse_mode = NOTHING;

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  int X0 = (glutGet(GLUT_SCREEN_WIDTH) - width) / 2;
  int Y0 = (glutGet(GLUT_SCREEN_HEIGHT) - height) / 2;
  glutInitWindowPosition(X0, Y0);
  glutInitWindowSize(width, height);

  main_window = glutCreateWindow("CONF VISUALIZER (MPMbox)");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // ==== Other initialisations
  glText::init();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  display();
  glutMainLoop();
  return 0;
}
