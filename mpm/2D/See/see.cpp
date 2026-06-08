#include "see.hpp"

#include <fstream>
#include <limits>
#include <typeinfo>

#include "Obstacles/Circle.hpp"
#include "Obstacles/Line.hpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next configuration file" << endl;
  cout << "-         load previous configuration file" << endl;
  cout << "=         fit the view" << endl;
  cout << "i         print current conf info" << endl;
  cout << "I         show Inertial Number (DEM/MPM)" << endl;
  cout << "M         show mu from principal stresses" << endl;
  cout << "R         show mu from deviatoric stress (rheology)" << endl;
  cout << "r         toggle raw/smoothed data for I and mu" << endl;
  cout << "p         save I, mu_stress, mu_rheology to file" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void printInfo() {
  using namespace std;

  cout << "\nCurrent Conf = " << confNum << "\n\n";
}

void printInertialNumbers() {
  char fileName[256];
  snprintf(fileName, 256, "conf%d_rheology%s.txt", confNum, use_raw_data ? "_raw" : "");
  std::ofstream outFile(fileName);

  if (!outFile.is_open()) {
    std::cout << "Error: Cannot open file " << fileName << " for writing.\n";
    return;
  }

  bool useDEM = !ADs.empty();
  outFile << "# MP_index  MP_x  MP_y  I  mu_stress  mu_rheology (" << (useDEM ? "DEM" : "MPM") << " mode, "
          << (use_raw_data ? "raw data" : "smoothed data") << ")\n";
  std::cout << "Exporting I, mu_stress (Mohr circle), mu_rheology (deviatoric norm) - "
            << (use_raw_data ? "RAW data" : "SMOOTHED data") << "\n";

  for (size_t i = 0; i < Conf.MP.size(); i++) {

    // --- stress components ---
    float sx, sy, sxy;
    if (use_raw_data) {
      sx  = (float)Conf.MP[i].stress.xx;
      sy  = (float)Conf.MP[i].stress.yy;
      sxy = 0.5f * ((float)Conf.MP[i].stress.xy + (float)Conf.MP[i].stress.yx);
    } else {
      sx  = (float)SmoothedData[i].stress.xx;
      sy  = (float)SmoothedData[i].stress.yy;
      sxy = 0.5f * ((float)SmoothedData[i].stress.xy + (float)SmoothedData[i].stress.yx);
    }

    // --- mu_stress: Mohr circle, même convention que gnuplot ---
    // sigma1 = -(sx+sy)/2 + R,  sigma3 = -(sx+sy)/2 - R
    // center = -(sx+sy)/2  (nén dương)
    // mu = R / |center|
    float center    = -(sx + sy) * 0.5f;
    float radius    = sqrtf(0.25f * (sx - sy) * (sx - sy) + sxy * sxy);
    float mu_stress = radius / (fabsf(center) + 1e-10f);

    // --- strain rate & pressure ---
    float d1, d2, epsdot, pres;
    if (use_raw_data) {
      d1   = (float)(Conf.MP[i].velGrad.xx - Conf.MP[i].velGrad.yy);
      d2   = (float)(Conf.MP[i].velGrad.xy + Conf.MP[i].velGrad.yx);
      pres = 0.5f * (float)(Conf.MP[i].stress.xx + Conf.MP[i].stress.yy);
    } else {
      d1   = (float)(SmoothedData[i].velGrad.xx - SmoothedData[i].velGrad.yy);
      d2   = (float)(SmoothedData[i].velGrad.xy + SmoothedData[i].velGrad.yx);
      pres = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
    }
    epsdot = sqrtf(d1 * d1 + d2 * d2);

    // --- Inertial Number ---
    float I;
    float rho = use_raw_data ? (float)Conf.MP[i].density : (float)SmoothedData[i].rho;
    if (useDEM) {
      I = (2.0f * ADs[i].Rmean * epsdot) / sqrtf((fabsf(pres) + 1e-6f) / rho);
    } else {
      float m = rho * (float)Conf.MP[i].vol;
      I       = epsdot * sqrtf(m / (fabsf(pres) + 1e-6f));
    }

    // --- mu_rheology: deviatoric stress norm / |P| ---
    float tau_xx      = sx - pres;
    float tau_yy      = sy - pres;
    float tau_xy      = sxy;
    float tau_norm    = sqrtf(tau_xx * tau_xx + tau_yy * tau_yy + 2.0f * tau_xy * tau_xy);
    float mu_rheology = tau_norm / (fabsf(pres) + 1e-10f);

    outFile << i << "  " << Conf.MP[i].pos.x << "  " << Conf.MP[i].pos.y << "  " << I << "  " << mu_stress << "  "
            << mu_rheology << "\n";
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

  case 'd': {
    show_node_dofs = 1 - show_node_dofs;
  } break;

  case 'b': {
    show_background = 1 - show_background;
  } break;

  case 'c': {
    MP_contour = 1 - MP_contour;
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

  case 'r': {
    use_raw_data = !use_raw_data;
    std::cout << "Data mode: " << (use_raw_data ? "RAW" : "SMOOTHED") << std::endl;
    // Recompute current visualization if showing I or mu
    if (color_option >= 10 && color_option <= 12) { precomputeColors(color_option); }
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
  };

  glutPostRedisplay();
}

void updateTextLine() {
  textZone.addLine("conf%d,  t = %0.4g s", confNum, Conf.t);
  // textZone.addLine("%s - time = %.3G", file_name, Conf.t);
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
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
  display();
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

  glutPostRedisplay();
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

  if (n > 0 && color_option != n) { color_option = n; }

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
    colorBar.setTitle("deps_q");
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
    std::cout << "MP colored by deps_q (deps_q_min = " << depsqmin << ", deps_q_max = " << depsqmax << ")\n";

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
    bool useDEM = !ADs.empty();

    if (useDEM) {
      colorBar.setTitle(use_raw_data ? "Inertial number (DEM, raw)" : "Inertial number (DEM)");
      std::cout << "Computing Inertial Number (DEM mode, " << (use_raw_data ? "raw data" : "smoothed data") << ")\n";
    } else {
      colorBar.setTitle(use_raw_data ? "Inertial number (MPM, raw)" : "Inertial number (MPM)");
      std::cout << "Computing Inertial Number (MPM mode, " << (use_raw_data ? "raw data" : "smoothed data") << ")\n";
      std::cout << "  Formula: I = epsdot * sqrt(vol * rho / P)\n";
    }

    float Imax = -std::numeric_limits<float>::max();
    float Imin = std::numeric_limits<float>::max();

    // Debug: track non-zero values
    int nonzero_count = 0;
    float epsdot_max  = 0.0f;

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      // Select data source based on flag
      float d1, d2, epsdot, pres;

      if (use_raw_data) {
        d1     = (float)(Conf.MP[i].velGrad.xx - Conf.MP[i].velGrad.yy);
        d2     = (float)(Conf.MP[i].velGrad.xy + Conf.MP[i].velGrad.yx);
        epsdot = sqrt(d1 * d1 + d2 * d2);
        pres   = 0.5f * (float)(Conf.MP[i].stress.xx + Conf.MP[i].stress.yy);
      } else {
        d1     = (float)(SmoothedData[i].velGrad.xx - SmoothedData[i].velGrad.yy);
        d2     = (float)(SmoothedData[i].velGrad.xy + SmoothedData[i].velGrad.yx);
        epsdot = sqrt(d1 * d1 + d2 * d2);
        pres   = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      }

      float I;
      if (useDEM) {
        // DEM: I = (2 * Rmean * epsdot) / sqrt(P/rho)
        float rho = use_raw_data ? (float)Conf.MP[i].density : (float)SmoothedData[i].rho;
        I         = (float)((2 * ADs[i].Rmean * epsdot) / sqrt((abs(pres) + 1e-6) / rho));
      } else {
        // MPM: I = epsdot * sqrt(m/P) with m = rho * vol
        float rho = use_raw_data ? (float)Conf.MP[i].density : (float)SmoothedData[i].rho;
        float m   = (float)(rho * Conf.MP[i].vol);
        I         = epsdot * sqrt(m / (abs(pres) + 1e-6));
      }

      if (I > Imax) Imax = I;
      if (I < Imin) Imin = I;
      if (I > 1e-10) nonzero_count++;
      if (epsdot > epsdot_max) epsdot_max = epsdot;
    }

    colorTable.setMinMax(Imin, Imax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by Inertial Number (I_min = " << Imin << ", I_max = " << Imax << ")\n";
    std::cout << "  Non-zero MPs: " << nonzero_count << "/" << Conf.MP.size() << ", epsdot_max = " << epsdot_max
              << "\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float d1, d2, epsdot, pres;

      if (use_raw_data) {
        d1     = (float)(Conf.MP[i].velGrad.xx - Conf.MP[i].velGrad.yy);
        d2     = (float)(Conf.MP[i].velGrad.xy + Conf.MP[i].velGrad.yx);
        epsdot = sqrt(d1 * d1 + d2 * d2);
        pres   = 0.5f * (float)(Conf.MP[i].stress.xx + Conf.MP[i].stress.yy);
      } else {
        d1     = (float)(SmoothedData[i].velGrad.xx - SmoothedData[i].velGrad.yy);
        d2     = (float)(SmoothedData[i].velGrad.xy + SmoothedData[i].velGrad.yx);
        epsdot = sqrt(d1 * d1 + d2 * d2);
        pres   = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      }

      float I;
      if (useDEM) {
        // DEM formula
        float rho = use_raw_data ? (float)Conf.MP[i].density : (float)SmoothedData[i].rho;
        I         = (float)((2 * ADs[i].Rmean * epsdot) / sqrt((abs(pres) + 1e-6) / rho));
      } else {
        // MPM formula
        float rho = use_raw_data ? (float)Conf.MP[i].density : (float)SmoothedData[i].rho;
        float m   = (float)(rho * Conf.MP[i].vol);
        I         = epsdot * sqrt(m / (abs(pres) + 1e-6));

        // Debug: print first 3 MPs with non-zero I
        if (i < 3 && I > 1e-10) {
          std::cout << "  MP[" << i << "]: epsdot=" << epsdot << ", pres=" << pres << ", rho=" << rho
                    << ", vol=" << Conf.MP[i].vol << ", I=" << I << "\n";
        }
      }

      colorTable.getRGB(I, &precompColors[i]);
    }
  } break;
  case 11: {
    colorBar.setTitle(use_raw_data ? "mu from principal stresses (raw)" : "mu from principal stresses");
    std::cout << "Computing mu from Mohr circle (" << (use_raw_data ? "raw data" : "smoothed data") << ")\n";

    float mumax = -std::numeric_limits<float>::max();
    float mumin = std::numeric_limits<float>::max();

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float sx, sy, sxy;
      if (use_raw_data) {
        sx  = (float)Conf.MP[i].stress.xx;
        sy  = (float)Conf.MP[i].stress.yy;
        sxy = (float)Conf.MP[i].stress.xy;
      } else {
        sx  = (float)SmoothedData[i].stress.xx;
        sy  = (float)SmoothedData[i].stress.yy;
        sxy = (float)SmoothedData[i].stress.xy;
      }

      // Gnuplot convention: nén dương
      // sigma1 = -(sx+sy)/2 + sqrt(((sx-sy)/2)^2 + sxy^2)
      // sigma3 = -(sx+sy)/2 - sqrt(...)
      float center = -(sx + sy) * 0.5f;
      float radius = sqrtf(0.25f * (sx - sy) * (sx - sy) + sxy * sxy);
      float mu     = radius / (fabsf(center) + 1e-10f);

      if (mu > mumax) mumax = mu;
      if (mu < mumin) mumin = mu;
    }

    colorTable.setMinMax(mumin, mumax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by mu from Mohr circle (mu_min = " << mumin << ", mu_max = " << mumax << ")\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float sx, sy, sxy;
      if (use_raw_data) {
        sx  = (float)Conf.MP[i].stress.xx;
        sy  = (float)Conf.MP[i].stress.yy;
        sxy = (float)Conf.MP[i].stress.xy;
      } else {
        sx  = (float)SmoothedData[i].stress.xx;
        sy  = (float)SmoothedData[i].stress.yy;
        sxy = (float)SmoothedData[i].stress.xy;
      }

      float center = -(sx + sy) * 0.5f;
      float radius = sqrtf(0.25f * (sx - sy) * (sx - sy) + sxy * sxy);
      float mu     = radius / (fabsf(center) + 1e-10f);

      colorTable.getRGB(mu, &precompColors[i]);
    }
  } break;

  case 12: {
    colorBar.setTitle(use_raw_data ? "mu from deviatoric stress (rheology, raw)"
                                   : "mu from deviatoric stress (rheology)");
    std::cout << "Computing mu from deviatoric stress (" << (use_raw_data ? "raw data" : "smoothed data") << ")\n";

    float mumax = -std::numeric_limits<float>::max();
    float mumin = std::numeric_limits<float>::max();

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float pres;
      if (use_raw_data) {
        pres = 0.5f * (float)(Conf.MP[i].stress.xx + Conf.MP[i].stress.yy);
      } else {
        pres = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      }

      // Deviatoric stress: tau = sigma - P*I
      float tau_xx, tau_yy, tau_xy;
      if (use_raw_data) {
        tau_xx = Conf.MP[i].stress.xx - pres;
        tau_yy = Conf.MP[i].stress.yy - pres;
        tau_xy = Conf.MP[i].stress.xy;
      } else {
        tau_xx = SmoothedData[i].stress.xx - pres;
        tau_yy = SmoothedData[i].stress.yy - pres;
        tau_xy = SmoothedData[i].stress.xy;
      }

      // Frobenius norm of deviatoric stress
      float tau_norm = sqrt(tau_xx * tau_xx + tau_yy * tau_yy + 2.0f * tau_xy * tau_xy);

      // mu = |tau| / P
      float mu = tau_norm / (abs(pres) + 1e-10);

      if (mu > mumax) mumax = mu;
      if (mu < mumin) mumin = mu;
    }

    colorTable.setMinMax(mumin, mumax);
    colorTable.setTableID(3);
    colorTable.Rebuild();
    std::cout << "MP colored by mu from deviatoric stress (mu_min = " << mumin << ", mu_max = " << mumax << ")\n";

    for (size_t i = 0; i < Conf.MP.size(); i++) {
      float pres;
      if (use_raw_data) {
        pres = 0.5f * (float)(Conf.MP[i].stress.xx + Conf.MP[i].stress.yy);
      } else {
        pres = 0.5f * (float)(SmoothedData[i].stress.xx + SmoothedData[i].stress.yy);
      }

      float tau_xx, tau_yy, tau_xy;
      if (use_raw_data) {
        tau_xx = Conf.MP[i].stress.xx - pres;
        tau_yy = Conf.MP[i].stress.yy - pres;
        tau_xy = Conf.MP[i].stress.xy;
      } else {
        tau_xx = SmoothedData[i].stress.xx - pres;
        tau_yy = SmoothedData[i].stress.yy - pres;
        tau_xy = SmoothedData[i].stress.xy;
      }

      float tau_norm = sqrt(tau_xx * tau_xx + tau_yy * tau_yy + 2.0f * tau_xy * tau_xy);
      float mu       = tau_norm / (abs(pres) + 1e-10);

      colorTable.getRGB(mu, &precompColors[i]);
    }
  } break;

  default: {
    for (size_t i = 0; i < SmoothedData.size(); i++) { precompColors[i].set(204, 204, 230, 255); }
  } break;
  }
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
      glVertex2d(L->pos.x + L->len * L->t.x, L->pos.y + L->len * L->t.y);
      glVertex2d(L->pos.x + L->len * L->t.x - w * L->n.x, L->pos.y + L->len * L->t.y - w * L->n.y);
      glVertex2d(L->pos.x - w * L->n.x, L->pos.y - w * L->n.y);
      glEnd();

      glColor4f(0.5f, 0.0f, 0.0f, 1.0f);
      glLineWidth(2.0f);
      glBegin(GL_LINES);
      glVertex2d(L->pos.x, L->pos.y);
      glVertex2d(L->pos.x + L->len * L->t.x, L->pos.y + L->len * L->t.y);
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

  }; // end switch

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
  glutAddMenuEntry("Delta Eps Deviator", 206);
  glutAddMenuEntry("Volume Variation", 207);
  glutAddMenuEntry("Fx (with Obstacle)", 208);
  glutAddMenuEntry("Fy (with Obstacle)", 209);
  glutAddMenuEntry("Inertial Number (DEM/MPM)", 210);
  glutAddMenuEntry("mu from Principal Stresses", 211);
  glutAddMenuEntry("mu from Deviatoric Stress (rheology)", 212);

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
