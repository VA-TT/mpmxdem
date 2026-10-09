#include "BeamDev.hpp"

#include "Core/MPMbox.hpp"
#include "Core/MaterialPoint.hpp"
#include "Obstacles/Obstacle.hpp"

#include "fileTool.hpp"

void BeamDev::read(std::istream& is) {
  std::string Filename;
  is >> nrec >> Filename;
  nstep = nrec;  

 
  filename = box->result_folder + fileTool::separator() + Filename;
  std::cout << "BeamDev: filename is " << filename << std::endl;
  // Only open file in computation mode to prevent overwriting during visualization
  if (box->computationMode) {
    file.open(filename.c_str());
  }
  file << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  
  double right = -1e10, top = -1e10, bottom = 1e10;
  for (size_t p = 0; p < box->MP.size(); p++) {
    double x = box->MP[p].pos.x; 
    double y = box->MP[p].pos.y; 
    if (x>right){
    	right = x;
    }
        if (y < bottom){
    	bottom = y;
    }
        if (y>top){
    	top = y;
    }
  }
  
  double ymid = 0.5*(bottom+top);
  id_MP_end = 0;
  double dst2Min = 1e10;
   for (size_t p = 0; p < box->MP.size(); p++) {
    double x = box->MP[p].pos.x; 
    double y = box->MP[p].pos.y; 
    double dx = x - right;
    double dy = y - ymid;
    double dst2 = dx*dx+dy*dy;
    if (dst2 < dst2Min) {
    dst2Min = dst2;
    id_MP_end = p;
    }
    
   }
  pos0 = box->MP[id_MP_end].pos;
  std::cout << "BeamDev: id MP is " << id_MP_end << std::endl;
   
}

void BeamDev::exec() {
KineticEnergy = 0.0;
  for (size_t p = 0; p < box->MP.size(); p++) {
    vec2r v = box->MP[p].vel;
    KineticEnergy += 0.5*box->MP[p].mass*v*v;
  }
}

void BeamDev::record() {
  // Only record if file is open (i.e., if we're in computation mode)
  if (!file.is_open()) return;
  file << box->t << " " << KineticEnergy << " " << box->MP[id_MP_end].pos - pos0 << std::endl;
}

void BeamDev::end() {

  // fileSlices.close();
  // file.close();
}
