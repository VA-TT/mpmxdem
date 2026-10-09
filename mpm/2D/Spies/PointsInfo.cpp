#include "PointsInfo.hpp"

#include "Core/MPMbox.hpp"
#include "Core/MaterialPoint.hpp"
#include "Obstacles/Obstacle.hpp"

#include "fileTool.hpp"


void PointsInfo::read(std::istream& is) {
  std::string Filename;
  is >> nrec >> Filename;
  nstep = nrec;  
	
  
  filename = box->result_folder + fileTool::separator() + Filename;
  std::cout << "PointsInfo: filename is " << filename << std::endl;
  std::cout << "Total Material Points: " << box->MP.size() << std::endl;
  // Only open file in computation mode to prevent overwriting during visualization
  if (box->computationMode) {
    file.open(filename.c_str());
  }
  file << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  pos0.resize(box->MP.size(), vec2r{0, 0});
  for (size_t p = 0; p < box->MP.size(); p++) {
    pos0[p] = box->MP[p].pos;
}
    
}

void PointsInfo::exec() {
  }



void PointsInfo::record() {
  // Only record if file is open (i.e., if we're in computation mode)
  if (!file.is_open()) return;
  for (size_t p = 0; p < box->MP.size(); p++) {
    file << box->t << " " << box->MP[p].pos << " " << box->MP[p].pos - pos0[p] << std::endl;
}
}
void PointsInfo::end() {
  pos0.clear();
  // fileSlices.close();
  // file.close();
}
