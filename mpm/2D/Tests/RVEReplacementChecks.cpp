#include "ConstitutiveModels/CHCL_DEM.hpp"
#include "ConstitutiveModels/CHCL_DEM_Multi.hpp"
#include "ConstitutiveModels/MohrCoulomb.hpp"
#include "Core/MPMbox.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
void near(double a,double b,const char* message) {require(std::fabs(a-b)<1e-10,message);}
mat9r cell(double y) {mat9r h;h.xx=10;h.yy=y;h.zz=10;return h;}
void point(MPMbox& box,ConstitutiveModel* dem,double y) {
 box.MP.emplace_back(0,.2,1603.72072810721,dem);
 auto& p=box.MP.back();p.PBC=new PBC3Dbox;p.PBC->Rmax=.5;p.PBC->Cell.h=cell(y);p.isDoubleScale=true;p.isTracked=true;
 p.stress.xx=-100;p.stress.yy=-100;p.outOfPlaneStress=-90;p.outOfPlaneShearXZ=2;p.outOfPlaneShearYZ=3;p.strain.xx=.3;p.plasticStrain.xx=.01;
 box.prevFstore.push_back(mat4r::unit());
 box.Elem.resize(1);box.nodes.resize(1);box.shapeGradN.resize(element::nbNodes*box.MP.size());box.dt=.001;
}
int main() {
 try {
 element::nbNodes=4;
 near(CHCL_DEM::shortestThickness(cell(2)),2,"orthogonal thickness");
 mat9r skew=cell(2);skew.xy=10;
 near(CHCL_DEM::shortestThickness(skew),20/std::sqrt(104.),"sheared-cell face thickness");
 mat9r rotated=cell(2);rotated.xx=0;rotated.xy=-2;rotated.yx=10;rotated.yy=0;
 near(CHCL_DEM::shortestThickness(rotated),2,"rotation invariance");
 auto invalid=cell(-1);bool rejected=false;try{CHCL_DEM::shortestThickness(invalid);}catch(const std::runtime_error&){rejected=true;}require(rejected,"reject inverted cell");
 MPMbox box;auto* dem=new CHCL_DEM;dem->box=&box;dem->key="DEM";box.models[dem->key]=dem;
 auto* mc=new MohrCoulomb(1e5,.2,.4,0,0);mc->box=&box;mc->key="MC";box.models[mc->key]=mc;
 box.CHCL.replacementModel="MC";point(box,dem,2);
 require(!dem->replaceIfNeeded(box,0),"disabled replacement");
 box.CHCL.replaceThicknessFactor=1.5;
 require(!dem->replaceIfNeeded(box,0),"thick RVE stays DEM");
 auto& p=box.MP[0];p.F.yy=.7;
 const double mass=p.mass,density=p.density,volume=p.vol;
 require(dem->replaceIfNeeded(box,0),"predicted threshold crossing");
 require(p.constitutiveModel==mc && !p.isDoubleScale && p.PBC==nullptr && !p.isTracked,"replacement flags and ownership");
 near(p.mass,mass,"mass preserved");near(p.density,density,"density preserved");near(p.vol,volume,"volume preserved");near(p.stress.xx,-100,"xx stress continuity");near(p.outOfPlaneStress,-90,"zz stress continuity");near(p.outOfPlaneShearXZ,2,"xz stress continuity");near(p.strain.xx,.3,"total strain continuity");near(p.plasticStrain.xx,.01,"plastic strain continuity");
 box.shapeGradN[0]=vec2r(1,0);box.nodes[0].vel=vec2r(-.01,0);
 mc->updateStrainAndStress(box,0);near(p.strain.xx,.3-.00001,"next MC increment exactly once");
 // Equality and inherited multi-VER path are covered independently.
 MPMbox multi;auto* md=new CHCL_DEM_Multi;md->box=&multi;md->key="MULTI";multi.models[md->key]=md;
 auto* mm=new MohrCoulomb(1e5,.2,.4,0,0);mm->box=&multi;mm->key="MC";multi.models[mm->key]=mm;
 multi.CHCL.replaceThicknessFactor=1.5;multi.CHCL.replacementModel="MC";point(multi,md,1.5);
 multi.nodes[0].vel=vec2r(-.01,0);multi.shapeGradN[0]=vec2r(1,0);
 require(md->replaceIfNeeded(multi,0),"threshold equality and multi-VER replacement");near(multi.MP[0].strain.xx,.3-.00001,"switch step MC increment exactly once");
 std::cout<<"PASS: geometry, disabled mode, threshold equality, predicted crossing, single/multi VER, continuity and single strain increment\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
