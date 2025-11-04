/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2021 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
*                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Drag coefficient of sphere.


// MechSys
#include <mechsys/lbmdem/Domain.h>
#include <vector>
#include <cmath>
#include <random>
#include <list>
#include <algorithm>
#include <iostream>
#include "json.hpp"
#include "getjsonpara.hpp"

struct UserData
{
    real3              * pacc;
    double                 nu;
    double                  R;
};

__global__ void Setup(real3 * BForce, real const * Rho, real3 const * acc, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    BForce[ic] = Rho[ic]*acc[0];
}
void Setup (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    Setup<<<dom.LBMDOM.Ncells/dom.Nthread+1,dom.Nthread>>>(dom.LBMDOM.pBForce,dom.LBMDOM.pRho,dat.pacc,dom.LBMDOM.plbmaux);
}

void Report (LBMDEM::Domain & dom, void * UD)
{
}

static struct  {
public:
    double nu = 1.0;
    size_t nx = 101;
    size_t ny = 101;
    size_t nz = 101;
    double dx = 0.2;
    double dt = 1.6e-2;
    double ParticleRadius = 0.9;
    double ParticleDensity  = 1.;
    double FluidDensity  = 1.;
    double ArtificialFluidAccelerationZ = -1.0e-4;
    double Tf = 1.0e5;
    int OutputStep = 200;
    double h  = 1.;
    double nlayer = 3;
    double d  = 3;
    double lambda  = 6;
    double lambdashift = 0.1 ;
    double dlayer  =1;
    double ObstaclePlaceCenterZPercent = 0.9;
    int SkipVisFactor = 2;

}  Parameter;

void LoadParameter(std::string jsonCFG_path) {  // should be rewritten into a template header
    std::cout << "Loading cfg..." << "\n";
    nlohmann::json j2;
    try{
        std::ifstream i(jsonCFG_path);
        i >> j2;
    }
    catch (std::exception& e){
        std::cout << "exception reading json file\n ";
        return ;
    }
    std::vector<std::string> strs;
    
    getJsonPara( j2, "nu", Parameter.nu);
    getJsonPara( j2, "nx", Parameter.nx);
    getJsonPara( j2, "ny", Parameter.ny);
    getJsonPara( j2, "nz", Parameter.nz);
    getJsonPara( j2, "dx", Parameter.dx);
    getJsonPara( j2, "dt", Parameter.dt);
   
    getJsonPara( j2, "ParticleRadius", Parameter.ParticleRadius);
    getJsonPara( j2, "ParticleDensity", Parameter.ParticleDensity);
    getJsonPara( j2, "FluidDensity", Parameter.FluidDensity);
    getJsonPara( j2, "ArtificialFluidAccelerationZ", Parameter.ArtificialFluidAccelerationZ);
    getJsonPara( j2, "Tf", Parameter.Tf);
    getJsonPara( j2, "OutputStep", Parameter.OutputStep);
    getJsonPara( j2, "nlayer", Parameter.nlayer);
    getJsonPara( j2, "h", Parameter.h);
    getJsonPara( j2, "d", Parameter.d);
    getJsonPara( j2, "lambda", Parameter.lambda);
    getJsonPara( j2, "lambdashift", Parameter.lambdashift);
    getJsonPara( j2, "dlayer", Parameter.dlayer);
    getJsonPara( j2, "PlaceCenterZPercent", Parameter.ObstaclePlaceCenterZPercent);
    getJsonPara( j2, "SkipVisFactor", Parameter.SkipVisFactor);
}

std::vector<Vec3_t> generate_obstacle_positions_lattice2D(
                                                        size_t nlayer,  //number of layers
                                                        double domainsize,
                                                        double lambda,
                                                        double lambdashift,
                                                        double lambdalayer // distance between layer,
                                                        ) {

    int nPoints = domainsize/lambda;
    
    std::vector<Vec3_t> positions;
    // Create approximate lattice
    for (size_t k = 0; k < nlayer ; ++k) {                                                        
        for (size_t j = 0; j < nPoints ; ++j) {
            double x = j * lambda + k * lambdashift;
            x = std::fmod(x, domainsize);
            if (x< 0){
                x+=domainsize;
            } 
            double z = k * (lambdalayer);
            std::cout << "Preliminary: particle pos: " << x << " "  << z << std::endl;
            positions.push_back(Vec3_t(x, 0, z));
        }
    }

    //finding center of positions
    double normalize_sum_x = 0.0;
    double normalize_sum_y = 0.0;
    double normalize_sum_z = 0.0;
    for (const auto &p : positions) {
        normalize_sum_x += p(0);
        normalize_sum_y += p(1);
        normalize_sum_z += p(2);
    }
    double normalize_avg_x = normalize_sum_x / positions.size();
    double normalize_avg_y = normalize_sum_y / positions.size();
    double normalize_avg_z = normalize_sum_z / positions.size();
    double avg_x = normalize_avg_x;
    double avg_y = normalize_avg_y;
    double avg_z = normalize_avg_z;
    std::cout << "avg x y: " << avg_x << " " << avg_y << " " << avg_z <<std::endl;
    // shifting the positions of the lattice to have an average center {0,0}
    for (auto &p : positions) {
        p(0) = p(0) - avg_x;
        p(1) = p(1) - avg_y;
        p(2) = p(2) - avg_z;
    }

    for (auto &p : positions) {
       std::cout << "Final: particle pos: " << p(0) << " " << p(1) << " " << p(2) << std::endl;
    }
    return positions;
}

int main(int argc, char **argv) try
{
    if (argc < 2) {
        std::cerr << "Usage: ./program <jsonpath>\n";
        return 1;
    }
    std::string jsonpath = argv[1];  // read first parameter
    std::cout << "jsonpath: " << jsonpath << std::endl;
    LoadParameter(jsonpath);

    size_t Nproc = 0.75*omp_get_max_threads();
    double nu = Parameter.nu;
    size_t nx =Parameter.nx;
    size_t ny = Parameter.ny;
    size_t nz = Parameter.nz;
    double dx = Parameter.dx;
    double dt = Parameter.dt;
    double R  = Parameter.ParticleRadius;
    double rho  = Parameter.ParticleDensity;
    double h  = Parameter.h;
    double nlayer = Parameter.nlayer;
    double d  = Parameter.d;
    double lambda  = Parameter.lambda;
    double lambdashift  = Parameter.lambdashift;
    double dlayer  = Parameter.dlayer;
    double rhof  = Parameter.FluidDensity;
    double accz = Parameter.ArtificialFluidAccelerationZ;
    double lx = nx * dx;
    double ly = ny * dx;
    double lz = nz * dx;
    double Tf = Parameter.Tf;
    double OutputStep =Parameter.OutputStep;
    double PlaceCenterZPercent = Parameter.ObstaclePlaceCenterZPercent;
    
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    dom.LBMDOM.Step = Parameter.SkipVisFactor; //it will reduce the save files by averagin every 2 cells
   
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    
    double spongeSize = 2*dx;

   // create Obstacle
    std::vector<Vec3_t> pos;
    pos = generate_obstacle_positions_lattice2D( nlayer,lx,lambda,lambdashift, h+dlayer);
    std::vector<Vec3_t> new_pos;
    for (auto &p : pos) {
        new_pos.push_back(Vec3_t(0.5*lx + p(0), 0.5*ly+p(1),PlaceCenterZPercent*lz + p(2)));
    }
    pos = new_pos;    
     Vec3_t axis0(OrthoSys::e0);
    for (size_t k = 0;k<pos.size();k++) {
        int id = -1 - k;
        dom.DEMDOM.AddRecBox(id, pos[k], Vec3_t(lambda-d,ly-spongeSize*2,h), spongeSize,rho,0,&axis0);
        dom.DEMDOM.GetParticle(id)->FixVeloc();
        dom.DEMDOM.GetParticle(id)->FixFree  = true;
    }

    // create Particle
    dom.DEMDOM.AddSphere(1,Vec3_t(0.5*lx, 0.5*ly,0.9*lz),R,rho);

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,rhof,v);
    }   

    real3 acc = make_real3(0.0,0.0,accz);
    cudaMalloc(&dat.pacc, sizeof(real3));
    cudaMemcpy(dat.pacc, &acc, sizeof(real3), cudaMemcpyHostToDevice);

    dom.Alpha = 2.0*dx;
    dom.PeriodicX= true;
    dom.PeriodicY= true;
    dom.PeriodicZ= true;

    dom.Solve(Tf,Tf/OutputStep,Setup,Report,"tlbmdem_chromatography",true,Nproc);
}
MECHSYS_CATCH

