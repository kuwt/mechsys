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
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include "json.hpp"
#include "getjsonpara.hpp"

struct UserData
{
    Vec3_t                acc;
    double                 nu;
    double                  R;
};


void Setup (LBMDEM::Domain & dom, void * UD)
{
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
    int ParticleNumber = 10;
    double ParticleRadius  = 0.9;
    double ParticleDist  = 4*dx;
    double ParticleDistVariance  = 0.5*dx;
    double ParticleDensity  = 1.;
    double FluidDensity  = 1.;
    double ArtificialParticleAccelerationZ = -1.0e-4;
    double Tf = 1.0e5;
    int OutputStep = 200;

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
    getJsonPara( j2, "ParticleNumber", Parameter.ParticleNumber);
    getJsonPara( j2, "ParticleRadius", Parameter.ParticleRadius);
    getJsonPara( j2, "ParticleDist", Parameter.ParticleDist);
    getJsonPara( j2, "ParticleDistVariance", Parameter.ParticleDistVariance);
    getJsonPara( j2, "ParticleDensity", Parameter.ParticleDensity);
    getJsonPara( j2, "FluidDensity", Parameter.FluidDensity);
    getJsonPara( j2, "ArtificialParticleAccelerationZ", Parameter.ArtificialParticleAccelerationZ);
    getJsonPara( j2, "Tf", Parameter.Tf);
    getJsonPara( j2, "OutputStep", Parameter.OutputStep);
}


std::vector<Vec3_t> generate_particle_positions_squarelattice(
                                                        size_t n,  //number of particles
                                                        double r,  //radius
                                                        double d, // distance
                                                        double s, //sd
                                                        unsigned seed = 42) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> jitter(0.0, s);

    // Estimate box size to roughly fit n particles with avg spacing d
    size_t nx = static_cast<size_t>(sqrt(n));
    size_t ny = static_cast<size_t>(ceil(static_cast<double>(n) / nx));
    double L = nx * d;

    std::vector<Vec3_t> positions;

    // Create approximate lattice
    for (size_t i = 0; i < ny && positions.size() < n; ++i) {
        for (size_t j = 0; j < nx && positions.size() < n; ++j) {
            double x = j * d + jitter(gen);
            double y = i * d + jitter(gen);
            double z = jitter(gen);
            std::cout << "Preliminary: particle pos: " << x << " " << y << " " << z << std::endl;
            positions.push_back(Vec3_t(x, y, z));
        }
    }

    //finding center of positions
    double normalize_sum_x = 0.0, normalize_sum_y = 0.0;
    for (const auto &p : positions) {
        normalize_sum_x += p(0)/L;
        normalize_sum_y += p(1)/L;
    }
    double normalize_avg_x = normalize_sum_x / positions.size();
    double normalize_avg_y = normalize_sum_y / positions.size();
    double avg_x = normalize_avg_x * L;
    double avg_y = normalize_avg_y * L;
    std::cout << "avg x y: " << avg_x << " " << avg_y << std::endl;
    // shifting the positions of the lattice to have an average center {0,0}
    for (auto &p : positions) {
        p(0) = p(0) - avg_x;
        p(1) = p(1) - avg_y;
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
    double ParticleDist  = Parameter.ParticleDist;
    double rho  = Parameter.ParticleDensity;
    double rhof  = Parameter.FluidDensity;
    double accz = Parameter.ArtificialParticleAccelerationZ;

    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    dom.LBMDOM.Step = 2; //it will reduce the save files by averagin every 2 cells
    
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    dat.acc = Vec3_t(0.0,0.0,accz);

    // add fixed box with z opening
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    int BoxInitialTag = -1;
    dom.DEMDOM.AddPlane (BoxInitialTag,   Vec3_t(0.*dx*nx,0.5*dx*ny,0.5*dx*nz),R* 0.2,dx*nz,dx*ny,rho, M_PI/2.0, &axis1);
    dom.DEMDOM.GetParticle(BoxInitialTag)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-1, Vec3_t(1*dx*nx,0.5*dx*ny,0.5*dx*nz),R* 0.2,dx*nz,dx*ny,rho, 3.0*M_PI/2.0, &axis1);
    dom.DEMDOM.GetParticle(BoxInitialTag-1)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-2, Vec3_t(0.5*dx*nx,0.*dx*ny,0.5*dx*nz),R* 0.2,dx*nx,dx*nz,rho, 3.0*M_PI/2.0, &axis0);
    dom.DEMDOM.GetParticle(BoxInitialTag-2)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-3, Vec3_t(0.5*dx*nx,1*dx*ny,0.5*dx*nz),R* 0.2,dx*nx,dx*nz,rho, M_PI/2.0, &axis0);
    dom.DEMDOM.GetParticle(BoxInitialTag-3)->FixVeloc();

    // add particles
    std::vector<Vec3_t> pos;
    pos = generate_particle_positions_squarelattice( Parameter.ParticleNumber,
                                               Parameter.ParticleRadius,
                                                Parameter.ParticleDist + 2 *Parameter.ParticleRadius ,
                                                Parameter.ParticleDistVariance  );

    for (size_t k = 0;k<Parameter.ParticleNumber;k++) {
        int id = k;
        dom.DEMDOM.AddSphere(id,Vec3_t(0.5*dx*nx + pos[k](0),0.5*dx*ny+pos[k](1),0.9*dx*nz + pos[k](2)),R,rho);
        dom.DEMDOM.GetParticle(id)->Ff = dom.DEMDOM.GetParticle(id)->Props.m*dat.acc;
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,rhof,v);
        if ((ix==0)||(ix==nx-1)||(iy==0)||(iy==ny-1)) dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
    }   

    dom.Alpha = 2.0*dx;
    //dom.PeriodicX= true;
    //dom.PeriodicY= true;
    dom.PeriodicZ= true;
    double Tf = Parameter.Tf;
    dom.Solve(Tf,Tf/Parameter.OutputStep,Setup,Report,"tlbmdem_cu_10",true,Nproc);
}
MECHSYS_CATCH

