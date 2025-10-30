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
    double ParticleRadius  = 0.9;
    double ParticleDist  = 4*dx;
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
    getJsonPara( j2, "ParticleRadius", Parameter.ParticleRadius);
    getJsonPara( j2, "ParticleDist", Parameter.ParticleDist);
    getJsonPara( j2, "ParticleDensity", Parameter.ParticleDensity);
    getJsonPara( j2, "FluidDensity", Parameter.FluidDensity);
    getJsonPara( j2, "ArtificialParticleAccelerationZ", Parameter.ArtificialParticleAccelerationZ);
    getJsonPara( j2, "Tf", Parameter.Tf);
    getJsonPara( j2, "OutputStep", Parameter.OutputStep);
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
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    dat.acc = Vec3_t(0.0,0.0,accz);
    int id = 1;
    dom.DEMDOM.AddSphere(id,Vec3_t(0.5*dx*nx + R + 0.5*ParticleDist,0.5*dx*ny,0.5*dx*nz),R,rho);
    dom.DEMDOM.GetParticle(id)->Ff = dom.DEMDOM.GetParticle(id)->Props.m*dat.acc;
    id = 2;
    dom.DEMDOM.AddSphere(id,Vec3_t(0.5*dx*nx - ( R + 0.5*ParticleDist),0.5*dx*ny,0.5*dx*nz),R,rho);
    dom.DEMDOM.GetParticle(id)->Ff = dom.DEMDOM.GetParticle(id)->Props.m*dat.acc;

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
    dom.Solve(Tf,Tf/Parameter.OutputStep,Setup,Report,"tlbmdem_cu_07",true,Nproc);
}
MECHSYS_CATCH

