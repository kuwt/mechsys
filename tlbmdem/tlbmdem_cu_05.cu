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

int main(int argc, char **argv) try
{
    size_t Nproc = 0.75*omp_get_max_threads();

        
    double nu = 1.0;
    size_t nx = 101;
    size_t ny = 101;
    size_t nz = 101;
    double dx = 0.2;
    double dt = 1.6e-2;
    double R  = 0.9;
    //double R  = 0.3;
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    double rho = 1.0;
    dat.acc = Vec3_t(0.0,0.0,-1.0e-2);
    int id = 1;
    dom.DEMDOM.AddSphere(id,Vec3_t(0.5*dx*nx + R + dx,0.5*dx*ny,0.5*dx*nz),R,rho);
    dom.DEMDOM.GetParticle(id)->Ff = dom.DEMDOM.GetParticle(id)->Props.m*dat.acc;
    id = 2;
    dom.DEMDOM.AddSphere(id,Vec3_t(0.5*dx*nx - ( R + dx),0.5*dx*ny,0.5*dx*nz),R,rho);
    dom.DEMDOM.GetParticle(id)->Ff = dom.DEMDOM.GetParticle(id)->Props.m*dat.acc;
    id = -1;
    dom.DEMDOM.AddPlane(id,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.2*dx*nz),R* 0.5,0.5*dx*nx,0.5*dx*ny,rho, 0,0 );
    dom.DEMDOM.GetParticle(id)->FixVeloc();

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1.0/*rho*/,v);
    }   

    dom.Alpha = 2.0*dx;

    double Tf = 1.0e4;
    dom.Solve(Tf,Tf/200,Setup,Report,"tlbmdem_cu_05",true,Nproc);
}
MECHSYS_CATCH

