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
    int ParticlePlacementMethod = 0;
    double ParticlePlaceCenterZPercent  = 0.9;
    int ParticleLatticeLayerNumber = 1;
    int ParticleNumber = 10;
    double ParticlePlaceRandomBoxLXpercent= 0.5;
    double ParticlePlaceRandomBoxLYpercent= 0.5;
    double ParticlePlaceRandomBoxLZpercent= 0.5;

    double ParticleRadius  = 0.9;
    double ParticleDist  = 4*dx;
    double ParticleDistVariance  = 0.5*dx;
    double ParticleHeightVariance = ParticleDistVariance;
    double ParticleDensity  = 1.;
    double FluidDensity  = 1.;
    double ArtificialParticleAccelerationZ = -1.0e-4;
    bool isPeriodic = false;
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
    getJsonPara( j2, "ParticlePlacementMethod", Parameter.ParticlePlacementMethod);
    getJsonPara( j2, "ParticlePlaceCenterZPercent", Parameter.ParticlePlaceCenterZPercent);
    getJsonPara( j2, "ParticleLatticeLayerNumber", Parameter.ParticleLatticeLayerNumber);
    getJsonPara( j2, "ParticleNumber", Parameter.ParticleNumber);
    getJsonPara( j2, "ParticlePlaceRandomBoxLXpercent", Parameter.ParticlePlaceRandomBoxLXpercent);
    getJsonPara( j2, "ParticlePlaceRandomBoxLYpercent", Parameter.ParticlePlaceRandomBoxLYpercent);
    getJsonPara( j2, "ParticlePlaceRandomBoxLZpercent", Parameter.ParticlePlaceRandomBoxLZpercent);
    getJsonPara( j2, "ParticleRadius", Parameter.ParticleRadius);
    getJsonPara( j2, "ParticleDist", Parameter.ParticleDist);
    getJsonPara( j2, "ParticleDistVariance", Parameter.ParticleDistVariance);
    getJsonPara( j2, "ParticleHeightVariance", Parameter.ParticleHeightVariance);
    getJsonPara( j2, "ParticleDensity", Parameter.ParticleDensity);
    getJsonPara( j2, "FluidDensity", Parameter.FluidDensity);
    getJsonPara( j2, "ArtificialParticleAccelerationZ", Parameter.ArtificialParticleAccelerationZ);
    getJsonPara( j2, "isPeriodic", Parameter.isPeriodic);
    getJsonPara( j2, "Tf", Parameter.Tf);
    getJsonPara( j2, "OutputStep", Parameter.OutputStep);
}

std::vector<Vec3_t> generate_particle_positions_lattice(
                                                        size_t nlayer,  //number of layers
                                                        size_t n,  //number of particles per layer
                                                        double d, // distance
                                                        double sxy, //sd
                                                        double sz, //sd
                                                        unsigned seed = 42) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> jitter(0.0, sxy);
    std::normal_distribution<double> jitterz(0.0, sz);
    // Estimate box size to roughly fit n particles with avg spacing d
    size_t nx = static_cast<size_t>(sqrt(n));
    size_t ny = static_cast<size_t>(ceil(static_cast<double>(n) / nx));
    double L = nx * d;

    std::vector<Vec3_t> positions;

    // Create approximate lattice
    for (size_t k = 0; k < nlayer && positions.size() < n*nlayer; ++k) {                                                        
        for (size_t i = 0; i < ny && positions.size() < n*nlayer; ++i) {
            for (size_t j = 0; j < nx && positions.size() < n*nlayer; ++j) {
                double x = j * d + jitter(gen);
                double y = i * d + jitter(gen);
                double z = k * d + jitterz(gen);
                std::cout << "Preliminary: particle pos: " << x << " " << y << " " << z << std::endl;
                positions.push_back(Vec3_t(x, y, z));
            }
        }
    }

    //finding center of positions
    double normalize_sum_x = 0.0;
    double normalize_sum_y = 0.0;
    double normalize_sum_z = 0.0;
    for (const auto &p : positions) {
        normalize_sum_x += p(0)/L;
        normalize_sum_y += p(1)/L;
        normalize_sum_z += p(2)/L;
    }
    double normalize_avg_x = normalize_sum_x / positions.size();
    double normalize_avg_y = normalize_sum_y / positions.size();
    double normalize_avg_z = normalize_sum_z / positions.size();
    double avg_x = normalize_avg_x * L;
    double avg_y = normalize_avg_y * L;
    double avg_z = normalize_avg_z * L;
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

struct Point3D {
    float x, y, z;
};
float distanceSq(const Point3D& p1, const Point3D& p2) {
    float dx = p1.x - p2.x;
    float dy = p1.y - p2.y;
    float dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}

double periodicSquaredDistance(const Point3D& p1, const Point3D& p2, double Lx, double Ly, double Lz) {
    auto dist_sq = [&](double d, double L) {
        // Minimum Image Convention
        d = std::fmod(d, L);
        if (d > L / 2.0) d -= L;
        if (d < -L / 2.0) d += L;
        return d * d;
    };

    double dx_sq = dist_sq(p1.x - p2.x, Lx);
    double dy_sq = dist_sq(p1.y - p2.y, Ly);
    double dz_sq = dist_sq(p1.z - p2.z, Lz);

    return dx_sq + dy_sq + dz_sq;
}

/***********
* To generate points that are truly periodic respecting the boundary, L should be multiple of r
*/
std::vector<Point3D> poissonDiskSampling3D_Periodic(double Lx, double Ly, double Lz, double r, int k=30) {
    // 1. Grid setup
    // Cell size S: remains isotropic (r / sqrt(3)) to ensure one disk per cell maximum.
    double S = r / std::sqrt(3.0); 
    
    // Grid dimensions are now independent
    int Nx = static_cast<int>(std::ceil(Lx / S));
    int Ny = static_cast<int>(std::ceil(Ly / S));
    int Nz = static_cast<int>(std::ceil(Lz / S));
    int N_total = Nx * Ny * Nz;
    
    // Grid stores the index of the point in the 'samples' vector, or -1 if empty.
    std::vector<int> grid(N_total, -1); 

    auto get_grid_index = [&](const Point3D& p) -> int {
        int gx = static_cast<int>(std::floor(p.x / S));
        int gy = static_cast<int>(std::floor(p.y / S));
        int gz = static_cast<int>(std::floor(p.z / S));

        // Use lambda for safe periodic wrapping, passing the correct dimension N
        auto safe_mod = [](int i, int N_dim) -> int {
            return (i % N_dim + N_dim) % N_dim;
        };

        gx = safe_mod(gx, Nx);
        gy = safe_mod(gy, Ny);
        gz = safe_mod(gz, Nz);
        
        // 1D indexing: k * (Ny * Nx) + j * Nx + i
        return gz * Ny * Nx + gy * Nx + gx;
    };
    
    auto insert_point = [&](const Point3D& p, int index) {
        grid[get_grid_index(p)] = index;
    };
    
    // 2. Initialization
    std::vector<Point3D> samples;
    std::vector<int> active_list;
    double r_sq = r * r;

    // C++11 Random setup
    std::mt19937 gen(static_cast<unsigned int>(std::time(0))); 
    std::uniform_real_distribution<> rand_x(0.0, Lx);
    std::uniform_real_distribution<> rand_y(0.0, Ly);
    std::uniform_real_distribution<> rand_z(0.0, Lz);
    
    // Initial random point
    Point3D p0 = {rand_x(gen), rand_y(gen), rand_z(gen)};
    samples.push_back(p0);
    active_list.push_back(0);
    insert_point(p0, 0);

    // Random angle generation setup
    std::uniform_real_distribution<> rand_dist_r(r, 2.0 * r);
    std::uniform_real_distribution<> rand_angle_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<> rand_z_polar(-1.0, 1.0); 

    // 3. Main loop
    while (!active_list.empty()) {
        std::uniform_int_distribution<> rand_active(0, active_list.size() - 1);
        int active_idx_in_list = rand_active(gen);
        int p_idx = active_list[active_idx_in_list];
        const Point3D& p = samples[p_idx];

        bool found_new_point = false;
        
        for (int i = 0; i < k; ++i) {
            // Generate candidate 'c' in the spherical annulus [r, 2r]
            double dist = rand_dist_r(gen);
            double phi = rand_angle_phi(gen);
            double cos_theta = rand_z_polar(gen); 
            double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
            
            // Cartesian coordinates (relative to p)
            double dx = dist * sin_theta * std::cos(phi);
            double dy = dist * sin_theta * std::sin(phi);
            double dz = dist * cos_theta;
            
            // New candidate point, applied with correct periodic wrapping (Lx, Ly, Lz)
            Point3D c = { 
                std::fmod(p.x + dx + Lx, Lx), 
                std::fmod(p.y + dy + Ly, Ly), 
                std::fmod(p.z + dz + Lz, Lz) 
            };

            // Check if the candidate 'c' is valid
            bool is_valid = true;
            
            // Determine the grid coordinates for the search (requires unwrapped floor)
            int gx = static_cast<int>(std::floor(c.x / S));
            int gy = static_cast<int>(std::floor(c.y / S));
            int gz = static_cast<int>(std::floor(c.z / S));
            
            int search_radius = 2; // Safe check radius

            // Nested search loops
            for (int sz = -search_radius; sz <= search_radius; ++sz) {
                for (int sy = -search_radius; sy <= search_radius; ++sy) {
                    for (int sx = -search_radius; sx <= search_radius; ++sx) {
                        
                        // Periodic wrapping for neighbor indices using Nx, Ny, Nz
                        auto safe_mod = [](int i, int N_dim) -> int {
                            return (i % N_dim + N_dim) % N_dim;
                        };
                        
                        int nx = safe_mod(gx + sx, Nx);
                        int ny = safe_mod(gy + sy, Ny);
                        int nz = safe_mod(gz + sz, Nz);
                        
                        int neighbor_grid_idx = nz * Ny * Nx + ny * Nx + nx;
                        int neighbor_point_idx = grid[neighbor_grid_idx];

                        if (neighbor_point_idx != -1) {
                            // Check distance using non-cubic periodic metric
                            if (periodicSquaredDistance(c, samples[neighbor_point_idx], Lx, Ly, Lz) < r_sq) {
                                is_valid = false;
                                goto next_candidate;
                            }
                        }
                    }
                }
            }
            
            // If the inner loops complete without breaking
            if (is_valid) {
                samples.push_back(c);
                int new_idx = samples.size() - 1;
                active_list.push_back(new_idx);
                insert_point(c, new_idx);
                found_new_point = true;
                break;
            }

            next_candidate:;
        }

        if (!found_new_point) {
            active_list.erase(active_list.begin() + active_idx_in_list);
        }
    }

    return samples;
}

std::vector<Vec3_t> generate_particle_positions_poissondisk3d(
                                                        double lx, 
                                                        double ly, 
                                                        double lz, 
                                                        double d) {
    std::vector<Point3D> points = poissonDiskSampling3D_Periodic(lx,ly,lz,d,true);
    
    std::vector<Vec3_t> positions;
    for (auto &p : points) {
        positions.push_back(Vec3_t(p.x,p.y,p.z));
    }

    //finding center of positions
    double normalize_sum_x = 0.0;
    double normalize_sum_y = 0.0;
    double normalize_sum_z = 0.0;
    for (const auto &p : positions) {
        normalize_sum_x += p(0)/lx;
        normalize_sum_y += p(1)/ly;
        normalize_sum_z += p(2)/lz;
    }
    double normalize_avg_x = normalize_sum_x / positions.size();
    double normalize_avg_y = normalize_sum_y / positions.size();
    double normalize_avg_z = normalize_sum_z / positions.size();
    double avg_x = normalize_avg_x * lx;
    double avg_y = normalize_avg_y * ly;
    double avg_z = normalize_avg_z * lz;
    std::cout << "avg x y z: " << avg_x << " " << avg_y << " " << avg_z <<std::endl;
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
    double ParticleDist  = Parameter.ParticleDist;
    double rho  = Parameter.ParticleDensity;
    double rhof  = Parameter.FluidDensity;
    double accz = Parameter.ArtificialParticleAccelerationZ;
    double lx = nx * dx;
    double ly = ny * dx;
    double lz = nz * dx;

    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    dom.LBMDOM.Step = 2; //it will reduce the save files by averagin every 2 cells
   
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    dat.acc = Vec3_t(0.0,0.0,accz);
    
    // add bottom plane
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    if (!Parameter.isPeriodic){
        dom.DEMDOM.AddPlane (-1, Vec3_t(0.5*lx,0.5*ly,0), R* 0.2, lx,ly, rho, M_PI, &axis0);
        dom.DEMDOM.GetParticle(-1)->FixVeloc();
    }

   // create particle position
    std::vector<Vec3_t> pos;
    if (Parameter.ParticlePlacementMethod == 0) {
        pos = generate_particle_positions_lattice( Parameter.ParticleLatticeLayerNumber,
                                            Parameter.ParticleNumber,
                                            Parameter.ParticleDist + 2 *Parameter.ParticleRadius ,
                                            Parameter.ParticleDistVariance,
                                            Parameter.ParticleHeightVariance);
    } else {
        pos = generate_particle_positions_poissondisk3d(
            Parameter.ParticlePlaceRandomBoxLXpercent * lx,
            Parameter.ParticlePlaceRandomBoxLYpercent * ly,
            Parameter.ParticlePlaceRandomBoxLZpercent * lz,
            Parameter.ParticleDist + 2 *Parameter.ParticleRadius);
    }
    // shift particle cloud center
    std::vector<Vec3_t> new_pos;
    for (auto &p : pos) {
        new_pos.push_back(Vec3_t(0.5*lx + p(0), 0.5*ly+p(1),Parameter.ParticlePlaceCenterZPercent*lz + p(2)));
    }
    pos = new_pos;    

    // elimate particles outside the fluid domain
    new_pos.clear();
    for (auto &p : pos) {
        if (p(0) > 0 && p(0) < lx && p(1) > 0 && p(1) < ly ){
            if ( Parameter.isPeriodic) {
                if (p(2) > 0 && p(2) < lz-2*R) {
                    new_pos.push_back(Vec3_t(p(0),p(1),p(2)));
                }
            } else {
                if (p(2) > 0+2*R && p(2) < lz-2*R) {
                    new_pos.push_back(Vec3_t(p(0),p(1),p(2)));
                }
            }
        }
    }
    pos = new_pos;

     // add particles
    for (size_t k = 0;k<pos.size();k++) {
        int id = k;
        dom.DEMDOM.AddSphere(id,pos[k],R,rho);
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
        if (!Parameter.isPeriodic){
            dom.LBMDOM.IsSolid[0][ix][iy][0] = true;
        }
    }   

    dom.Alpha = 2.0*dx;
    dom.PeriodicX= true;
    dom.PeriodicY= true;
    if (Parameter.isPeriodic){
        dom.PeriodicZ= true;
    }
    double Tf = Parameter.Tf;
    dom.Solve(Tf,Tf/Parameter.OutputStep,Setup,Report,"tlbmdemParticleSettling",true,Nproc);
}
MECHSYS_CATCH

