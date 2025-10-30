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

struct Point {
    float x, y, z;
};
float distanceSq(const Point& p1, const Point& p2) {
    float dx = p1.x - p2.x;
    float dy = p1.y - p2.y;
    float dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}
/**
 * @brief Generates a Poisson disk sample set within a 3D box domain using Bridson's algorithm.
 *
 * @param Lx Domain size along the X-axis.
 * @param Ly Domain size along the Y-axis.
 * @param Lz Domain size along the Z-axis.
 * @param R Minimum distance between any two points (radius).
 * @param K Max attempts to generate a new point from an active point (e.g., 30 is common).
 * @return std::vector<Point> The generated Poisson disk sample points.
 */
std::vector<Point> poissonDiskSampling3D(float Lx, float Ly, float Lz, float R, int K = 30) {
    // --- Initialization and Setup ---
    const float R_SQR = R * R;
    
    // Cell size s = R / sqrt(dimension) to ensure at most one point per cell (3D: sqrt(3))
    const float CELL_SIZE = R / std::sqrt(3.0f); 
    
    // Calculate grid dimensions
    const int GRID_WIDTH_X = static_cast<int>(std::ceil(Lx / CELL_SIZE));
    const int GRID_WIDTH_Y = static_cast<int>(std::ceil(Ly / CELL_SIZE));
    const int GRID_WIDTH_Z = static_cast<int>(std::ceil(Lz / CELL_SIZE));
    const int GRID_SIZE = GRID_WIDTH_X * GRID_WIDTH_Y * GRID_WIDTH_Z;

    // Random number generation setup
    std::random_device rd;
    std::mt19937 gen(rd());
    // For generating coordinates and distances
    std::uniform_real_distribution<float> distrib_01(0.0f, 1.0f); 

    std::vector<Point> samples;
    // Grid stores 1-based index into 'samples' (0 means empty) for fast lookups
    std::vector<int> grid(GRID_SIZE, 0); 
    
    // Active list stores 0-based indices of points in 'samples' to generate new candidates from
    std::list<int> active_list; 

    // Lambda to calculate grid index from coordinates
    auto getGridIndex = [&](const Point& p) -> int {
        int x = static_cast<int>(p.x / CELL_SIZE);
        int y = static_cast<int>(p.y / CELL_SIZE);
        int z = static_cast<int>(p.z / CELL_SIZE);
        
        // Check if the point is outside the bounds of the grid
        if (x < 0 || x >= GRID_WIDTH_X || 
            y < 0 || y >= GRID_WIDTH_Y || 
            z < 0 || z >= GRID_WIDTH_Z) {
            return -1;
        }
        
        // Flatten 3D index to 1D
        return x + y * GRID_WIDTH_X + z * GRID_WIDTH_X * GRID_WIDTH_Y;
    };
    
    // --- Step 1: Place Initial Sample ---
    
    // Choose a random point in the domain for the first sample
    Point p0 = {
        distrib_01(gen) * Lx, 
        distrib_01(gen) * Ly, 
        distrib_01(gen) * Lz
    };
    
    samples.push_back(p0);
    active_list.push_back(0); // 0-based index
    grid[getGridIndex(p0)] = 1; // 1-based index (to distinguish from 0/empty)

    // --- Step 2: Sample Generation Loop ---
    
    while (!active_list.empty()) {
        // Pick the first point 'p_i' from the active list
        int p_i_index = active_list.front();
        active_list.pop_front();
        
        const Point& p_i = samples[p_i_index];
        bool found_new_sample = false;

        // Generate K candidate points
        for (int i = 0; i < K; ++i) {
            // Generate a candidate point 'p_k' in the annulus [R, 2R] around p_i
            
            // 1. Random distance d in [R, 2R]
            float d = R + R * distrib_01(gen); 
            
            // 2. Generate random spherical coordinates for uniform sampling on a sphere
            float phi = 2.0f * M_PI * distrib_01(gen); // Azimuthal angle [0, 2*PI]
            // Polar angle calculation for uniform distribution over sphere surface area
            float theta = std::acos(1.0f - 2.0f * distrib_01(gen)); 

            Point p_k;
            // 3. Convert spherical to Cartesian coordinates and offset by p_i
            p_k.x = p_i.x + d * std::cos(phi) * std::sin(theta);
            p_k.y = p_i.y + d * std::sin(phi) * std::sin(theta);
            p_k.z = p_i.z + d * std::cos(theta);

            // Check if p_k is inside the domain [0, L]
            if (p_k.x < 0 || p_k.x >= Lx || 
                p_k.y < 0 || p_k.y >= Ly || 
                p_k.z < 0 || p_k.z >= Lz) {
                continue;
            }

            // Check for minimum distance constraint
            int grid_idx = getGridIndex(p_k);
            if (grid_idx == -1) continue; 

            bool is_valid = true;
            int x_cell = static_cast<int>(p_k.x / CELL_SIZE);
            int y_cell = static_cast<int>(p_k.y / CELL_SIZE);
            int z_cell = static_cast<int>(p_k.z / CELL_SIZE);

            // Check neighbors in a local 5x5x5 box of cells
            // A 3x3x3 neighborhood is generally sufficient for the R to 2R sampling property,
            // but 5x5x5 ensures all points within R distance are considered safely.
            for (int z = std::max(0, z_cell - 2); z <= std::min(GRID_WIDTH_Z - 1, z_cell + 2); ++z) {
                for (int y = std::max(0, y_cell - 2); y <= std::min(GRID_WIDTH_Y - 1, y_cell + 2); ++y) {
                    for (int x = std::max(0, x_cell - 2); x <= std::min(GRID_WIDTH_X - 1, x_cell + 2); ++x) {
                        
                        int neighbor_idx = x + y * GRID_WIDTH_X + z * GRID_WIDTH_X * GRID_WIDTH_Y;
                        int point_index_in_samples = grid[neighbor_idx];

                        if (point_index_in_samples > 0) { // Found an existing point
                            const Point& p_neighbor = samples[point_index_in_samples - 1]; // Convert back to 0-based index

                            if (distanceSq(p_k, p_neighbor) < R_SQR) {
                                is_valid = false;
                                goto check_end; // Exit all three loops
                            }
                        }
                    }
                }
            }
            
            check_end:; // Label for breaking out of multiple loops

            // 4. If valid, add it to the samples, grid, and active list
            if (is_valid) {
                samples.push_back(p_k);
                int new_index = samples.size() - 1; // 0-based index
                
                // The parent point p_i is pushed back to the active list
                active_list.push_back(p_i_index);
                // The new point p_k is added to the active list
                active_list.push_back(new_index);
                
                grid[grid_idx] = new_index + 1; // Store 1-based index
                found_new_sample = true;
                break; // Stop generating candidates for this p_i and continue the main loop
            }
        }

        // If no new sample was found after K attempts, the point p_i is removed
        // (it was removed by active_list.pop_front() earlier and not re-added).
    }

    return samples;
}


std::vector<Vec3_t> generate_particle_positions_poissondisk3d(
                                                        double lx, 
                                                        double ly, 
                                                        double lz, 
                                                        double d) {
    std::vector<Point> points = poissonDiskSampling3D(lx,ly,lz,d);
    
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
    
    // add fixed box with z opening
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    int BoxInitialTag = -1;
    dom.DEMDOM.AddPlane (BoxInitialTag,   Vec3_t(0.*lx,0.5*ly,0.5*lz),R* 0.2,lz,ly,rho, M_PI/2.0, &axis1);
    dom.DEMDOM.GetParticle(BoxInitialTag)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-1, Vec3_t(1*lx,0.5*ly,0.5*lz),R* 0.2,lz,ly,rho, 3.0*M_PI/2.0, &axis1);
    dom.DEMDOM.GetParticle(BoxInitialTag-1)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-2, Vec3_t(0.5*lx,0.*ly,0.5*lz),R* 0.2,lx,lz,rho, 3.0*M_PI/2.0, &axis0);
    dom.DEMDOM.GetParticle(BoxInitialTag-2)->FixVeloc();
    dom.DEMDOM.AddPlane (BoxInitialTag-3, Vec3_t(0.5*lx,1*ly,0.5*lz),R* 0.2,lx,lz,rho, M_PI/2.0, &axis0);
    dom.DEMDOM.GetParticle(BoxInitialTag-3)->FixVeloc();
    if (!Parameter.isPeriodic){
        dom.DEMDOM.AddPlane (BoxInitialTag-4, Vec3_t(0.5*lx,0.5*ly,0), R* 0.2, lx, ly, rho, M_PI, &axis0);
        dom.DEMDOM.GetParticle(BoxInitialTag-4)->FixVeloc();
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
            Parameter.ParticlePlaceRandomBoxLXpercent *lx,
            Parameter.ParticlePlaceRandomBoxLYpercent *ly,
            Parameter.ParticlePlaceRandomBoxLZpercent *lz,
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
        if (p(0) > 0+R && p(0) < lx-R && p(1) > 0+R && p(1) < ly-R && p(2) > 0+2*R && p(2) < lz-2*R ){
            new_pos.push_back(Vec3_t(p(0),p(1),p(2)));
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
        if ((ix==0)||(ix==nx-1)||(iy==0)||(iy==ny-1)) {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
        if (!Parameter.isPeriodic){
            dom.LBMDOM.IsSolid[0][ix][iy][0] = true;
        }
    }   

    dom.Alpha = 2.0*dx;
    if (Parameter.isPeriodic){
        dom.PeriodicZ= true;
    }
    double Tf = Parameter.Tf;
    dom.Solve(Tf,Tf/Parameter.OutputStep,Setup,Report,"tlbmdemParticleSettling",true,Nproc);
}
MECHSYS_CATCH

