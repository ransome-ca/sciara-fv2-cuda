/*                                                                      
 * GPL3 License                                                                         
 *                                                                      
 *                                                                      
 * Copyright (c) 2022 Ransome CA                              
 *                                                                      
 * This file is part of SCIARA-fv2-CUDA.                                          
 *                                                                      
 * SCIARA-fv2-CUDA is free software: you can redistribute it and/or modify        
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or    
 * (at your option) any later version.                                  
 *                                                                      
 * SCIARA-fv2-CUDA is distributed in the hope that it will be useful,             
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        
 * GNU General Public License for more details.                         
 *                                                                      
 * You should have received a copy of the GNU General Public License    
 * along with SCIARA-fv2-CUDA.  If not, see <http://www.gnu.org/licenses/>.       
 */       

#include "cal2DBuffer.cuh"
#include "configurationPathLib.cuh"
#include "GISInfo.cuh"
#include "io.cuh"
#include "vent.cuh"
#include "Sciara.cuh"
#include "util.cu"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <new>
#include <cassert>



#define VENTS_COUNT             2


#define GET(M, c, i, j)         ((M)[(c) * (i) + (j)])
#define SET(M, c, i, j, val)    ((M)[(c) * (i) + (j)] = (val))

#define BUF_SET(M, r, c, n, i, j, val)    \
    ((M)[((n) * (r) * (c)) + ((i) * (c)) + (j)] = (val))

#define BUF_GET(M, r, c, n, i, j)         \
    ((M)[((n) * (r) * (c)) + ((i) * (c)) + (j)])


#if defined(PREC_DIV)
#   define DIV(a, b)            ((a) / (b))
#else
#   define DIV(a, b)            ((a) * __drcp_rn(b))
#endif





__constant__ int Xi[MOORE_NEIGHBORS];
__constant__ int Xj[MOORE_NEIGHBORS];
__constant__ TVent vents[VENTS_COUNT];



__global__
void emit_lava (
    size_t r, 
    size_t c, 
    double elapsed_time, 
    double Pclock, 
    double emission_time,
    double Pac,
    double PTvent,
    double* Sh,
    double* Sh_next,
    double* ST_next,
    float* total_emitted_lava
) {

    constexpr const size_t vents_size = sizeof(vents) / sizeof(TVent);


    const size_t k = blockIdx.x * blockDim.x + threadIdx.x;

    if(k < vents_size) {

        const size_t i = vents[k].y();
        const size_t j = vents[k].x();

        double v = vents[k].thickness(elapsed_time, Pclock, emission_time, Pac);

        SET(Sh_next, c, i, j, GET(Sh, c, i, j) + v);
        SET(ST_next, c, i, j, PTvent);

        if(v != 0.0) {
            atomicAdd(total_emitted_lava, v);
        }

    }

}


__global__
void compute_outflows(
    size_t r,
    size_t c,
    double* Sz,
    double* Sh,
    double* ST,
    double* Mf,
    double Pc,
    double _a,
    double _b,
    double _c,
    double _d
) {

    bool eliminated[MOORE_NEIGHBORS] = {};
    double z[MOORE_NEIGHBORS] = {};
    double h[MOORE_NEIGHBORS] = {};
    double H[MOORE_NEIGHBORS] = {};
    double theta[MOORE_NEIGHBORS] = {};
    
    double sz0;
    double sz;
    double T;
    double avg;
    double rr;
    double hc;

    bool loop;
    size_t counter;

    
    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i > 0 && j > 0 && i < r - 1 && j < c - 1) {

        if(GET(Sh, c, i, j) <= 0)
            return;

        T  = GET(ST, c, i, j);
        rr = pow(10, _a + _b * T);
        hc = pow(10, _c + _d * T);

        #pragma unroll
        for(size_t k = 0; k < MOORE_NEIGHBORS; k++) {

            sz0   = GET(Sz, c, i, j);
            sz    = GET(Sz, c, i + Xi[k], j + Xj[k]);
            h[k]  = GET(Sh, c, i + Xi[k], j + Xj[k]);


            if(k < VON_NEUMANN_NEIGHBORS) {
                z[k] = sz;
            } else {
                z[k] = sz0 - DIV((sz0 - sz), sqrt(2.0));
            }

        }


        H[0] = z[0];
        theta[0] = 0;
        eliminated[0] = false;


        #pragma unroll
        for(size_t k = 1; k < MOORE_NEIGHBORS; k++) {

            if(z[0] + h[0] > z[k] + h[k]) {

                H[k]          = z[k] + h[k];
                theta[k]      = atan(DIV(((z[0] + h[0]) - (z[k] + h[k])), Pc));
                eliminated[k] = false; 
            
            } else {
                
                eliminated[k] = true;

            }

        }



        do {

            loop    = false;
            avg     = h[0];
            counter = 0;


            #pragma unroll
            for(size_t k = 0; k < MOORE_NEIGHBORS; k++) {

                if(!eliminated[k]) {

                    avg += H[k];
                    counter++;

                }

            }

            if(counter != 0) {

                avg = DIV(avg, double(counter));

            }


            #pragma unroll
            for(size_t k = 0; k < MOORE_NEIGHBORS; k++) {

                if(!eliminated[k] && avg <= H[k]) {

                    eliminated[k] = true;
                    loop = true;

                }

            }

        } while(loop);


        #pragma unroll
        for(size_t k = 1; k < MOORE_NEIGHBORS; k++) {

            if(!eliminated[k] && h[0] > hc * cos(theta[k])) {

                BUF_SET(Mf, r, c, k - 1, i, j, rr * (avg - H[k]));

            } else {

                BUF_SET(Mf, r, c, k - 1, i, j, 0.0);

            }

        }

    }

}



__global__
void mass_balance (
    size_t r,
    size_t c,
    double* Sh,
    double* Sh_next,
    double* ST,
    double* ST_next,
    double* Mf
) {

    const uint8_t inflows_indices[NUMBER_OF_OUTFLOWS] = { 3, 2, 1, 0, 6, 7, 4, 5 };


    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;


    if(i > 0 && j > 0 && i < r - 1 && j < c - 1) {

        double inflow;
        double outflow;
        double t_neigh;

        double h_initial = GET(Sh, c, i, j);
        double t_initial = GET(ST, c, i, j);

        double h_next = h_initial;
        double t_next = h_initial * t_initial;


        #pragma unroll
        for(size_t n = 1; n < MOORE_NEIGHBORS; n++) {

            t_neigh = GET(ST, c, i + Xi[n], j + Xj[n]);

            inflow  = BUF_GET(Mf, r, c, inflows_indices[n - 1], i + Xi[n], j + Xj[n]);
            outflow = BUF_GET(Mf, r, c, n - 1, i, j);

            h_next += inflow - outflow;
            t_next += inflow * t_neigh - outflow * t_initial;

        }

        if(h_next > 0.0) {

            t_next = DIV(t_next, h_next);

            SET(ST_next, c, i, j, t_next);
            SET(Sh_next, c, i, j, h_next);

        }


    }

}


__global__
void compute_new_temperature_and_solidification (
    size_t r,
    size_t c,
    double Pepsilon,
    double Psigma,
    double Pclock,
    double Pcool,
    double Prho,
    double Pcv,
    double Pac,
    double PTsol,
    double *Sz,
    double *Sz_next,
    double *Sh, 
    double *Sh_next, 
    double *ST,
    double *ST_next,
    double *Mf,
    double *Mhs,
    bool   *Mb
) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i < r && j < c) {

        double nT;
        double aus;

        double z = GET(Sz, c, i, j);
        double h = GET(Sh, c, i, j);
        double T = GET(ST, c, i, j);

        if(h > 0 && GET(Mb, c, i, j) == false) {

            aus = 1.0 + DIV((3 * pow(T, 3.0) * Pepsilon * Psigma * Pclock * Pcool), Prho * Pcv * h * Pac);
            nT  = DIV(T, pow(aus, DIV(1.0, 3.0)));

            if(nT > PTsol) {

                SET(ST_next, c, i, j, nT);
            
            } else {

                SET(Sz_next, c, i, j, z + h);
                SET(Sh_next, c, i, j, 0.0);
                SET(ST_next, c, i, j, PTsol);

                SET(Mhs, c, i, j, GET(Mhs, c, i, j) + h);

            }

        }
        
    }

}



__global__
void reduce_add(size_t r, size_t c, double* buffer, float* acc) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i < r && j < c) {

        double v = buffer[i * c + j];

        if(v != 0.0) {
            atomicAdd(acc, v);
        }

    }

}


__global__
void memcpy_gpu(double *dst, double *src, int r, int c) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < r && j < c) {

        dst[i * c + j] = src[i * c + j];

    }

}



int main(int argc, char** argv) {


#if defined(DEBUG)

    const char* __argv[] = {
        argv[0],
        "8",
        "./data/2006/2006_000000000000.cfg",
        "./data/2006/output_2006",
        "1000",
        "1000",
        "1.0",
        NULL
    };

    argc = 7;
    argv = const_cast<char**>(&__argv[0]);

#endif


    if(argc < 7) {
        std::cout << "Usage: " << argv[0] << " <8|16|32> <INPUT> <OUTPUT> <STEPS> <REDUCE> <THICKNESS>" << std::endl;
        return 1;
    }


    size_t THREADS_PER_BLOCK = atol(argv[1]);
    const char* INPUT_CFG = argv[2];
    const char* OUTPUT_PATH = argv[3];
    size_t STEPS = atol(argv[4]);
    size_t REDUCE = atol(argv[5]);
    float THICKNESS = atof(argv[6]);

    
    
    cudaDeviceProp props;
    if(cudaGetDeviceProperties(&props, 0) != cudaSuccess) {
        std::cout << "Error: cudaGetDeviceProperties" << std::endl;
        return 1;
    }

    if(THREADS_PER_BLOCK != 8 && THREADS_PER_BLOCK != 16 && THREADS_PER_BLOCK != 32) {
        std::cout << "Error: THREADS_PER_BLOCK must be 8, 16 or 32" << std::endl;
        return 1;
    }

    if(THREADS_PER_BLOCK * THREADS_PER_BLOCK > props.maxThreadsPerBlock) {
        std::cout << "Error: THREADS_PER_BLOCK must be <= " << props.maxThreadsPerBlock << std::endl;
        return 1;
    }

    if(THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double) * 2 > props.sharedMemPerBlock) {
        std::cout << "Error: THREADS_PER_BLOCK^2 must be <= " << props.sharedMemPerBlock / (sizeof(double) * 2) << std::endl;
        return 1;
    }

    if(THREADS_PER_BLOCK * THREADS_PER_BLOCK % props.warpSize != 0) {
        std::cout << "Error: THREADS_PER_BLOCK must be divisible by " << props.warpSize << std::endl;
        return 1;
    }

    if(STEPS < 1) {
        std::cout << "Error: STEPS must be >= 1" << std::endl;
        return 1;
    }

    if(REDUCE < 1) {
        std::cout << "Error: REDUCE must be >= 1" << std::endl;
        return 1;
    }

    if(THICKNESS < 0.0f) {
        std::cout << "Error: THICKNESS must be >= 0.0" << std::endl;
        return 1;
    }


#if defined(DEBUG)
    fprintf(stdout, "Device '%s' props:\n", props.name);
    fprintf(stdout, " - props.warpSize: %lu\n", props.warpSize);
    fprintf(stdout, " - props.maxThreadsPerBlock: %lu\n", props.maxThreadsPerBlock);
    fprintf(stdout, " - props.sharedMemPerBlock: %lu\n", props.sharedMemPerBlock);
    fprintf(stdout, " - props.totalConstMem: %lu\n", props.totalConstMem);
#endif



    Sciara* sciara;
    init(sciara);
    
    if(loadConfiguration(INPUT_CFG, sciara) != 1) {
        std::cout << "Error: loadConfiguration '" << INPUT_CFG << "'" << std::endl;
        return 1;
    }


    size_t M = sciara->domain->rows;
    size_t N = sciara->domain->cols;


    dim3 threads_1d (
        THREADS_PER_BLOCK * THREADS_PER_BLOCK
    );

    dim3 threads_2d (
        THREADS_PER_BLOCK,
        THREADS_PER_BLOCK
    );

    dim3 wgrid {
        max(1U, (unsigned int) ceil(double(M) / THREADS_PER_BLOCK)),
        max(1U, (unsigned int) ceil(double(N) / THREADS_PER_BLOCK))
    };

    dim3 vgrid {
        max(1U, (unsigned int) ceil(double(sizeof(vents) / sizeof(vents[0])) / THREADS_PER_BLOCK)),
    };



    simulationInitialize(sciara);



    TVent h_vent[sciara->simulation->vent.size()];

    for(size_t i = 0; i < sciara->simulation->vent.size(); i++) {
        memcpy(&h_vent[i], &sciara->simulation->vent[i], sizeof(TVent));
    }


    assert(sciara->simulation->vent.size() == VENTS_COUNT);

    if(cudaMemcpyToSymbol(vents, h_vent, sizeof(TVent) * sciara->simulation->vent.size()) != cudaSuccess) {
        std::cout << "Error: cudaMemcpyToSymbol 'h_vent'" << std::endl;
        return 1;
    }

    if(cudaMemcpyToSymbol(Xi, sciara->X->Xi, sizeof(int) * MOORE_NEIGHBORS) != cudaSuccess) {
        std::cout << "Error: cudaMemcpyToSymbol 'Xi'" << std::endl;
        return 1;
    }

    if(cudaMemcpyToSymbol(Xj, sciara->X->Xj, sizeof(int) * MOORE_NEIGHBORS) != cudaSuccess) {
        std::cout << "Error: cudaMemcpyToSymbol 'Xj'" << std::endl;
        return 1;
    }



    float* total_emitted_lava;
    float* total_current_lava;

    if(cudaMallocManaged(&total_emitted_lava, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMallocManaged 'total_emitted_lava'" << std::endl;
        return 1;
    }

    if(cudaMallocManaged(&total_current_lava, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMallocManaged 'total_current_lava'" << std::endl;
        return 1;
    }


    *total_emitted_lava =  0.0f;
    *total_current_lava = -1.0f;



    util::Timer cl_timer;

    while((*total_current_lava == -1.0f || *total_current_lava > THICKNESS) || (sciara->simulation->elapsed_time <= sciara->simulation->effusion_duration)) {

        sciara->simulation->elapsed_time += sciara->parameters->Pclock;
        sciara->simulation->step++;


        emit_lava<<<vgrid, threads_1d>>> (
            M, N,
            sciara->simulation->elapsed_time,
            sciara->parameters->Pclock,
            sciara->simulation->emission_time,
            sciara->parameters->Pac,
            sciara->parameters->PTvent,
            sciara->substates->Sh,
            sciara->substates->Sh_next,
            sciara->substates->ST_next,
            total_emitted_lava
        );

        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);


        compute_outflows<<<wgrid, threads_2d>>> (
            M, N,
            sciara->substates->Sz,
            sciara->substates->Sh,
            sciara->substates->ST,
            sciara->substates->Mf,
            sciara->parameters->Pc,
            sciara->parameters->a,
            sciara->parameters->b,
            sciara->parameters->c,
            sciara->parameters->d
        );


        mass_balance<<<wgrid, threads_2d>>> (
            M, N,
            sciara->substates->Sh,
            sciara->substates->Sh_next,
            sciara->substates->ST,
            sciara->substates->ST_next,
            sciara->substates->Mf
        );

        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);


        compute_new_temperature_and_solidification<<<wgrid, threads_2d>>> (
            M, N,
            sciara->parameters->Pepsilon,
            sciara->parameters->Psigma,
            sciara->parameters->Pclock,
            sciara->parameters->Pcool,
            sciara->parameters->Prho,
            sciara->parameters->Pcv,
            sciara->parameters->Pac,
            sciara->parameters->PTsol,
            sciara->substates->Sz,
            sciara->substates->Sz_next,
            sciara->substates->Sh,
            sciara->substates->Sh_next,
            sciara->substates->ST,
            sciara->substates->ST_next,
            sciara->substates->Mf,
            sciara->substates->Mhs,
            sciara->substates->Mb
        );

        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->Sz, sciara->substates->Sz_next, M, N);
        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<wgrid, threads_2d>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);

        cudaDeviceSynchronize();

        
        if(sciara->simulation->step % REDUCE == 0) {

            *total_current_lava = 0.0f;

            reduce_add<<<wgrid, threads_2d>>>(
                M, N,
                sciara->substates->Sh,
                total_current_lava
            );

            cudaDeviceSynchronize();

        }


#if defined(DEBUG)
        if(sciara->simulation->step % 100 == 0) {
            fprintf(stdout, "\r[%08d]: %3d%% (%.2f s) [%f]", 
                sciara->simulation->step,
                (int)(double(sciara->simulation->step) / 16000.0 * 100.0),
                cl_timer.getTimeMilliseconds() / 1000.0f,
                *total_current_lava);
            fflush(stdout);
        }
#endif

    }

    cudaDeviceSynchronize();


    double cl_time = static_cast<double>(cl_timer.getTimeMilliseconds()) / 1000.0;

    fprintf(stdout, "\n");
    fprintf(stdout, "Step %d\n", sciara->simulation->step);
    fprintf(stdout, "Elapsed time [s]: %lf\n", cl_time);
    fprintf(stdout, "Emitted lava [m]: %lf\n", *total_emitted_lava);
    fprintf(stdout, "Current lava [m]: %lf\n", *total_current_lava);


    fprintf(stdout, "Saving output to %s...\n", OUTPUT_PATH);

    if(saveConfiguration(OUTPUT_PATH, sciara) != 1) {
        std::cout << "Error: saveConfiguration" << std::endl;
        return 1;
    }


    fprintf(stdout, "Releasing memory...\n");
    finalize(sciara);

    return 0;

}