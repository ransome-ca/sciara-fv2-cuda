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
__constant__ TVent vents[2];



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

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i < r && j < c) {

        constexpr const size_t vents_size = sizeof(vents) / sizeof(TVent);

        for(size_t k = 0; k < vents_size; k++) {

            if(vents[k].y() == i && vents[k].x() == j) {

                double v = vents[k].thickness(elapsed_time, Pclock, emission_time, Pac);

                SET(Sh_next, c, i, j, GET(Sh, c, i, j) + v);
                SET(ST_next, c, i, j, PTvent);

                if(v != 0.0) {
                    atomicAdd(total_emitted_lava, v);
                }

            }

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
    double w[MOORE_NEIGHBORS] = {};
    double Pr[MOORE_NEIGHBORS] = {};
    
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
            w[k]  = Pc;
            Pr[k] = rr;

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
                theta[k]      = atan(DIV(((z[0] + h[0]) - (z[k] + h[k])), w[k]));
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

                BUF_SET(Mf, r, c, k - 1, i, j, Pr[k] * (avg - H[k]));

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




    Sciara* sciara;
    init(sciara);
    
    if(loadConfiguration(INPUT_CFG, sciara) != 1) {
        std::cout << "Error: loadConfiguration '" << INPUT_CFG << "'" << std::endl;
        return 1;
    }


    size_t M = sciara->domain->rows;
    size_t N = sciara->domain->cols;

    dim3 threads (
        THREADS_PER_BLOCK,
        THREADS_PER_BLOCK
    );

    dim3 grid {
        max(1U, (unsigned int) ceil(double(M) / THREADS_PER_BLOCK)),
        max(1U, (unsigned int) ceil(double(N) / THREADS_PER_BLOCK))
    };



    simulationInitialize(sciara);


    assert(sciara->simulation->vent.size() == 2);


    TVent h_vent[2];

    for(size_t i = 0; i < 2; i++) {
        memcpy(&h_vent[i], &sciara->simulation->vent[i], sizeof(TVent));
    }


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


        emit_lava<<<grid, threads>>> (
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

        memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);


        compute_outflows<<<grid, threads>>> (
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


        mass_balance<<<grid, threads>>> (
            M, N,
            sciara->substates->Sh,
            sciara->substates->Sh_next,
            sciara->substates->ST,
            sciara->substates->ST_next,
            sciara->substates->Mf
        );

        memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);


        compute_new_temperature_and_solidification<<<grid, threads>>> (
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

        memcpy_gpu<<<grid, threads>>>(sciara->substates->Sz, sciara->substates->Sz_next, M, N);
        memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
        memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);

        cudaDeviceSynchronize();

        
        if(sciara->simulation->step % REDUCE == 0) {

            *total_current_lava = 0.0f;

            reduce_add<<<grid, threads>>>(
                M, N,
                sciara->substates->Sh,
                total_current_lava
            );

            cudaDeviceSynchronize();

        }


#if defined(DEBUG)
        if(sciara->simulation->step % 100 == 0) {
            fprintf(stderr, "\r[%08d]: %3d%% (%.2f s) [%f / %f]", 
                sciara->simulation->step,
                (int)(double(sciara->simulation->step) / 16000.0 * 100.0),
                cl_timer.getTimeMilliseconds() / 1000.0f,
                *total_current_lava,
                *total_emitted_lava);
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