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

#include <omp.h>



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


#define OUTSIDE(i, o, s)  \
    (i < o || i >= o + s)




__constant__ int Xi[MOORE_NEIGHBORS];
__constant__ int Xj[MOORE_NEIGHBORS];
__constant__ TVent vents[VENTS_COUNT];



__global__
void emit_lava (
    size_t r, 
    size_t c,
    size_t offset,
    size_t offset_size,
    double elapsed_time, 
    double Pclock, 
    double emission_time,
    double Pac,
    double PTvent,
    double* Sh,
    double* ST,
    double* ST_next,
    float* total_emitted_lava
) {

    constexpr const size_t vents_size = sizeof(vents) / sizeof(TVent);


    const size_t k = blockIdx.x * blockDim.x + threadIdx.x;

    if(k < vents_size) {

        const size_t i = vents[k].y();
        const size_t j = vents[k].x();


        if(i >= offset && i < offset_size) {

            double v = vents[k].thickness(elapsed_time, Pclock, emission_time, Pac);

            SET(Sh, c, i, j, GET(Sh, c, i, j) + v);
            SET(ST, c, i, j, PTvent);
            SET(ST_next, c, i, j, PTvent);

            if(v != 0.0) {
                atomicAdd(total_emitted_lava, v);
            }

        }

    }

}



__global__
void compute_outflows(
    size_t r,
    size_t c,
    size_t offset,
    size_t offset_size,
    double* Sz,
    double* Sz_halo,
    double* Sh,
    double* Sh_halo,
    double* ST,
    double* Mf,
    double Pc,
    double _a,
    double _b,
    double _c,
    double _d,
    size_t TILE_PITCH
) {

    extern __shared__ double shared[];

    double* sh_shared = &shared[0];
    double* sz_shared = &shared[TILE_PITCH * TILE_PITCH];



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

    
    const size_t i = blockIdx.x * blockDim.x + threadIdx.x + offset;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;


    if(OUTSIDE(i, offset, offset_size))
        return;



    if(i < r && j < c) {

        sh_shared[threadIdx.y * TILE_PITCH + threadIdx.x] = GET(Sh, c, i, j);
        sz_shared[threadIdx.y * TILE_PITCH + threadIdx.x] = GET(Sz, c, i, j);

    } else {

        sh_shared[threadIdx.y * TILE_PITCH + threadIdx.x] = 0.0;
        sz_shared[threadIdx.y * TILE_PITCH + threadIdx.x] = 0.0;

    }

    __syncthreads();


    if(i > 0 && j > 0 && i < r - 1 && j < c - 1) {


        if(sh_shared[threadIdx.y * TILE_PITCH + threadIdx.x] <= 0)
            return;


        T  = GET(ST, c, i, j);
        rr = pow(10, _a + _b * T);
        hc = pow(10, _c + _d * T);


        #pragma unroll
        for(size_t k = 0; k < MOORE_NEIGHBORS; k++) {


            if(OUTSIDE(i + Xi[k], offset, offset_size)) {

                sz    = GET(Sz_halo, c, i + Xi[k], j + Xj[k]);
                h[k]  = GET(Sh_halo, c, i + Xi[k], j + Xj[k]);   
            
            } else {
            
                if(threadIdx.x + Xi[k] >= blockDim.x || threadIdx.y + Xj[k] >= blockDim.y) {

                    sz    = GET(Sz, c, i + Xi[k], j + Xj[k]);
                    h[k]  = GET(Sh, c, i + Xi[k], j + Xj[k]);        

                } else {

                    sz    = sz_shared[TILE_PITCH * (threadIdx.y + Xj[k]) + threadIdx.x + Xi[k]];
                    h[k]  = sh_shared[TILE_PITCH * (threadIdx.y + Xj[k]) + threadIdx.x + Xi[k]];

                }

            }

            sz0   = sz_shared[threadIdx.y * TILE_PITCH + threadIdx.x];


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
    size_t offset,
    size_t offset_size,
    double* Sh,
    double* ST,
    double* ST_halo,
    double* ST_next,
    double* Mf,
    double* Mf_halo,
    size_t TILE_PITCH
) {

    extern __shared__ double st_shared[];


    const uint8_t inflows_indices[NUMBER_OF_OUTFLOWS] = { 3, 2, 1, 0, 6, 7, 4, 5 };

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x + offset;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(OUTSIDE(i, offset, offset_size))
        return;


    if(i > 0 && j > 0 && i < r - 1 && j < c - 1) {

        double t_initial = GET(ST, c, i, j);
        double h_initial = GET(Sh, c, i, j);

        st_shared[TILE_PITCH * threadIdx.y + threadIdx.x] = t_initial;

        __syncthreads();


        double inflow;
        double outflow;
        double t_neigh;

        double h_next = h_initial;
        double t_next = h_initial * t_initial;


        #pragma unroll
        for(size_t n = 1; n < MOORE_NEIGHBORS; n++) {


            if(OUTSIDE(i + Xi[n], offset, offset_size)) {

                t_neigh = GET(ST_halo, c, i + Xi[n], j + Xj[n]);
                inflow  = BUF_GET(Mf_halo, r, c, inflows_indices[n - 1], i + Xi[n], j + Xj[n]);
            
            } else {

                if(threadIdx.x + Xi[n] >= blockDim.x || threadIdx.y + Xj[n] >= blockDim.y) {

                    t_neigh = GET(ST, c, i + Xi[n], j + Xj[n]);

                } else {

                    t_neigh = st_shared[TILE_PITCH * (threadIdx.y + Xj[n]) + threadIdx.x + Xi[n]];

                }


                inflow  = BUF_GET(Mf, r, c, inflows_indices[n - 1], i + Xi[n], j + Xj[n]);

            }

            
            outflow = BUF_GET(Mf, r, c, n - 1, i, j);

            h_next += inflow - outflow;
            t_next += inflow * t_neigh - outflow * t_initial;

        }

        if(h_next > 0.0) {

            t_next = DIV(t_next, h_next);

            SET(ST_next, c, i, j, t_next);
            SET(Sh, c, i, j, h_next);

        }


    }

}


__global__
void compute_new_temperature_and_solidification (
    size_t r,
    size_t c,
    size_t offset,
    size_t offset_size,
    double Pepsilon,
    double Psigma,
    double Pclock,
    double Pcool,
    double Prho,
    double Pcv,
    double Pac,
    double PTsol,
    double *Sz,
    double *Sh,
    double *ST,
    double *ST_next,
    double *Mf,
    double *Mhs,
    bool   *Mb
) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x + offset;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(OUTSIDE(i, offset, offset_size))
        return;


    double nT;
    double aus;

    double h = GET(Sh, c, i, j);
    double T = GET(ST, c, i, j);


    if(h > 0 && GET(Mb, c, i, j) == false) {

        aus = 1.0 + DIV((3 * pow(T, 3.0) * Pepsilon * Psigma * Pclock * Pcool), Prho * Pcv * h * Pac);
        nT  = DIV(T, pow(aus, DIV(1.0, 3.0)));

        if(nT > PTsol) {

            SET(ST, c, i, j, nT);
            SET(ST_next, c, i, j, nT);
        
        } else {

            SET(Sz, c, i, j, GET(Sz, c, i, j) + h);
            SET(Sh, c, i, j, 0.0);

            SET(ST, c, i, j, PTsol);
            SET(ST_next, c, i, j, PTsol);

            SET(Mhs, c, i, j, GET(Mhs, c, i, j) + h);

        }

    }
        

}



__global__
void reduce_add(size_t size, double* buffer, float* acc) {

    extern __shared__ double shared[];


    shared[threadIdx.x] = 0.0;

    for(size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < size; i += gridDim.x * blockDim.x) {

        shared[threadIdx.x] += buffer[i];

    }

    __syncthreads();


    const size_t blocksize = blockDim.x * blockDim.y;

    if(blocksize >= 512 && threadIdx.x < 256) {
        shared[threadIdx.x] += shared[threadIdx.x + 256];
        __syncthreads();
    }

    if(blocksize >= 256 && threadIdx.x < 128) {
        shared[threadIdx.x] += shared[threadIdx.x + 128];
        __syncthreads();
    }

    if(blocksize >= 128 && threadIdx.x < 64) {
        shared[threadIdx.x] += shared[threadIdx.x + 64];
        __syncthreads();
    }


    double v = shared[threadIdx.x];

    if(threadIdx.x < 32) {
        v += __shfl_down_sync(0xFFFFFFFF, v, 32);
        v += __shfl_down_sync(0xFFFFFFFF, v, 16);
        v += __shfl_down_sync(0xFFFFFFFF, v,  8);
        v += __shfl_down_sync(0xFFFFFFFF, v,  4); 
        v += __shfl_down_sync(0xFFFFFFFF, v,  2); 
        v += __shfl_down_sync(0xFFFFFFFF, v,  1); 
    }

    if(threadIdx.x == 0 && v != 0.0) {
        atomicAdd(acc, v);
    }

}



__global__
void memcpy_gpu(double *dst, double *src, size_t size) {

    const size_t GRIDE_STRIDE_SIZE = gridDim.x * blockDim.x;

    for(size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < size; i += GRIDE_STRIDE_SIZE) {

        dst[i] = src[i];

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


    size_t THREADS_PER_BLOCK    = atol(argv[1]);
    const char* INPUT_CFG       = argv[2];
    const char* OUTPUT_PATH     = argv[3];
    size_t STEPS                = atol(argv[4]);
    size_t REDUCE               = atol(argv[5]);
    float THICKNESS             = atof(argv[6]);

    
    
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



    int num_procs;
    if(cudaGetDeviceCount(&num_procs) != cudaSuccess) {
        std::cout << "Error: cudaGetDeviceCount" << std::endl;
        return 1;
    }


    for(size_t rank = 0; rank < num_procs; rank++) {

        if(cudaSetDevice(rank) != cudaSuccess) {
            std::cout << "Error: cudaSetDevice #" << rank << std::endl;
            return 1;
        }

        int is_able;
        if(cudaDeviceCanAccessPeer(&is_able, rank, !rank) != cudaSuccess) {
            std::cout << "Error: cudaDeviceCanAccessPeer #" << rank << std::endl;
            return 1;
        }

        if(is_able) {
            if(cudaDeviceEnablePeerAccess(!rank, 0) != cudaSuccess) {
                std::cout << "Error: cudaDeviceEnablePeerAccess #" << rank << std::endl;
                return 1;
            }

#if defined(DEBUG)
            std::cout << "cudaDeviceEnablePeerAccess #" << rank << " to #" << !rank << std::endl;
#endif
        }


#if defined(DEBUG)
        {
            cudaDeviceProp props;
            cudaGetDeviceProperties(&props, rank);
            fprintf(stdout, "Device '%s' (%d) props:\n", props.name, props.multiGpuBoardGroupID);
            fprintf(stdout, " - props.warpSize: %lu\n", props.warpSize);
            fprintf(stdout, " - props.maxThreadsPerBlock: %lu\n", props.maxThreadsPerBlock);
            fprintf(stdout, " - props.sharedMemPerBlock: %lu\n", props.sharedMemPerBlock);
            fprintf(stdout, " - props.totalConstMem: %lu\n", props.totalConstMem);
            fprintf(stdout, " - props.unifiedAddressing: %d\n", props.unifiedAddressing);
        }
#endif


    }


    cudaSetDevice(0);


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
        max(1U, (unsigned int) ceil(double(M) / THREADS_PER_BLOCK / num_procs)),
        max(1U, (unsigned int) ceil(double(N) / THREADS_PER_BLOCK))
    };
    
    dim3 mgrid_stride_8 {
        max(1U, (unsigned int) ceil(double(M * N) / THREADS_PER_BLOCK / 8 / num_procs)),
    };

    dim3 rgrid_stride_8 {
        max(1U, (unsigned int) ceil(double(M * N) / THREADS_PER_BLOCK / 8 / num_procs))
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


    for(size_t rank = 0; rank < num_procs; rank++) {

        if(cudaSetDevice(rank) != cudaSuccess) {
            std::cout << "Error: cudaSetDevice #" << rank << std::endl;
            return 1;
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

    }

    cudaSetDevice(0);


    float* d_total_emitted_lava;
    float* d_total_current_lava[num_procs];


    for(size_t rank = 0; rank < num_procs; rank++) {

        cudaSetDevice(rank);
        
        if(cudaMalloc(&d_total_current_lava[rank], sizeof(float)) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'total_current_lava'" << std::endl;
            return 1;
        }

    }


    if(cudaMalloc(&d_total_emitted_lava, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMallocManaged 'total_emitted_lava'" << std::endl;
        return 1;
    }

    if(cudaMemset(d_total_emitted_lava, 0, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMemset 'total_emitted_lava'" << std::endl;
        return 1;
    }


    float total_emitted_lava =  0.0f;
    float total_current_lava = -1.0f;
    float local_current_lava = -1.0f;


    double* Sh[num_procs];
    double* Sz[num_procs];
    double* ST[num_procs];
    double* ST_next[num_procs];
    double* ST_temp[num_procs];
    double* Mf[num_procs];
    double* Mhs[num_procs];
    bool*   Mb[num_procs];


    #pragma omp parallel for num_threads(num_procs)
    for(size_t rank = 0; rank < num_procs; rank++) {

        cudaSetDevice(rank);

        cudaMalloc(&Sz[rank], sizeof(double) * M * N);
        cudaMalloc(&Sh[rank], sizeof(double) * M * N);
        cudaMalloc(&ST[rank], sizeof(double) * M * N);
        cudaMalloc(&ST_next[rank], sizeof(double) * M * N);
        cudaMalloc(&ST_temp[rank], sizeof(double) * M * N);
        cudaMalloc(&Mf[rank], sizeof(double) * M * N * NUMBER_OF_OUTFLOWS);
        cudaMalloc(&Mhs[rank], sizeof(double) * M * N);
        cudaMalloc(&Mb[rank], sizeof(bool) * M * N);

        cudaMemcpy(Sh[rank], sciara->substates->Sh, sizeof(double) * M * N, cudaMemcpyDefault);
        cudaMemcpy(Sz[rank], sciara->substates->Sz, sizeof(double) * M * N, cudaMemcpyDefault);
        cudaMemcpy(ST[rank], sciara->substates->ST, sizeof(double) * M * N, cudaMemcpyDefault);
        cudaMemcpy(ST_next[rank], sciara->substates->ST_next, sizeof(double) * M * N, cudaMemcpyDefault);
        cudaMemcpy(Mf[rank], sciara->substates->Mf, sizeof(double) * M * N * NUMBER_OF_OUTFLOWS, cudaMemcpyDefault);
        cudaMemcpy(Mhs[rank], sciara->substates->Mhs, sizeof(double) * M * N, cudaMemcpyDefault);
        cudaMemcpy(Mb[rank], sciara->substates->Mb, sizeof(bool) * M * N, cudaMemcpyDefault);

    }

    cudaSetDevice(0);




    #define WORLD_RANK       !rank
    #define LOCAL_RANK      !!rank

    uintptr_t WORLD_SIZE = M * N;
    uintptr_t LOCAL_SIZE = M * N / num_procs;

    uintptr_t LOCAL_ROWS = M / num_procs;
    uintptr_t LOCAL_COLS = N / num_procs;




    util::Timer cl_timer;

    while((total_current_lava == -1.0f || total_current_lava > THICKNESS) || (sciara->simulation->elapsed_time <= sciara->simulation->effusion_duration)) {

        sciara->simulation->elapsed_time += sciara->parameters->Pclock;
        sciara->simulation->step++;


        #pragma omp paraller for num_threads(num_procs) schedule(static)
        for(size_t rank = 0; rank < num_procs; rank++) {

            cudaSetDevice(rank);

            emit_lava<<<vgrid, threads_1d>>> (
                M, N,
                LOCAL_ROWS * LOCAL_RANK,
                LOCAL_ROWS,
                sciara->simulation->elapsed_time,
                sciara->parameters->Pclock,
                sciara->simulation->emission_time,
                sciara->parameters->Pac,
                sciara->parameters->PTvent,
                Sh[LOCAL_RANK],
                ST[LOCAL_RANK],
                ST_next[LOCAL_RANK],
                d_total_emitted_lava
            );

        }



        #pragma omp paraller for num_threads(num_procs) schedule(static)
        for(size_t rank = 0; rank < num_procs; rank++) {

            cudaSetDevice(rank);

            compute_outflows<<<wgrid, threads_2d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double) * 2>>> (
                M, N,
                LOCAL_ROWS * LOCAL_RANK,
                LOCAL_ROWS,
                Sz[LOCAL_RANK], Sz[WORLD_RANK],
                Sh[LOCAL_RANK], Sh[WORLD_RANK],
                ST[LOCAL_RANK],
                Mf[LOCAL_RANK],
                sciara->parameters->Pc,
                sciara->parameters->a,
                sciara->parameters->b,
                sciara->parameters->c,
                sciara->parameters->d,
                THREADS_PER_BLOCK
            );

        }



        cudaMemcpyPeerAsync(&ST_temp[0][LOCAL_SIZE - 0], 0, &ST[1][LOCAL_SIZE - 0], 1, sizeof(double) * N);
        cudaMemcpyPeerAsync(&ST_temp[1][LOCAL_SIZE - N], 1, &ST[0][LOCAL_SIZE - N], 0, sizeof(double) * N);

        #pragma omp paraller for num_threads(num_procs) schedule(static)
        for(size_t rank = 0; rank < num_procs; rank++) {

            cudaSetDevice(rank);

            mass_balance<<<wgrid, threads_2d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double)>>> (
                M, N,
                LOCAL_ROWS * LOCAL_RANK,
                LOCAL_ROWS,
                Sh[LOCAL_RANK],
                ST[LOCAL_RANK], 
                ST_temp[LOCAL_RANK],
                ST_next[LOCAL_RANK],
                Mf[LOCAL_RANK], Mf[WORLD_RANK],
                THREADS_PER_BLOCK
            );

        }

        #pragma omp paraller for num_threads(num_procs) schedule(static)
        for(size_t rank = 0; rank < num_procs; rank++) {

            cudaSetDevice(rank);

            memcpy_gpu<<<mgrid_stride_8, threads_1d>>>(&ST[rank][LOCAL_SIZE * LOCAL_RANK], &ST_next[rank][LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE);

        }



        #pragma omp paraller for num_threads(num_procs) schedule(static)
        for(size_t rank = 0; rank < num_procs; rank++) {

            cudaSetDevice(rank);

            compute_new_temperature_and_solidification<<<wgrid, threads_2d>>> (
                M, N,
                LOCAL_ROWS * LOCAL_RANK,
                LOCAL_ROWS,
                sciara->parameters->Pepsilon,
                sciara->parameters->Psigma,
                sciara->parameters->Pclock,
                sciara->parameters->Pcool,
                sciara->parameters->Prho,
                sciara->parameters->Pcv,
                sciara->parameters->Pac,
                sciara->parameters->PTsol,
                Sz[LOCAL_RANK],
                Sh[LOCAL_RANK],
                ST[LOCAL_RANK],
                ST_next[LOCAL_RANK],
                Mf[LOCAL_RANK],
                Mhs[LOCAL_RANK],
                Mb[LOCAL_RANK]
            );

        }


        if(sciara->simulation->step % REDUCE == 0) {

            total_current_lava = 0.0f;

            #pragma omp paraller for num_threads(num_procs) reduction(+:total_current_lava) schedule(static)
            for(size_t rank = 0; rank < num_procs; rank++) {

                cudaSetDevice(rank);
                cudaMemset(d_total_current_lava[rank], 0, sizeof(float));

                reduce_add<<<rgrid_stride_8, threads_1d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double)>>>(
                    LOCAL_SIZE,
                    &Sh[rank][LOCAL_SIZE * rank],
                    d_total_current_lava[rank]
                );


                cudaDeviceSynchronize();
                cudaMemcpy(&local_current_lava, d_total_current_lava[rank], sizeof(float), cudaMemcpyDeviceToHost);


                total_current_lava += local_current_lava;
            
            }

        }


#if defined(DEBUG)
        if(sciara->simulation->step % 100 == 0) {
            fprintf(stderr, "\r[%08d]: %3d%% (%.2f s) [%f]",
                sciara->simulation->step,
                (int)(double(sciara->simulation->step) / 16000.0 * 100.0),
                cl_timer.getTimeMilliseconds() / 1000.0f,
                total_current_lava);
        }
#endif

    }

    double cl_time = static_cast<double>(cl_timer.getTimeMilliseconds()) / 1000.0;


    cudaSetDevice(0);
    cudaDeviceSynchronize();
    cudaMemcpy(&total_emitted_lava, d_total_emitted_lava, sizeof(float), cudaMemcpyDeviceToHost);


    #pragma omp paraller for num_threads(num_procs)
    for(size_t rank = 0; rank < num_procs; rank++) {

        cudaSetDevice(rank);

        cudaMemcpy(&sciara->substates->Sz[LOCAL_SIZE * LOCAL_RANK],      &Sz[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK],        sizeof(double) * LOCAL_SIZE, cudaMemcpyDefault);
        cudaMemcpy(&sciara->substates->Sh[LOCAL_SIZE * LOCAL_RANK],      &Sh[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK],        sizeof(double) * LOCAL_SIZE, cudaMemcpyDefault);
        cudaMemcpy(&sciara->substates->ST[LOCAL_SIZE * LOCAL_RANK],      &ST[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK],        sizeof(double) * LOCAL_SIZE, cudaMemcpyDefault);
        cudaMemcpy(&sciara->substates->Mhs[LOCAL_SIZE * LOCAL_RANK],     &Mhs[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK],       sizeof(double) * LOCAL_SIZE, cudaMemcpyDefault);
        cudaMemcpy(&sciara->substates->Mb[LOCAL_SIZE * LOCAL_RANK],      &Mb[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK],        sizeof(bool)   * LOCAL_SIZE, cudaMemcpyDefault);

        for(size_t k = 0; k < NUMBER_OF_OUTFLOWS; k++) {

            cudaMemcpy(&sciara->substates->Mf[LOCAL_SIZE * LOCAL_RANK + WORLD_SIZE * k], &Mf[LOCAL_RANK][LOCAL_SIZE * LOCAL_RANK + WORLD_SIZE * k], sizeof(double) * LOCAL_SIZE, cudaMemcpyDefault);

        }

    }



    fprintf(stdout, "\n");
    fprintf(stdout, "Step %d\n", sciara->simulation->step);
    fprintf(stdout, "Elapsed time [s]: %lf\n", cl_time);
    fprintf(stdout, "Emitted lava [m]: %lf\n", total_emitted_lava);
    fprintf(stdout, "Current lava [m]: %lf\n", total_current_lava);


    fprintf(stdout, "Saving output to %s...\n", OUTPUT_PATH);

    if(saveConfiguration(OUTPUT_PATH, sciara) != 1) {
        std::cout << "Error: saveConfiguration" << std::endl;
        return 1;
    }



    fprintf(stdout, "Releasing memory...\n");
    finalize(sciara);


    return 0;

}