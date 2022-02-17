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

#include <mpi.h>



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
    size_t rank,
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


        if((i >= ((r * c) >> 1) * rank) && (i < ((r * c) >> 1) * (rank + 1))) {

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
    double* Sz,
    double* Sh,
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
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y + offset;


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

            if(threadIdx.x + Xi[k] >= blockDim.x || threadIdx.y + Xj[k] >= blockDim.y) {
                
                sz    = GET(Sz, c, i + Xi[k], j + Xj[k]);
                h[k]  = GET(Sh, c, i + Xi[k], j + Xj[k]);

            } else {

                sz    = sz_shared[TILE_PITCH * (threadIdx.y + Xj[k]) + threadIdx.x + Xi[k]];
                h[k]  = sh_shared[TILE_PITCH * (threadIdx.y + Xj[k]) + threadIdx.x + Xi[k]];

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
    double* Sh,
    double* ST,
    double* ST_next,
    double* Mf,
    size_t TILE_PITCH
) {

    extern __shared__ double st_shared[];


    const uint8_t inflows_indices[NUMBER_OF_OUTFLOWS] = { 3, 2, 1, 0, 6, 7, 4, 5 };

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x + offset;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y + offset;


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

            if(threadIdx.x + Xi[n] >= blockDim.x || threadIdx.y + Xj[n] >= blockDim.y) {

                t_neigh = GET(ST, c, i + Xi[n], j + Xj[n]);

            } else {

                t_neigh = st_shared[TILE_PITCH * (threadIdx.y + Xj[n]) + threadIdx.x + Xi[n]];

            }


            inflow  = BUF_GET(Mf, r, c, inflows_indices[n - 1], i + Xi[n], j + Xj[n]);
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
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y + offset;


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


    int rank;
    int num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


#if defined(DEBUG)
    std::cout << "Running node #" << rank << " of " << num_procs << std::endl;
#endif


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
    if(cudaGetDeviceProperties(&props, rank) != cudaSuccess) {
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

    fprintf(stdout, "Device '%s' (%d) props:\n", props.name, props.multiGpuBoardGroupID);
    fprintf(stdout, " - props.warpSize: %lu\n", props.warpSize);
    fprintf(stdout, " - props.maxThreadsPerBlock: %lu\n", props.maxThreadsPerBlock);
    fprintf(stdout, " - props.sharedMemPerBlock: %lu\n", props.sharedMemPerBlock);
    fprintf(stdout, " - props.totalConstMem: %lu\n", props.totalConstMem);
    fprintf(stdout, " - props.unifiedAddressing: %d\n", props.unifiedAddressing);

#endif


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
    }




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
    
    dim3 wgrid_stride_2 {
        max(1U, (unsigned int) ceil(double(M) / THREADS_PER_BLOCK / 2 / num_procs)),
        max(1U, (unsigned int) ceil(double(N) / THREADS_PER_BLOCK / 2))
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




    float* d_total_emitted_lava;
    float* d_total_current_lava;

    if(cudaMalloc(&d_total_emitted_lava, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMallocManaged 'total_emitted_lava'" << std::endl;
        return 1;
    }

    if(cudaMalloc(&d_total_current_lava, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMallocManaged 'total_current_lava'" << std::endl;
        return 1;
    }

    if(cudaMemset(d_total_emitted_lava, 0, sizeof(float)) != cudaSuccess) {
        std::cout << "Error: cudaMemset 'total_emitted_lava'" << std::endl;
        return 1;
    }


    float total_emitted_lava =  0.0f;
    float total_current_lava = -1.0f;

    float local_current_lava = 0.0f;
    float world_current_lava = 0.0f;


    double* Sh;
    double* Sz;
    double* ST;
    double* ST_next;
    double* Mf;
    double* Mhs;
    bool*   Mb;

    cudaMalloc(&Sh, sizeof(double) * M * N);
    cudaMalloc(&Sz, sizeof(double) * M * N);
    cudaMalloc(&ST, sizeof(double) * M * N);
    cudaMalloc(&ST_next, sizeof(double) * M * N);
    cudaMalloc(&Mf, sizeof(double) * M * N * NUMBER_OF_OUTFLOWS);
    cudaMalloc(&Mhs, sizeof(double) * M * N);
    cudaMalloc(&Mb, sizeof(bool) * M * N);

    cudaMemcpy(Sh, sciara->substates->Sh, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(Sz, sciara->substates->Sz, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(ST, sciara->substates->ST, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(ST_next, sciara->substates->ST_next, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(Mf, sciara->substates->Mf, sizeof(double) * M * N * NUMBER_OF_OUTFLOWS, cudaMemcpyDefault);
    cudaMemcpy(Mhs, sciara->substates->Mhs, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(Mb, sciara->substates->Mb, sizeof(bool) * M * N, cudaMemcpyDefault);



    uintptr_t WORLD_RANK =  !rank;
    uintptr_t LOCAL_RANK = !!rank;
    uintptr_t WORLD_SIZE = M * N;
    uintptr_t LOCAL_SIZE = M * N / num_procs;


    MPI_Barrier(MPI_COMM_WORLD);


    util::Timer cl_timer;

    while((total_current_lava == -1.0f || total_current_lava > THICKNESS) || (sciara->simulation->elapsed_time <= sciara->simulation->effusion_duration)) {

        sciara->simulation->elapsed_time += sciara->parameters->Pclock;
        sciara->simulation->step++;


        emit_lava<<<vgrid, threads_1d>>> (
            M, N, LOCAL_RANK,
            sciara->simulation->elapsed_time,
            sciara->parameters->Pclock,
            sciara->simulation->emission_time,
            sciara->parameters->Pac,
            sciara->parameters->PTvent,
            Sh,
            ST,
            ST_next,
            d_total_emitted_lava
        );

        cudaDeviceSynchronize();

        cudaMemcpy(&sciara->substates->Sh[LOCAL_SIZE - WORLD_RANK], &Sh[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);
        cudaMemcpy(&sciara->substates->Sz[LOCAL_SIZE - WORLD_RANK], &Sz[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);

        MPI_Sendrecv(&sciara->substates->Sh[LOCAL_SIZE - WORLD_RANK], N, MPI_DOUBLE, WORLD_RANK, 0,
                     &sciara->substates->Sh[LOCAL_SIZE - LOCAL_RANK], N, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&sciara->substates->Sz[LOCAL_SIZE - WORLD_RANK], N, MPI_DOUBLE, WORLD_RANK, 0,
                     &sciara->substates->Sz[LOCAL_SIZE - LOCAL_RANK], N, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        cudaMemcpy(&Sz[LOCAL_SIZE - WORLD_RANK], &sciara->substates->Sz[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);
        cudaMemcpy(&Sh[LOCAL_SIZE - WORLD_RANK], &sciara->substates->Sh[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);



        compute_outflows<<<wgrid, threads_2d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double) * 2>>> (
            M, N, LOCAL_SIZE * LOCAL_RANK,
            Sz,
            Sh,
            ST,
            Mf,
            sciara->parameters->Pc,
            sciara->parameters->a,
            sciara->parameters->b,
            sciara->parameters->c,
            sciara->parameters->d,
            THREADS_PER_BLOCK
        );

        cudaDeviceSynchronize();

        cudaMemcpy(&sciara->substates->ST[LOCAL_SIZE - WORLD_RANK], &ST[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);

        MPI_Sendrecv(&sciara->substates->ST[LOCAL_SIZE - WORLD_RANK], N, MPI_DOUBLE, WORLD_RANK, 0,
                     &sciara->substates->ST[LOCAL_SIZE - LOCAL_RANK], N, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        cudaMemcpy(&ST[LOCAL_SIZE - WORLD_RANK], &sciara->substates->ST[LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);


        for(size_t k = 0; k < NUMBER_OF_OUTFLOWS; k++) {

            cudaMemcpy(&sciara->substates->Mf[(WORLD_SIZE * k) + LOCAL_SIZE - WORLD_RANK], &Mf[(WORLD_SIZE * k) + LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);

            MPI_Sendrecv(&sciara->substates->Mf[(WORLD_SIZE * k) + LOCAL_SIZE - WORLD_RANK], N, MPI_DOUBLE, WORLD_RANK, 0,
                         &sciara->substates->Mf[(WORLD_SIZE * k) + LOCAL_SIZE - LOCAL_RANK], N, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            cudaMemcpy(&Mf[(WORLD_SIZE * k) + LOCAL_SIZE - WORLD_RANK], &sciara->substates->Mf[(WORLD_SIZE * k) + LOCAL_SIZE - WORLD_RANK], N * sizeof(double), cudaMemcpyDefault);

        }


        mass_balance<<<wgrid, threads_2d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double)>>> (
            M, N, LOCAL_SIZE * LOCAL_RANK,
            Sh,
            ST,
            ST_next,
            Mf,
            THREADS_PER_BLOCK
        );

        memcpy_gpu<<<mgrid_stride_8, threads_1d>>>(&ST[LOCAL_SIZE * LOCAL_RANK], &ST_next[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE);



        compute_new_temperature_and_solidification<<<wgrid_stride_2, threads_2d>>> (
            M, N, LOCAL_SIZE * LOCAL_RANK,
            sciara->parameters->Pepsilon,
            sciara->parameters->Psigma,
            sciara->parameters->Pclock,
            sciara->parameters->Pcool,
            sciara->parameters->Prho,
            sciara->parameters->Pcv,
            sciara->parameters->Pac,
            sciara->parameters->PTsol,
            Sz,
            Sh,
            ST,
            ST_next,
            Mf,
            Mhs,
            Mb
        );


        if(sciara->simulation->step % REDUCE == 0) {

            cudaMemset(d_total_current_lava, 0, sizeof(float));

            reduce_add<<<rgrid_stride_8, threads_1d, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double)>>>(
                LOCAL_SIZE,
                &Sh[LOCAL_SIZE * LOCAL_RANK],
                d_total_current_lava
            );

            cudaDeviceSynchronize();
            cudaMemcpy(&local_current_lava, d_total_current_lava, sizeof(float), cudaMemcpyDeviceToHost);


            MPI_Sendrecv(&local_current_lava, 1, MPI_FLOAT, WORLD_RANK, 0,
                         &world_current_lava, 1, MPI_FLOAT, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            total_current_lava = local_current_lava + world_current_lava;

        }


#if defined(DEBUG)
        if(sciara->simulation->step % 100 == 0) {
            fprintf(stderr, "<%d> [%08d]: %3d%% (%.2f s) [%f]\n",
                rank, 
                sciara->simulation->step,
                (int)(double(sciara->simulation->step) / 16000.0 * 100.0),
                cl_timer.getTimeMilliseconds() / 1000.0f,
                total_current_lava);
        }
#endif

    }


    double cl_time = static_cast<double>(cl_timer.getTimeMilliseconds()) / 1000.0;


    cudaDeviceSynchronize();
    cudaMemcpy(&total_emitted_lava, d_total_emitted_lava, sizeof(float), cudaMemcpyDeviceToHost);


    MPI_Sendrecv(&Sh[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                 &Sh[LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&Sz[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                 &Sz[LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&ST[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                 &ST[LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&ST_next[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                 &ST_next[LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&Mhs[LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                 &Mhs[LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(size_t k = 0; k < NUMBER_OF_OUTFLOWS; k++) {

        MPI_Sendrecv(&Mf[(WORLD_SIZE * k) + LOCAL_SIZE * LOCAL_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0,
                     &Mf[(WORLD_SIZE * k) + LOCAL_SIZE * WORLD_RANK], LOCAL_SIZE, MPI_DOUBLE, WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }



    cudaMemcpy(sciara->substates->Sh,      Sh, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->Sz,      Sz, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->ST,      ST, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->ST_next, ST_next, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->Mf,      Mf, sizeof(double) * M * N * NUMBER_OF_OUTFLOWS, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->Mhs,     Mhs, sizeof(double) * M * N, cudaMemcpyDefault);
    cudaMemcpy(sciara->substates->Mb,      Mb, sizeof(bool) * M * N, cudaMemcpyDefault);







    fprintf(stdout, "\n");
    fprintf(stdout, "Step %d\n", sciara->simulation->step);
    fprintf(stdout, "Elapsed time [s]: %lf\n", cl_time);
    fprintf(stdout, "Emitted lava [m]: %lf\n", total_emitted_lava);
    fprintf(stdout, "Current lava [m]: %lf\n", total_current_lava);


    if(rank == 0) {

        fprintf(stdout, "Saving output to %s...\n", OUTPUT_PATH);

        if(saveConfiguration(OUTPUT_PATH, sciara) != 1) {
            std::cout << "Error: saveConfiguration" << std::endl;
            return 1;
        }

    }


    fprintf(stdout, "Releasing memory...\n");
    finalize(sciara);


    MPI_Finalize();

    return 0;

}