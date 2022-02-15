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
    double* ST,
    float* total_emitted_lava
) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;


    if(i < r && j < c) {


        double sh_local = GET(Sh, c, i, j);

        __syncthreads();    


        constexpr const size_t vents_size = sizeof(vents) / sizeof(TVent);

        for(size_t k = 0; k < vents_size; k++) {

            if(vents[k].y() == i && vents[k].x() == j) {

                double v = vents[k].thickness(elapsed_time, Pclock, emission_time, Pac);

                SET(Sh, c, i, j, sh_local + v);
                SET(ST, c, i, j, PTvent);

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
    double* ST,
    double* ST_halo,
    double* Mf,
    size_t TILE_PITCH
) {

    extern __shared__ double st_shared[];


    const uint8_t inflows_indices[NUMBER_OF_OUTFLOWS] = { 3, 2, 1, 0, 6, 7, 4, 5 };

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;


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

                t_neigh = GET(ST_halo, c, i + Xi[n], j + Xj[n]);

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

            SET(ST, c, i, j, t_next);
            SET(Sh, c, i, j, h_next);

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
    double *Sh,
    double *ST,
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

        __syncthreads();


        if(h > 0 && GET(Mb, c, i, j) == false) {

            aus = 1.0 + DIV((3 * pow(T, 3.0) * Pepsilon * Psigma * Pclock * Pcool), Prho * Pcv * h * Pac);
            nT  = DIV(T, pow(aus, DIV(1.0, 3.0)));

            if(nT > PTsol) {

                SET(ST, c, i, j, nT);
            
            } else {

                SET(Sz, c, i, j, z + h);
                SET(Sh, c, i, j, 0.0);
                SET(ST, c, i, j, PTsol);

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
void prepare_halo(double *dst, double *src, int r, int c) {

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < r && j < c) {

        if(threadIdx.x == 0 || threadIdx.x == blockDim.x - 1 || threadIdx.y == 0 || threadIdx.y == blockDim.y - 1) {

            dst[i * c + j] = src[i * c + j];

        }

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

    

    int cuda_devices = 1;
    if(cudaGetDeviceCount(&cuda_devices) != cudaSuccess) {
        std::cout << "Error: cudaGetDeviceCount" << std::endl;
        return 1;
    }

    
    

    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);


        cudaDeviceProp props;
        if(cudaGetDeviceProperties(&props, d) != cudaSuccess) {
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

    }





#if defined(DEBUG)

    int count;
    if(cudaGetDeviceCount(&count) != cudaSuccess) {
        std::cout << "Error: cudaGetDeviceCount" << std::endl;
        return 1;
    }

    for(size_t i = 0; i < count; i++) {

        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);

        fprintf(stdout, "Device '%s' (%d) props:\n", props.name, props.multiGpuBoardGroupID);
        fprintf(stdout, " - props.warpSize: %lu\n", props.warpSize);
        fprintf(stdout, " - props.maxThreadsPerBlock: %lu\n", props.maxThreadsPerBlock);
        fprintf(stdout, " - props.sharedMemPerBlock: %lu\n", props.sharedMemPerBlock);
        fprintf(stdout, " - props.totalConstMem: %lu\n", props.totalConstMem);

    }

#endif



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


    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

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




    float* d_total_emitted_lava[cuda_devices];
    float* d_total_current_lava[cuda_devices];

    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

        if(cudaMalloc(&d_total_emitted_lava[d], sizeof(float)) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'total_emitted_lava'" << std::endl;
            return 1;
        }

        if(cudaMalloc(&d_total_current_lava[d], sizeof(float)) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'total_current_lava'" << std::endl;
            return 1;
        }

        if(cudaMemset(d_total_emitted_lava[d], 0, sizeof(float)) != cudaSuccess) {
            std::cout << "Error: cudaMemset 'total_emitted_lava'" << std::endl;
            return 1;
        }

    }


    float total_emitted_lava =  0.0f;
    float total_current_lava = -1.0f;






    double* Sh[cuda_devices];
    double* Sz[cuda_devices];
    double* ST[cuda_devices];
    double* Mf[cuda_devices];
    double* Mhs[cuda_devices];
    bool*   Mb[cuda_devices];
    

    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

        if(cudaMallocManaged(&Sh[d], sizeof(double) * (M + 2) * N / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'Sh'" << std::endl;
            return 1;
        }

        if(cudaMallocManaged(&Sz[d], sizeof(double) * (M + 2) * N / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'Sz'" << std::endl;
            return 1;
        }

        if(cudaMallocManaged(&ST[d], sizeof(double) * (M + 2) * N / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'ST'" << std::endl;
            return 1;
        }

        if(cudaMallocManaged(&Mf[d], sizeof(double) * (M + 2) * N * (NUMBER_OF_OUTFLOWS + 2) / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'Mf'" << std::endl;
            return 1;
        }

        if(cudaMallocManaged(&Mhs[d], sizeof(double) * (M + 2) * N / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'Mhs'" << std::endl;
            return 1;
        }

        if(cudaMallocManaged(&Mb[d], sizeof(bool) * (M + 2) * N / cuda_devices) != cudaSuccess) {
            std::cout << "Error: cudaMallocManaged 'Mb'" << std::endl;
            return 1;
        }



        const size_t stride = M * N / cuda_devices;

        if(d == 0) {

            memcpy(&Sz [d][N],  &sciara->substates->Sz [stride * d],  sizeof(double) * stride + N);
            memcpy(&ST [d][N],  &sciara->substates->ST [stride * d],  sizeof(double) * stride + N);
            memcpy(&Sh [d][N],  &sciara->substates->Sh [stride * d],  sizeof(double) * stride + N);
            memcpy(&Mhs[d][N], &sciara->substates->Mhs [stride * d],  sizeof(double) * stride + N);
            memcpy(&Mb [d][N],  &sciara->substates->Mb [stride * d],  sizeof(bool)   * stride + N);

        } else if(d == cuda_devices - 1) {

            memcpy(&Sz [d][0],  &sciara->substates->Sz [stride * d - N],  sizeof(double) * stride);
            memcpy(&ST [d][0],  &sciara->substates->ST [stride * d - N],  sizeof(double) * stride);
            memcpy(&Sh [d][0],  &sciara->substates->Sh [stride * d - N],  sizeof(double) * stride);
            memcpy(&Mhs[d][0], &sciara->substates->Mhs [stride * d - N],  sizeof(double) * stride);
            memcpy(&Mb [d][0],  &sciara->substates->Mb [stride * d - N],  sizeof(bool)   * stride);

        } else {

            memcpy(&Sz [d][0],  &sciara->substates->Sz [stride * d - N],  sizeof(double) * stride + N);
            memcpy(&ST [d][0],  &sciara->substates->ST [stride * d - N],  sizeof(double) * stride + N);
            memcpy(&Sh [d][0],  &sciara->substates->Sh [stride * d - N],  sizeof(double) * stride + N);
            memcpy(&Mhs[d][0], &sciara->substates->Mhs [stride * d - N],  sizeof(double) * stride + N);
            memcpy(&Mb [d][0],  &sciara->substates->Mb [stride * d - N],  sizeof(bool)   * stride + N);

        }

        memset(Mf[d], 0, sizeof(double) * (M + 2) * N * (NUMBER_OF_OUTFLOWS + 2) / cuda_devices);


    }



    double* ST_halo[cuda_devices];

    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

        if(cudaMallocManaged(&ST_halo[d], sizeof(double) * M * N) != cudaSuccess) {
            std::cout << "Error: cudaMalloc 'ST_halo'" << std::endl;
            return 1;
        }

    }


    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

        for(int k = 0; k < cuda_devices; k++) {

            if(d != k) {

                int can;

                if(cudaDeviceCanAccessPeer(&can, d, k) != cudaSuccess) {
                    std::cout << "Error: cudaDeviceCanAccessPeer" << std::endl;
                    return 1;
                }

                if(can) {
                    if(cudaDeviceEnablePeerAccess(k, 0) != cudaSuccess) {
                        std::cout << "Error: cudaDeviceEnablePeerAccess" << std::endl;
                        return 1;
                    }
                }

            }

        }

    }



    util::Timer cl_timer;

    while((total_current_lava == -1.0f || total_current_lava > THICKNESS) || (sciara->simulation->elapsed_time <= sciara->simulation->effusion_duration)) {

        sciara->simulation->elapsed_time += sciara->parameters->Pclock;
        sciara->simulation->step++;


        cudaSetDevice(0);

        #pragma omp parallel for num_threads(cuda_devices)
        for(int d = 0; d < cuda_devices; d++) {

            cudaSetDevice(d);

            emit_lava<<<grid, threads>>> (
                M / cuda_devices, N,
                sciara->simulation->elapsed_time,
                sciara->parameters->Pclock,
                sciara->simulation->emission_time,
                sciara->parameters->Pac,
                sciara->parameters->PTvent,
                &Sh[d][N],
                &ST[d][N],
                d_total_emitted_lava[d]
            );

        }

        #pragma omp parallel for num_threads(cuda_devices)
        for(int d = 0; d < cuda_devices; d++) {

            cudaSetDevice(d);

            if(d == 0) {

                cudaMemcpyPeerAsync(&Sz[d][M * N - N], d, &Sz[d + 1][N], d + 1, sizeof(double) * N);
                cudaMemcpyPeerAsync(&Sh[d][M * N - N], d, &Sh[d + 1][N], d + 1, sizeof(double) * N);

            } else if(d == cuda_devices - 1) {

                cudaMemcpyPeerAsync(&Sz[d][0], d, &Sz[d - 1][M * N - N], d - 1, sizeof(double) * N);
                cudaMemcpyPeerAsync(&Sh[d][0], d, &Sh[d - 1][M * N - N], d - 1, sizeof(double) * N);

            } else {

                cudaMemcpyPeerAsync(&Sz[d][0],         d, &Sz[d - 1][M * N - N], d - 1, sizeof(double) * N);
                cudaMemcpyPeerAsync(&Sz[d][M * N - N], d, &Sz[d + 1][N],         d + 1, sizeof(double) * N);

                cudaMemcpyPeerAsync(&Sh[d][0],         d, &Sh[d - 1][M * N - N], d - 1, sizeof(double) * N);
                cudaMemcpyPeerAsync(&Sh[d][M * N - N], d, &Sh[d + 1][N],         d + 1, sizeof(double) * N);

            }

            cudaStreamSynchronize(0);

        }

        #pragma omp parallel for num_threads(cuda_devices)
        for(int d = 0; d < cuda_devices; d++) {

            compute_outflows<<<grid, threads, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double) * 2>>> (
                M / cuda_devices, N,
                Sz[d],
                Sh[d],
                ST[d],
                Mf[d],
                sciara->parameters->Pc,
                sciara->parameters->a,
                sciara->parameters->b,
                sciara->parameters->c,
                sciara->parameters->d,
                THREADS_PER_BLOCK
            );

        }

        

        #pragma omp parallel for num_threads(cuda_devices)
        for(int d = 0; d < cuda_devices; d++) {

            cudaSetDevice(d);

            prepare_halo<<<grid, threads>>> (ST_halo[d], ST[d], M / cuda_devices, N);

            mass_balance<<<grid, threads, THREADS_PER_BLOCK * THREADS_PER_BLOCK * sizeof(double)>>> (
                M / cuda_devices, N,
                &Sh[d][N],
                &ST[d][N],
                ST_halo[d],
                Mf[d][N],
                THREADS_PER_BLOCK
            );

        }


        #pragma omp parallel for num_threads(cuda_devices)
        for(int d = 0; d < cuda_devices; d++) {

            cudaSetDevice(d);

            compute_new_temperature_and_solidification<<<grid, threads>>> (
                M / cuda_devices, N,
                sciara->parameters->Pepsilon,
                sciara->parameters->Psigma,
                sciara->parameters->Pclock,
                sciara->parameters->Pcool,
                sciara->parameters->Prho,
                sciara->parameters->Pcv,
                sciara->parameters->Pac,
                sciara->parameters->PTsol,
                Sz[d],
                Sh[d],
                ST[d],
                Mf[d],
                Mhs[d],
                Mb[d]
            );

        }


        if(sciara->simulation->step % REDUCE == 0) {

            float __total_current_lava[cuda_devices] = {};

            #pragma omp parallel for num_threads(cuda_devices)
            for(int d = 0; d < cuda_devices; d++) {

                cudaSetDevice(d);
                cudaMemset(d_total_current_lava[d], 0, sizeof(float));

                reduce_add<<<grid, threads>>>(
                    M / cuda_devices, N,
                    Sh[d],
                    d_total_current_lava[d]
                );

                cudaDeviceSynchronize();
                cudaMemcpy(&__total_current_lava[d], d_total_current_lava[d], sizeof(float), cudaMemcpyDeviceToHost);

            }

            for(int d = 0; d < cuda_devices; d++) {
                total_current_lava += __total_current_lava[d];
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

    #pragma omp parallel for num_threads(cuda_devices)
    for(int d = 0; d < cuda_devices; d++) {

        cudaSetDevice(d);

        cudaDeviceSynchronize();
        cudaMemcpy(&total_emitted_lava, d_total_emitted_lava[d], sizeof(float), cudaMemcpyDeviceToHost);

    }


    double cl_time = static_cast<double>(cl_timer.getTimeMilliseconds()) / 1000.0;

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