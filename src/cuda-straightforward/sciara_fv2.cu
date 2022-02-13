#include "cal2DBuffer.cuh"
#include "configurationPathLib.cuh"
#include "GISInfo.cuh"
#include "io.cuh"
#include "vent.cuh"
#include <omp.h>
#include <new>
#include "Sciara.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.cu"


#define BLKSIZE     8

// ----------------------------------------------------------------------------
// I/O parameters used to index argv[]
// ----------------------------------------------------------------------------

#define INPUT_PATH_ID          1
#define OUTPUT_PATH_ID         2
#define MAX_STEPS_ID           3
#define REDUCE_INTERVL_ID      4
#define THICKNESS_THRESHOLD_ID 5

// ----------------------------------------------------------------------------
// Read/Write access macros linearizing single/multy layer buffer 2D indices
// ----------------------------------------------------------------------------

#define SET(M, columns, i, j, value) ((M)[(((i) * (columns)) + (j))] = (value))
#define GET(M, columns, i, j) (M[(((i) * (columns)) + (j))])
#define BUF_SET(M, rows, columns, n, i, j, value) ( (M)[( ((n)*(rows)*(columns)) + ((i)*(columns)) + (j) )] = (value) )
#define BUF_GET(M, rows, columns, n, i, j) ( M[( ((n)*(rows)*(columns)) + ((i)*(columns)) + (j) )] )

// ----------------------------------------------------------------------------
// computing kernels, aka elementary processes in the XCA terminology
// ----------------------------------------------------------------------------

__device__
double emitLava(
    int i,
    int j,
    int r,
    int c,
    TVent* vent_1,
    TVent* vent_2,
    double elapsed_time,
    double Pclock,
    double emission_time,
    double Pac,
    double PTvent,
    double* Sh,
    double* Sh_next,
    double* ST_next)
{

    TVent* vent = NULL;

    if(i == vent_1->y() && j == vent_1->x())
      vent = vent_1;

    else if(i == vent_2->y() && j == vent_2->x())
      vent = vent_2;


    if (vent)
    {
      double t = vent->thickness(elapsed_time, Pclock, emission_time, Pac);

      SET(Sh_next, c, i, j, GET(Sh, c, i, j) + t);
      SET(ST_next, c, i, j, PTvent); 

      return t;

    }

    return .0;

}


__device__
void computeOutflows (
    int i, 
    int j, 
    int r, 
    int c, 
    int* Xi, 
    int* Xj, 
    double *Sz, 
    double *Sh, 
    double *ST,
    double *Mf, 
    double  Pc, 
    double  _a,
    double  _b,
    double  _c,
    double  _d)
{
  bool   eliminated[MOORE_NEIGHBORS];
  double z[MOORE_NEIGHBORS];
  double h[MOORE_NEIGHBORS];
  double H[MOORE_NEIGHBORS];
  double theta[MOORE_NEIGHBORS];
  double w[MOORE_NEIGHBORS];		//Distances between central and adjecent cells
  double Pr[MOORE_NEIGHBORS];		//Relaiation rate arraj
  bool loop;
  int counter;
  double sz0, sz, T, avg, rr, hc;


  if (GET(Sh,c,i,j) <=0)
    return;

  T  = GET(ST, c, i, j);
  rr = pow(10, _a+_b*T);
  hc = pow(10, _c+_d*T);

  for (int k = 0; k < MOORE_NEIGHBORS; k++)
  {
    sz0      = GET(Sz, c, i,       j      );
    sz       = GET(Sz, c, i+Xi[k], j+Xj[k]);
    h[k]     = GET(Sh, c, i+Xi[k], j+Xj[k]);
    w[k]     = Pc;
    Pr[k]    = rr;

    if (k < VON_NEUMANN_NEIGHBORS)
      z[k] = sz;
    else
      z[k] = sz0 - (sz0 - sz) * __drcp_rn(sqrt(2.0));
  }

  H[0] = z[0];
  theta[0] = 0;
  eliminated[0] = false;
  for (int k = 1; k < MOORE_NEIGHBORS; k++)
    if (z[0] + h[0] > z[k] + h[k])
    {
      H[k] = z[k] + h[k];
      theta[k] = atan(((z[0] + h[0]) - (z[k] + h[k])) * __drcp_rn(w[k]));
      eliminated[k] = false;
    } 
    else
    {
      //H[k] = 0;
      //theta[k] = 0;
      eliminated[k] = true;
    }

  do {
    loop = false;
    avg = h[0];
    counter = 0;
    for (int k = 0; k < MOORE_NEIGHBORS; k++)
      if (!eliminated[k])
      {
        avg += H[k];
        counter++;
      }
    if (counter != 0)
      avg = avg * __drcp_rn(double(counter));
    for (int k = 0; k < MOORE_NEIGHBORS; k++)
      if (!eliminated[k] && avg <= H[k])
      {
        eliminated[k] = true;
        loop = true;
      }
  } while (loop);

  for (int k = 1; k < MOORE_NEIGHBORS; k++) 
    if (!eliminated[k] && h[0] > hc * cos(theta[k]))
      BUF_SET(Mf,r,c,k-1,i,j, Pr[k]*(avg - H[k]));
    else
      BUF_SET(Mf,r,c,k-1,i,j,0.0);
}

__device__
void massBalance(
    int i, 
    int j, 
    int r, 
    int c, 
    int* Xi, 
    int* Xj, 
    double *Sh, 
    double *Sh_next, 
    double *ST,
    double *ST_next,
    double *Mf)
{
  const int inflowsIndices[NUMBER_OF_OUTFLOWS] = { 3, 2, 1, 0, 6, 7, 4, 5 };
  double inFlow;
  double outFlow;
  double neigh_t;
  double initial_h = GET(Sh,c,i,j);
  double initial_t = GET(ST,c,i,j);
  double h_next = initial_h;
  double t_next = initial_h * initial_t;

  for (int n = 1; n < MOORE_NEIGHBORS; n++)
  {
    neigh_t = GET(ST,c,i+Xi[n],j+Xj[n]);
    inFlow  = BUF_GET(Mf,r,c,inflowsIndices[n-1],i+Xi[n],j+Xj[n]);

    outFlow = BUF_GET(Mf,r,c,n-1,i,j);

    h_next +=  inFlow - outFlow;
    t_next += (inFlow * neigh_t - outFlow * initial_t);
  }

  if (h_next > 0)
  {
    t_next *= __drcp_rn(h_next);
    SET(ST_next,c,i,j,t_next);
    SET(Sh_next,c,i,j,h_next);
  }

}

__device__
void computeNewTemperatureAndSolidification(
    int i, 
    int j, 
    int r, 
    int c,
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
    bool   *Mb)
{
  double nT, aus;
  double z = GET(Sz,c,i,j);
  double h = GET(Sh,c,i,j);
  double T = GET(ST,c,i,j);

  if (h > 0 && GET(Mb,c,i,j) == false ) 
  {
    aus = 1.0 + (3 * pow(T, 3.0) * Pepsilon * Psigma * Pclock * Pcool) * __drcp_rn((Prho * Pcv * h * Pac));
    nT = T * __drcp_rn(pow(aus, 1.0 * __drcp_rn(3.0)));

    if (nT > PTsol) // no solidification
      SET(ST_next,c,i,j, nT);
    else            // solidification
    {
      SET(Sz_next,c,i,j,z+h);
      SET(Sh_next,c,i,j,0.0);
      SET(ST_next,c,i,j,PTsol);
      SET(Mhs,c,i,j, GET(Mhs,c,i,j)+h);
    }

  }
}


__global__
void emitLava_kernel(
    int r,
    int c,
    TVent* vent_1,
    TVent* vent_2,
    double elapsed_time,
    double Pclock,
    double emission_time,
    float* total_emitted_lava,
    double Pac,
    double PTvent,
    double* Sh,
    double* Sh_next,
    double* ST_next) {

  __shared__ double s_acc[BLKSIZE][BLKSIZE];

  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;


  if (i > 0 && j > 0 && i < r - 1 && j < c - 1) {

    s_acc[threadIdx.y][threadIdx.x] = emitLava(i, j, r, c, vent_1, vent_2, elapsed_time, Pclock, emission_time, Pac, PTvent, Sh, Sh_next, ST_next);

    __syncthreads();

    if(threadIdx.x == 0 && threadIdx.y == 0) {
      
      double v = 0.0;

      for(size_t y = 0; y < BLKSIZE; y++) {
        for(size_t x = 0; x < BLKSIZE; x++) {
          v += s_acc[y][x];
        }
      }

      atomicAdd(total_emitted_lava, (float) v);
      
    }

  } else {

    s_acc[threadIdx.y][threadIdx.x] = 0.0f;

  }

}


__global__
void massBalance_kernel(
    int r, 
    int c, 
    int* Xi, 
    int* Xj, 
    double *Sh, 
    double *Sh_next, 
    double *ST,
    double *ST_next,
    double *Mf) {
  
  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && j > 0 && i < r - 1 && j < c - 1)
    massBalance(i, j, r, c, Xi, Xj, Sh, Sh_next, ST, ST_next, Mf);

}


__global__
void computeOutflows_kernel(
    int r, 
    int c, 
    int* Xi, 
    int* Xj, 
    double *Sz, 
    double *Sh, 
    double *ST,
    double *Mf, 
    double  Pc, 
    double  _a,
    double  _b,
    double  _c,
    double  _d) {
  
  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && j > 0 && i < r - 1 && j < c - 1)
    computeOutflows(i, j, r, c, Xi, Xj, Sz, Sh, ST, Mf, Pc, _a, _b, _c, _d);

}



__global__
void computeNewTemperatureAndSolidification_kernel(
    int r, 
    int c,
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
    bool   *Mb) {
  
  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i > 0 && j > 0 && i < r - 1 && j < c - 1)
    computeNewTemperatureAndSolidification(i, j, r, c, Pepsilon, Psigma, Pclock, Pcool, Prho, Pcv, Pac, PTsol, Sz, Sz_next, Sh, Sh_next, ST, ST_next, Mf, Mhs, Mb);

}



__global__
void memcpy_gpu(double *dst, double *src, int r, int c) {

  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < r && j < c)
    dst[i * c + j] = src[i * c + j];

}

__global__
void reduceAdd_kernel(size_t r, size_t c, double* buffer, float* acc) {

  __shared__ double s_acc[BLKSIZE][BLKSIZE];

  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < r && j < c) {

    s_acc[threadIdx.x][threadIdx.y] = buffer[i * c + j];
    
    __syncthreads();

    if(threadIdx.x == 0 && threadIdx.y == 0) {
      
      double v = 0.0;

      for(size_t y = 0; y < BLKSIZE; y++) {
        for(size_t x = 0; x < BLKSIZE; x++) {
          v += s_acc[y][x];
        }
      }

      atomicAdd(acc, (float) v);
      
    }

  } else {

    s_acc[threadIdx.y][threadIdx.x] = 0.0;

  }

}



// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  Sciara *sciara;
  init(sciara);

  // Input data 
  int max_steps = 1000;
  loadConfiguration("./data/2006/2006_000000000000.cfg", sciara);


  size_t M = sciara->domain->rows;
  size_t N = sciara->domain->cols;
  size_t THREADS_PER_BLOCK = BLKSIZE;
  size_t GRID_STRIDE_SIZE = 1;

  dim3 threads(THREADS_PER_BLOCK, THREADS_PER_BLOCK);

  dim3 grid (
      fmax(1, ceil(double(M) / THREADS_PER_BLOCK / GRID_STRIDE_SIZE)),
      fmax(1, ceil(double(N) / THREADS_PER_BLOCK / GRID_STRIDE_SIZE))
  );




  // simulation initialization and loop
  float total_current_lava = -1.0f;
  simulationInitialize(sciara);


  TVent* h_vent = new TVent[sciara->simulation->vent.size()];

  for (int i = 0; i < sciara->simulation->vent.size(); i++)
    memcpy(&h_vent[i], &sciara->simulation->vent[i], sizeof(TVent));


  TVent* vent_1;
  TVent* vent_2;

  cudaMalloc(&vent_1, sizeof(TVent));
  cudaMalloc(&vent_2, sizeof(TVent));

  cudaMemcpy(vent_1, &h_vent[0], sizeof(TVent), cudaMemcpyHostToDevice);
  cudaMemcpy(vent_2, &h_vent[1], sizeof(TVent), cudaMemcpyHostToDevice);

  size_t vent_size = sciara->simulation->vent.size();


  float* total_emitted_lava;
  cudaMalloc(&total_emitted_lava, sizeof(float));

  float* d_total_current_lava;
  cudaMalloc(&d_total_current_lava, sizeof(float));



  util::Timer cl_timer;

  int reduceInterval = 1000;
  double thickness_threshold = 1.0;



  while (    (max_steps > 0 && sciara->simulation->step < max_steps) 
          || (sciara->simulation->elapsed_time <= sciara->simulation->effusion_duration)
          || (total_current_lava == -1 ||  total_current_lava > thickness_threshold) 
        )
  {

    sciara->simulation->elapsed_time += sciara->parameters->Pclock;
    sciara->simulation->step++;

    if(sciara->simulation->step % 100 == 0) {
      fprintf(stderr, "\rStep: %d (%d %%) [%f / %f m]", 
        sciara->simulation->step, 
        int(100.0 * sciara->simulation->elapsed_time / sciara->simulation->effusion_duration), 
        total_current_lava, 
        thickness_threshold);
    }

    // Apply the emitLava kernel to the whole domain and update the Sh and ST state variables

    cudaMemcpy(total_emitted_lava, &sciara->simulation->total_emitted_lava, sizeof(float), cudaMemcpyHostToDevice);

    emitLava_kernel<<<grid, threads>>>( 
        sciara->domain->rows, 
        sciara->domain->cols, 
        vent_1,
        vent_2, 
        sciara->simulation->elapsed_time, 
        sciara->parameters->Pclock, 
        sciara->simulation->emission_time, 
        total_emitted_lava,
        sciara->parameters->Pac, 
        sciara->parameters->PTvent, 
        sciara->substates->Sh, 
        sciara->substates->Sh_next,
        sciara->substates->ST_next);

    cudaDeviceSynchronize();

    memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, M, N);
    memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, M, N);

    cudaMemcpy(&sciara->simulation->total_emitted_lava, total_emitted_lava, sizeof(float), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();



    computeOutflows_kernel<<<grid, threads>>>(
        sciara->domain->rows, 
        sciara->domain->cols, 
        sciara->X->Xi, 
        sciara->X->Xj, 
        sciara->substates->Sz, 
        sciara->substates->Sh, 
        sciara->substates->ST, 
        sciara->substates->Mf, 
        sciara->parameters->Pc, 
        sciara->parameters->a, 
        sciara->parameters->b, 
        sciara->parameters->c, 
        sciara->parameters->d);


    cudaDeviceSynchronize();

    // Apply the massBalance mass balance kernel to the whole domain and update the Sh and ST state variables
    massBalance_kernel<<<grid, threads>>>(
        sciara->domain->rows, 
        sciara->domain->cols, 
        sciara->X->Xi, 
        sciara->X->Xj, 
        sciara->substates->Sh, 
        sciara->substates->Sh_next, 
        sciara->substates->ST, 
        sciara->substates->ST_next, 
        sciara->substates->Mf);


    cudaDeviceSynchronize();

    memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, sciara->domain->rows, sciara->domain->cols);
    memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, sciara->domain->rows, sciara->domain->cols);

    cudaDeviceSynchronize();


    computeNewTemperatureAndSolidification_kernel<<<grid, threads>>> (
      sciara->domain->rows, 
      sciara->domain->cols, 
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
      sciara->substates->Mb);

    cudaDeviceSynchronize();


    memcpy_gpu<<<grid, threads>>>(sciara->substates->Sz, sciara->substates->Sz_next, sciara->domain->rows, sciara->domain->cols);
    memcpy_gpu<<<grid, threads>>>(sciara->substates->Sh, sciara->substates->Sh_next, sciara->domain->rows, sciara->domain->cols);
    memcpy_gpu<<<grid, threads>>>(sciara->substates->ST, sciara->substates->ST_next, sciara->domain->rows, sciara->domain->cols);

    cudaDeviceSynchronize();

    // Global reduction
    if (sciara->simulation->step % reduceInterval == 0) {
     
      cudaMemset(d_total_current_lava, 0, sizeof(float));
      reduceAdd_kernel<<<grid, threads>>>(sciara->domain->rows, sciara->domain->cols, sciara->substates->Sh, d_total_current_lava);
      cudaDeviceSynchronize();
      cudaMemcpy(&total_current_lava, d_total_current_lava, sizeof(float), cudaMemcpyDeviceToHost);

    }

    // if(cudaGetLastError() != cudaSuccess) {
    //   fprintf(stderr, "Error: %s in step %d\n", cudaGetErrorString(cudaGetLastError()), sciara->simulation->step);
    //   exit(1);
    // }

  }

  double cl_time = static_cast<double>(cl_timer.getTimeMilliseconds()) / 1000.0;
  printf("Step %d\n", sciara->simulation->step);
  printf("Elapsed time [s]: %lf\n", cl_time);
  printf("Emitted lava [m]: %lf\n", sciara->simulation->total_emitted_lava);
  printf("Current lava [m]: %lf\n", total_current_lava);

  printf("Saving output to %s...\n", "./data/2006/output_2006");
  saveConfiguration("./data/2006/output_2006", sciara);

  printf("Releasing memory...\n");
  finalize(sciara);

  return 0;
}
