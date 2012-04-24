#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 25


extern double size;
//
//  benchmarking program
//


__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

 __device__ inline void myAtomicAdd(double *address, double value)  //See CUDA official forum
{
    unsigned long long oldval, newval, readback;

    oldval = __double_as_longlong(*address);
    newval = __double_as_longlong(__longlong_as_double(oldval) + value);
    while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
    {
        oldval = readback;
        newval = __double_as_longlong(__longlong_as_double(oldval) + value);
    }
}

__global__ void compute_forces_gpu(particle_t ** bins, int * binParticleNum, double blks_len)
{

  int tId = threadIdx.x;
  int bIdx = blockIdx.x;
  int bIdy = blockIdx.y;
  int blockId = gridDim.x * bIdx + bIdy;
  int blks_num = gridDim.x;



  if (tId >= binParticleNum[blockId]) return;


  // Check all particles in bth subBlock

  for (int j=0; j<binParticleNum[blockId]; ++j) { // The jth particle in bth bin
       apply_force_gpu(*bins[blockId*NUM_THREADS + tId], *bins[blockId*NUM_THREADS + j]);
  }

  // Compute Forces
  double leftBnd, rightBnd, topBnd, botBnd;
  double leftDist, rightDist, topDist, botDist;
  int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;

  //printf("botBnd is %f \n", botBnd);
  leftBnd = bIdx*blks_len;
  rightBnd = (bIdx*blks_len) + blks_len;
  topBnd = bIdy*blks_len;
  botBnd = bIdy*blks_len + blks_len;

  //printf("botDist is %f \n", botDist);
  leftDist = fabs((bins[blockId*NUM_THREADS + tId]->x) - leftBnd);
  rightDist = fabs((bins[blockId*NUM_THREADS + tId]->x) - rightBnd);
  topDist = fabs((bins[blockId*NUM_THREADS + tId]->y) - topBnd);
  botDist = fabs((bins[blockId*NUM_THREADS + tId]->y) - botBnd);

  // Consider 8 different adjacent subBlocks
  if (leftDist<=cutoff) {
      if (bIdx != 0) { // Left subBlock index is valid
	  bLeft = blockId - blks_num;
	  for (int k=0; k<binParticleNum[bLeft]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bLeft*NUM_THREADS + k]);
          }
      }
  }

  if (rightDist<=cutoff) {
      if (bIdx != blks_num-1) { 
          bRight = blockId + blks_num; 
	  for (int k=0; k<binParticleNum[bRight]; ++k) { 
	       apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bRight*NUM_THREADS + k]);
          }
      }
  }

  if (topDist<=cutoff) {
      if (bIdy != 0) { 
          bTop = blockId - 1; 
          for (int k=0; k<binParticleNum[bTop]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bTop*NUM_THREADS + k]);
          }
      }
  }

  if (botDist<=cutoff) {
      if (bIdy != blks_num-1) {
          bBottom = blockId + 1;
          for (int k=0; k<binParticleNum[bBottom]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bBottom*NUM_THREADS + k]);
          }
      }
  }

  if (topDist<=cutoff && leftDist<=cutoff) { 
      if (bIdy != 0 && bIdx !=0) {
          bTopLeft = blockId-blks_num-1;     
          for (int k=0; k<binParticleNum[bTopLeft]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bTopLeft*NUM_THREADS + k]);
          }
      }
  }

  if (botDist<=cutoff && leftDist<=cutoff) { 
      if (bIdy != blks_num-1 && bIdx != 0) { 
          bBotLeft = blockId-blks_num+1;     
          for (int k=0; k<binParticleNum[bBotLeft]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bBotLeft*NUM_THREADS + k]);
          }
      }
  }

  if (topDist<=cutoff && rightDist<=cutoff) { 
      if (bIdy != 0 && bIdx != blks_num-1) {
          bTopRight = blockId+blks_num-1;
          for (int k=0; k<binParticleNum[bTopRight]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bTopRight*NUM_THREADS + k]);
          }
      }
  }

	    
  if (botDist<=cutoff && rightDist<=cutoff) { 
      if (bIdy!=blks_num-1 && bIdx!=blks_num-1) {
          bBotRight = blockId+blks_num+1;     
          for (int k=0; k<binParticleNum[bBotRight]; ++k) { 
               apply_force_gpu(*bins[blockId*NUM_THREADS + tId],*bins[bBotRight*NUM_THREADS + k]);
          }
      }
  }


}

__global__ void move_gpu (particle_t ** bins, int * part_num, double size)
{

  int tId = threadIdx.x;
  int bIdx = blockIdx.x;
  int bIdy = blockIdx.y;
  int blockId = gridDim.x * bIdx + bIdy;


  if (tId >= part_num[blockId]) return;

  particle_t * p = bins[blockId*NUM_THREADS + tId];

  //
  //  slightly simplified Velocity Verlet integration
  //  conserves energy better than explicit Euler method
  //
  p->vx += p->ax * dt;
  p->vy += p->ay * dt;
  p->x  += p->vx * dt;
  p->y  += p->vy * dt;

  //
  //  bounce from walls
  //
  while( p->x < 0 || p->x > size )
  {
      p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
      p->vx = -(p->vx);
  }
  while( p->y < 0 || p->y > size )
  {
      p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
      p->vy = -(p->vy);
  }

}

__global__ void rebin_gpu (particle_t *particles, particle_t **bins, int *part_nums, int n, double blks_size, int blks_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (tid >= n) return;

	int x_bin = particles[tid].x/blks_size;
	int y_bin = particles[tid].y/blks_size;
	int blockId = blks_num * x_bin + y_bin;
	int pos = atomicAdd( part_nums + blockId, 1);

	bins[blockId * NUM_THREADS + pos] = particles + tid;
	particles[tid].ax = particles[tid].ay = 0;

}


int main( int argc, char **argv )
{    

    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 


    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", "gpu.txt" );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    set_size( n );

    //creates the correct number of blocks given fixed thread size
    int blks_num = (int) ceil(size / sqrt((NUM_THREADS)*(3.14159)*cutoff*cutoff));
    double blks_size =  size / ((double) blks_num);
    dim3 blks (blks_num, blks_num);
 
    int bin_blks = (n + NUM_THREADS - 1) / NUM_THREADS;

    // GPU particle data structure
    particle_t * d_particles;
    particle_t ** d_blocks;
    int * d_blk_part_num;

    // cudaMalloc shit
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));
    cudaMalloc((void ***) &d_blocks, blks_num * blks_num * NUM_THREADS * sizeof(particle_t*));
    cudaMalloc((void **) &d_blk_part_num, blks_num * blks_num * sizeof(int));


    //initialize shit
    cudaMemset(d_blk_part_num, 0, blks_num * blks_num * sizeof(int));
    init_particles( n, particles );
    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();

    //bin the particles
    rebin_gpu <<< bin_blks, NUM_THREADS >>> (d_particles, d_blocks, d_blk_part_num, n, blks_size, blks_num);

    //calculate time to bin particles
    cudaThreadSynchronize();
    copy_time = read_timer() - copy_time;
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);

    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    

    for( int step = 0; step < NSTEPS; step++ )
    {

        //
        //rebin particles
        //
	cudaMemset(d_blk_part_num, 0, blks_num * blks_num * sizeof(int));
        rebin_gpu <<< bin_blks, NUM_THREADS >>> (d_particles, d_blocks, d_blk_part_num, n, blks_size, blks_num);

	//cudaDeviceSynchronize();

        //
        //  compute forces
        //
	compute_forces_gpu <<< blks, NUM_THREADS >>> (d_blocks, d_blk_part_num, blks_size);
	//cudaThreadSynchronize();
	//compute_border_forces_gpu <<< blks, NUM_THREADS >>> (d_blocks, d_blk_part_num, blks_size);
	
        
    	//cudaThreadSynchronize();

        //
        //  move particles
        //
	move_gpu <<< blks, NUM_THREADS >>> (d_blocks, d_blk_part_num, size);

    	//cudaThreadSynchronize();
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	    // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
	}
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    

    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);


    cudaFree(d_blocks);
    cudaFree(d_blk_part_num);


    if( fsave )
        fclose( fsave );
    
    return 0;
}
