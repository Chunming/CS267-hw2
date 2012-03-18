#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define EMPTY (-1)

int nthreads = 256; 
int numBinsX = 10;
int numBinsY = 10;
int numBins = 100; // Should be perfect square! TODO: change this
int n = 500; 
int maxParticlesPerBin;

particle_t *particles;
particle_t *d_particles = NULL;
int *binParticlesIds = NULL;
int *d_binParticlesIds = NULL;
int *freeLocationPerBin = NULL;
int *d_freeLocationPerBin = NULL;
int *particlesPerBin = NULL;
int *d_particlesPerBin = NULL;



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

__global__ void compute_forces_gpu(int n, particle_t *d_particles, int *d_binParticlesIds, int *d_particlesPerBin, int *d_freeLocationPerBin, int numBinsX, int numBinsY, int maxParticlesPerBin)
{
  // Get thread (particle) ID
	int bid = blockIdx.x; // block idx
  	int tid = threadIdx.x ; // thread idx

	if(bid > numBinsX*numBinsY) return;

  if(tid >= d_particlesPerBin[bid]) return;

	int pid = d_binParticlesIds[bid*maxParticlesPerBin + tid];

  d_particles[pid].ax = d_particles[pid].ay = 0;

	// apply force on each particle within the bin
  for(int j = 0 ; j < d_particlesPerBin[bid] ; j++) 
	{
		int jpid = d_binParticlesIds[bid*maxParticlesPerBin + j];
    	apply_force_gpu(d_particles[pid], particles[jpid]);
	}
	
	// apply force on each particle within the bin to the WEST
	int bdx = bid - 1;
	if(bid%numBinsX != 0) // if not on left edge of space
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin to the EAST
	int bdx = bid + 1;
	if((bid+1)%numBinsX != 0) // if not on right edge
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin NORTH
	int bdx = bid - numBinsX;
	if(bid >= numBinsX) // if bid is not in first row
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin SOUTH
	int bdx = bid + numBinsX;
	if(bdx < (numBins - numBinsX))
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin NE
	int bdx = bid + numBinsX;
	if(bid >= numBinsX && (bid+1)%numBinsX != 0)
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin NW
	int bdx = bid + numBinsX;
	if(bid >= numBinsX && bid%numBinsX != 0)
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin SE
	int bdx = bid + numBinsX;
	if(bdx < (numBins - numBinsX) && (bid+1)%numBinsX != 0)
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}
	
	// apply force on each particle within the bin SW
	int bdx = bid + numBinsX;
	if(bdx < (numBins - numBinsX) && bid%numBinsX != 0)
	{
	  for(int j = 0 ; j < d_particlesPerBin[bdx] ; j++) 
		{
			int jpid = d_binParticlesIds[bdx*maxParticlesPerBin + j];
	    	apply_force_gpu(d_particles[pid], particles[jpid]);
		}
	}		
}



__global__ void move_gpu (int n, double size, particle_t *d_particles, int *d_binParticlesIds, int *d_particlesPerBin, int *d_freeLocationPerBin, int numBinsX, int numBinsY, int maxParticlesPerBin)
{

  // Get thread (particle) ID within blocks
	int bid = blockIdx.x;
	int tid = threadIdx.x;
	
	if(bid > numBinsX*numBinsY) return;

  	if(tid >= d_particlesPerBin[bid]) return;

	int pid = d_binParticlesIds[bid*maxParticlesPerBin + tid];

  particle_t *p = &d_particles[pid];
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

__device__ transferParticle(int pid, int fromBdx, int toBdx)
{
	
}

void doBinning()
{
	numBins = numBinsX * numBinsY;

	double binWidth = sqrt(density * n);
	
	// find maximum no of particles per bin
	double bin_area = (binWidth*binWidth) / numBins; // area of space = size * size 
	maxParticlesPerBin = 3 * (int)( bin_area / (3.14 * (cutoff/2) * (cutoff/2)) ); // radius of particle = cutoff/2
	
	binParticlesIds = (int*) malloc(numBins * sizeof(int) * maxParticlesPerBin);
	if(binParticlesIds == NULL) printf("\ndoBinning(): malloc failed\n");
	
	cudaMalloc((void **) &d_binParticlesIds, numBins * sizeof(int) * maxParticlesPerBin);
	if(d_binParticlesIds == NULL) printf("\ndoBinning(): cudamalloc() failed\n");
	
	freeLocationPerBin = (int*) malloc(numBins * sizeof(int));
	if(freeLocationPerBin == NULL) printf("\ndoBinning(): malloc failed\n");
	
	cudaMalloc((void **) &d_freeLocationPerBin, numBins * sizeof(int));
	if(d_freeLocationPerBin == NULL) printf("\ndoBinning(): cudamalloc() failed\n");
	
	particlesPerBin = (int*) malloc(numBins * sizeof(int));
	if(particlesPerBin == NULL) printf("\ndoBinning(): malloc failed\n");
	
	cudaMalloc((void **) &d_particlesPerBin, numBins * sizeof(int));
	if(d_particlesPerBin == NULL) printf("\ndoBinning(): cudamalloc() failed\n");
	
	for(int i=0; i<maxParticlesPerBin; i++)
	{
		// set all particles ids to EMPTY
		binParticles[i] = EMPTY;
	}

	for(int i=0; i<numBins; i++)
	{
		freeLocationPerBin[i] = 0;
		particlesPerBin[i] = 0;
	}	
	
	for(int ndx=0; ndx<n; ndx++) // for each particle
	{
		int binx = (particles[ndx].x) / binWidth;
		int biny = (particles[ndx].y) / binWidth;
		int bdx = biny*numBinsX + binx;
		
		// add particle ndx to bin bdx
		int loc = freeLocationPerBin[bdx];
		freeLocationPerBin[bdx]++;
		particlesPerBin[bdx]++;
		binParticles[bdx*maxParticlesPerBin + loc] = particles[ndx];
	}
}

int copyToDevice()
{
	cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_binParticlesIds, binParticlesIds, numBins * sizeof(int) * maxParticlesPerBin, cudaMemcpyHostToDevice);
	cudaMemcpy(d_freeLocationPerBin, freeLocationPerBin, numBins * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_particlesPerBin, particlesPerBin, numBins * sizeof(int), cudaMemcpyHostToDevice);
}

void freeBins()
{
	free(binParticlesIds);
	cudaFree(d_binParticlesIds);
	
	free(freeLocationPerBin);
	cudaFree(d_freeLocationPerBin);
	
	free(particlesPerBin);
	cudaFree(d_particlesPerBin);
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
    
    n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );

    init_particles( n, particles );

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

	int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);
        
        //
        //  move particles
        //
	move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        
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
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
