#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;

//VIRAJ

particle_t **binParticles = NULL;
unsigned char **binParticlesFlag = NULL;
int *particlesPerBin = NULL;
int *freeLocationPerBin = NULL;

int numBins = -1;
int maxParticlesPerBin = -1;
double binLength = -1;

int isCloseToEdge(particle_t &particle, double binEdge)
{
	double dx = binEdge - particle.x;
	double r2 = dx * dx;
	if( r2 > cutoff*cutoff )
		return 0;
	return 1;
}

void setupBins()
{
	// setup bins
	numBins = n_threads;
	binLength = size / numBins;

	// find maximum no of particles per bin
	double bin_area = (size*size) / numBins; // area of space = size * size 
	maxParticlesPerBin = (int)( bin_area / (3.14 * (cutoff/2) * (cutoff/2)) ); // radius of particle = cutoff/2

	// allocate memory
	particlesPerBin = (int*) malloc(sizeof(int) * numBins);
	freeLocationPerBin = (int*) malloc(sizeof(int) * numBins);
	binParticles = (particle_t**) malloc(sizeof(particle_t*) * numBins);
	binParticlesFlag = (unsigned char**) malloc(sizeof(unsigned char*) * numBins);
	for(int i=0; i<numBins; i++)
	{
		particlesPerBin[i] = 0;
		freeLocationPerBin[i] = 0;
		binParticles[i] = (particle_t*) malloc(sizeof(particle_t) * maxParticlesPerBin);
		memset(binParticles[i], 0, sizeof(particle_t) * maxParticlesPerBin);
		binParticlesFlag[i] = (unsigned char*) malloc(sizeof(unsigned char) * maxParticlesPerBin);
		memset(binParticlesFlag[i], 0, sizeof(unsigned char) * maxParticlesPerBin);
	}
}

void freeBins()
{
	for(int i=0; i<numBins; i++)
	{
		free(binParticles[i]);
		free(binParticlesFlag[i]);
	}
	free(binParticles);
	free(binParticlesFlag);
	
	free(particlesPerBin);
	free(freeLocationPerBin);
}

int getFreeLocation(int bin)
{
	int loc = freeLocationPerBin[bin];
	if(binParticlesFlag[bin][loc] != 0)
	{
		// Error check
		printf("\ngetFreeLocation(): RED FLAG! @ binParticlesFlag[bin][loc] != 0")
	}	
}

void copyParticleToBin(particle_t particle, int bin)
{
	if(particlesPerBin[bin] > maxParticlesPerBin)
		printf("\nputParticleData(): RED FLAG! THIS SHOULD NEVER HAPPEN");
	
	// copy particle data into binParticles[][] and update everything
	
	int loc = getFreeLocation(bin); // get index of first free location

	binParticles[bin][loc] = particle; // copy data
	binParticlesFlag[bin][loc] = 1; // set presence flag
	particlesPerBin[bin]++; // increment particles per bin
	
	// update freeLocationPerBin[bin]
	int i = freeLocationPerBin[bin];
	
	while(binParticlesFlag[bin][i] != 0)
	{
		i++; // search for next free location
		if(i == maxParticlesPerBin) 
		{
			i=0; // reset if pointer reaches end of bin
			printf("\putParticleData(): Pointer reached end of bin! That's surprising!\n");
		}
		
		if(i == loc)
		{
			printf("\nputParticleData(): RED FLAG! THIS SHOULD NEVER HAPPEN");
		}
	}
}

void doBinning()
{
	for(int ndx=0; ndx<n; ndx++) // for each particle
	{
		int bin = (particles[ndx].x) / binLength;

		// copy particle data into bin
		copyParticleToBin(particles[ndx], bin);
	}
	
	// Error check loop
	// Check if number of particles in bin sums up
	// REMOVE this code when timing
	int sum = 0;
	for(int bin=0; bin<numBins; bin++)
	{
		sum += particlesPerBin[bin];
	}
	if(sum != n) printf("\ndoBinning(): RED FLAG! Sums don't match\n");
}

//!VIRAJ

//
//  check that pthreads routine call was successful
//

#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        for( int i = first; i < last; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j] );
        }
        
        pthread_barrier_wait( &barrier );
        
        //
        //  move particles
        //
        for( int i = first; i < last; i++ ) 
            move( particles[i] );
        
        pthread_barrier_wait( &barrier );
        
        //
        //  save if necessary
        //
        if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

	// VIRAJ
	doBinning();
	// !VIRAJ

	// create threads
    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
