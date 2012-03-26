#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include <string.h>

#define __VIRAJ_DEBUG
					
//
//  global variables
//

int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;

//VIRAJ

// Duplicating from common.cpp
#define cutoff  0.01
#define density 0.0005

double spaceDim = 0.0;
particle_t **binParticles = NULL;
unsigned char **binParticlesFlag = NULL;
int *particlesPerBin = NULL;
int *freeLocationPerBin = NULL;
pthread_mutex_t *binMutexes = NULL;

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
	spaceDim = sqrt(density * n);
	
	// setup bins
	numBins = n_threads;
	binLength = spaceDim / numBins;
	
#ifdef __VIRAJ_DEBUG
	printf("\nsetupBins(): Entered\n");
	printf("\nBox Size = %f\n", binLength);
	
	if(binLength <= cutoff)
	{
		printf("\nsetupBins(): RED FLAG! Box size is smaller than cutoff!\n");
	}
#endif

	// find maximum no of particles per bin
	double bin_area = (spaceDim*spaceDim) / numBins; // area of space = size * size 
	maxParticlesPerBin = 3 * (int)( bin_area / (3.14 * (cutoff/2) * (cutoff/2)) ); // radius of particle = cutoff/2

	// allocate memory
	particlesPerBin = (int*) malloc(sizeof(int) * numBins);
	if(particlesPerBin == NULL) printf("\nsetupBins(): RED FLAG! ALLOC Failure at particlesPerBin\n");
	
	freeLocationPerBin = (int*) malloc(sizeof(int) * numBins);
	if(freeLocationPerBin == NULL) printf("\nsetupBins(): RED FLAG! ALLOC Failure at freeLocationPerBin\n");
		
	binParticles = (particle_t**) malloc(sizeof(particle_t*) * numBins);
	if(binParticles == NULL) printf("\nsetupBins(): RED FLAG! ALLOC Failure at binParticles\n");
	
	binParticlesFlag = (unsigned char**) malloc(sizeof(unsigned char*) * numBins);
	if(binParticlesFlag == NULL) printf("\nsetupBins(): RED FLAG! ALLOC Failure at binParticlesFlag\n");

	for(int i=0; i<numBins; i++) // for each bin
	{
		particlesPerBin[i] = 0;
		freeLocationPerBin[i] = 0;
		
		binParticles[i] = (particle_t*) malloc(sizeof(particle_t) * maxParticlesPerBin);
		if(binParticles[i] == NULL) printf("\nsetupBins(): RED FLAG!  ALLOC Failure at binParticles[i]\n");
		memset(binParticles[i], 0, sizeof(particle_t) * maxParticlesPerBin);
		
		binParticlesFlag[i] = (unsigned char*) malloc(sizeof(unsigned char) * maxParticlesPerBin);
		if(binParticlesFlag[i] == NULL) printf("\nsetupBins(): RED FLAG! ALLOC Failure at binParticlesFlag[i]\n");
		memset(binParticlesFlag[i], 0, sizeof(unsigned char) * maxParticlesPerBin);
	}
	
	// setup mutexes
	binMutexes = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * numBins);
	if(binMutexes == NULL) printf("\nsetupBins():RED FLAG! Mutex ALLOC Failure\n");
	for(int i=0; i<numBins; i++) // for each bin, initialize its mutex
	{
		if(pthread_mutex_init(&binMutexes[i], NULL))
	    {
	        printf("\nsetupBins(): RED FLAG! Unable to initialize a mutex\n");
	    }	
	}	
}

void freeBins()
{
	for(int i=0; i<numBins; i++)
	{
		free(binParticles[i]);
		free(binParticlesFlag[i]);
		pthread_mutex_destroy(&binMutexes[i]);
	}
	free(binParticles);
	free(binParticlesFlag);
	free(binMutexes);
	
	free(particlesPerBin);
	free(freeLocationPerBin);
}

int getFreeLocation(int bin)
{
	int loc = freeLocationPerBin[bin];
	
#ifdef __VIRAJ_DEBUG
		if(binParticlesFlag[bin][loc] != 0)
			printf("\ngetFreeLocation(): RED FLAG! @ binParticlesFlag[bin=%d][loc=%d] = %d, should be 0", bin, loc, binParticlesFlag[bin][loc]);		
#endif
		
	return loc;
}

void copyParticleToBin(particle_t& particle, int bdx)
{	
	// MUST BE CALLED FROM A THREAD SAFE CONTEXT
	
#ifdef __VIRAJ_DEBUG
	if(particlesPerBin[bdx] > maxParticlesPerBin)
		printf("\ncopyParticleToBin(): RED FLAG! More particlesPerBin[bdx=%d], maxParticlesPerBin=%d\n", bdx, maxParticlesPerBin);
#endif
	
	// copy particle data into binParticles[][] and update everything
	
	int loc = getFreeLocation(bdx);			// get index of first free location

	binParticles[bdx][loc] = particle; 		// copy data
	binParticlesFlag[bdx][loc] = 1; 		// set presence flag
	particlesPerBin[bdx]++; 				// increment particles per bin
	
#ifdef __VIRAJ_DEBUG
	if(particlesPerBin[bdx] == maxParticlesPerBin)
	{
		printf("\ncopyParticleToBin(): RED FLAG! BIN FULL! particlesPerBin[bdx=%d] = maxParticlesPerBin = %d\n", bdx, maxParticlesPerBin);
	}
#endif
		
	// update freeLocationPerBin[bin]
	int i = loc;
	
	while(binParticlesFlag[bdx][i] != 0)
	{	
		i++; // search for next free location
		if(i == maxParticlesPerBin) i = 0;
		
#ifdef __VIRAJ_DEBUG
		if(i == loc)
			printf("\ncopyParticleToBin(): RED FLAG! No space left in bin!!\n");
#endif

	}
	
	freeLocationPerBin[bdx] = i; // update freeLocationPerBin[bin]
}

void removeParticleFromBin(int bdx, int loc)
{
//	binParticles[bin][loc] = 0;			// reset particle data -> we dont need to do this
	binParticlesFlag[bdx][loc] = 0; 	// set presence flag to false
	particlesPerBin[bdx]--; 			// decrement counter
}

void doBinning()
{
#ifdef __VIRAJ_DEBUG
	printf("\ndoBinning(): Entered\n");
#endif

	int ndx = 0;
	for(ndx=0; ndx<n; ndx++) // for each particle
	{
		int bdx = (particles[ndx].x) / binLength;
		copyParticleToBin(particles[ndx], bdx);
	}
	
#ifdef __VIRAJ_DEBUG
	printf("\ndoBinning(): Binned %d particles\n", ndx);
#endif
}

void collectParticlesFromBins()
{
#ifdef __VIRAJ_DEBUG
	printf("\ncollectParticlesFromBins(): Entered\n");
#endif

	// Iterate over all bins and copy each particle to its correct location in the original particles array
	int particlesCollected = 0; // check variable
	for(int bdx=0; bdx<numBins; bdx++) // for each bin
	{
		particle_t *bin = binParticles[bdx];
		unsigned char *flags = binParticlesFlag[bdx];
		int numParticles = particlesPerBin[bdx];
				
		int loc_p = 0;
		for(int p=0; p<numParticles; p++) // for each particle in the bin
		{
			while(flags[loc_p]==0) loc_p++; // find next valid particle in the bin
		
#ifdef __VIRAJ_DEBUG	
			if(loc_p >= maxParticlesPerBin)
				printf("\ncollectParticlesFromBins(): RED FLAG! THIS SHOULD NEVER HAPPEN");
#endif
		
			particle_t particle = bin[loc_p];
			int id = particle.globalID;
			
			// Write particle to original particles array
			particles[id] = particle;
			
			// Increment check variable
			particlesCollected++;
			
			// Increment pointer
			loc_p++;
		}
	}
	
#ifdef __VIRAJ_DEBUG
	if(particlesCollected != n)
	{
		// We collected less particles than n - we missed some
		printf("\ncollectParticlesFromBins(): RED FLAG! We missed some! num of particles collected = %d, total = %d\n", particlesCollected, n);
	}
#endif
}

//
// Mark for move PARTICLE (particleIdx) to NEW BIN (targetBinIdx)
// Update: PARTICLES_TO_MOVE (particleArray) to NEW_BINS (targetBinsArray), NUMBER (moveParticleCounter)
//

void markForMove(int particleIdx, int targetBinIdx, int *particleArray, int *targetBinsArray, int& moveParticleCounter)
{
	particleArray[moveParticleCounter] = particleIdx;
	targetBinsArray[moveParticleCounter] = targetBinIdx;
	moveParticleCounter++;
	
#ifdef __VIRAJ_DEBUG
	if(moveParticleCounter >= maxParticlesPerBin)
		printf("\nmarkForMove(): RED FLAG! moveParticleCounter = %d exceeds bounds\n", moveParticleCounter);
#endif
}

//
// From BIN (fromBinIdx), move (moveParticleCounter) NUMBER of PARTICLES (particleArray) to NEW BINS (targetBinsArray). 
//

void moveBetweenBins(int fromBinIdx, int *particleArray, int *targetBinsArray, int& moveParticleCounter)
{
	// Must be threadsafe
	// Perform the operation of moving particles between bins
	
	particle_t *bin = binParticles[fromBinIdx];
	
	for(int i=0; i<moveParticleCounter; i++) // for each particle to be moved
	{
		int pIdx = particleArray[i];
		int toBinIdx = targetBinsArray[i];
		
		// THREAD_SAFE: Acquire exclusive access to bin (toBinIdx)
		pthread_mutex_lock(&binMutexes[toBinIdx]);
		copyParticleToBin(bin[pIdx], toBinIdx); // COPY particle to bin
		pthread_mutex_unlock(&binMutexes[toBinIdx]);
		// THREAD_SAFE: Release mutex

		// remove particle from current bin
		removeParticleFromBin(fromBinIdx, pIdx);
	}
	
	// Reset moveParticleCounter once particles have been moved
	moveParticleCounter = 0;
}

void compactBin(int bdx)
{
	// MUST BE THREADSAFE
	
	pthread_mutex_lock(&binMutexes[bdx]);
	particle_t *bin = binParticles[bdx];
	unsigned char *flags = binParticlesFlag[bdx];
	int numParticles = particlesPerBin[bdx];
	
	for(int loc=0; loc<numParticles; loc++)
	{		
		if(flags[loc] == 0) // we need to shift a particle into this slot
		{
			int nextLoc = loc + 1;
			while(flags[nextLoc]==0) nextLoc++; // find the next valid particle
			
			bin[loc] = bin[nextLoc]; // shift particle to loc
			flags[loc] = 1; // set this slot to 1
			flags[nextLoc] = 0; // mark nextLoc as invalid
		}
	}
	pthread_mutex_unlock(&binMutexes[bdx]);
}

void printBin(int bdx)
{
	// DEBUG CODE - Function should not be called while measuring performance
	
	printf("\nthreadRoutine(): # particlesPerBin[bdx=%d] = %d, maxParticlesPerBin=%d\n", bdx, particlesPerBin[bdx], maxParticlesPerBin);
	for(int z=0; z<maxParticlesPerBin; z++)
	{
		printf("\nthreadRoutine(): #printbin: binParticlesFlag[bdx=%d][z=%d] = %d\n", bdx, z, binParticlesFlag[bdx][z]);
	}
	
}

// 
// The REAL action
//

void* thread_routine(void *pthread_id)
{
	int threadIdx = *(int*) pthread_id;
	int bdx = threadIdx;

	// bin represents the bin array for this particular thread
	particle_t *bin = binParticles[bdx];
	unsigned char *flags = binParticlesFlag[bdx];
	int numParticles = particlesPerBin[bdx]; // !must be updated on every iteration of STEP loop
	
	// THESE ARE LOCAL THREAD VARIABLES
	int *moveParticleIdx = (int*) malloc(sizeof(int) * maxParticlesPerBin);
	int *moveTargetBinIdx = (int*) malloc(sizeof(int) * maxParticlesPerBin);
	int numMoveParticles = 0;

	for( int step = 0; step < NSTEPS; step++ )	// for each step
	{
		// printf("\nThread: %d, step: %d, ENTERING STEP LOOP\n", threadIdx, step);
		
#ifdef __VIRAJ_DEBUG
		// Check if we lost some particles
		
		if(threadIdx == 0)
		{
			// Find sum of number of particles in each bin
			int sum = 0;
			for(int i=0; i<numBins; i++)
			{
				sum += particlesPerBin[i];
			}
			
			if((n-sum) != 0)
				printf("\nthreadRoutine(): (TID=0) RED FLAG! Sum does not match!\n");
		}

#endif
		
		// update number of particles in THIS bin
		numParticles = particlesPerBin[bdx];
		
		//
		// DEFRAGMENT bins after a every NUM iterations
		// Defragmentation occurs every SAVEFREQ iterations anyways
		// For good perf, SAVEFREQ should be a multiple of NUM
		//
		if(step % (1) == 0) compactBin(bdx); //todo
		
		// THREAD_SYNC
		pthread_barrier_wait(&barrier);	

		//
		// SELF_LOOP - Compute O(n2) particle-particle interactions within the bin
		//

		int loc_i = 0;
		int loc_j = 0;
		for(int i=0; i<numParticles; i++) // for each particle in bin
		{
			// Find next valid particle
			while(flags[loc_i]==0) loc_i++;
			
#ifdef __VIRAJ_DEBUG
			if(loc_i >= maxParticlesPerBin) printf("\nthreadRoutine(): RED FLAG! (SELF) SegFault: loc_i, threadId= %d, loc_i=%d, num_particles=%d\n", threadIdx, loc_i, numParticles);
#endif		
			// for each particle, reset forces first time
			bin[loc_i].ax = bin[loc_i].ay = 0;
			
			loc_j = 0;
			for(int j=0; j<numParticles; j++) // again, for each particle
			{	
				while(flags[loc_j]==0) loc_j++;
				
#ifdef __VIRAJ_DEBUG
				if(loc_j >= maxParticlesPerBin) printf("\nthreadRoutine(): RED FLAG! (SELF) SegFault: loc_j, threadId= %d, loc_i=%d, loc_j=%d, j=%d, num_particles=%d\n", threadIdx, loc_i, loc_j, j, numParticles);
#endif
				// if(loc_i == loc_j) continue; // trivial case // todo
				
				apply_force(bin[loc_i], bin[loc_j]);
				loc_j++;
			}
			
			loc_i++;
		}
		
		
		//
		// Handle particles close to the edge of the bin
		//

		// TODO

		pthread_barrier_wait(&barrier);	

		//
		// PREV LOOP: Compute interactions with particles from previous bin
		//
		if(bdx != 0)
		{		
			particle_t *prevBin = binParticles[bdx-1];
			int prevNumParticles = particlesPerBin[bdx-1];
			unsigned char *prevFlags = binParticlesFlag[bdx-1];

			int loc_i = 0;
			int loc_j = 0;
			for(int i=0; i<numParticles; i++) // for each particle in THIS bin
			{
				while(flags[loc_i]==0) loc_i++;
				
#ifdef __VIRAJ_DEBUG
				if(loc_i >= maxParticlesPerBin) printf("\nthreadRoutine(): RED flag! (PREV) SEG_FAULT loc_i\n");
#endif			
				loc_j = 0;
				for(int j=0; j<prevNumParticles; j++) // for each particle in PREV bin
				{
					while(prevFlags[loc_j]==0) loc_j++;
					
#ifdef __VIRAJ_DEBUG
					if(loc_j >= maxParticlesPerBin) printf("\nthreadRoutine(): RED flag! SEG_FAULT for PREV bin at loc_j\n");
#endif
					apply_force(bin[loc_i], prevBin[loc_j]);
					loc_j++;
				}
				
				loc_i++;
			}
		}
		

		pthread_barrier_wait(&barrier);	
		
		//
		// NEXT LOOP: Compute interactions with particles from next bin
		//
		if(bdx != (numBins-1))
		{
			particle_t *nextBin = binParticles[bdx+1];
			int nextNumParticles = particlesPerBin[bdx+1];
			unsigned char *nextFlags = binParticlesFlag[bdx+1];

			int loc_i = 0;
			int loc_j = 0;
			for(int i=0; i<numParticles; i++) // for each particle in THIS bin
			{
				while(flags[loc_i]==0) loc_i++;
				
#ifdef __VIRAJ_DEBUG
				if(loc_i >= maxParticlesPerBin) printf("\nthreadRoutine(): RED flag! SEG_FAULT for NEXT bin at loc_i\n");
#endif
				loc_j = 0;
				for(int j=0; j<nextNumParticles; j++) // for each particle in NEXT bin
				{
					while(nextFlags[loc_j]==0) loc_j++;
					
#ifdef __VIRAJ_DEBUG
					if(loc_j >= maxParticlesPerBin) printf("\nthreadRoutine(): (TID=%d) RED flag! (NEXT) SEG_FAULT at loc_j = %d, nextNumParticles = %d\n", threadIdx, loc_j, nextNumParticles);
#endif
					apply_force(bin[loc_i], nextBin[loc_j]);
					loc_j++;
				}
				
				loc_i++;
			}
		}
		
		pthread_barrier_wait(&barrier);		
		
		
		//
		// MOVE_MARK LOOP: Find which particles should be moved
		//
		numMoveParticles = 0; // keeps a track of how many particles have been marked for move
		
		loc_i = 0;
		for(int i=0; i<numParticles; i++) // for each particle in this bin
		{
			while(flags[loc_i]==0) loc_i++;
			
			move(bin[loc_i]); // move particle

			int newBdx = (bin[loc_i].x) / binLength; // find if particle crossed bin boundries

			if(bdx != newBdx) // if it did cross bin boundaries
			{	
				// Mark for Move PARTICLE (loc_i) to BIN (newBinIdx)
				markForMove(loc_i, newBdx, moveParticleIdx, moveTargetBinIdx, numMoveParticles);
			}
			
			loc_i++;
		}
	
		pthread_barrier_wait(&barrier);
		
		//
		// TRANSFER particles between bins
		//
		moveBetweenBins(bdx, moveParticleIdx, moveTargetBinIdx, numMoveParticles);
		
		pthread_barrier_wait(&barrier);

		//
		// Collect particles and their states to disk
		//
		if(fsave && (step%SAVEFREQ) == 0)
		{
			compactBin(bdx);
		}
		if(threadIdx == 0 && fsave && (step%SAVEFREQ) == 0 )
		{
			printf("\nthreadRoutine(): SAVE FILES\n");
			collectParticlesFromBins();
			save(fsave, n, particles);
		}
		
		pthread_barrier_wait(&barrier);
	}	
	
	// cleanup
	free(moveParticleIdx);
	free(moveTargetBinIdx);
	return NULL;
}

//!VIRAJ

//
//  check that pthreads routine call was successful
//

#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

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
	setupBins();
	doBinning();
	memset(particles, 0, sizeof(particle_t) * n);
	// !VIRAJ

	// setup threads
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
