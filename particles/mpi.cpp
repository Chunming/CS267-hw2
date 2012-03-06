#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <math.h>
#include <string.h>

// Global variables
#define cutoff 0.01
#define density 0.0005
double spaceDim = 0.0;
particle_t **binParticles = NULL;
unsigned char **binParticlesFlag = NULL;
int *particlesPerBin = NULL;
int *freeLocationPerBin = NULL;
int numBins = -1;
int maxParticlesPerBin = -1;
double binLength = -1;

//
//  benchmarking program
//

// Initialize each local bin/processor
// nlocal       : No. of particles in specific bin/processor
// nlocalMax    : Max no. of particles allowed in bin/processor 
// freeIdx : Local free location

void compactBin(particle_t* local, unsigned char* localFlags, int &nlocal) {

   for (int loc=0; loc<nlocal; ++loc) {
      if (0 == localFlags[loc]) {
         int nextLoc = loc + 1;
	 while (0 == localFlags[nextLoc]) nextLoc++;
	 local[loc] = local[nextLoc];
	 localFlags[loc] = 1;
	 localFlags[nextLoc] = 0
      }
   }
}


void copyParticleToBin(particle_t *src, particle_t *dst, unsigned char *localFlags, int bdx, int &nlocal, int &nlocalMax, int &freeIdx) {
   dst[freeIdx] = *src;
   localFlags[freeIdx] = 1;
   nlocal++;

   int i = freeIdx;
   while (0 != localFlags[i]) {
      i++; // search for next free location
      if (i == nlocalMax) {i = 0;}
   }
   freeIdx = i;
}

int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    printf("Outside the MPI loop \n");

    //
    //  Set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc ); // n_proc = 24
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );


    //
    //  Set up the data partitioning across processors
    //



    //
    // Allocate storage for local partition/ Set up bins
    //
    
    spaceDim = sqrt(density * n); // 0.5 default
    numBins = n_proc; // No. of bins 
    binLength = spaceDim / numBins; // 0.5 / 24 default

    double bin_area = (spaceDim*spaceDim) / numBins; // Find max no. of particles per bin
    int particle_per_proc = 3* (int)( bin_area / (3.14*(cutoff/2)*(cutoff/2)) ); // Max particle num per processor
    // int particle_per_proc = (n + n_proc - 1) / n_proc;
    
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];

    int nlocalMax= partition_sizes[rank]; // Same as maxParticlesPerBin
    int localFreeLoc = 0; // Same as freeLocationPerBin
    int nlocal = 0; // Same as particlesPerBin

    particle_t *local = (particle_t*) malloc( nlocalMax * sizeof(particle_t) ); // Same as binParticles
    memset(local, 0, nlocalMax*sizeof(particle_t));

    particle_t *prevBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
    memset(prevBin, 0, nlocalMax*sizeof(particle_t));

    particle_t *nextBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
    memset(nextBin, 0, nlocalMax*sizeof(particle_t));

    unsigned char *localFlags = (unsigned char *) malloc( nlocal * sizeof(unsigned char)  ); // Same as binPariclesFlag
    memset(localFlags, 0, nlocal*sizeof(unsigned char));


    //
    //  Initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
//    if( rank == 0 )
    init_particles( n, particles );


/*
    // Scatters a buffer in parts to all processes in a communicator
    void* sSendBuf = particles;// address of the send buffer
    int* sSendCount = partition_sizes;// int array specifies no. of elems to send to each processor
    int* sDispls = partition_offsets; // int array specifies the disp. relative to the send buffer
    void* sRecvBuf = local; // Receiving buffer
    int sRecvCount = nlocal; // No. of elems in receive buffer
    MPI_Scatterv( sSendBuf, sSendCount, sDispls, PARTICLE, sRecvBuf, sRecvCount, PARTICLE, 0, MPI_COMM_WORLD );

    // Gathers data from all tasks and deliver combined data to all tasks
    // Blk of data sent frm jth process is received by every process
    // and placed in the jth block of the buffer gRecvBuf
    void* gSendBuf = local;// address of the send buffer
    int gSendCount = nlocal;// int array specifies no. of elems to send to each processor
    int* gDispls = partition_offsets; // int array specifies the disp. relative to the send buffer
    void* gRecvBuf = particles; // Receiving buffer
    int* gRecvCount = partition_sizes; // No. of elems in receive buffer
    MPI_Allgatherv( gSendBuf, gSendCount, PARTICLE, gRecvBuf, gRecCount, gRecBuf, PARTICLE, MPI_COMM_WORLD );
*/


   //
   // Do Binning onto local array
   //
   int ndx = 0;
   for (ndx=0; ndx<n; ndx++) {
      int bdx = (particles[ndx].y / binLength);
      if (bdx == rank) {
         copyParticleToBin(&particles[ndx], local, localFlags, bdx, nlocal, nlocalMax, localFreeLoc);
      }
   }





    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {   
        // 
	// 1. MPI send/receive particles to/from adjacent bins
	// MPI_Send(void* start, int numElem, DATA_TYPE, int source, int tag, COMM)
	// MPI_Recv(void* start, int numElem, DATA_TYPE, int source, int tag, COMM, status)
        
	int sTag1 = 100; // Message ID
	int recvdCount;
        int nPrevBin;
	int nNextBin;	
        MPI_Status status; 
	if (rank-1 >= 0) { // Check if top bin exists 
	   MPI_Send(local, nlocal, PARTICLE, rank-1, sTag1, MPI_COMM_WORLD); // Send to top bin	  
	   MPI_Recv(prevBin, nlocalMax, PARTICLE, rank-1, sTag1, MPI_COMM_WORLD, &status); //Recv from top bin
	   MPI_Get_count(&status, PARTICLE, &recvdCount); // Get received count
           nPrevBin = recvdCount; 

	}

	if (rank+1 <= 23) { // Check if bottom bin exists
	   MPI_Send(local, nlocal, PARTICLE, rank+1, sTag1, MPI_COMM_WORLD); // Send to bot bin
	   MPI_Recv(nextBin, nlocalMax, PARTICLE, rank+1, rank+1, MPI_COMM_WORLD, &status); // Recv from bot bin
	   MPI_Get_count(&status, PARTICLE, &recvdCount); // Get received count
           nNextBin = recvdCount; 
	}

	//
	// 2. Apply Force
	//
	for (int i=0; i<nlocal; i++) { // For each particle in local bin
           local[i].ax = local[i].ay = 0;

	   // SELF LOOP - Compute interactions with particles within the bin
	   for (int j=0; j<nlocal; ++j) {
	      apply_force (local[i], local[j]);
	   }

	   // PREV LOOP Compute interactions with particles from top bin
	   for (int j=0; j<nPrevBin; ++j) {
	      apply_force (local[i], prevBin[j]);
	   }

	   // NEXT LOOP Compute interactions with particles from bot bin
	   for (int j=0; j<nNextBin; ++j) {
	      apply_force (local[i], nextBin[j]);
	   }

	}

	// 
	// 3. Move Particles
	//
	for (int i=0; i<nlocal; i++) {
	   move( local[i] );
	}

	//
	// 4. Re-bin Particles
	//
	int tag4 = 40;
	int idx = 0;
	for (int i=0; i<nlocal; i++) {
           int bdx = (local[idx].y / binLength);
           if (bdx == rank) { // If particle is still in same bin, do nothing
	   }

	   else if (bdx == rank-1) { // Particle moved to top bin
	      MPI_Send(&local[idx], 1, PARTICLE, rank-1, tag4, MPI_COMM_WORLD); //Send to top bin
	      localFlags[idx] = 0;
	      nlocal--;
	   }

	   else { // (bdx == rank+1) Particle moved to bot bin
	      MPI_Send(&local[idx], 1, PARTICLE, rank+1, tag4, MPI_COMM_WORLD); //Send to bot bin
	      localFlags[idx] = 0;
	      nlocal--;
	   }
	   idx++;
	}
	
	// Beyond current idx value, all other array elements are invalid
	int recvdCount;
	if (rank-1 >= 0) { // Check if top bin exists	
	   MPI_Recv(local[idx], nlocalMax, PARTICLE, rank-1, tag4, MPI_COMM_WORLD, &status); //Recv from top bin
	   MPI_Get_count(&status, PARTICLE, &recvdCount); // Get received count
	   for (int j=idx; j<recvdCount; ++j) localFlags[j]=1;
	   idx = idx + recvdCount;
	   nlocal += recvdCount;
	}

	if (rank+1 <= 23) { // Check if bot bin exists
	   MPI_Recv(local[i], nlocalMax, PARTICLE, rank+1, tag4, MPI_COMM_WORLD, &status); //Recv from bot bin
	   MPI_Get_count(&status, PARTICLE, &recvd_count);
	   for (int j=idx; j<recvdCount; ++j) localFlags[j]=1;
	   idx = idx + recvdCount;
	   nlocal += recvdCount;
	}

	//
	// 5. Compact Particles
	//
	

        compactBin(local, localFlags, nlocal);
	 
	// 
        //  Collect all global data locally (not good idea to do)
        //  gathers data from all tasks & deliver combined data to all tasks
	//
	MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );


	//
        //  save current step if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( prevBin );
    free( nextBin );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
