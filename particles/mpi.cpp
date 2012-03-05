#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

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
void copyParticleToBin(particle_t *src, particle_t *dst, unsigned char *localFlag, int bdx, int &nlocal, int &nlocalMax, int &freeIdx) {
   dst[freeIdx] = *src;
   localFlag[freeIdx] = 1;
   nlocal++;

   int i = freeIdx;
   while (0 != localFlag[i]) {
      i++ // search for next free location
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
    memset(local, 0, sizeof(particle_t) * nlocal);
    unsigned char *localFlag = (unsigned char *) malloc( nlocal * sizeof(unsigned char)  ); // Same as binPariclesFlag
    memset(localFlag, 0, sizeof(unsigned char) * nlocal);


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
         copyParticleToBin(particles[ndx], local, localFlag, bdx, nlocal, nlocalMax, localFreeLoc);
      }
   }





    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute all forces
        //
        // nlocal is the no. of particles in a specific bin/processor per time step
        for( int i = 0; i < nlocal; i++ ) 
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( local[i], particles[j] );
        }
        
        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );

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
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
