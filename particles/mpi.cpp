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

//
//  benchmarking program
//

// Initialize each local bin/processor
// nlocal       : No. of particles in specific bin/processor
// nlocalMax    : Max no. of particles allowed in bin/processor 
// freeIdx : Local free location

int isCloseToEdge(particle_t &particle, double binEdge) {
   double dy = binEdge - particle.y;
   double r2 = dy * dy;
   if (r2 > cutoff * cutoff) 
      return 0;
   else 
      return 1;
}

void copyParticleToBin(particle_t *src, particle_t *dst, unsigned char *localFlags, int bdx, int *nlocal, int &nlocalMax, int &freeIdx) {
   dst[freeIdx] = *src;
   localFlags[freeIdx] = 1;
   (*nlocal)++; // Increment no. of elems in localBin

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
    if (NULL == particles) {
       printf("ERR allocating *particles \n");
       return -1;
    }
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

  
    //
    // Allocate storage for local partition/ Set up bins
    //
    
    double spaceDim = sqrt(density * n); // 0.5 default
    int numBins = n_proc; // No. of bins 
    double binLength = spaceDim / numBins; // 0.5/24 = 0.020833 by default
    double bin_area = (spaceDim*spaceDim) / numBins; // Find max no. of particles per bin

    // 
    //  Set up data partitioning across processors
    //
    //int particle_per_proc = (n + n_proc - 1) / n_proc;
    int particle_per_proc = (int)( bin_area / (3.14*(cutoff/2)*(cutoff/2)) ); // Max particle num per proc

    printf("particle_per_proc is %d \n", particle_per_proc);

    int nlocalMax = particle_per_proc;

    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    if (NULL == partition_offsets) {
       printf("ERR allocating *partition_offsets \n");
       return -1;
    }

    for( int i = 0; i < n_proc+1; i++ )
       partition_offsets[i] = min( i * particle_per_proc, n );
 
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    if (NULL == partition_sizes) {
       printf("ERR allocating *partition_sizes \n");
       return -1;
    }

    for( int i = 0; i < n_proc; i++ )
       partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    	 
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    printf("nlocal is %d \n", nlocal);

    particle_t *localBin = (particle_t*) malloc( nlocal * sizeof(particle_t) ); // Replace nlocal with nlocalMax
    if (NULL == localBin) {
       printf("ERR allocating *localBin \n");
       return -1;
    }
    //memset(localBin, 0, nlocalMax*sizeof(particle_t));

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
       init_particles( n, particles );



    // Scatters a buffer in parts to all processes in a communicator

    // MPI_Scatterv(void* sendBuf, int* sendCount, int* Displs, void* recvBuf, int recvCount, data_type, int rank, handle)
    // void* sSendBuf = particles;// address of the send buffer
    // int* sSendCount = partition_sizes;// int array specifies no. of elems to send to each processor
    // int* sDispls = partition_offsets; // int array specifies the disp. relative to the send buffer
    // void* sRecvBuf = local; // Receiving buffer
    // int sRecvCount = nlocal; // No. of elems in receive buffer
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, localBin, nlocal, PARTICLE, 0, MPI_COMM_WORLD );


    // Gathers data from all tasks and deliver combined data to all tasks
    // Blk of data sent frm jth process is received by every process
    // and placed in the jth block of the buffer gRecvBuf

    // MPI_Allgatherv(void* sendBuf, int* sendCount, data_type, void* recvBuf, int* recvCount, int* Displs, data_type, handle);
    // void* sendBuf  : address of the send buffer
    // int sendCount  : int array specifies no. of elems to send to each processor
    // int* Displs    : int array specifies the disp. relative to the send buffer
    // void* recvBuf  : Receiving buffer
    // int* recvCount : No. of elems in receive buffer
    MPI_Allgatherv( localBin, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );


    // totalN is only used by rank = 1
    int* totalN = (int*) malloc(sizeof(int));
    if (NULL == totalN) {
       printf("ERR allocating *totalN \n");
       return -1;
    }
    (*totalN) = 0;

    // Tmp placeholder for particles from previous bin
    particle_t *prevBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
    if (NULL == prevBin) {
       printf("ERR allocating *prevBin \n");
       return -1;
    }
    memset(prevBin, 0, nlocalMax*sizeof(particle_t));

    // Tmp placeholder for particles from next bin
    particle_t *nextBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
    if (NULL == nextBin) {
       printf("ERR allocating *nextBin \n");
       return -1;
    }
    memset(nextBin, 0, nlocalMax*sizeof(particle_t));

    // Use array of flags to refer to particles in localBin
    unsigned char *localFlags = (unsigned char *) malloc( nlocalMax * sizeof(unsigned char)  ); // Same as binPariclesFlag
    if (NULL == localFlags) {
       printf("ERR allocating *localFlags \n");
       return -1;
    }
    memset(localFlags, 0, nlocalMax*sizeof(unsigned char));



   //
   // Do initial binning onto localBin array
   //

   memset(localBin, 0, nlocal*sizeof(particle_t));
   nlocal = 0;
   int localFreeLoc = 0; // Same as freeLocationPerBin
   for (int ndx=0; ndx<n; ++ndx) {
      int bdx = (particles[ndx].y / binLength);
      if (bdx == rank) {
         copyParticleToBin(&particles[ndx], localBin, localFlags, bdx, &nlocal, nlocalMax, localFreeLoc);
      }
   }


    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    //int recvSig = 0;
    //int prevSig = 0;
    //int nextSig = 0;
    for( int step = 0; step < NSTEPS; step++ )
    {  
	//printf("Time step is %d from rank %d \n", step, rank);
	//printf("AFT nlocal from rank %d is %d \n", rank, *nlocal);	    


	// 
	// 1. MPI send/receive particles to/from adjacent bins
	// MPI_Send(void* start, int numElem, DATA_TYPE, int source, int tag, COMM)
	// MPI_Recv(void* start, int numElem, DATA_TYPE, int source, int tag, COMM, status)
        // Each MPI_Send must be paired with a MPI_Recv
        // Cannot call successive MPI_Sends and call just one MPI_Recv

	int tag1 = 100; // Message ID
	int adjCount; // Received count from adjacent bins
        int nPrevBin=0; // No. of elems from prevBin
	int nNextBin=0;
        double binEdge = rank*binLength;
	int idx = 0;
	bool fPrevCheck = 0;	
	bool fNextCheck = 0;	
	MPI_Status status; 

        // If top bin exists, then can send
	if (rank-1 >= 0) {            
	   idx = 0;
	   fPrevCheck = 0; // Check if comms with prev bin is req
    	   for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
	      while (localFlags[idx]==0) idx++; 
  	      if (isCloseToEdge(localBin[idx], binEdge)) { // Comms with prev bin is req
		    fPrevCheck = 1;
		    break;
	      }
	      idx++;
	   }

	   if (1==fPrevCheck) MPI_Send(localBin, nlocal, PARTICLE, rank-1, tag1, MPI_COMM_WORLD);
	   else MPI_Send(localBin, 0, PARTICLE, rank-1, tag1, MPI_COMM_WORLD); 
	  
	}

        // If bottom bin exists, then can receive
	if (rank+1 <= 23) {
	   MPI_Recv(nextBin, nlocalMax, PARTICLE, rank+1, tag1, MPI_COMM_WORLD, &status); // Recv from bot bin
	   MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
	   nNextBin = adjCount; 
	}

        // If bottom bin exists, then can send
	if (rank+1 <= 23) { 
     	   idx = 0;
	   fNextCheck = 0;
	   for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
     	      while (localFlags[idx]==0) idx++;
	      if (isCloseToEdge(localBin[idx], binEdge)) {
	         fNextCheck = 1;
		 break;
	      }
	      idx++;
	   }

	   if (1==fNextCheck) MPI_Send(localBin, nlocal, PARTICLE, rank+1, tag1, MPI_COMM_WORLD); // Send to bot bin
	   else MPI_Send(localBin, 0, PARTICLE, rank+1, tag1, MPI_COMM_WORLD); // Send to bot bin
	}

        // If top bin exists, then can receive
	if (rank-1 >= 0) { 
	   // Check receive signal from prevBin
	   MPI_Recv(prevBin, nlocalMax, PARTICLE, rank-1, tag1, MPI_COMM_WORLD, &status); //Recv from top bin
	   MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
	   nPrevBin = adjCount; 
	}


	//
	// 2. Apply Force
	//
	int loc_i = 0;
	int loc_j = 0;
	for (int i=0; i<(nlocal); ++i) { // For each particle in local bin

           while(localFlags[loc_i]==0) loc_i++;
	   localBin[i].ax = localBin[i].ay = 0;

	   
	   // SELF LOOP - Compute interactions with particles within the bin
	   loc_j = 0;
	   for (int j=0; j<(nlocal); ++j) {
	      while(localFlags[loc_j]==0) loc_j++;
	      apply_force (localBin[loc_i], localBin[loc_j]);
	      loc_j++;
	   }

	   // PREV LOOP Compute interactions with particles from top bin
	   for (int j=0; j<nPrevBin; ++j) {
	      apply_force (localBin[loc_i], prevBin[j]);
	   }

	   // NEXT LOOP Compute interactions with particles from bot bin
	   for (int j=0; j<nNextBin; ++j) {
	      apply_force (localBin[loc_i], nextBin[j]);
	   }

	   loc_i++;
	}
	
	memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
	memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
	nPrevBin = 0;
	nNextBin = 0;

	//printf("COMPLETED Apply Force in rank %d \n", rank);

	// 
	// 3. Move Particles
	//
	loc_i = 0;
	for (int i=0; i<(nlocal); ++i) {
	   while(localFlags[loc_i]==0) loc_i++; // Make sure particle_t has valid data 
	   move( localBin[loc_i] );
	   loc_i++;
	}

	//
	// 4. Re-bin Particles
	//
	int bdx;
	int jdx = 0; int kdx = 0;
	idx = 0; 
	for (int i=0; i<(nlocal); ++i) { // Analyze each particle in localBin
	   while(localFlags[idx]==0) idx++;
           bdx = (localBin[idx].y / binLength);

	   if (bdx == rank-1) { // Particle moved to top bin
	      localFlags[idx] = 0; // Remove from localBin
	      prevBin[jdx] = localBin[idx];
	      jdx++;
	      nlocal--;
	   }

	   if (bdx == rank+1) { // Particle moved to top bin
	      localFlags[idx] = 0; // Remove from localBin
	      nextBin[kdx] = localBin[idx];
	      kdx++;
	      nlocal--;
	   }
	   idx++;
	}
	nPrevBin = jdx; // No. of elems to shift from prevBin to localBin
	nNextBin = kdx; // No. of elems to shift from nextBin to localBin


	// Get particles from adjacent bins
	// Using synchronous blocking send/receive
	int rebinCount;
	int tag4 = 400;


        // Have EVEN processors send particles first
	// and ODD processors receive
	if (0 == rank%2) { // Even Processors  
	   if (rank-1 >= 0) { // If top bin exists, send
     	      MPI_Send(prevBin, nPrevBin, PARTICLE, rank-1, tag4, MPI_COMM_WORLD); //Send to bot bin
	   }

	}

	if (1 == rank%2) { // Odd Processors
	   if (rank+1 <= 23) { // If bottom bin exists, receive
	      MPI_Recv(&localBin[idx], nlocalMax, PARTICLE, rank+1, tag4, MPI_COMM_WORLD, &status); //Recv from top bin
	      MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
	      for (int j=idx; j<rebinCount; ++j) localFlags[j]=1;
	      idx = idx + rebinCount;
	      nlocal += rebinCount;
	   }
	}

	if (0 == rank%2) { // Even Processors  
	   if (rank+1 <= 23) { // If bottom bin exists, send 
     	      MPI_Send(nextBin, nNextBin, PARTICLE, rank+1, tag4, MPI_COMM_WORLD); //Send to bot bin
   	   }
	}

	if (1 == rank%2) { // Odd Processors
	   if (rank-1 >=0) { // If top bin exists, receive
	      MPI_Recv(&localBin[idx], nlocalMax, PARTICLE, rank-1, tag4, MPI_COMM_WORLD, &status); //Recv from bot bin
	      MPI_Get_count(&status, PARTICLE, &rebinCount);
	      for (int j=idx; j<rebinCount; ++j) localFlags[j]=1;
	      idx = idx + rebinCount;
	      nlocal += rebinCount;
	   }
	}

	// Have ODD processors send next
	// and EVEN processors receive
	if (1 == rank%2) { // Odd Processors  
	   if (rank-1 >= 0) { // If top bin exists, send
     	      MPI_Send(prevBin, nPrevBin, PARTICLE, rank-1, tag4, MPI_COMM_WORLD); //Send to bot bin
	   }

	}

	if (0 == rank%2) { // Even Processors
	   if (rank+1 <= 23) { // If bottom bin exists, receive
	      MPI_Recv(&localBin[idx], nlocalMax, PARTICLE, rank+1, tag4, MPI_COMM_WORLD, &status); //Recv from top bin
	      MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
	      for (int j=idx; j<rebinCount; ++j) localFlags[j]=1;
	      idx = idx + rebinCount;
	      nlocal += rebinCount;
	   }
	}

	if (1 == rank%2) { // Odd Processors  
	   if (rank+1 <= 23) { // If bottom bin exists, send 
     	      MPI_Send(nextBin, nNextBin, PARTICLE, rank+1, tag4, MPI_COMM_WORLD); //Send to bot bin
   	   }
	}

	if (0 == rank%2) { // Even Processors
	   if (rank-1 >=0) { // If top bin exists, receive
	      MPI_Recv(&localBin[idx], nlocalMax, PARTICLE, rank-1, tag4, MPI_COMM_WORLD, &status); //Recv from bot bin
	      MPI_Get_count(&status, PARTICLE, &rebinCount);
	      for (int j=idx; j<rebinCount; ++j) localFlags[j]=1;
	      idx = idx + rebinCount;
	      nlocal += rebinCount;
	   }
	}





















	memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
	memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
	nPrevBin = 0; // No. of elems to shift from prevBin to localBin
	nNextBin = 0; // No. of elems to shift from nextBin to localBin
	
        MPI_Reduce(&nlocal, totalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );        
        if (rank == 0) printf("Total N is %d \n", *totalN);

	//
	// 5. Compact Particles
	//
	int nextLoc;
         for (int loc=0; loc<nlocal; ++loc) {
            if (0 == localFlags[loc]) {
               nextLoc = loc + 1;
	       while (0 == localFlags[nextLoc]) nextLoc++;
	       localBin[loc] = localBin[nextLoc];
	       localFlags[loc] = 1;
	       localFlags[nextLoc] = 0;
            }
         }









	// 
        //  Collect all global data locally (not good idea to do)
        //  gathers data from all tasks & deliver combined data to all tasks
	//
	//MPI_Allgatherv( localBin, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );


	//
        //  save current step if necessary
        //
        //if( fsave && (step%SAVEFREQ) == 0 )
            //save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    //free( nlocal );
    free( totalN );
    free( partition_offsets );
    free( partition_sizes );
    free( localBin );
    free( prevBin );
    free( nextBin );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
