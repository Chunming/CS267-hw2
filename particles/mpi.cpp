#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <math.h>
#include <string.h>
#include <vector>

int isCloseToEdge(particle_t &particle, double binEdge, double cutoff) {
    double dy = binEdge - particle.y;
    double r2 = dy * dy;
    if (r2 > cutoff * cutoff) 
       return 0;
    else 
       return 1;
}

int compare(const void *a, const void* b) {
   particle_t *a0 = *(particle_t**) a;
   particle_t *b0 = *(particle_t**) b;

   if (a0->globalID > b0->globalID) return 1;
   else return -1;
}

//
//  benchmarking program
//
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
    //  set up MPI
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
    particle_t *particles_mpi = (particle_t*) malloc( n * sizeof(particle_t) );
 
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    double cutoff = 0.01;
    double density = 0.0005;
    double spaceDim = sqrt(density * n); // 0.5 default
    int numBins = n_proc; // No. of bins 
    double binLength = spaceDim / numBins; // 0.5/24 = 0.020833 by default

    double bin_area = (spaceDim*spaceDim) / numBins; // Find max no. of particles per bin
    int nlocalMax = 3 * (int)( bin_area / (3.14*(cutoff/2)*(cutoff/2)) ); // Max particle num per proc
    vector<double> xVect;
    vector<double> yVect;
    vector<int> globalIDVect;

   particle_t** compactVect;
   compactVect = (particle_t**)malloc( n * sizeof(particle_t*));

   particle_t* particleVect;
   particleVect = (particle_t*)malloc( n_proc*nlocalMax*sizeof(particle_t));
   if (NULL == particleVect) {
      printf("ERR allocating *particleVect \n");
      return -1;
   }

   int* flagVect;
   flagVect = (int*)malloc( n_proc*nlocalMax*sizeof(int));
   if (NULL == flagVect) {
      printf("ERR allocating *flagVect \n");
      return -1;
   }

   // No. of local particles from each proc
   int* nlocalVect;
   nlocalVect = (int*)malloc( n_proc*sizeof(int));
   if (NULL == nlocalVect) {
      printf("ERR allocating *nlocalVect \n");
      return -1;
   }
   memset(flagVect, 0, n_proc*nlocalMax*sizeof(int));
   memset(particleVect, 0, n_proc*nlocalMax*sizeof(particle_t));
   memset(nlocalVect, 0, n_proc*sizeof(int));



    int *displs = (int*) malloc( (n_proc) * sizeof(int) );
    int *rcounts = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ ) {
        displs[i] = i * nlocalMax;
        rcounts[i] = nlocalMax; 
    }

    int *displsC = (int*) malloc( (n_proc) * sizeof(int) );
    int *rcountsC = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ ) {
        displsC[i] = i;
        rcountsC[i] = 1; 
    }
 
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *localBin = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, localBin, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    MPI_Allgatherv( localBin, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
  

    // Do initial binning onto localBin array
    if (binLength<cutoff) {
       printf("ERROR, subBlock width cannot be smaller than cutoff value \n");
       return -1;
    }

    free (localBin);
    localBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) ); // Replace nlocal with nlocalMax
    if (NULL == localBin) {
       printf("ERR allocating *localBin \n");
       return -1;
    }
    memset(localBin, 0, nlocalMax*sizeof(particle_t));
    nlocal = 0;

    int* totalN = (int*) malloc(sizeof(int));
     (*totalN) = 0;
 
     // Tmp placeholder to receive particles from previous bin
     particle_t *prevBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
     if (NULL == prevBin) {
        printf("ERR allocating *prevBin \n");
        return -1;
     }
     memset(prevBin, 0, nlocalMax*sizeof(particle_t));
 
     // Tmp placeholder to receive particles from next bin
     particle_t *nextBin = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
     if (NULL == nextBin) {
        printf("ERR allocating *nextBin \n");
        return -1;
     }
     memset(nextBin, 0, nlocalMax*sizeof(particle_t));

     // Tmp placeholder to send particles to previous bin 
     particle_t *prevBinSend = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
     if (NULL == prevBinSend) {
        printf("ERR allocating *prevBinSend \n");
        return -1;
     }
     memset(prevBinSend, 0, nlocalMax*sizeof(particle_t));

     // Tmp placeholder to send particles to next bin 
     particle_t *nextBinSend = (particle_t*) malloc( nlocalMax * sizeof(particle_t) );
     if (NULL == nextBinSend) {
        printf("ERR allocating *nextBinSend \n");
        return -1;
     }
     memset(nextBinSend, 0, nlocalMax*sizeof(particle_t));


    int *localFlags = (int *) malloc( nlocalMax * sizeof(int)  ); // Same as binPariclesFlag
    if (NULL == localFlags) {
       printf("ERR allocating *localFlags \n");
       return -1;
    }
    memset(localFlags, 0, nlocalMax*sizeof(int));


    // 
    // Start binning
    //
    int freeIdx = 0; // Same as freeLocationPerBin
    for (int ndx=0; ndx<n; ++ndx) {
       int bdx = (particles[ndx].y / binLength);
       if (bdx == rank) {
          localBin[freeIdx] = particles[ndx];
          localFlags[freeIdx] = 1;
          nlocal++; // Increment no. of elems in localBin
 
          int i = freeIdx + 1;
          while (0 != localFlags[i]) {
             i++; // search for next free location
             if (i == nlocalMax) {i = 0;}
          }
          freeIdx = i;
       }
    }
   //free( particles ); // Particle array is not used after binning


    //
    //  Simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {

        int tag1 = 100; // Message ID
        int adjCount; // Received count from adjacent bins
        int nPrevBin=0; // No. of elems from prevBin
        int nNextBin=0;
        double binEdge = rank*binLength;
        int idx = 0;
	int sIdx = 0; // Send index
        MPI_Status status;

        // First, EVEN will send, ODD will receive
        if (0 == rank%2) { // EVEN
           if (rank-1 >= 0) { // If top bin exists, then can send
              idx = 0;
	      sIdx = 0;
              //fPrevCheck = 1; // Check if comms with prev bin is req
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge, cutoff)) { // Comms with prev bin is req
		       prevBinSend[sIdx] = localBin[idx]; // Set send index
		       sIdx++;
                 }
                 idx++;
              }
              MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1, MPI_COMM_WORLD);
           }
        }

        if (1 == rank%2) { // ODD
           if (rank+1 <= n_proc-1) { // If bottom bin exists, then can receive
              MPI_Recv(&nextBin[nNextBin], nlocalMax, PARTICLE, rank+1, tag1, MPI_COMM_WORLD, &status); // Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
              nNextBin += adjCount;
           }
        }

        if (0 == rank%2) { // EVEN 
           if (rank+1 <= n_proc-1) { // If bottom bin exists, then can send
              idx = 0;
	      sIdx = 0;
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge+binLength, cutoff)) {
		    nextBinSend[sIdx] = localBin[idx];
		    sIdx++;
                 }
                 idx++;
              }
	      MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+1, MPI_COMM_WORLD); // Send to bot bin
           }
        }

        if (1 == rank%2) { // ODD
           if (rank-1 >= 0) { // If top bin exists, then can receive
              // Check receive signal from prevBin
              MPI_Recv(&prevBin[nPrevBin], nlocalMax, PARTICLE, rank-1, tag1+1, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
	      nPrevBin += adjCount;
           }
        }


        // Now ODD will send, EVEN will receive
        if (1 == rank%2) { // ODD
           if (rank-1 >= 0) { // If top bin exists, then can send
              idx = 0;
	      sIdx = 0;
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge, cutoff)) { // Comms with prev bin is req
		       prevBinSend[sIdx] = localBin[idx];
		       sIdx++;
                 }
                 idx++;
              }
	      MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1+2, MPI_COMM_WORLD);
           }
        }

        if (0 == rank%2) { // EVEN
           if (rank+1 <= n_proc-1) { // If bottom bin exists, then can receive
              MPI_Recv(&nextBin[nNextBin], nlocalMax, PARTICLE, rank+1, tag1+2, MPI_COMM_WORLD, &status); // Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
              nNextBin += adjCount;
           }
        }

        if (1 == rank%2) { // ODD 
           if (rank+1 <= n_proc-1) { // If bottom bin exists, then can send
              idx = 0;
	      sIdx = 0;
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge+binLength, cutoff)) {
		    nextBinSend[sIdx] = localBin[idx];
		    sIdx++;
                 }
                 idx++;
              }
	      MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+3, MPI_COMM_WORLD); // Send to bot bin 
           }
        }

        if (0 == rank%2) { // EVEN
           if (rank-1 >= 0) { // If top bin exists, then can receive
              // Check receive signal from prevBin
              MPI_Recv(&prevBin[nPrevBin], nlocalMax, PARTICLE, rank-1, tag1+3, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
              nPrevBin += adjCount;
           }
        }
        memset(prevBinSend, 0, nPrevBin*sizeof(particle_t));
        memset(nextBinSend, 0, nNextBin*sizeof(particle_t));


	// 
        // 2. Apply Force
        //
        for (int i=0; i< nlocal; ++i) { // For each particle in local bin

           localBin[i].ax = localBin[i].ay = 0;

           // SELF LOOP - Compute interactions with particles within the bin
           for (int j=0; j< nlocal; ++j) {
              apply_force (localBin[i], localBin[j]);
           }

           // PREV LOOP Compute interactions with particles from top bin
           for (int j=0; j<nPrevBin; ++j) {
              apply_force (localBin[i], prevBin[j]);
           }

           // NEXT LOOP Compute interactions with particles from bot bin
           for (int j=0; j<nNextBin; ++j) {
              apply_force (localBin[i], nextBin[j]);
	   }
        }
        memset(prevBin, 0, nPrevBin*sizeof(particle_t));
        memset(nextBin, 0, nNextBin*sizeof(particle_t)); 
        nPrevBin = 0;
        nNextBin = 0;


        // 
        // 3. Move Particles
        //
        for (int i=0; i< nlocal; ++i) {move( localBin[i] );}


        //
        // 4. Re-bin Particles
        //
        int bdx;
        int jdx = 0; int kdx = 0;
        idx = 0;
        for (int i=0; i<nlocal; ++i) { // Analyze each particle in localBin
           while(0==localFlags[idx]) { 
	      idx++; 
	      if (idx == nlocalMax) {idx = 0;}
	   }
           bdx = (localBin[idx].y / binLength);

	   if (bdx == rank) {
	      // Do nothing
	   }
           else if (bdx == rank-1) { // Particle moved to top bin
              localFlags[idx] = 0; // Remove from localBin
              prevBin[jdx] = localBin[idx];
              jdx++;
           }

           else if (bdx == rank+1) { // Particle moved to top bin
              localFlags[idx] = 0; // Remove from localBin
              nextBin[kdx] = localBin[idx];
              kdx++;
           }
	   else {
	      printf("ERROR outside blocks! \n");
	   }

	
           idx++;

        }
        nPrevBin = jdx; // No. of elems to shift from prevBin to localBin
        nNextBin = kdx; // No. of elems to shift from nextBin to localBin
	nlocal = nlocal - nPrevBin - nNextBin;


        //
        // 5. Compact Particles
        //
        int nextLoc;
         for (int loc=0; loc<nlocal; ++loc) {
            if (0 == localFlags[loc]) {
               nextLoc = loc + 1;
               while (0 == localFlags[nextLoc]) { 
	          nextLoc++;
		  if (nextLoc == nlocalMax) {nextLoc = 0;}
	       }
               localBin[loc] = localBin[nextLoc];
               localFlags[loc] = 1;
               localFlags[nextLoc] = 0;
            }  
         }  


	//
        // 6. Get particles from adj bins. Use sync blocking send/receive
	//
        int rebinCount;
        int tag4 = 400;

        // Have EVEN procs send particles first and ODD procs receive
        if (0 == rank%2) { // Even Processors  
           if (rank-1 >= 0) { // If top bin exists, send
              MPI_Send(prevBin, nPrevBin, PARTICLE, rank-1, tag4, MPI_COMM_WORLD); //Send to top bin
           }
         }

        if (1 == rank%2) { // Odd Processors
           if (rank+1 <= n_proc-1) { // If bottom bin exists, receive
              MPI_Recv(&localBin[nlocal], nlocalMax, PARTICLE, rank+1, tag4, MPI_COMM_WORLD, &status); //Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
              for (int j=nlocal; j<(nlocal+rebinCount); ++j) {localFlags[j]=1;}
              nlocal += rebinCount;
           }
        }

        if (0 == rank%2) { // Even Processors  
           if (rank+1 <= n_proc-1) { // If bottom bin exists, send 
              MPI_Send(nextBin, nNextBin, PARTICLE, rank+1, tag4+1, MPI_COMM_WORLD); //Send to bot bin
           }
        }

        if (1 == rank%2) { // Odd Processors
           if (rank-1 >=0) { // If top bin exists, receive
              MPI_Recv(&localBin[nlocal], nlocalMax, PARTICLE, rank-1, tag4+1, MPI_COMM_WORLD, &status); //Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &rebinCount);
              for (int j=nlocal; j<(nlocal+rebinCount); ++j) {localFlags[j]=1;}
              nlocal += rebinCount;
           }
        }

         // Have ODD processors send next
        // and EVEN processors receive
        if (1 == rank%2) { // Odd Processors  
           if (rank-1 >= 0) { // If top bin exists, send
              MPI_Send(prevBin, nPrevBin, PARTICLE, rank-1, tag4+2, MPI_COMM_WORLD); //Send to top bin
           }
        }

        if (0 == rank%2) { // Even Processors
           if (rank+1 <= n_proc-1) { // If bottom bin exists, receive
              MPI_Recv(&localBin[nlocal], nlocalMax, PARTICLE, rank+1, tag4+2, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
              for (int j=nlocal; j<(nlocal+rebinCount); ++j) {localFlags[j]=1;}
              nlocal += rebinCount;
           }
        }

        if (1 == rank%2) { // Odd Processors  
           if (rank+1 <= n_proc-1) { // If bottom bin exists, send 
              MPI_Send(nextBin, nNextBin, PARTICLE, rank+1, tag4+3, MPI_COMM_WORLD); //Send to bot bin
           }
        }

        if (0 == rank%2) { // Even Processors
           if (rank-1 >=0) { // If top bin exists, receive
              MPI_Recv(&localBin[nlocal], nlocalMax, PARTICLE, rank-1, tag4+3, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &rebinCount);
              for (int j=nlocal; j<(nlocal+rebinCount); ++j) {localFlags[j]=1;}
	      nlocal += rebinCount;
	   }
	}	

        memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
        memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
        nPrevBin = 0; // No. of elems to shift from prevBin to localBin
        nNextBin = 0; // No. of elems to shift from nextBin to localBin

   	// Check total num of particles
   	//MPI_Reduce(&nlocal, totalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );        
   	//if (rank == 0) printf("Total N is %d \n", *totalN);

/*
	//
        //  Save current step if necessary
        //
	MPI_Gatherv( localBin, nlocalMax, PARTICLE, particleVect, rcounts, displs, PARTICLE, 0, MPI_COMM_WORLD );
        MPI_Gatherv(&nlocal, 1, MPI_INT, nlocalVect, rcountsC, displsC, MPI_INT, 0, MPI_COMM_WORLD);

	//
	// At Master Processor, gather all particles, reorder based on globalID, and save.
	//
        if (rank == 0) {
           if( fsave && (step%SAVEFREQ) == 0 ) {

	   int cdx=0;
           for (int rdx=0; rdx<n_proc; rdx++) {
	      for (int sdx=0; sdx<nlocalVect[rdx]; sdx++) {
	         compactVect[cdx] = &particleVect[displs[rdx] + sdx]; // displs[rdx] is displac of one proc
	         cdx++;
	      }
           }

	   for (int rdx=0; rdx<cdx; rdx++){
              for (int sdx=0; sdx<cdx; sdx++) {
	         if (compactVect[rdx]->globalID == compactVect[sdx]->globalID && rdx!=sdx) {
	            printf("ERROR: Same globalID found \n");
	         }
	      }
	   }

	   qsort(compactVect, n, sizeof(particle_t*), compare);

	   for (int rdx=0; rdx<cdx; rdx++){
	      particles_mpi[rdx] = *compactVect[rdx];
	   }

           save( fsave, n, particles_mpi );
	   }
	}

	//}
*/

	MPI_Barrier(MPI_COMM_WORLD);

    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //

    free ( displs );
    free ( rcounts );
    free ( displsC );
    free ( rcountsC );
    free ( particleVect );
    free ( flagVect );
    free ( nlocalVect );
    free( partition_offsets );
    free( partition_sizes );
    free( localBin );
    free( particles );
    free( particles_mpi);
    free( prevBin );
    free( nextBin );
    free( compactVect );
    free( prevBinSend );
    free( nextBinSend );
    free( localFlags );
    free( totalN );

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
