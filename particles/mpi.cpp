#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <math.h>
#include <string.h>

int isCloseToEdge(particle_t &particle, double binEdge, double cutoff) {
    double dy = binEdge - particle.y;
    double r2 = dy * dy;
    if (r2 > cutoff * cutoff) 
       return 0;
    else 
       return 1;
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
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
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

       double cutoff = 0.01;
       double density = 0.0005;
       double spaceDim = sqrt(density * n); // 0.5 default
       int numBins = n_proc; // No. of bins 
       double binLength = spaceDim / numBins; // 0.5/24 = 0.020833 by default
       double bin_area = (spaceDim*spaceDim) / numBins; // Find max no. of particles per bin
       int nlocalMax = 3 * (int)( bin_area / (3.14*(cutoff/2)*(cutoff/2)) ); // Max particle num per proc


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








 
    // Start binning
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
   free( particles ); // Particle array is not used after binning


   particle_t* gParticles;
   int* gFlags;
   gParticles = (particle_t*)malloc( n_proc*30*sizeof(particle_t));
   if (NULL == gParticles) {
      printf("ERR allocating *localFlags \n");
      return -1;
   }

   gFlags = (int*)malloc( n_proc*30*sizeof(int));
   if (NULL == gFlags) {
      printf("ERR allocating *localFlags \n");
      return -1;
   }
   memset(gFlags, 0, n_proc*30*sizeof(int));
   memset(gParticles, 0, n_proc*30*sizeof(particle_t));


   int* displs = (int *)malloc(n_proc*sizeof(int)); 
   int* rcounts = (int *)malloc(n_proc*sizeof(int)); 
   for (int i=0; i<n_proc; ++i) { 
       displs[i] = i*30; 
       rcounts[i] = 30; 
   } 






  //
    //  simulate a number of time steps
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
        //bool fPrevCheck = 0;
        //bool fNextCheck = 0;
        MPI_Status status;

        // First, EVEN will send, ODD will receive
        if (0 == rank%2) { // EVEN
           if (rank-1 >= 0) { // If top bin exists, then can send
              idx = 0;
	      sIdx = 0;
              //fPrevCheck = 0; // Check if comms with prev bin is req
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge, cutoff)) { // Comms with prev bin is req
		       prevBinSend[sIdx] = localBin[idx]; // Set send index
		       sIdx++;
                       //fPrevCheck = 1;
                       //break;
                 }
                 idx++;
              }

              //if (1==fPrevCheck) {
                 MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1, MPI_COMM_WORLD);
                 //printf(" %d particles from rank %d is sent to rank %d \n", nlocal, rank, rank-1);
              //}
              //else {
                 //MPI_Send(localBin, 0, PARTICLE, rank-1, tag1, MPI_COMM_WORLD);
                 //printf(" 0 particles from rank %d is sent to rank %d \n", rank, rank-1);
              //}
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
              //fNextCheck = 0;
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge+binLength, cutoff)) {
		    nextBinSend[sIdx] = localBin[idx];
		    sIdx++;
                    //fNextCheck = 1;
                    //break;
                 }
                 idx++;
              }

              //if (1==fNextCheck) 
		MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+1, MPI_COMM_WORLD); // Send to bot bin
              //else MPI_Send(localBin, 0, PARTICLE, rank+1, tag1+1, MPI_COMM_WORLD); // Send to bot bin
           }
        }

        if (1 == rank%2) { // ODD
           if (rank-1 >= 0) { // If top bin exists, then can receive
              // Check receive signal from prevBin
              MPI_Recv(&prevBin[nPrevBin], nlocalMax, PARTICLE, rank-1, tag1+1, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
           }
        }



        // Now ODD will send, EVEN will receive
        if (1 == rank%2) { // ODD
           if (rank-1 >= 0) { // If top bin exists, then can send
              idx = 0;
	      sIdx = 0;
              //fPrevCheck = 0; // Check if comms with prev bin is req
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge, cutoff)) { // Comms with prev bin is req
		       prevBinSend[sIdx] = localBin[idx];
		       sIdx++;
                       //fPrevCheck = 1;
                       //break;
                 }
                 idx++;
              }

              //if (1==fPrevCheck) 
	      MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1+2, MPI_COMM_WORLD);
              //else MPI_Send(localBin, 0, PARTICLE, rank-1, tag1+2, MPI_COMM_WORLD);
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
              //fNextCheck = 0;
              for (int i=0; i< (nlocal); ++i) { // Only send particles that are close to edge
                 while (0==localFlags[idx]) idx++;
                 if (isCloseToEdge(localBin[idx], binEdge+binLength, cutoff)) {
		    nextBinSend[sIdx] = localBin[idx];
		    sIdx++;
                    //fNextCheck = 1;
                    //break;
                 }
                 idx++;
              }

              //if (1==fNextCheck) 
	      MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+3, MPI_COMM_WORLD); // Send to bot bin
              //else MPI_Send(localBin, 0, PARTICLE, rank+1, tag1+3, MPI_COMM_WORLD); // Send to bot bin
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
        //memset(prevBinSend, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
        //memset(nextBinSend, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation


	// 
        // 2. Apply Force
        //
        int loc_i = 0;
        int loc_j = 0;
        for (int i=0; i<(nlocal); ++i) { // For each particle in local bin
           while(0==localFlags[loc_i]) {
              loc_i++;
              if (loc_i == nlocalMax) {loc_i = 0;}
	   }	
           localBin[loc_i].ax = localBin[loc_i].ay = 0;

           // SELF LOOP - Compute interactions with particles within the bin
           loc_j = 0;
           for (int j=0; j<(nlocal); ++j) {
              while(0==localFlags[loc_j]) { 
	         loc_j++;
		 if (loc_j == nlocalMax) {loc_j = 0;}
	      }
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


        // 
        // 3. Move Particles
        //
        loc_i = 0;
        for (int i=0; i<(nlocal); ++i) {
           while(0==localFlags[loc_i]) {
	      loc_i++;
	      if (loc_i == nlocalMax) {loc_i=0;}
	   }
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
           while(0==localFlags[idx]) { 
	      idx++; 
	      if (idx == nlocalMax) {idx = 0;}
	   }
           bdx = (localBin[idx].y / binLength);

           if (bdx == rank-1) { // Particle moved to top bin
              localFlags[idx] = 0; // Remove from localBin
              prevBin[jdx] = localBin[idx];
              jdx++;
           }

           if (bdx == rank+1) { // Particle moved to top bin
              localFlags[idx] = 0; // Remove from localBin
              nextBin[kdx] = localBin[idx];
              kdx++;
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

        // Get particles from adjacent bins
        // Using synchronous blocking send/receive
        int rebinCount;
        int tag4 = 400;

/*
   if (rank == 1) {
      int tdx = 0;
      int udx = 0;
      printf("BEF nlocal is %d \n", nlocal);
      while (tdx < nlocal) {
         if (0!=localFlags[udx]) {
            printf("BEF udx is %d \n", udx);
            tdx++;
         }
	 printf("BEF udx is incremented \n ");
	 udx++;
      }
   }
*/



        // Have EVEN processors send particles first
        // and ODD processors receive
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
              MPI_Send(prevBin, nPrevBin, PARTICLE, rank-1, tag4+2, MPI_COMM_WORLD); //Send to bot bin
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
              MPI_Recv(&localBin[nlocal], nlocalMax, PARTICLE, rank-1, tag4+3, MPI_COMM_WORLD, &status); //Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &rebinCount);
              for (int j=nlocal; j<(nlocal+rebinCount); ++j) {localFlags[j]=1;}
	      nlocal += rebinCount;
	   }
	}	

        memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
        memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
        nPrevBin = 0; // No. of elems to shift from prevBin to localBin
        nNextBin = 0; // No. of elems to shift from nextBin to localBin


/*
   if (rank == 1) {
      int tdx = 0;
      int udx = 0;
      printf("AFT nlocal is %d \n", nlocal);
      while (tdx < nlocal) {
         if (0!=localFlags[udx]) {
            printf("AFT udx is %d \n", udx);
            tdx++;
         }
	 printf("AFT udx is incremented \n ");
	 udx++;
      }
   }
*/

   printf("Timestep %d nlocal at rank %d is %d \n", step, rank, nlocal);
   // Check total num of particles
   MPI_Reduce(&nlocal, totalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );        
   if (rank == 0) printf("Total N is %d \n", *totalN);


	// 
        //  collect all global data locally (not good idea to do)
        //  gathers data from all tasks & deliver combined data to all tasks
	//
	//MPI_Allgatherv( localBin, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );


	//
        //  save current step if necessary
        //

//MPI_Barrier(MPI_COMM_WORLD);

        if( fsave && (step%SAVEFREQ) == 0 ) {
            //save( fsave, n, particles );

//	      if (rank != 0) {
//                 MPI_Send(localBin, nlocal, PARTICLE, 0, 500, MPI_COMM_WORLD); //Send to bot bin
//	      }

//	      if (rank == 0) {
//		 int pIdx = 0;
//	         for (int p=1; p<n_proc; p++) {
  //                  MPI_Recv(&gParticles[pIdx], nlocalMax, PARTICLE, p, 500, MPI_COMM_WORLD, &status); //Recv from bot bin
//                    MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
                    //for (int j=pIdx; j<(pIdx+rebinCount); ++j) {gFlags[j]=1;}
  //                  pIdx += rebinCount;
//	         }
//	      }

           // Check answer
           // MPI_Gatherv( sendarray, 100, MPI_INT, rbuf, rcounts, displs, MPI_INT, root, comm); 
          //MPI_Gatherv(localBin, nlocal, PARTICLE, gParticles, rcounts, displs, PARTICLE, 0, MPI_COMM_WORLD);
         // MPI_Gatherv(localFlags, nlocal, MPI_INT, gFlags, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

/*
        int nextLoc;
         for (int loc=0; loc<n; ++loc) {
            if (0 == gFlags[loc]) {
               nextLoc = loc + 1;
               while (0 == gFlags[nextLoc]) { 
	          nextLoc++;
	       }
               gParticles[loc] = gParticles[nextLoc];
               gFlags[loc] = 1;
               gFlags[nextLoc] = 0;
            }  
         }  

       for (int zdx=0; zdx< (n_proc*nlocalMax); ++zdx) {
          if (gFlags!=0) printf("gFlag index is %d \n", zdx);     
       }

*/

//MPI_Barrier(MPI_COMM_WORLD);





	}

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
    free ( gParticles );
    free ( gFlags );
    free( partition_offsets );
    free( partition_sizes );
    free( localBin );
    // free( particles );
    free( prevBin );
    free( nextBin );

    free( prevBinSend );
    free( nextBinSend );


    free( localFlags );
    free( totalN );

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
