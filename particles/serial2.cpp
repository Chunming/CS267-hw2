#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <string.h>
#include <vector>

#define subBlockLen 0.025

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
 
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
    // Set up MPI
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
    double boxSize = sqrt(density*n); // Length on each side of box 
    int numBins = n_proc; // No. of bins 
    double binLength = boxSize / numBins; // 0.5/24 = 0.020833 by default
    int subBlockNum = ceil(boxSize/subBlockLen); // No. of sub blocks along a row/column (Default: 5)
    int binNum = subBlockNum*subBlockNum; // No. of bins


    //
    // Set up serial code
    //  
    int* binParticleNum = new int [binNum]; // No. of particles in each bin
    if (binParticleNum == NULL) {
       printf("ERROR binParticleNum mem alloc failed \n");
       return -1;
    }

    particle_t*** binArray = new particle_t** [binNum]; // Array of ptrs to particle*
    if (NULL == binArray) {
       printf("ERROR binArray mem alloc failed \n");
       return -1;
    }

   for (int i = 0; i<binNum; ++i) {binParticleNum[i] = 0;}
   memset(binParticleNum, 0, sizeof(int)*binNum); // sizeof(binParticleNum[b] is 4 bytes
   int max_per_bin = floor( 2 * ((2*subBlockLen)/(sqrt(3)*cutoff)+1) * (subBlockLen/(2*cutoff) + 1) ) + 2;

   for (int b=0; b<binNum; ++b) {
            binArray[b] = new particle_t* [max_per_bin]; // Set it to max 100 particles first
            if (NULL == binArray[b]) {
               printf("ERROR binArray index  mem alloc failed \n");
               return -1;
            }
            memset(binArray[b], NULL, max_per_bin*sizeof(particle_t*));
    }

    double bin_area = (boxSize*boxSize) / numBins; // Find max no. of particles per bin
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
        rcounts[i] = 100;
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
    if( rank == 0 ) {init_particles( n, particles );}
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






    // Initialize particle binning
    int xIdx, yIdx; // x/y index in 1D array
    int bdx; 

    for (int ndx=0; ndx<n; ++ndx) { // For each particle
       xIdx = (particles[ndx].x)/subBlockLen; // Takes values 0-4
       yIdx = (particles[ndx].y)/subBlockLen; // Takes values 0-4
//       bdx = 0;
//       while (binArray[yIdx+(xIdx*subBlockNum)][bdx]!=NULL) {
//          bdx++;
//       }
//       printf("Bin index is %d \n", yIdx+xIdx*subBlockNum);
       bdx = binParticleNum[yIdx+(xIdx*subBlockNum)];
       binArray[yIdx+(xIdx*subBlockNum)][bdx] = &particles[ndx];
       binParticleNum[yIdx+(xIdx*subBlockNum)]++; // Increment bin count
    }

    //
    //  Simulate a number of time steps
    //
    double simulation_time = read_timer( );
    int count;
    int idx;
    int jdx;
    int kdx;
    for( int step = 0; step < NSTEPS; step++ )
    {
	    

	//
        //  Compute forces
        //
	/*
	// Orginal apply_force function
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j] );
        }
	*/

        // Set particle accelerations to 0 before using apply_force
//        for( int i = 0; i < 500; i++ ) 
//		particles[i].ax = particles[i].ay = 0;

	// Compute Forces
	double leftBnd, rightBnd, topBnd, botBnd;
	double leftDist, rightDist, topDist, botDist;
	int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;
        int idxCount = 0;
	for (int b=0; b<binNum; ++b) { // The bth bin

	   idx=0;
	   for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
	     

	       while ((NULL == binArray[b][idx])) { idx++; } // Index skips NULL values
      	       (*binArray[b][idx]).ax = (*binArray[b][idx]).ay = 0;

	      // Check all particles in bth subBlock
	      jdx=0;
     	      for (int j=0; j<binParticleNum[b]; ++j) { // The jth particle in bth bin
     	         while (NULL == (binArray[b][jdx])) { jdx++; }
		 apply_force(*binArray[b][idx],*binArray[b][jdx]);
		 jdx++;
	      }

	      // Check particles to the left/right/top/bottom subBlocks of bth subBlock	 
	      //printf("xIdx is %d \n", xIdx);
	      xIdx = (*binArray[b][idx]).x/subBlockLen;
	      yIdx = (*binArray[b][idx]).y/subBlockLen;

	      //printf("botBnd is %f \n", botBnd);
	      leftBnd = xIdx*subBlockLen;
	      rightBnd = (xIdx*subBlockLen) + subBlockLen;
	      topBnd = yIdx*subBlockLen;
	      botBnd = yIdx*subBlockLen + subBlockLen;


	      //printf("botDist is %f \n", botDist);
	      leftDist = fabs((*binArray[b][idx]).x - leftBnd);
	      rightDist = fabs((*binArray[b][idx]).x - rightBnd);
	      topDist = fabs((*binArray[b][idx]).y - topBnd);
	      botDist = fabs((*binArray[b][idx]).y - botBnd);

	      // Consider 8 different adjacent subBlocks
	      if (leftDist<=cutoff) {
		 if (xIdx != 0) { // Left subBlock index is valid
		    bLeft = b - subBlockNum;
		    kdx=0;
		    for (int k=0; k<binParticleNum[bLeft]; ++k) { 
		       while (NULL == (binArray[bLeft][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bLeft][kdx]);
		       kdx++;
		    }
		 }
	      }


	      if (rightDist<=cutoff) {
		 if (xIdx != subBlockNum-1) { 
  		    bRight = b + subBlockNum; 
		    kdx=0;
		    for (int k=0; k<binParticleNum[bRight]; ++k) { 
		       while (NULL == (binArray[bRight][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bRight][kdx]);
		       kdx++;
		    }
		 }

	      }


	      if (topDist<=cutoff) {
		if (yIdx != 0) { 
		  bTop = b - 1; 
		  kdx=0;
		  for (int k=0; k<binParticleNum[bTop]; ++k) { 
		    while (NULL == (binArray[bTop][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bTop][kdx]);
		    kdx++;
		  }
		}
	      }


	      if (botDist<=cutoff) {
		if (yIdx != subBlockNum-1) {
		  bBottom = b + 1;
		  kdx=0;
		  for (int k=0; k<binParticleNum[bBottom]; ++k) { 
		    while (NULL == (binArray[bBottom][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bBottom][kdx]);
		    kdx++;
		  }
		}
	      }




	      if (topDist<=cutoff && leftDist<=cutoff) { 
		if (yIdx != 0 && xIdx !=0) {
		  bTopLeft = b-subBlockNum-1;     
		  kdx=0;
		  for (int k=0; k<binParticleNum[bTopLeft]; ++k) { 
		    while (NULL == (binArray[bTopLeft][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bTopLeft][kdx]);
		    kdx++;
		  }
		}
	      }



	      if (botDist<=cutoff && leftDist<=cutoff) { 
		if (yIdx!=subBlockNum-1 && xIdx != 0) { 
		  bBotLeft = b-subBlockNum+1;     
		  kdx=0;
		  for (int k=0; k<binParticleNum[bBotLeft]; ++k) { 
		    while (NULL == (binArray[bBotLeft][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bBotLeft][kdx]);
		    kdx++;
		  }
		}
	      }



	      if (topDist<=cutoff && rightDist<=cutoff) { 
		if (yIdx != 0 && xIdx != subBlockNum-1) {
		  bTopRight = b+subBlockNum-1;     
		  kdx=0;
		  for (int k=0; k<binParticleNum[bTopRight]; ++k) { 
		    while (NULL == (binArray[bTopRight][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bTopRight][kdx]);
		    kdx++;
		  }
		}
	      }

	      if (botDist<=cutoff && rightDist<=cutoff) { 
		if (yIdx!=subBlockNum-1 && xIdx!=subBlockNum-1) {
		  bBotRight = b+subBlockNum+1;     
		  kdx=0;

		  for (int k=0; k<binParticleNum[bBotRight]; ++k) { 
		    while (NULL == (binArray[bBotRight][kdx])) { kdx++; }
		    apply_force(*binArray[b][idx],*binArray[bBotRight][kdx]);
		    kdx++;
		  }
		}
	      }

               idx++;
            }

	   idxCount += binParticleNum[b];
	}

        //printf("The idxCount is %d \n", idxCount);
        
	
	//
        //  Move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        
	//
	// Re-bin particles
	//
	int index;
	for (int b=0; b<binNum; b++) { // The bth bin
	   idx = 0; // Skips across
	   for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin
	      
	      // Check if particle shld be moved to another bin	  
	     
     	      while (NULL == (binArray[b][idx])) { idx++; }
	      xIdx = (*binArray[b][idx]).x/subBlockLen;
              yIdx = (*binArray[b][idx]).y/subBlockLen;
	      index = yIdx+(xIdx*subBlockNum); // Map 2D to 1D index
	      if (index != b) { // Particle has moved out of of the bin

                 // Store into first non-NULL index in array
		 bdx = 0;
                 while (NULL != binArray[index][bdx]) {
	            bdx++;
                    if (bdx > n) {
                       printf("ERROR: Overflow \n");
                       return -1;
		    }	 
		 }
                 binArray[index][bdx] = binArray[b][idx]; // Add element into first non-NULL index   	
		 binArray[b][idx] = NULL;
		 binParticleNum[b]--;
		 i--;
		 binParticleNum[index]++;

	      }
	      idx++;
	   }
	}

	// Check count
	count = 0;
	for (int b=0; b<binNum; ++b) {
	   count += binParticleNum[b];
	}

        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {

            save( fsave, n, particles );

         }
    }
    simulation_time = read_timer( ) - simulation_time;


    if(rank==0)
       printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );

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

   
    delete [] binParticleNum;
    for (int b=0; b<binNum; ++b) {
       delete [] binArray[b];
    }
    delete binArray;


    free( particles );
    if( fsave )
        fclose( fsave );
  
    MPI_Finalize();
  
    return 0;
}
