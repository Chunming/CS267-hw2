#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <math.h>
#include <string.h>
#include <vector>

#define subBlockLen 0.025

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
    double boxSize = sqrt(density * n); // 0.5 default
    int numBins = n_proc; // No. of bins set to 20
    double binLength = boxSize / numBins; // 0.5/24 = 0.020833 by default
    int subBlockNum = ceil(boxSize/subBlockLen); // 0.5/0.025=20  No. of sub blocks along a row/column 
    int binNum = subBlockNum*subBlockNum; // No. of bins (400)
    double bin_area = (boxSize*boxSize) / numBins; // Find max no. of particles per bin
    int nlocalMax = 3 * (int)( bin_area / (3.14*(cutoff/2)*(cutoff/2)) ); // Max particle num per proc
    int max_per_bin = floor( 2 * ((2*subBlockLen)/(sqrt(3)*cutoff)+1) * (subBlockLen/(2*cutoff) + 1) ) + 2;


    vector<double> xVect;
    vector<double> yVect;
    vector<int> globalIDVect;


    int* binParticleNum = new int [binNum]; // No. of particles in each bin
    if (binParticleNum == NULL) {
       printf("ERROR binParticleNum mem alloc failed \n");
       return -1;
    }

    particle_t* binArray = new particle_t [binNum*max_per_bin]; // Array of ptrs to particle*
    if (NULL == binArray) {
       printf("ERROR binArray mem alloc failed \n");
       return -1;
    }

    int* binFlags = new int [binNum*max_per_bin]; // Array of ptrs to particle*
    if (NULL == binFlags) {
       printf("ERROR binFlags mem alloc failed \n");
       return -1;
    }

   memset(binParticleNum, 0, sizeof(int)*binNum); // sizeof(binParticleNum[b] is 4 bytes
   memset(binArray, 0, binNum*max_per_bin*sizeof(particle_t));
   memset(binFlags, 0, binNum*max_per_bin*sizeof(int));

   particle_t** compactVect = (particle_t**)malloc( n * sizeof(particle_t*));
   particle_t* particleVect = (particle_t*)malloc( n_proc*nlocalMax*sizeof(particle_t));
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
    int xIdx, yIdx; // x/y index in 1D array
    int bdx;
    int offsetIdx;
    for (int ndx=0; ndx<n; ++ndx) { // For each particle
       xIdx = (particles[ndx].x)/subBlockLen; // Takes values 0-4
       yIdx = (particles[ndx].y)/subBlockLen; // Takes values 0-4
       if (xIdx == rank) {
          offsetIdx =  yIdx+(xIdx*subBlockNum);
    	  bdx = binParticleNum[offsetIdx];
          binArray[(offsetIdx*max_per_bin)+bdx] = particles[ndx]; // Offset from binArray
          binFlags[(offsetIdx*max_per_bin)+bdx] = 1; // Offset from binArray
          binParticleNum[offsetIdx]++; // Increment bin count
       }
    }
    int actual_n_proc = subBlockNum; // 20











    //
    //  Simulate a number of time steps
    //
    double simulation_time = read_timer( );
   // int count;
    int idx;
    int jdx;
    int kdx;

//printf("Check 1\n");




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
        int leftmostBins = rank;
        int rightmostBins = rank;
        double leftBnd, rightBnd, topBnd, botBnd;
        double leftDist, rightDist, topDist, botDist;
        int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;


    // First, EVEN will send, ODD will receive
     if (0 == rank%2 && rank<actual_n_proc) { // EVEN
           if (rank-1 >= 0) { // If left bin exists, then can send

           sIdx = 0; // Send index
           // Only traverse bins that are local to the bin
           for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin

              idx = 0;
              for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
                 while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values

                 xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
                 yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
                 leftBnd = xIdx*subBlockLen;
                 rightBnd = (xIdx*subBlockLen) + subBlockLen;
                 leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
                 rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);

                 if (leftDist<=cutoff && xIdx==leftmostBins) {
                    if(xIdx != 0 ) {
                       prevBinSend[sIdx] = binArray[(b*max_per_bin)+idx];
                       sIdx++;
                    }
                 }
              }
           }

              MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1, MPI_COMM_WORLD);
           }
     }


     if (1 == rank%2 && rank<actual_n_proc) { // ODD
        if (rank+1 <= actual_n_proc-1) { // If right bin exists, then can receive
           MPI_Recv(&nextBin[nNextBin], nlocalMax, PARTICLE, rank+1, tag1, MPI_COMM_WORLD, &status); // Recv from bot bin
           MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
           nNextBin += adjCount;
        }
     }


     if (0 == rank%2 && rank<actual_n_proc) { // EVEN
        if (rank+1 <= actual_n_proc-1) { // If bottom bin exists, then can send

           sIdx = 0;
           for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin

              idx = 0;
              for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
                 while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values

                 xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
                 yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
                 leftBnd = xIdx*subBlockLen;
                 rightBnd = (xIdx*subBlockLen) + subBlockLen;
                 leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
                 rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);

                 if (rightDist<=cutoff && xIdx==rightmostBins) {
                    if (xIdx != subBlockNum-1) { // Right subBlock is valid
                       nextBinSend[sIdx] = binArray[(b*max_per_bin)+idx];
                       sIdx++;
                    }
                 }
              }
           }

           MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+1, MPI_COMM_WORLD);
        }
     }


     if (1 == rank%2 && rank<actual_n_proc) { // ODD
        if (rank-1 >= 0) { // If top bin exists, then can receive
           MPI_Recv(&prevBin[nPrevBin], nlocalMax, PARTICLE, rank-1, tag1+1, MPI_COMM_WORLD, &status); //Recv from top bin
           MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
           nPrevBin += adjCount;
        }
     }


     // Now ODD will send, EVEN will receive
     if (1 == rank%2 && rank<actual_n_proc) { // ODD
        if (rank-1 >= 0) { // If top bin exists, then can send

           sIdx = 0; // Send index
           for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin

              idx = 0;
              for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
                 while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values

                 xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
                 yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
                 leftBnd = xIdx*subBlockLen;
                 rightBnd = (xIdx*subBlockLen) + subBlockLen;
                 leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
                 rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);

                 if (leftDist<=cutoff && xIdx==leftmostBins) {
                    if(xIdx != 0 ) {
                       prevBinSend[sIdx] = binArray[(b*max_per_bin)+idx];
                       sIdx++;
                    }
                 }
              }
           }


           MPI_Send(prevBinSend, sIdx, PARTICLE, rank-1, tag1+2, MPI_COMM_WORLD);
        }
     }


     if (0 == rank%2 && rank< actual_n_proc) { // EVEN
        if (rank+1 <= actual_n_proc-1) { // If bottom bin exists, then can receive
           MPI_Recv(&nextBin[nNextBin], nlocalMax, PARTICLE, rank+1, tag1+2, MPI_COMM_WORLD, &status); // Recv from bot bin
           MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
           nNextBin += adjCount;
        }
     }


     if (1 == rank%2 && rank<actual_n_proc) { // ODD 
        if (rank+1 <= actual_n_proc-1) { // If bottom bin exists, then can send

           sIdx = 0;
           for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin

              idx = 0;
              for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
                 while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values

                 xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
                 yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
                 leftBnd = xIdx*subBlockLen;
                 rightBnd = (xIdx*subBlockLen) + subBlockLen;
                 leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
                 rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);

                 if (rightDist<=cutoff && xIdx==rightmostBins) {
                    if (xIdx != subBlockNum-1) { // Right subBlock is valid
                       nextBinSend[sIdx] = binArray[(b*max_per_bin)+idx];
                       sIdx++;
                    }
                 }
              }
           }

           MPI_Send(nextBinSend, sIdx, PARTICLE, rank+1, tag1+3, MPI_COMM_WORLD); // Send to bot bin 
        }
     }


     if (0 == rank%2 && rank<actual_n_proc) { // EVEN
        if (rank-1 >= 0) { // If top bin exists, then can receive
           MPI_Recv(&prevBin[nPrevBin], nlocalMax, PARTICLE, rank-1, tag1+3, MPI_COMM_WORLD, &status); //Recv from top bin
           MPI_Get_count(&status, PARTICLE, &adjCount); // Get received count
           nPrevBin += adjCount;
        }
     }
     memset(prevBinSend, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
     memset(nextBinSend, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation





//printf("Check 2\n");

     // 
     // Bin adjacent particles for comparison
     //
    int xIdx, yIdx;
    int bdx;
    int offsetIdx;
    for (int i=0; i<nPrevBin; ++i) { // For each particle
       xIdx = (prevBin[i].x)/subBlockLen; // Takes values 0-4
       if (xIdx==rank) {printf("ERROR \n");}
       yIdx = (prevBin[i].y)/subBlockLen;

       //offsetIdx = yIdx+(xIdx*subBlockNum); 
       //bdx = binParticleNum[offsetIdx];
       //binArray[(offsetIdx*max_per_bin)+bdx] = prevBin[i];
       //binFlags[(offsetIdx*max_per_bin)+bdx] = 1;
       //binParticleNum[offsetIdx]++; // Increment bin count
    }

    for (int i=0; i<nNextBin; ++i) {
       xIdx = (nextBin[i].x)/subBlockLen;
       if (xIdx==rank) {printf("ERROR \n");}
       yIdx = (nextBin[i].y)/subBlockLen;

       //offsetIdx = yIdx+(xIdx*subBlockNum); 
       //bdx = binParticleNum[offsetIdx];
       //binArray[(offsetIdx*max_per_bin)+bdx] = nextBin[i];
       //binFlags[(offsetIdx*max_per_bin)+bdx] = 1;
       //binParticleNum[offsetIdx]++; // Increment bin count
    }
    //memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
    //memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
    //nPrevBin = 0;
    //nNextBin = 0;

MPI_Barrier(MPI_COMM_WORLD);
//printf("Check 3\n");




        //
        //  Compute forces
        //
        int idxCount = 0;
        for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin
           idx=0;
           for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin

               while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values
               (binArray[(b*max_per_bin)+idx]).ax = (binArray[(b*max_per_bin)+idx]).ay = 0;

              // Check all particles in bth subBlock
              jdx=0;
              for (int j=0; j<binParticleNum[b]; ++j) { // The jth particle in bth bin
                 while (0 == (binFlags[(b*max_per_bin)+jdx])) { jdx++; }
                 apply_force(binArray[(b*max_per_bin)+idx], binArray[(b*max_per_bin)+jdx]);
                 jdx++;
              }

              // Tmp operation
              for (int i=0; i<nPrevBin; ++i) { // For each particle
                 apply_force(binArray[(b*max_per_bin)+idx], prevBin[i]);
              }

              // Tmp operation
              for (int i=0; i<nNextBin; ++i) {
                 apply_force(binArray[(b*max_per_bin)+idx], nextBin[i]);
              }
/*

              // Check particles to the left/right/top/bottom subBlocks of bth subBlock  
              xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
              yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;

              leftBnd = xIdx*subBlockLen;
              rightBnd = (xIdx*subBlockLen) + subBlockLen;
              topBnd = yIdx*subBlockLen;
              botBnd = yIdx*subBlockLen + subBlockLen;

              leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
              rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);
              topDist = fabs((binArray[(b*max_per_bin)+idx]).y - topBnd);
              botDist = fabs((binArray[(b*max_per_bin)+idx]).y - botBnd);


              if (leftDist<=cutoff) {
                 if (xIdx != 0) { // Left subBlock index is valid
                    bLeft = b - subBlockNum;
                    kdx=0;
                    for (int k=0; k<binParticleNum[bLeft]; ++k) {
                       while (0 == (binFlags[(bLeft*max_per_bin)+kdx])) { kdx++; }
                       apply_force(binArray[(b*max_per_bin)+idx],binArray[(bLeft*max_per_bin)+kdx]);
                       kdx++;
                    }
                 }
              }

              if (rightDist<=cutoff) {
                 if (xIdx != subBlockNum-1) {
                    bRight = b + subBlockNum;
                    kdx=0;
                    for (int k=0; k<binParticleNum[bRight]; ++k) {
                       while (0 == (binFlags[(bRight*max_per_bin)+kdx])) { kdx++; }
                       apply_force(binArray[(b*max_per_bin)+idx],binArray[(bRight*max_per_bin)+kdx]);
                       kdx++;
                    }
                 }
              }

              if (topDist<=cutoff) {
                if (yIdx != 0) {
                  bTop = b - 1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bTop]; ++k) {
                    while (0 == (binFlags[(bTop*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bTop*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

              if (botDist<=cutoff) {
                if (yIdx != subBlockNum-1) {
                  bBottom = b + 1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bBottom]; ++k) {
                    while (0 == (binFlags[(bBottom*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bBottom*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

             if (topDist<=cutoff && leftDist<=cutoff) {
                if (yIdx != 0 && xIdx !=0) {
                  bTopLeft = b-subBlockNum-1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bTopLeft]; ++k) {
                    while (0 == (binFlags[(bTopLeft*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bTopLeft*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

              if (botDist<=cutoff && leftDist<=cutoff) {
                if (yIdx!=subBlockNum-1 && xIdx != 0) {
                  bBotLeft = b-subBlockNum+1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bBotLeft]; ++k) {
                    while (0 == (binFlags[(bBotLeft*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bBotLeft*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

              if (topDist<=cutoff && rightDist<=cutoff) {
                if (yIdx != 0 && xIdx != subBlockNum-1) {
                  bTopRight = b+subBlockNum-1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bTopRight]; ++k) {
                    while (0 == (binFlags[(bTopRight*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bTopRight*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

              if (botDist<=cutoff && rightDist<=cutoff) {
                if (yIdx!=subBlockNum-1 && xIdx!=subBlockNum-1) {
                  bBotRight = b+subBlockNum+1;
                  kdx=0;
                  for (int k=0; k<binParticleNum[bBotRight]; ++k) {
                    while (0 == (binFlags[(bBotRight*max_per_bin)+kdx])) { kdx++; }
                    apply_force(binArray[(b*max_per_bin)+idx],binArray[(bBotRight*max_per_bin)+kdx]);
                    kdx++;
                  }
                }
              }

*/
               idx++;
            }

           idxCount += binParticleNum[b];
        }
        memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
        memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
        nPrevBin = 0;
        nNextBin = 0;

//printf("Check 4\n");
MPI_Barrier(MPI_COMM_WORLD);


        //   
        //  Move particles
        //
       for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); ++b) { // The bth bin
           for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin
              move( binArray[(b*max_per_bin)+i] );
           } 
        }

//printf("Check 5\n");
MPI_Barrier(MPI_COMM_WORLD);

        //
        // Re-bin particles
        //
        int index;
        int jdx = 0;
        int kdx = 0;
        for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); b++) { // The bth bin
           idx = 0; // Skips across
           for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin

              // Check if particle shld be moved to another bin   
              while (0 == (binFlags[(b*max_per_bin)+idx])) { idx++; }
              xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
              yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
              index = yIdx+(xIdx*subBlockNum); // Map 2D to 1D index
              if (index != b) { // Particle has moved out of of the bin

                 // Store into first ZERO index in array
                 bdx = 0;
                 while (0 != binFlags[(index*max_per_bin)+bdx]) {
                    bdx++;
                    if (bdx > n) {
                       printf("ERROR: Overflow \n");
                       return -1;
                    }
                 }

                 // Outside the bins bounded by the processor
                 if (index < rank*subBlockNum && index>=(rank-1)*subBlockNum ) { // Left side
                    prevBinSend[jdx] = binArray[(b*max_per_bin)+idx];
                    jdx++;
                    binFlags[(b*max_per_bin)+idx] = 0;
                    binParticleNum[b]--;
                    i--;
                 }
                 else if (index >= (rank*subBlockNum+subBlockNum) && index<((rank+1)*subBlockNum+subBlockNum)) { // Right side
                    nextBinSend[kdx] = binArray[(b*max_per_bin)+idx];
                    kdx++;
                    binFlags[(b*max_per_bin)+idx] = 0;
                    binParticleNum[b]--;
                    i--;
                 }
                 else if (index>=rank*subBlockNum && index<(rank*subBlockNum+subBlockNum)) { // Same side
                    binArray[(index*max_per_bin)+bdx] = binArray[(b*max_per_bin)+idx]; // Add element into first non-NULL index           
                    binFlags[(b*max_per_bin)+idx] = 0;
                    binFlags[(index*max_per_bin)+bdx] = 1;
                    binParticleNum[b]--;
                    i--;
                    binParticleNum[index]++;
                 }
		 else {
		    printf("BIG ERROR, question assumptions\n");

		 }
              }
              idx++;
           }
        }
//        nPrevBin = jdx; // No. of elems to shift from prevBin to localBin
//        nNextBin = kdx; // No. of elems to shift from nextBin to localBin

MPI_Barrier(MPI_COMM_WORLD);
//printf("Check 6\n");

        //
        // 5. Compact Particles
        //
        int nextLoc;
        for (int b=rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); b++) { // The bth bin
           for (int i=0; i<binParticleNum[b]; ++i) {
              if (0 == binFlags[(b*max_per_bin)+i]) {
                 nextLoc = i + 1;
                 while (0 == binFlags[(b*max_per_bin)+nextLoc]) {
                    nextLoc++;
                    if (nextLoc == max_per_bin) {nextLoc = 0;}
                 }
                 binArray[(b*max_per_bin)+i] = binArray[(b*max_per_bin)+nextLoc];
                 binFlags[(b*max_per_bin)+i] = 1;
                 binFlags[(b*max_per_bin)+nextLoc] = 0;
              }
           }
        }


MPI_Barrier(MPI_COMM_WORLD);
//printf("Check 7\n");

	//
        // 6. Get particles from adj bins. Use sync blocking send/receive
	//
        int rebinCount;
        int tag4 = 400;
	nPrevBin = 0;
	nNextBin = 0;
        // Have EVEN processors send particles first
        // and ODD processors receive
        if (0 == rank%2 && rank < actual_n_proc) { // Even Processors  
           if (rank-1 >= 0) { // If left bin exists, send
              MPI_Send(prevBinSend, jdx, PARTICLE, rank-1, tag4, MPI_COMM_WORLD); //Send to left bin
           }
         }

        if (1 == rank%2 && rank < actual_n_proc) { // Odd Processors
           if (rank+1 <= actual_n_proc-1) { // If right bin exists, receive
              MPI_Recv(nextBin, nlocalMax, PARTICLE, rank+1, tag4, MPI_COMM_WORLD, &status); //Recv from right bin
              MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
              nNextBin = rebinCount;
           }
        }

        if (0 == rank%2 && rank < actual_n_proc) { // Even Processors  
           if (rank+1 <= actual_n_proc-1) { // If right bin exists, send 
              MPI_Send(nextBinSend, kdx, PARTICLE, rank+1, tag4+1, MPI_COMM_WORLD); //Send to bot bin
           }
        }

        if (1 == rank%2 && rank < actual_n_proc) { // Odd Processors
           if (rank-1 >=0) { // If left bin exists, receive
              MPI_Recv(prevBin, nlocalMax, PARTICLE, rank-1, tag4+1, MPI_COMM_WORLD, &status);
              MPI_Get_count(&status, PARTICLE, &rebinCount);
              nPrevBin = rebinCount;
           }
        }

        // Have ODD processors send next
        // and EVEN processors receive
        if (1 == rank%2  && rank < actual_n_proc) { // Odd Processors  
           if (rank-1 >= 0) { // If top bin exists, send
              MPI_Send(prevBinSend, jdx, PARTICLE, rank-1, tag4+2, MPI_COMM_WORLD); //Send to bot bin
           }
        }

        if (0 == rank%2  && rank < actual_n_proc) { // Even Processors
           if (rank+1 <= actual_n_proc-1) { // If bottom bin exists, receive
              MPI_Recv(nextBin, nlocalMax, PARTICLE, rank+1, tag4+2, MPI_COMM_WORLD, &status); //Recv from top bin
              MPI_Get_count(&status, PARTICLE, &rebinCount); // Get received count
              nNextBin = rebinCount;
           }
        }

        if (1 == rank%2 && rank < actual_n_proc) { // Odd Processors  
           if (rank+1 <= actual_n_proc-1) { // If bottom bin exists, send 
              MPI_Send(nextBinSend, kdx, PARTICLE, rank+1, tag4+3, MPI_COMM_WORLD); //Send to bot bin
           }
        }

        if (0 == rank%2 && rank < actual_n_proc)  { // Even Processors
           if (rank-1 >=0) { // If top bin exists, receive
              MPI_Recv(prevBin, nlocalMax, PARTICLE, rank-1, tag4+3, MPI_COMM_WORLD, &status); //Recv from bot bin
              MPI_Get_count(&status, PARTICLE, &rebinCount);
              nPrevBin = rebinCount;
           }
        }

MPI_Barrier(MPI_COMM_WORLD);
//printf("Check 8\n");

     // 
     // Bin adjacent particles for comparison
     //
     for (int i=0; i<nPrevBin; ++i) { // For each particle
        xIdx = (prevBin[i].x)/subBlockLen; // Takes values 0-4
	printf("Prev Test xIdx: %d, rank: %d \n", xIdx, rank);
        if (xIdx!=rank) {printf("ERROR 3 transfered elems shld be in right rank \n");}
        yIdx = (prevBin[i].y)/subBlockLen;

        offsetIdx = yIdx+(xIdx*subBlockNum); 
        bdx = binParticleNum[offsetIdx];
        binArray[(offsetIdx*max_per_bin)+bdx] = prevBin[i];
        binFlags[(offsetIdx*max_per_bin)+bdx] = 1;
        binParticleNum[offsetIdx]++; // Increment bin count
     }

     for (int i=0; i<nNextBin; ++i) {
        xIdx = (nextBin[i].x)/subBlockLen;
	printf("Next Test xIdx: %d, rank: %d \n", xIdx, rank);
        if (xIdx!=rank) {printf("ERROR 4 transfered elems shld be in right rank \n");}
        yIdx = (nextBin[i].y)/subBlockLen;

        offsetIdx = yIdx+(xIdx*subBlockNum); 
        bdx = binParticleNum[offsetIdx];
        binArray[(offsetIdx*max_per_bin)+bdx] = nextBin[i];
        binFlags[(offsetIdx*max_per_bin)+bdx] = 1;
        binParticleNum[offsetIdx]++; // Increment bin count
     }

        memset(prevBin, 0, nPrevBin*sizeof(particle_t)); // Reset prevBin ptr for next itereation
        memset(nextBin, 0, nNextBin*sizeof(particle_t)); // Reset nextBin ptr for next itereation
        nPrevBin = 0; // No. of elems to shift from prevBin to localBin
        nNextBin = 0; // No. of elems to shift from nextBin to localBin


MPI_Barrier(MPI_COMM_WORLD);


	//
   	// Check total num of particles
	//
	int count = 0;
        for (int b = rank*subBlockNum; b<(rank*subBlockNum+subBlockNum); b++) { // Bins 0 to 399 
	   count += binParticleNum[b];	
	}
   	MPI_Reduce(&count, totalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );        
   	if (rank == 0) printf("Total N is %d \n", *totalN);


MPI_Barrier(MPI_COMM_WORLD);


//printf("Check 9\n");

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



	MPI_Barrier(MPI_COMM_WORLD);
*/
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

    delete [] binArray;
    delete [] binFlags;
    delete [] binParticleNum;

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
