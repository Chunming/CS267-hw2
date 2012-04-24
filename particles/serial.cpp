#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <string.h>


#define cutoff 0.01 // Was defined in common.cpp
#define density 0.0005
#define subBlockLen 0.025


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
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );

    double boxSize = sqrt(density*n); // Length on each side of box 
    int subBlockNum = ceil(boxSize/subBlockLen); // No. of sub blocks along a row/column (Default: 5)
    int binNum = subBlockNum*subBlockNum; // No. of bins


   int max_per_bin = floor( 2 * ((2*subBlockLen)/(sqrt(3)*cutoff)+1) * (subBlockLen/(2*cutoff) + 1) ) + 2;

    // Determine number of sublocks, represent as a vector;
    // Set no. of blocks as 5 x 5
    int* binParticleNum = new int [binNum]; // No. of particles in each bin
    if (binParticleNum == NULL) {
       printf("ERROR binParticleNum mem alloc failed \n");
       return -1;
    }

    // binNum * max_per_bin
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

   particle_t** compactVect = (particle_t**)malloc( n * sizeof(particle_t*));
   particle_t *particles_mpi = (particle_t*) malloc( n * sizeof(particle_t) );


   memset(binParticleNum, 0, sizeof(int)*binNum); // sizeof(binParticleNum[b] is 4 bytes
   memset(binArray, 0, binNum*max_per_bin*sizeof(particle_t));
   memset(binFlags, 0, binNum*max_per_bin*sizeof(int));

    // Initialize particles 
    init_particles( n, particles );

    // Initialize particle binning
    int xIdx, yIdx; // x/y index in 1D array
    int bdx; 
    int offsetIdx;
    for (int ndx=0; ndx<n; ++ndx) { // For each particle
       xIdx = (particles[ndx].x)/subBlockLen; // Takes values 0-4
       yIdx = (particles[ndx].y)/subBlockLen; // Takes values 0-4
//       bdx = 0;
//       while (binArray[yIdx+(xIdx*subBlockNum)][bdx]!=NULL) {
//          bdx++;
//       }
//       printf("Bin index is %d \n", yIdx+xIdx*subBlockNum);
       bdx = binParticleNum[yIdx+(xIdx*subBlockNum)];
       offsetIdx =  yIdx+(xIdx*subBlockNum);
       binArray[(offsetIdx*max_per_bin)+bdx] = particles[ndx]; // Offset from binArray
       binFlags[(offsetIdx*max_per_bin)+bdx] = 1; // Offset from binArray
       binParticleNum[offsetIdx]++; // Increment bin count
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
	    
	// Compute Forces
	double leftBnd, rightBnd, topBnd, botBnd;
	double leftDist, rightDist, topDist, botDist;
	int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;
        int idxCount = 0;
	for (int b=0; b<binNum; ++b) { // The bth bin

	   idx=0;
	   for (int i=0; i<binParticleNum[b]; ++i) { // The ith particle in bth bin
	     

	       while ((0 == binFlags[(b*max_per_bin)+idx])) { idx++; } // Index skips NULL values
      	       (binArray[(b*max_per_bin)+idx]).ax = (binArray[(b*max_per_bin)+idx]).ay = 0;

	      // Check all particles in bth subBlock
	      jdx=0;
     	      for (int j=0; j<binParticleNum[b]; ++j) { // The jth particle in bth bin
     	         while (0 == (binFlags[(b*max_per_bin)+jdx])) { jdx++; }
		 apply_force(binArray[(b*max_per_bin)+idx],binArray[(b*max_per_bin)+jdx]);
		 jdx++;
	      }

	      // Check particles to the left/right/top/bottom subBlocks of bth subBlock	 
	      //printf("xIdx is %d \n", xIdx);
	      xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
	      yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;

	      //printf("botBnd is %f \n", botBnd);
	      leftBnd = xIdx*subBlockLen;
	      rightBnd = (xIdx*subBlockLen) + subBlockLen;
	      topBnd = yIdx*subBlockLen;
	      botBnd = yIdx*subBlockLen + subBlockLen;


	      //printf("botDist is %f \n", botDist);
	      leftDist = fabs((binArray[(b*max_per_bin)+idx]).x - leftBnd);
	      rightDist = fabs((binArray[(b*max_per_bin)+idx]).x - rightBnd);
	      topDist = fabs((binArray[(b*max_per_bin)+idx]).y - topBnd);
	      botDist = fabs((binArray[(b*max_per_bin)+idx]).y - botBnd);

	      // Consider 8 different adjacent subBlocks
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

               idx++;
            }

	   idxCount += binParticleNum[b];
	}

	
	//
        //  Move particles
        //
       for (int b=0; b<binNum; b++) { // The bth bin
           for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin
              move( binArray[(b*max_per_bin)+i] );
           }
        }

 
	//
	// Re-bin particles
	//
	int index;
	for (int b=0; b<binNum; b++) { // The bth bin
	   idx = 0; // Skips across
	   for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin
	      
	      // Check if particle shld be moved to another bin	  
	     
     	      while (0 == (binFlags[(b*max_per_bin)+idx])) { idx++; }
	      xIdx = (binArray[(b*max_per_bin)+idx]).x/subBlockLen;
              yIdx = (binArray[(b*max_per_bin)+idx]).y/subBlockLen;
	      index = yIdx+(xIdx*subBlockNum); // Map 2D to 1D index
	      if (index != b) { // Particle has moved out of of the bin

                 // Store into first non-NULL index in array
		 bdx = 0;
                 while (0 != binFlags[(index*max_per_bin)+bdx]) {
	            bdx++;
                    if (bdx > n) {
                       printf("ERROR: Overflow \n");
                       return -1;
		    }	 
		 }
                 binArray[(index*max_per_bin)+bdx] = binArray[(b*max_per_bin)+idx]; // Add element into first non-NULL index   	
                 binFlags[(index*max_per_bin)+bdx] = 1;
		 binFlags[(b*max_per_bin)+idx] = 0; // Add element into first non-NULL index   	

		 binParticleNum[b]--;
		 i--;
		 binParticleNum[index]++;

	      }
	      idx++;
	   }
	}


        //
        // 5. Compact Particles
        //
        int nextLoc;
        for (int b=0; b<binNum; b++) { // The bth bin
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



	// Check count
	count = 0;
	for (int b=0; b<binNum; ++b) {
	   count += binParticleNum[b];
	}

//        if( fsave && (step%SAVEFREQ) == 0 ) {save( fsave, n, particles );}

           if( fsave && (step%SAVEFREQ) == 0 ) {

           int cdx=0;
           for (int rdx=0; rdx<binNum; rdx++) {
              for (int sdx=0; sdx<binParticleNum[rdx]; sdx++) {
                 compactVect[cdx] = &binArray[(rdx*max_per_bin)+sdx]; // displs[rdx] is displac of one proc
                 //printf("PRE GlobalID is %d \n", compactVect[cdx]->globalID);
                 cdx++;
              }
           }
           //printf("CDX is %d \n", cdx);

           for (int rdx=0; rdx<cdx; rdx++){
              for (int sdx=0; sdx<cdx; sdx++) {
                 if (compactVect[rdx]->globalID == compactVect[sdx]->globalID && rdx!=sdx) {
                    printf("ERROR GLOBAL ID CONFLICT, %d \n", compactVect[rdx]->globalID);
                 }
              }
           }

           qsort(compactVect, n, sizeof(particle_t*), compare);

           int currID = 0;
           for (int rdx=0; rdx<cdx; rdx++){
              particles_mpi[rdx] = *compactVect[rdx];
              //if (particles_mpi[rdx].globalID != currID  ) { 
                 //printf("POST GlobalID is %d \n", particles_mpi[rdx].globalID);
              //}
              currID++; 
           }

           save( fsave, n, particles_mpi );

           }


    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
   
    delete [] binParticleNum;
    delete [] binArray;
    delete [] binFlags;

    free( compactVect );
    free( particles );
    free ( particles_mpi );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
