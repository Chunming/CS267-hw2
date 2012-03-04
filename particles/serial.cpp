#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <string.h>

#define cutoff 0.01 // Was defined in common.cpp
#define density 0.0005
#define subBlockLen 0.01

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
    int subBlockNum = boxSize/subBlockLen; // No. of sub blocks along a row/column (Default: 5)
    int binNum = subBlockNum*subBlockNum; // No. of bins

    // Determine number of sublocks, represent as a vector;
    // Set no. of blocks as 5 x 5
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


   for (int i = 0; i<binNum; ++i) {
      binParticleNum[i] = 0; // sizeof(binParticleNum[b] is 4 bytes
   }	   

   memset(binParticleNum, 0, sizeof(int)*binNum); // sizeof(binParticleNum[b] is 4 bytes


   for (int b=0; b<binNum; ++b) { 
	    binArray[b] = new particle_t* [100]; // Set it to max 100 particles first
            if (NULL == binArray[b]) {
               printf("ERROR binArray index  mem alloc failed \n");
               return -1;
            }
          
            //for (int i=0; i<100; i++) {
	    //   binArray[b][i] = NULL;
	    //}
	    memset(binArray[b], NULL, 100*sizeof(particle_t*));
    }

    // Initialize particles 
    init_particles( n, particles );

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
		 if (leftBnd!=0) { // Left subBlock index is valid
		    bLeft = b - subBlockNum; 
		    //printf("bLeft is %d \n", bLeft);
		    kdx=0;
		    for (int k=0; k<binParticleNum[bLeft]; ++k) { 
		       while (NULL == (binArray[bLeft][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bLeft][kdx]);
		       kdx++;
		    }
		 }
	      }

	      if (rightDist<=cutoff) {
		 if (rightBnd!=boxSize) { 
  		    bRight = b + subBlockNum; 
		    //printf("bRight is %d \n", bRight);
		    kdx=0;
		    for (int k=0; k<binParticleNum[bRight]; ++k) { 
		       while (NULL == (binArray[bRight][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bRight][kdx]);
		       kdx++;
		    }
		 }
	      }

	      if (topDist<=cutoff) {
		 if (topBnd!=0) { 
		    bTop = b - 1; 
		    //printf("bTop is %d \n", bTop);
  		    kdx=0;
		    for (int k=0; k<binParticleNum[bTop]; ++k) { 
		       while (NULL == (binArray[bTop][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bTop][kdx]);
		       kdx++;
		    }
		 }
	      }

	      if (botDist<=cutoff) {
		 if (botBnd!=boxSize) {
       		    bBottom = b + 1;
		    //printf("bBottom is %d \n", bBottom);
     		    kdx=0;
		    for (int k=0; k<binParticleNum[bBottom]; ++k) { 
		       while (NULL == (binArray[bBottom][kdx])) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bBottom][kdx]);
		       kdx++;
		    }
		 }
	      }

	      if (topDist<=cutoff && leftDist<=cutoff) { 
		 if (topBnd!=0 && leftBnd !=0) {
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
		 if (botBnd!=boxSize && leftBnd!=0) { 
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
		 if (topBnd!=0 && rightBnd!=boxSize) {
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
		 if (botBnd!=boxSize && rightBnd!=boxSize) {
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
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
   
    delete [] binParticleNum;
    for (int b=0; b<binNum; ++b) {
       delete [] binArray[b];
    }
    delete binArray;


    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
