#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

#define cutoff 0.01 // Was defined in common.cpp

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

    double boxSize = sqrt(0.0005*n); // Density is set at 0.0005
    //printf("The value of n is  %d \n", n);
    //printf("The box size is %f \n", boxSize); // Box size is 0.5
     
	      
    // Determine number of sublocks, represent as a vector;
    // Set no. of blocks as 5 x 5
    int binNum = 25; // No. of bins
    int* binParticleNum = new int [binNum]; // No. of particles in each bin
    if (binParticleNum == NULL) {
       printf("ERROR binParticleNum mem alloc failed \n");
       return -1;
    }

    particle_t*** binArray = new particle_t** [binNum]; // Array of ptrs to particle*
    if (binArray == NULL) {
       printf("ERROR binArray mem alloc failed \n");
       return -1;
    }

    for (int b=0; b<binNum; b++) { 
	    binParticleNum[b] = 0;
	    binArray[b] = new particle_t* [100]; // Set it to max 100 particles first
            if (binArray[b] == NULL) {
               printf("ERROR binArray index  mem alloc failed \n");
               return -1;
            }
    
	    for (int i=0; i<100; i++) {
	       binArray[b][i] = NULL; // Set all ptrs to NULL
	    }
    }
    printf("Bins are allocated, error check done \n");

    // Initialize particles 
    init_particles( n, particles );

    // Initialize particle binning
    int xIdx, yIdx; // x/y index in 1D array
    int bdx; 
    double subBlockLen = 0.1;
    int subBlockNum = 5; // No. of sub blocks along a row/column

    for (int ndx=0; ndx<n; ndx++) { // For each particle
       xIdx = (particles[ndx].x)/subBlockLen; 
//       printf("xIdx is %d \n", xIdx);
       yIdx = (particles[ndx].y)/subBlockLen;
//       printf("yIdx is %d \n", yIdx);
       bdx = 0;
       // Store into first non-NULL index in array
       while (binArray[yIdx+(xIdx*subBlockNum)][bdx]!=NULL) {
//          printf("bdx is %d \n", bdx)
          bdx++;
       }
//       printf("Bin index is %d \n", yIdx+xIdx*subBlockNum);
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

	double leftBnd, rightBnd, topBnd, botBnd;
	double leftDist, rightDist, topDist, botDist;
	int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;

	for (int b=0; b<binNum; b++) { // The bth bin

	   idx=0;
	   for (int i=0; i<binParticleNum[b]; i++) { // The ith particle in bth bin
     	      while ((binArray[b][idx])==NULL) { idx++; } // Index skips NULL values
      	      (*binArray[b][idx]).ax = (*binArray[b][idx]).ay = 0;


	      // Check all particles in bth subBlock
	      jdx=0;
     	      for (int j=0; j<binParticleNum[b]; j++) { // The jth particle in bth bin
     	         while ((binArray[b][jdx])==NULL) { jdx++; }
		 apply_force(*binArray[b][idx],*binArray[b][jdx]);
	      }

	      // Check particles to the left/right/top/bottom subBlocks of bth subBlock	 
	      xIdx = (*binArray[b][idx]).x/subBlockLen;
	      //printf("xIdx is %d \n", xIdx);
	      yIdx = (*binArray[b][idx]).y/subBlockLen;
	      //printf("yIdx is %d \n", yIdx);

	      leftBnd = xIdx*subBlockLen;
	      //printf("leftBnd is %f \n", leftBnd);

	      rightBnd = (xIdx*subBlockLen) + subBlockLen;
	      //printf("rightBnd is %f \n", rightBnd);

	      topBnd = yIdx*subBlockLen;
	      //printf("topBnd is %f \n", topBnd);

	      botBnd = yIdx*subBlockLen + subBlockLen;
	      //printf("botBnd is %f \n", botBnd);

	      leftDist = fabs((*binArray[b][idx]).x - leftBnd);
	      //printf("leftDist is %f \n", leftDist);

	      rightDist = fabs((*binArray[b][idx]).x - rightBnd);
	      //printf("rightDist is %f \n", rightDist);

	      topDist = fabs((*binArray[b][idx]).y - topBnd);
	      //printf("topDist is %f \n", topDist);

	      botDist = fabs((*binArray[b][idx]).y - botBnd);
	      //printf("botDist is %f \n", botDist);


	      // Consider 8 different adjacent subBlocks
	      if (leftDist<cutoff) {
		 if (leftBnd!=0) { // Left subBlock index is valid
		    bLeft = b - subBlockNum; // subBlockLen=0.1, subBlockNum = 5 
		    //printf("bLeft is %d \n", bLeft);
		    kdx=0;
		    for (int k=0; k<binParticleNum[bLeft]; k++) { 
		       while ((binArray[bLeft][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bLeft][kdx]);
		    }
		 }
	      }

	      if (rightDist<cutoff) {
		 if (rightBnd!=boxSize) { 
  		    bRight = b + subBlockNum; 
		    //printf("bRight is %d \n", bRight);
		    kdx=0;
		    for (int k=0; k<binParticleNum[bRight]; k++) { 
		       while ((binArray[bRight][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bRight][kdx]);
		    }
		 }
	      }

	      if (topDist<cutoff) {
		 if (topBnd!=0) { 
		    bTop = b - 1; 
		    //printf("bTop is %d \n", bTop);
  		    kdx=0;
		    for (int k=0; k<binParticleNum[bTop]; k++) { 
		       while ((binArray[bTop][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bTop][kdx]);
		    }
		 }
	      }

	      if (botDist<cutoff) {
		 if (botBnd!=boxSize) {
       		    bBottom = b + 1;
		    //printf("bBottom is %d \n", bBottom);
     		    kdx=0;
		    for (int k=0; k<binParticleNum[bBottom]; k++) { 
		       while ((binArray[bBottom][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bBottom][kdx]);
		    }
		 }
	      }

	      if (topDist<cutoff && leftDist<cutoff) { 
		 if (topBnd!=0 && leftBnd !=0) {
		    bTopLeft = b-subBlockNum-1;     
		    kdx=0;
		    for (int k=0; k<binParticleNum[bTopLeft]; k++) { 
		       while ((binArray[bTopLeft][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bTopLeft][kdx]);
		    }
		 }
	      }

	      if (botDist<cutoff && leftDist<cutoff) { 
		 if (botBnd!=boxSize && leftBnd!=0) { 
		    bBotLeft = b-subBlockNum+1;     
     		    kdx=0;
		    for (int k=0; k<binParticleNum[bBotLeft]; k++) { 
		       while ((binArray[bBotLeft][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bBotLeft][kdx]);
		    }
		 }
	      }

	      if (topDist<cutoff && rightDist<cutoff) { 
		 if (topBnd!=0 && rightBnd!=boxSize) {
		    bTopRight = b+subBlockNum-1;     
     		    kdx=0;
		    for (int k=0; k<binParticleNum[bTopRight]; k++) { 
		       while ((binArray[bTopRight][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bTopRight][kdx]);
		    }
		 }
	      }

	      if (botDist<cutoff && rightDist<cutoff) { 
		 if (botBnd!=boxSize && rightBnd!=boxSize) {
   		    bBotRight = b+subBlockNum+1;     
     		    kdx=0;
		    for (int k=0; k<binParticleNum[bBotRight]; k++) { 
		       while ((binArray[bBotRight][kdx])==NULL) { kdx++; }
		       apply_force(*binArray[b][idx],*binArray[bBotRight][kdx]);
		    }
		 }
	      }

	      

            }
	}

        
	
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

     	      while ((binArray[b][idx])==NULL) { idx++; }
	      xIdx = (*binArray[b][idx]).x/subBlockLen;
              yIdx = (*binArray[b][idx]).y/subBlockLen;
	      index = yIdx+(xIdx*subBlockNum); // Map 2D to 1D index
	      if (index != b) { // Particle has moved out of of the bin

                 // Store into first non-NULL index in array
		 bdx = 0;
                 while (binArray[index][bdx]!=NULL) {
	            bdx++;
                    if (bdx > 500) {
                       printf("ERROR: Overflow \n");
                       return -1;
		    }	 
		 }
                 binArray[index][bdx] = binArray[b][idx]; // Add element into first non-NULL index   	
		 binArray[b][idx] = NULL;
		 binParticleNum[b]--;
		 binParticleNum[index]++;

	      }
	   }
	}

	// Check count
	count = 0;
	for (int b=0; b<binNum; b++) {
	   count += binParticleNum[b];
	}
	//printf("The count is %d \n", count); // Check count, should be 500


        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {

	    //printf("The time step is %d \n", step);
            save( fsave, n, particles );

         }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
   
    delete [] binParticleNum;
    for (int b=0; b<binNum; b++) {
       delete [] binArray[b];
    }
    delete binArray;


    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
