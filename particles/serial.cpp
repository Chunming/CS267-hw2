#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

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
    printf("The value of n is  %d \n", n);
    printf("The box size is %f \n", boxSize); // Box size is 0.5
     
	      
    // Determine number of sublocks, represent as a vector;
    // Set no. of blocks as 5 x 5
    int binNum = 25; // No. of bins
    int* binParticleNum = new int [binNum]; // No. of particles in each bin
    particle_t*** binArray = new particle_t** [binNum]; // Array of ptrs to particle*
   
    for (int idx=0; idx<binNum; idx++) { 
	    binParticleNum[idx] = 0;
	    binArray[idx] = new particle_t* [100]; // Set it to max 100 particles first
    }
    


    // Initialize particles 
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j] );
        }
        
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {

	    printf("The time step is %d \n", step);
            save( fsave, n, particles );

         }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
   
    delete [] binParticleNum;
    for (int idx=0; idx<binNum; idx++) {
       delete [] binArray[idx];
    }
    delete binArray;


    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
