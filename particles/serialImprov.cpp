#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>
#include <vector>

// Improve version of serial implementation
// Size of walls in x/y direction is size
// Declare new class called Node
class Node {
   public:
   particle_t* pAddress; // Value of curr particle's address
   Node* nextNode;
};

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
   
    vector<vector<int> > subBlock[25];
    vector <int> A;

    A.push_back(0);
    subBlock.push_back(A); 

  //vector<particle_t*> subBlock; // Each subBlock is represented by a vector


    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {


        // Do binning first for all particles
        // Initialize linked list of pointers to &particles[i] 
        // (i.e. Address corressponding to index in particles array)
        // Use dynamic array since bins will have particles << n
        // Node* head = new Node;
        // head->pAddress = &particles[0];
        // head->nextNode = NULL;
        // delete head;


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

	    //printf("The time step is %d \n", step);
            save( fsave, n, particles );

         }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
