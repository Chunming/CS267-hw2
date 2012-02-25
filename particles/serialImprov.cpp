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
class ListNode {
   public:
   particle_t* pAddress; // Value of curr particle's address
   ListNode* nextNode;
};

bool insert(ListNode** root, ListNode* targetNode) {

   if (targetNode==NULL) return 0; // Error
	
   else if ((*root)==NULL) *root = targetNode; // Initialize array

   else {
      ListNode* tmp;
      tmp = *root      
      while (tmp->nextNode != NULL) {
         tmp = tmp->next;
      }
      tmp->next = targetNode;
   }

   return 1; 


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

    double boxSize = sqrt(0.0005*n); // Density is set at 0.0005
    printf("The value of n is  %d \n", n);
    printf("The box size is %f \n", boxSize); // Box size is 0.5


    // Determine number of sublocks, represent as a vector;
    // Set no. of blocks as 5 x 5
    int arrSize = 25;
    ListNode** subBlocks = new ListNode* [arrSize]; // Array of ptrs ListNode*
    for (int idx=0; idx<arrSize; idx++) { // Ini all entries to NULL
       subBlocks[idx] = NULL;
    }



    init_particles( n, particles );
   

    //
    //  simulate a number of time steps
    //

    Node* head;
    int xIdx;
    int yIdx;
    double subBlockLen = 0.1;
    int subBlockNum = 5; // No. of sub blocks along a row/column
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {


        // Do binning first for all particles
        // Initialize linked list of pointers to &particles[i] 
        // (i.e. Address corressponding to index in particles array)
        // Use dynamic array since bins will have particles << n

	for (int ndx=0; ndx < n; n++) { // For each particle
           head = new Node; // Remember to delete later
           head->pAddress = &particles[ndx];
           head->nextNode = NULL;
           
	   xIdx = (particles[ndx].x)/subBlockLen;
	   yIdx = (particles[ndx].y)/subBlockLen;

	   insert(&subBlocks[yIdx+(xIdx*subBlockNum)], head); // Insert to end of Linked List


	}	




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
