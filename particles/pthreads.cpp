#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>
#include "common.h"

extern double cutoff; //make sure this shit works.
extern double size;

//
//  global variables
//
int n, n_threads, blks_num, subblks_num, totblks_num, part_int;
double blks_size, subblks_size, simulation_time;
particle_t *particles;
particle_t **bins;
int *bin_part_num;
FILE *fsave;
pthread_barrier_t barrier;
pthread_mutex_t *bin_mutexes;
pthread_mutex_t set_count_zero;
pthread_cond_t count_zero_cond;

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

#define MAX_PART_THREADS 5000
#define MAX_PART_SUBB 5


//
// This is where we put particles into bins
//
void bin_particles ( int block_id, int step ) {
	
  int loop_limit = min(n, (block_id + 1)*part_int);

  for (int i = block_id * part_int; i < loop_limit; ++i) {
    int x_bin = particles[i].x/subblks_size;
    int y_bin = particles[i].y/subblks_size;
    int blockId = totblks_num * x_bin + y_bin;

    pthread_mutex_lock(&bin_mutexes[blockId]);
    int part_loc = bin_part_num[blockId]++;
    pthread_mutex_unlock(&bin_mutexes[blockId]);

    particles[i].ax = particles[i].ay = 0;
    bins[blockId * MAX_PART_SUBB + part_loc] = particles + i;
  }
   
}

//
//  This is the subblock force computing routine
//                  
void compute_sub_forces( int blockId, int bIdx, int bIdy )
{
  // Compute Forces 
  double leftBnd, rightBnd, topBnd, botBnd;
  double leftDist, rightDist, topDist, botDist;
  int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;

  for (int i = 0; i<bin_part_num[blockId]; ++i) {
    for (int j = 0; j<bin_part_num[blockId]; ++j) { // The jth particle in bth                          
      apply_force(*bins[blockId*MAX_PART_SUBB + i], *bins[blockId*MAX_PART_SUBB + j]);
    }

    leftBnd = bIdx*subblks_size;
    rightBnd = (bIdx*subblks_size) + subblks_size;
    topBnd = bIdy*subblks_size;
    botBnd = bIdy*subblks_size + subblks_size;

    leftDist = bins[blockId*MAX_PART_SUBB + i]->x - leftBnd;
    rightDist = rightBnd - bins[blockId*MAX_PART_SUBB + i]->x;
    topDist = bins[blockId*MAX_PART_SUBB + i]->y - topBnd;
    botDist = botBnd - bins[blockId*MAX_PART_SUBB + i]->y;

    // Consider 8 different adjacent subBlocks
    if (leftDist<=cutoff && bIdx != 0) {
      bLeft = blockId - totblks_num;
      for (int k=0; k<bin_part_num[bLeft]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bLeft*MAX_PART_SUBB + k]);
      }
    }
    //2
    if (rightDist<=cutoff && bIdx != totblks_num-1) {
      bRight = blockId + totblks_num;
      for (int k=0; k<bin_part_num[bRight]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bRight*MAX_PART_SUBB + k]);
      }
    }
    //3
    if (topDist<=cutoff && bIdy != 0) {
      bTop = blockId - 1;
      for (int k=0; k<bin_part_num[bTop]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTop*MAX_PART_SUBB + k]);
      }
    }
    //4
    if (botDist<=cutoff && bIdy != totblks_num-1) {
      bBottom = blockId + 1;
      for (int k=0; k<bin_part_num[bBottom]; ++k) {
	apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBottom*MAX_PART_SUBB + k]);
      }
    }
    //5
    if (topDist<=cutoff && leftDist<=cutoff && bIdy != 0 && bIdx !=0) {
      bTopLeft = blockId - totblks_num - 1;
      for (int k=0; k<bin_part_num[bTopLeft]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTopLeft*MAX_PART_SUBB + k]);
      }
    }
    //6
    if (botDist<=cutoff && leftDist<=cutoff && bIdy != totblks_num-1 && bIdx != 0) {
      bBotLeft = blockId-totblks_num+1;
      for (int k=0; k<bin_part_num[bBotLeft]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBotLeft*MAX_PART_SUBB + k]);
      }
    }
    //7
    if (topDist<=cutoff && rightDist<=cutoff && bIdy != 0 && bIdx != totblks_num-1) {
      bTopRight = blockId+totblks_num-1;
      for (int k=0; k<bin_part_num[bTopRight]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTopRight*MAX_PART_SUBB + k]);
      }
    }
    //8
    if (botDist<=cutoff && rightDist<=cutoff && bIdy!=totblks_num-1 && bIdx!=totblks_num-1) {
      bBotRight = blockId+totblks_num+1;
      for (int k=0; k<bin_part_num[bBotRight]; ++k) {
        apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBotRight*MAX_PART_SUBB + k]);
      }
    }

  }

}


//
//  This is where we compute forces for a given block
//
void compute_sup_forces( int sup_blockId, int step )
{

  int xIdx = (sup_blockId / blks_num) * subblks_num;
  int xIdy = (sup_blockId % blks_num) * subblks_num;
  int bIdx, bIdy, blockId;

  double leftBnd, rightBnd, topBnd, botBnd;
  double leftDist, rightDist, topDist, botDist;
  int bLeft, bRight, bBottom, bTop, bTopLeft, bTopRight, bBotLeft, bBotRight;

  
  //handle the middle
  for (int xsub = 0; xsub < subblks_num; ++xsub) {
    bIdx = xIdx + xsub;
    for (int ysub = 0; ysub < subblks_num; ++ysub) {
      bIdy = xIdy + ysub;
      blockId = bIdx * totblks_num + bIdy;
      //compute_sub_forces(blockId, bIdx, bIdy);
      for (int i = 0; i<bin_part_num[blockId]; ++i) {
	for (int j = 0; j<bin_part_num[blockId]; ++j) { 
	  apply_force(*bins[blockId*MAX_PART_SUBB + i], *bins[blockId*MAX_PART_SUBB + j]);
	}

	leftBnd = bIdx*subblks_size;
	rightBnd = (bIdx*subblks_size) + subblks_size;
	topBnd = bIdy*subblks_size;
	botBnd = bIdy*subblks_size + subblks_size;

	leftDist = bins[blockId*MAX_PART_SUBB + i]->x - leftBnd;
	rightDist = rightBnd - bins[blockId*MAX_PART_SUBB + i]->x;
	topDist = bins[blockId*MAX_PART_SUBB + i]->y - topBnd;
	botDist = botBnd - bins[blockId*MAX_PART_SUBB + i]->y;

	if (leftDist < 0) printf("leftDist <0 \n");
	if (rightDist < 0) printf("rightDist <0 \n");
	if (botDist < 0) printf("botDist <0 \n");
	if (topDist < 0) printf("topDist <0 \n");

	// Consider 8 different adjacent subBlocks   
	if (leftDist<=cutoff && bIdx != 0) {
	  bLeft = blockId - totblks_num;
	  for (int k=0; k<bin_part_num[bLeft]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bLeft*MAX_PART_SUBB + k]);
	  }
	}
	//2                                                                     
	if (rightDist<=cutoff && bIdx != totblks_num-1) {
	  bRight = blockId + totblks_num;
	  for (int k=0; k<bin_part_num[bRight]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bRight*MAX_PART_SUBB + k]);
	  }
	}
	//3    
	if (topDist<=cutoff && bIdy != 0) {
	  bTop = blockId - 1;
	  for (int k=0; k<bin_part_num[bTop]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTop*MAX_PART_SUBB + k]);
	  }
	}
	//4    
	if (botDist<=cutoff && bIdy != totblks_num-1) {
	  bBottom = blockId + 1;
	  for (int k=0; k<bin_part_num[bBottom]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBottom*MAX_PART_SUBB + k]);
	  }
	}
	//5
	if (topDist<=cutoff && leftDist<=cutoff && bIdy != 0 && bIdx !=0) {
	  bTopLeft = blockId - totblks_num - 1;
	  for (int k=0; k<bin_part_num[bTopLeft]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTopLeft*MAX_PART_SUBB + k]);
	  }
	}
	//6
	if (botDist<=cutoff && leftDist<=cutoff && bIdy != totblks_num-1 && bIdx != 0) {
	  bBotLeft = blockId-totblks_num+1;
	  for (int k=0; k<bin_part_num[bBotLeft]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBotLeft*MAX_PART_SUBB + k]);
	  }
	}
	//7        
	if (topDist<=cutoff && rightDist<=cutoff && bIdy != 0 && bIdx != totblks_num-1) {
	  bTopRight = blockId+totblks_num-1;
	  for (int k=0; k<bin_part_num[bTopRight]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bTopRight*MAX_PART_SUBB + k]);
	  }
	}
	//8            
	if (botDist<=cutoff && rightDist<=cutoff && bIdy!=totblks_num-1 && bIdx!=totblks_num-1) {
	  bBotRight = blockId+totblks_num+1;
	  for (int k=0; k<bin_part_num[bBotRight]; ++k) {
	    apply_force(*bins[blockId*MAX_PART_SUBB + i],*bins[bBotRight*MAX_PART_SUBB + k]);
	  }
	}

      }

    }
  }
    
}


//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{

 

  int thread_id = *(int*)pthread_id;
  int sup_x = (thread_id / blks_num) * subblks_num;
  int sup_y = (thread_id % blks_num) * subblks_num;
  int block_id;


  //
  //  simulate a number of time steps
  //

  for( int step = 0; step < NSTEPS; step++ )
  {

      //
      //  bin particles
      //
      bin_particles( thread_id, step ); 

      pthread_barrier_wait( &barrier );
      //
      //  compute forces in your bin
      //
      compute_sup_forces( thread_id, step );
      //  printf("extied compute_sup_forces %d \n", step);      
      //play with this barrier!!!
      pthread_barrier_wait( &barrier );

      //
      //  move particles in your bin
      //
     
      for (int sub_x = 0; sub_x < subblks_num; ++sub_x) {
	for (int sub_y = 0; sub_y < subblks_num; ++sub_y) {
	  block_id = (sup_x + sub_x) * totblks_num + (sup_y + sub_y);
	  for (int i = 0; i < bin_part_num[block_id]; ++i) {
	    move( *bins[block_id * MAX_PART_SUBB + i] );
	  }
	  bin_part_num[block_id] = 0;
	}
      } 
      
	pthread_barrier_wait( &barrier );

      //
      //  save if necessary
      //
      if(fsave && (step%SAVEFREQ) == 0) {

	if( thread_id == 0 ) {
	  save( fsave, n, particles );
	}
      }
     
    }
    
  return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
  //
  //  process command line
  //
  if( find_option( argc, argv, "-h" ) >= 0 )
    {
      printf( "Options:\n" );
      printf( "-h to see this help\n" );
      printf( "-n <int> to set the number of particles\n" );
      printf( "-p <int> to set the number of threads\n" );
      printf( "-o <filename> to specify the output file name\n" );
      return 0;
    }
    
  n = read_int( argc, argv, "-n", 1000 );
  n_threads = read_int( argc, argv, "-p", 2 );
  char *savename = read_string( argc, argv, "-o", "pthreads.txt" );
  
  //
  //  allocate resources
  //
  fsave = savename ? fopen( savename, "w" ) : NULL;

  particles = (particle_t*) malloc( n * sizeof(particle_t) );
  
  set_size( n );
  init_particles( n, particles );

  blks_num = max(2 * size / sqrt((MAX_PART_THREADS)*(3.14159)*cutoff*cutoff), 1);
  blks_size =  size / ((double) blks_num);
  subblks_num = max(2 * blks_size / sqrt((MAX_PART_SUBB)*(3.14159)*cutoff*cutoff), 1);
  subblks_size = blks_size / ((double) subblks_num);
  totblks_num = blks_num * subblks_num;
  part_int = n/( blks_num * blks_num ) + 1;

  printf("blks_num %d, blks_size %f, subblks_num %d, subblks_size %f, totblks_num %d \n", blks_num, blks_size, subblks_num, subblks_size, totblks_num);

  printf("size %f, totblks_num * subblks_size %f \n", size, totblks_num * subblks_size);

  //printf("n %d, part_int %d, blks_num %d, cutoff %f \n", n, part_int, blks_num, cutoff);
  
  bins = (particle_t **) malloc( MAX_PART_SUBB * totblks_num * totblks_num * sizeof(particle_t*));
  bin_part_num = (int *) malloc( totblks_num * totblks_num * sizeof(int));
  bin_mutexes = (pthread_mutex_t*) malloc( totblks_num * totblks_num * sizeof(pthread_mutex_t));
  memset(bin_part_num, 0, totblks_num * totblks_num * sizeof(int));
  

  pthread_attr_t attr;
  P( pthread_attr_init( &attr ) );
  P( pthread_barrier_init( &barrier, NULL, blks_num * blks_num ) );
  pthread_mutex_init( &set_count_zero, NULL );
  pthread_cond_init( &count_zero_cond, NULL );

  for(int i=0; i<totblks_num*totblks_num; ++i)
  {
      if(pthread_mutex_init(&bin_mutexes[i], NULL))
	{
	  printf("\nsetupBins(): Unable to initialize the %i th mutex\n", i);
	}
  }

  int *thread_ids = (int *) malloc( blks_num * blks_num * sizeof( int ) );
  for( int i = 0; i < blks_num*blks_num; i++ ) 
    thread_ids[i] = i;

  pthread_t *threads = (pthread_t *) malloc( blks_num * blks_num * sizeof( pthread_t ) );
    
  printf("shit allocated \n");

  //
  //  do the parallel work
  //
  simulation_time = read_timer( );
  for( int i = 1; i < blks_num * blks_num; i++ ) 
    P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
  printf("pthreads created \n");

  thread_routine( &thread_ids[0] );
    
  printf("first thread runs \n");

  for( int i = 1; i < blks_num * blks_num; i++ ) 
    P( pthread_join( threads[i], NULL ) );

  simulation_time = read_timer( ) - simulation_time;
    
  printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, blks_num*blks_num, simulation_time );
    
  //
  //  release resources
  //
  P( pthread_barrier_destroy( &barrier ) );
  P( pthread_attr_destroy( &attr ) );
  pthread_mutex_destroy(&set_count_zero);
  pthread_cond_destroy(&count_zero_cond);

  for(int i=0; i<totblks_num*totblks_num; i++)
  {
      pthread_mutex_destroy(&bin_mutexes[i]);
  }

  free( thread_ids );
  free( threads );
  free( bin_mutexes );
  free( bin_part_num );
  free( bins );
  free( particles );
  if( fsave )
    fclose( fsave );
    
  return 0;
}

