
#include "simulation.h"
#include "readconfig.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <cinttypes>
#include <ctime>

#ifdef _MPI_
#include <mpi.h>
#include "helper.h"
#endif

using namespace LBM;

int main(int argc, char** argv) {

	int verbose 	= 0;
	std::unique_ptr<char[]> problem_config( new char[1024] );
	std::unique_ptr<char[]> problem_geometry( new char[1024] );

	if( argc == 2 ) {
		sprintf( problem_config.get(), "dat/%s.dat", argv[1] );
	} else if( argc == 3 ) {
		if( strcmp( argv[1], "-v" ) == 0 ) {
				printf( "Running Verbose mode\n" );
				verbose 	= 1;
		}
		sprintf( problem_config.get(), "dat/%s.dat", argv[2] );
	} else {
#ifdef _SEQ_
		sprintf( problem_config.get(), "%s", "dat/movingwall_20x20x20.dat" );
#elif _MPI_
		sprintf( problem_config.get(), "%s", "dat/movingwall_50x50x50_Parallel.dat" );
#endif
	}

#ifdef _MPI_
	int mpi_size, rank;

	// MPI Init:
	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status status;
	MPI_Comm cart_comm, widthComm, heightComm, depthComm;
	int dim[DIM];
	int periods[DIM] = { 0 };

#endif
	/* READ PARAMETERS */
	int size[3] = { 0 };
	int timesteps;
	double tau;
	int t_vis;

#ifdef _MPI_
	std::unique_ptr<int[]> GLOBAL_FLAG;
	std::unique_ptr<int[]> FLAG;

	std::unique_ptr<int[]>flag_disp( new int[ mpi_size ] );
	std::unique_ptr<int[]>flag_count( new int[ mpi_size ] );

	int x_proc, y_proc, z_proc;
	std::unique_ptr<int[]>scatter_size( new int[ DIM * mpi_size ] );
	if (rank==0) {
		
		int remainder[DIM];

		int GLOBAL_size[DIM] = { 0 };
		int GLOBAL_x;
		int GLOBAL_y;
		int GLOBAL_z;

		readConfig( problem_config.get(), problem_geometry.get()
					, &GLOBAL_x, &GLOBAL_y, &GLOBAL_z
					, &timesteps, &tau, &t_vis
					, &x_proc, &y_proc, &z_proc
				);

		int domCheck = x_proc * y_proc * z_proc;

		if( domCheck != mpi_size ) {
			printf( "INVALID DOMAIN SIZE: number of processsors does not match domain decomposition\n%i != %i"
						, domCheck, mpi_size  );
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		dim[0] = x_proc;
		dim[1] = y_proc;
		dim[2] = z_proc;

		GLOBAL_size[0] = GLOBAL_x;
		GLOBAL_size[1] = GLOBAL_y;
		GLOBAL_size[2] = GLOBAL_z;

		// DOMAIN DECOMPOSITION:
		int uneven_decomp_x = 0;
		int uneven_decomp_y = 0;
		int uneven_decomp_z = 0;
		if( ( GLOBAL_size[0] % x_proc != 0 ) ) {
			remainder[0] = GLOBAL_size[0] % x_proc;
			uneven_decomp_x = 1;
		} 
		if ( GLOBAL_size[1] % y_proc != 0 ) { 
			remainder[1] = GLOBAL_size[1] % y_proc;
			uneven_decomp_z = 1;
		}
		if ( GLOBAL_size[2] % z_proc != 0 ) {
			remainder[2] = GLOBAL_size[2] % z_proc;
			uneven_decomp_z = 1;
		}
		int idx;
		int GLOBAL_totLen = 0;
		for( int i = 0; i < x_proc; i++ ) {
			for( int j = 0; j < y_proc; j++ ) {
				for( int k = 0; k < z_proc; k++ ) {

					idx = k *(x_proc*y_proc) + j * x_proc + i;
					scatter_size[ DIM * idx + 0] = GLOBAL_size[0] / x_proc;
					scatter_size[ DIM * idx + 1] = GLOBAL_size[1] / y_proc;
					scatter_size[ DIM * idx + 2] = GLOBAL_size[2] / z_proc;

					if( uneven_decomp_x && i == x_proc - 1 ) {
						scatter_size[ DIM * idx + 0] += remainder[0];
					}
					if( uneven_decomp_y && j == y_proc - 1 ) {
						scatter_size[ DIM * idx + 1] += remainder[1];
					}
					if( uneven_decomp_z && j == z_proc - 1 ) {
						scatter_size[ DIM * idx + 2] += remainder[2];
					}

					flag_count[ idx ] = ( scatter_size[ DIM * idx + 0] + 2 ) 
									* ( scatter_size[ DIM * idx + 1] + 2 ) 
									* ( scatter_size[ DIM * idx + 2] + 2 );

					flag_disp[ idx ] = GLOBAL_totLen;
					GLOBAL_totLen += flag_count[ idx ];
				}
			}
		}
		// Global Flag:
		GLOBAL_FLAG = std::unique_ptr<int[]>( new int[ GLOBAL_totLen ] );

		buildFlag( problem_geometry.get(), GLOBAL_FLAG.get(), GLOBAL_size
					, x_proc, y_proc, z_proc, scatter_size.get() );	
	}
	// Scatterv/Bcast information to all other processors:
	MPI_Bcast( &x_proc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &y_proc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &z_proc, 1, MPI_INT, 0, MPI_COMM_WORLD );

	MPI_Bcast( &tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( dim, 3, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &timesteps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &t_vis, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	std::unique_ptr<int[]>disp( new int[ mpi_size ] );
	std::unique_ptr<int[]>count( new int[ mpi_size ] );
	
	for( int i = 0; i < mpi_size; i++ ) {
		disp[i] = i * DIM;
		count[i] = DIM;
	}

	MPI_Scatterv( scatter_size.get(), count.get(), disp.get(), MPI_INT, size, 3, MPI_INT, 0, MPI_COMM_WORLD );

	MPI_Cart_create(MPI_COMM_WORLD, 3, dim, periods, 1, &cart_comm);
	int freecoords[] = { 1, 0, 0 };
	// width comms:
	int widthRank, depthRank, heightRank;
	MPI_Cart_sub(cart_comm, freecoords, &widthComm);
	MPI_Comm_rank(widthComm, &widthRank);
	// depth comms:
	freecoords[0] = 0;
	freecoords[1] = 1;
	MPI_Cart_sub( cart_comm, freecoords, &depthComm );
	MPI_Comm_rank( depthComm, &depthRank );
	// height comms:
	freecoords[1] = 0;
	freecoords[2] = 1;
	MPI_Cart_sub( cart_comm, freecoords, &heightComm );
	MPI_Comm_rank( heightComm, &heightRank );

	int totLen = ( size[0] + 2 ) * ( size[1] + 2 ) * ( size[2] + 2 );
	FLAG = std::unique_ptr<int[]>( new int[totLen] );

	MPI_Scatterv( GLOBAL_FLAG.get(), flag_count.get(), flag_disp.get(), MPI_INT, FLAG.get(), totLen, MPI_INT, 0, MPI_COMM_WORLD );

#elif _SEQ_
	int x, y, z;
	readConfig( problem_config.get(), problem_geometry.get()
					, &x, &y, &z, &timesteps, &tau, &t_vis );
	size[0] = x;
	size[1] = y;
	size[2] = z;

	int totLen = (size[0] + 2)*(size[1] + 2)*(size[2] + 2);
	std::unique_ptr<int[]> FLAG( new int[totLen] );

	buildFlag( problem_geometry.get(), FLAG.get(), size );

#endif
	/* SIM */
	std::unique_ptr<SimDomain> sim ( new SimDomain( size, tau, FLAG.get() ) );
	/*************/

	/* SOLVER */
	printf("SOLVER START\n");

	clock_t Begin = clock();

	for( int t = 0; t <= timesteps; t++ ) {
#ifdef _MPI_
		sim->exchangePressureGhostCells( 0, 0, 0, widthComm, status );
		sim->exchangePressureGhostCells( 1, 2, 2, depthComm, status );
		sim->exchangePressureGhostCells( 2, 4, 4, heightComm, status );
#endif
		// STREAM:
		if( verbose ) { printf( "Streaming\n" ); }
        sim->stream();
        // SWAP:
		if( verbose ) { printf( "Swapping\n" ); }
        sim->swapFields();
        // COLLIDE:
		if( verbose ) { printf( "Colliding\n" ); }
        sim->collide();
        // TREAT BOUNDARIES:
		if( verbose ) { printf( "Treating boundaries\n" ); }
        sim->treatBoundaries();

		// VISUALIZATION:
		if( t % t_vis == 0 ) {
			if( verbose ) { printf( "Saving VTK at time %i\n", t ); }
#ifdef _MPI_
			sim->visualize( t, rank, widthRank, depthRank, heightRank );
#elif _SEQ_
			sim->visualize( t );
#endif
		}

		if( verbose ) {
			printf( "Time step %i completed\n", t );
		} else {
#ifdef _SEQ_
			printf("\rRunning Solver: %3.2f%%, time step: %i", ( 100.0 * t / timesteps ), t );
#elif _MPI_
			if(rank == 0 )
				printf("\rRunning Solver: %3.2f%%, time step: %i", ( 100.0 * t / timesteps ), t );
#endif
		}
	}
	/**********/
	clock_t End = clock();


	/* PERFORMANCE MEASURE */
	double ConsumedTime = (double)( End - Begin ) / CLOCKS_PER_SEC;

#ifdef _SEQ_
	
	double mLUPs = ( Q * totLen * timesteps ) / ( ConsumedTime * 1e6 );

#elif _MPI_
	double mLUPs;
	double Metrics[ 2 ] = { 0.0, 0.0 };

    // Metrics[ 0 ] contains MLUPS
    Metrics[ 0 ] = ( Q * totLen * timesteps ) / ( ConsumedTime * 1e6 );
    Metrics[ 1 ] = totLen;

    double TotalMetrics[ 2 ] = { 0.0, 0.0 };
    MPI_Reduce( Metrics,
                TotalMetrics,
                2,
                MPI_DOUBLE,
                MPI_SUM,
                0,
                MPI_COMM_WORLD );

    if( rank == 0 ) {
    	mLUPs = Metrics[ 0 ];
    	totLen = Metrics[ 1 ];
#endif
    
    printf("\nMLUPS: %f\nTotal Lattices: %i\n", mLUPs, Q * totLen );

#ifdef _MPI_
	}
#endif
	/**********/

#ifdef _MPI_
    MPI_Finalize();
#endif

    return 0;
}
