
#include "simulation.h"
#include "LB_consts.h"
#include "visual.h"
#include <cstdlib>
#include <cstdio>
#include <cassert>

#ifdef _MPI_
#include <mpi.h>
#include "helper.h"
#include <vector>
#endif

using namespace LBM;

SimDomain::SimDomain( int *domSize, double inTau, int *flagField )
						: tau(inTau)
{
	// Init problem domain:
	size[0] = domSize[0];
	size[1] = domSize[1];
	size[2] = domSize[2];

	xlen = size[0] + 2;
	ylen = size[1] + 2;
	zlen = size[2] + 2;

	xylen = xlen * ylen;

	int totLen = xlen * ylen * zlen;

    visDensity = std::unique_ptr<double[]>( new double[ xylen * zlen ] );
    visVelocity = std::unique_ptr<double[]>( new double[  DIM * xylen * zlen ] );

	collisionField = std::unique_ptr<double[]>( new double[ Q * totLen ] );
	streamField = std::unique_ptr<double[]>( new double[ Q * totLen ] );
	FLAG = flagField;

	int idx;
	for( int x = 0; x < xlen; x++ ) {
		for( int y = 0; y < ylen; y++ ) {
			for( int z = 0; z < zlen; z++ ) {

				for( int i = 0; i < Q; i++ ){

					idx = Q * ( z * xylen + y * xlen + x ) + i ;

					collisionField[ idx ] 	= LATTICEWEIGHTS[i];
					streamField[ idx ] 		= LATTICEWEIGHTS[i];
				}
			}
		}
	}
}

void SimDomain::stream() {

	double fi;
	int dx, dy, dz, zd, yd, xd, collisionIDX, flagIDX;

	for( int x = 1; x <= size[0]; x++) {
		for( int y = 1; y <= size[1]; y++ ) {
			for( int z = 1; z <= size[2]; z++ ) {
				flagIDX = z * xylen + y * xlen + x;

				if( FLAG[ flagIDX ] == OBS )
					continue;

				collisionIDX = Q * flagIDX;

				for (int i = 0; i < Q; ++i) {

					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					xd = x - dx;
					yd = y - dy;
					zd = z - dz;

					fi = collisionField[ Q * ( zd * xylen + yd * xlen + xd ) + i ];

					streamField[ collisionIDX + i ] = fi;
				}

			}
		}
	}
}

void SimDomain::swapFields() {
    streamField.swap(collisionField);
}

void SimDomain::collide() {

	int cellidx;
	double density = 0.0;
	double velocity[DIM] = { 0.0 };
	double f_eq[Q] = { 0.0 };

    for( int x = 1; x <= size[0]; x++) {
        for( int y = 1; y <= size[1]; y++ ) {
            for( int z = 1; z <= size[2]; z++ ) {

                cellidx = Q * (z * xylen + y * xlen + x);

                computeDensity( cellidx, &density );
                computeVelocity( cellidx, &density, velocity );
                computeFEQ( velocity, density, f_eq );

                computePostCollisionDistributions ( cellidx, f_eq );
            }
        }
    }
}

void SimDomain::computeDensity( int cellidx, double *density ) {
    *( density ) = 0.0;

    for (int i=0; i<Q; i++){
        *( density ) += collisionField[ cellidx + i ];
    }

    assert( *( density ) >= 0.0);
}

void SimDomain::computeVelocity( int cellidx, double *density, double *velocity ) {

    double temp[3] = { 0.0 };

	for( int i = 0; i < Q; i++) {
		temp[0] += collisionField[ cellidx + i ] * LATTICEVELOCITIES[i][0];
		temp[1] += collisionField[ cellidx + i ] * LATTICEVELOCITIES[i][1];
		temp[2] += collisionField[ cellidx + i ] * LATTICEVELOCITIES[i][2];
	}

	velocity[0] = temp[0] / *(density);
	velocity[1] = temp[1] / *(density);
	velocity[2] = temp[2] / *(density);
}

void SimDomain::computeFEQ( double *velocity, double density, double *f_eq ) {

	double dot_c_u;
	double cssqr = 1.0 / ( CS * CS) ;
	double dot_u_u = 0.0;
	double constFactor;

	dot_u_u += velocity[0] * velocity[0]
				+ velocity[1] * velocity[1]
				+ velocity[2] * velocity[2];

	constFactor = 1 - 0.5 * dot_u_u * ( cssqr );

	for( int i = 0; i < Q; i++ ) {

		// Compute dot product:
		dot_c_u = LATTICEVELOCITIES[ i ][ 0 ] * velocity[ 0 ]
					+ LATTICEVELOCITIES[ i ][ 1 ] * velocity[ 1 ]
					+ LATTICEVELOCITIES[ i ][ 2 ] * velocity[ 2 ];

		f_eq[ i ] = LATTICEWEIGHTS[ i ] * density * ( constFactor + dot_c_u * cssqr
														+ 0.5 * ( dot_c_u ) * ( dot_c_u ) * cssqr * cssqr
													);
	}
}

void SimDomain::computePostCollisionDistributions( int cellidx, double *f_eq ) {
	double _invTau = 1.0 / tau;

	for( int i = 0; i < Q; i++ ) {
		collisionField[ cellidx + i ] -= _invTau * ( collisionField[ cellidx + i ] - f_eq[ i ] );
	}
}

void SimDomain::treatBoundaries() {

	int cellIDX;
	int flagIDX;
	int neighbors[ NEIGHBORS ];

	for( int x = 0; x < xlen; x++ ) {
        for( int y = 0; y < ylen; y++ ) {
			for( int z = 0; z < zlen; z++ ) {

				flagIDX = z * xylen + y * xlen + x;
				cellIDX = Q * flagIDX;

				if( FLAG[ flagIDX ] == FLUID )
					continue;

				getNeighbors( x, y, z, neighbors );

				switch( FLAG[ flagIDX ] ) {
					case NOSLIP:
						noSlip( cellIDX, neighbors );
						break;
					case FREESLIP:
						freeSlip( cellIDX, x, y, z );
						break;
					case INFLOW:
						inflow( cellIDX, neighbors );
						break;
					case OUTFLOW:
						outflow( cellIDX, neighbors );
						break;
					case MOVINGWALL:
						movingWall( cellIDX, neighbors );
						break;
			    default:
		        break;
				}
			}
		}
	}
}

void SimDomain::getNeighbors( int x, int y, int z, int *neighbors ) {
	int dx, dy, dz, xd, yd, zd;
	int neighborFLAGIDX;

	for (int i = 0; i < Q; ++i) {

		neighbors[i] = -1;

		dx = LATTICEVELOCITIES[i][0];
		dy = LATTICEVELOCITIES[i][1];
		dz = LATTICEVELOCITIES[i][2];
		if( dx == 0 && dy == 0 && dz == 0)
			continue;
		xd = x + dx;
		yd = y + dy;
		zd = z + dz;

		if ( xd > 0 && xd <= size[0] && yd > 0 && yd <= size[1] && zd > 0 && zd <= size[2] ) {
			neighborFLAGIDX =  zd * xylen + yd * xlen + xd;

			if( FLAG[ neighborFLAGIDX ] != FLUID )
				continue;

			neighbors[i] = Q * neighborFLAGIDX;
		}
	}
}

void SimDomain::noSlip( int cellidx, int *neighbors ) {
	int neighborCellIDX;

	for( int i = 0; i < NEIGHBORS; i++ ) {
		neighborCellIDX = neighbors[i];

		if( neighborCellIDX == -1 )
			continue;

		collisionField[ cellidx + i ] = collisionField[ neighborCellIDX + ( Q - 1 ) - i ];
	}
}

void SimDomain::freeSlip( int cellidx, int x, int y, int z ) {

	int cellCenter, neighborFLAGIDX, neighborCellIDX, face_vel, midx, mirror_vel;
	int nx, ny, nz;
	for( int i = 0; i < FACES; i++ ) {
		cellCenter = CELLCENTERS[i];

		nx = x - LATTICEVELOCITIES[ cellCenter ][0];
		ny = y - LATTICEVELOCITIES[ cellCenter ][1];
		nz = z - LATTICEVELOCITIES[ cellCenter ][2];

		if ( nx > 0 && nx <= size[0] && ny > 0 && ny <= size[1] && nz > 0 && nz <= size[2] ) {
			neighborFLAGIDX = nz * xylen + ny * xlen + nx;

			if( FLAG[ neighborFLAGIDX ] != FLUID )
				continue;

			neighborCellIDX = Q * neighborFLAGIDX;
			// face velocities:
			midx = ( FACES - 1 ) - i;
			for( int j = 0; j < SHAREDVEL; j++ ) {
				face_vel = CELLFACES[ i ][ j ];
				mirror_vel = CELLFACES[ midx ][ j ];

				collisionField[ cellidx + face_vel ] = collisionField[ neighborCellIDX + mirror_vel ];
			}
		}
	}
}

void SimDomain::inflow( int cellidx, int *neighbors ) {
	computeFEQ( inflowVelocity, p_ref_in, &collisionField[ cellidx ] );
}

void SimDomain::outflow( int cellidx, int *neighbors ) {

	int inv_i, neighborCellIDX;
	double density;
	double feq[Q]= { 0.0 };
	double velocity[DIM] = { 0.0 };

	for( int i = 0; i < NEIGHBORS; i++ ) {
		neighborCellIDX = neighbors[i];

		if( neighborCellIDX == -1 )
			continue;

		computeDensity( neighborCellIDX, &density );
		computeVelocity( neighborCellIDX, &density, velocity);
		computeFEQ( velocity, p_ref_out, feq );

		inv_i = ( Q-1 ) - i;
		collisionField[ cellidx + i] = feq[inv_i] + feq[i] - collisionField[ neighborCellIDX + inv_i ];
	}
}

void SimDomain::movingWall( int cellidx, int *neighbors ) {

	int neighborCellIDX;

	for( int i = 0; i < NEIGHBORS; i++ ) {
		neighborCellIDX = neighbors[i];

		if( neighborCellIDX == -1 )
			continue;

		double dot_cu = 0;
		double neighborDensity;
		double cssqr = 1.0 / ( CS * CS );

		dot_cu =      wallVelocity[0] * LATTICEVELOCITIES[ i ][0]
					+ wallVelocity[1] * LATTICEVELOCITIES[ i ][1]
					+ wallVelocity[2] * LATTICEVELOCITIES[ i ][2];

		computeDensity( neighborCellIDX, &neighborDensity );

		collisionField[ cellidx + i ] = collisionField[ neighborCellIDX + ( Q - 1 ) - i ]
											+ 2.0 * LATTICEWEIGHTS[ i ] * neighborDensity * dot_cu * cssqr;
	}
}

#ifdef _SEQ_
void SimDomain::visualize( int t ) {

    double density = 0.0;
    double velocity[DIM] = { 0.0 };
    int idx, cellidx;

	for(int x = 1; x <= size[0]; x++ ) {
		for(int y = 1; y <= size[1]; y++ ) {
			for(int z = 1; z <= size[2]; z++ ) {

				idx = z * xylen + y * xlen + x;
				cellidx = Q * idx;

				computeDensity( cellidx, &density );
				computeVelocity( cellidx, &density, velocity );

				visDensity[ idx ] = density;
                visVelocity[ DIM * idx + 0 ] = velocity[ 0 ];
                visVelocity[ DIM * idx + 1 ] = velocity[ 1 ];
                visVelocity[ DIM * idx + 2 ] = velocity[ 2 ];
			}
		}
	}

	writeVtkOutput( "VTK/SEQ/simulation", t, size, visDensity.get(), visVelocity.get() );
}

#elif _MPI_
void SimDomain::visualize( int t, int rank, int widthRank, int depthRank, int heightRank ) {

    double density = 0.0;
    double velocity[DIM] = { 0.0 };
    int idx, cellidx;

	for(int x = 1; x <= size[0]; x++ ) {
		for(int y = 1; y <= size[1]; y++ ) {
			for(int z = 1; z <= size[2]; z++ ) {

				idx = z * xylen + y * xlen + x;
				cellidx = Q * idx;

				computeDensity( cellidx, &density );
				computeVelocity( cellidx, &density, velocity );

				visDensity[ idx ] = density;
                visVelocity[ DIM * idx + 0 ] = velocity[ 0 ];
                visVelocity[ DIM * idx + 1 ] = velocity[ 1 ];
                visVelocity[ DIM * idx + 2 ] = velocity[ 2 ];
			}
		}
	}

	writeVtkOutput( "VTK/MPI/simulation", t, rank, widthRank, depthRank, heightRank, size, visDensity.get(), visVelocity.get() );
}

void SimDomain::exchangePressureGhostCells( int axis, int sendtag, int recvtag, MPI_Comm comm, MPI_Status status )
{
	int i, j, k, idx_f, idx_b, array_size, cnt;
	std::unique_ptr<double[]> send_f;
	std::unique_ptr<double[]> send_b;

	std::unique_ptr<double[]> recv_f;
	std::unique_ptr<double[]> recv_b;

	if( axis == 0 ) {			// left / right send
		cnt = 0;
		array_size = 5 * ( size[1] ) * ( size[2] );
		send_f = std::unique_ptr<double[]>(new double[array_size]);
		send_b = std::unique_ptr<double[]>(new double[array_size]);

		recv_f = std::unique_ptr<double[]>(new double[array_size]);
		recv_b = std::unique_ptr<double[]>(new double[array_size]);

		for( j = 1; j <= size[1]; j++ ) {
			for( k = 1 ; k <= size[2]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( k * xylen + j * xlen + 1 );
					idx_b = Q * ( k * xylen + j * xlen + ( size[0] ) );

					send_f[ cnt ] = collisionField[ idx_f + CELLFACES[2][i] ];
					send_b[ cnt ] = collisionField[ idx_b + CELLFACES[3][i] ];
					cnt++;
				}
			}
		}
	} else if ( axis == 1 ) {  	// front / back send
		cnt = 0;
		array_size = 5 * ( size[0] ) * ( size[2] );
		send_f = std::unique_ptr<double[]>(new double[array_size]);
		send_b = std::unique_ptr<double[]>(new double[array_size]);

		recv_f = std::unique_ptr<double[]>(new double[array_size]);
		recv_b = std::unique_ptr<double[]>(new double[array_size]);

		for( j = 1; j <= size[0]; j++ ) {
			for( k = 1 ; k <= size[2]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( k * xylen + 1 * xlen + j );
					idx_b = Q * ( k * xylen + ( size[1] ) * xlen + j );

					send_f[ cnt ] = collisionField[ idx_f + CELLFACES[1][i] ];
					send_b[ cnt ] = collisionField[ idx_b + CELLFACES[4][i] ];
					cnt++;
				}
			}
		}
	} else {					// top / bottom send
		cnt = 0;
		array_size = 5 * ( size[0] ) * ( size[1] );
		send_f = std::unique_ptr<double[]>(new double[array_size]);
		send_b = std::unique_ptr<double[]>(new double[array_size]);

		recv_f = std::unique_ptr<double[]>(new double[array_size]);
		recv_b = std::unique_ptr<double[]>(new double[array_size]);

		for( j = 1; j <= size[0]; j++ ) {
			for( k = 1 ; k <= size[1]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( 1 * xylen + k * xlen + j );
					idx_b = Q * ( ( size[2] ) * xylen + k * xlen + j );

					send_f[ cnt ] = collisionField[ idx_f + CELLFACES[0][i] ];
					send_b[ cnt ] = collisionField[ idx_b + CELLFACES[5][i] ];
					cnt++;
				}
			}
		}
	}

	int lsrc, ldest;
	// send/recv data from top cell to bottom, right to left, back to front
	MPI_Cart_shift( comm, 0, -1, &lsrc, &ldest );
	MPI_Sendrecv( send_f.get(), array_size, MPI_DOUBLE, ldest, sendtag, recv_f.get(), array_size, MPI_DOUBLE, lsrc, recvtag, comm, &status );
	int rsrc, rdest;
	// send/recv data from bottom cell to top, left to right, front to back
	MPI_Cart_shift( comm, 0, 1, &rsrc, &rdest );
	MPI_Sendrecv( send_b.get(), array_size, MPI_DOUBLE, rdest, sendtag+1, recv_b.get(), array_size, MPI_DOUBLE, rsrc, recvtag+1, comm, &status );

	// Reconstruct:
	if( axis == 0 ) {			// left / right send
		cnt = 0;
		for( j = 1; j <= size[1]; j++ ) {
			for( k = 1 ; k <= size[2]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( k * xylen + j * xlen + 0 );
					idx_b = Q * ( k * xylen + j * xlen + ( size[0] + 1 ) );

					if( !(rsrc == -1) )
						collisionField[ idx_f + CELLFACES[2][i] ] = recv_f[ cnt ];
					if( !(lsrc == -1) )
						collisionField[ idx_b + CELLFACES[3][i] ] = recv_b[ cnt ];
					cnt++;
				}
			}
		}
	} else if ( axis == 1 ) {  	// front / back send
		cnt = 0;
		for( j = 1; j <= size[0]; j++ ) {
			for( k = 1 ; k <= size[2]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( k * xylen + 0 * xlen + j );
					idx_b = Q * ( k * xylen + ( size[1] + 1 ) * xlen + j );

					if( !(rsrc == -1) )
						collisionField[ idx_f + CELLFACES[1][i] ] = recv_f[ cnt ];
					if( !(lsrc == -1) )
						collisionField[ idx_b + CELLFACES[4][i] ] = recv_b[ cnt ];
					cnt++;
				}
			}
		}
	} else {					// top / bottom send
		cnt = 0;
		for( j = 1; j <= size[0]; j++ ) {
			for( k = 1 ; k <= size[1]; k++ ) {
				for( i = 0; i < 5; i++ ) {
					idx_f = Q * ( 0 * xylen + k * xlen + j );
					idx_b = Q * ( ( size[2] + 1 ) * xylen + k * xlen + j );

					if( !(rsrc == -1) )
						collisionField[ idx_f + CELLFACES[0][i] ] = recv_f[ cnt ];
					if( !(lsrc == -1) )
						collisionField[ idx_b + CELLFACES[5][i] ] = recv_b[ cnt ];
					cnt++;
				}
			}
		}
	}
}
#endif
