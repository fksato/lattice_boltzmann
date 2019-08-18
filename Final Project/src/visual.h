#ifndef __VISUAL_H__
#define __VISUAL_H__

namespace LBM {
#ifdef _SEQ_
    void writeVtkOutput( const char * filename, int t, int const *size, double *density, const double *velocity );
#elif _MPI_
    void writeVtkOutput( const char * filename, int t, int rank, int widthRank, int depthRank, int heightRank
    					, int const *size, double *density, const double *velocity );
#endif
}
#endif
