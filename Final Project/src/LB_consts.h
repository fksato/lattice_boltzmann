
#ifndef LB_CONST_H
#define LB_CONST_H

namespace LBM {

    static const int LATTICEVELOCITIES[19][3] = {{0,  -1, -1},
                                                 {-1, 0,  -1},
                                                 {0,  0,  -1},
                                                 {1,  0,  -1},
                                                 {0,  1,  -1},
                                                 {-1, -1, 0},
                                                 {0,  -1, 0},
                                                 {1,  -1, 0},
                                                 {-1, 0,  0},
                                                 {0,  0,  0},
                                                 {1,  0,  0},
                                                 {-1, 1,  0},
                                                 {0,  1,  0},
                                                 {1,  1,  0},
                                                 {0,  -1, 1},
                                                 {-1, 0,  1},
                                                 {0,  0,  1},
                                                 {1,  0,  1},
                                                 {0,  1,  1}};

    static const double LATTICEWEIGHTS[19] = {1.0 / 36.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                              1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 3.0,
                                              1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                              1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0};

    // BOTTOM, FRONT, LEFT, RIGHT, BACK, TOP
    static const int CELLCENTERS[6] = { 2, 6, 8, 10, 12, 16 };
    static const int CELLFACES[6][5] = { 
                                            { 0, 1, 2, 3, 4 }, { 0, 5, 6, 7, 14 }, { 1, 5, 8, 11, 15 }
                                            , { 3, 7, 10, 13, 17 }, { 4, 11, 12, 13, 18 }, { 14, 15, 16, 17, 18 }
                                            };

    static const int Q = 19;
    static const int DIM = 3;
    static const double CS = 0.577350269189626;

}
#endif