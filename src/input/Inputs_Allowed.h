#ifndef Inputs_Allowed_h
#define Inputs_Allowed_h

/// SPATIAL_DISCRETIZATION block

enum QUADRATURE_TYPE { GAUSS , LOBATTO , EQUAL_SPACED , INVALID_QUADRATURE_TYPE};

enum MATRIX_INTEGRATION { SELF_LUMPING , TRAD_LUMPING , EXACT , INVALID_MATRIX_INTEGRATION};

enum OPACITY_TREATMENT { MOMENT_PRESERVING , SLXS, INTERPOLATING, INVALID_OPACITY_TREATMENT};

/// ANGULAR_DISCRETIZATION block

enum ANGULAR_QUADRATURE_TYPE { GAUSS_ANGLE, LOBATTO_ANGLE, INVALID_ANGULAR_QUADRATURE_TYPE};

/// Introduced in the time TIME block

enum TIME_SOLVER {IMPLICIT_EULER , INVALID_TIME_SOLVER};

enum STARTING_METHOD { EXPONENTIAL , VECTOR , INVALID_STARTING_METHOD };

/// Introduced in the REGIONS block

enum GRID_SPACING {EQUAL , LOG, INVALID_GRID_SPACING};

/// Introduced by the MATERIALS block

enum OPACITY_TYPE { CONSTANT_XS, RATIONAL, TABLE_LOOKUP ,  INVALID_OPACITY_TYPE};

enum CV_TYPE { CONSTANT_CV , INVALID_CV_TYPE} ;

enum FIXED_SOURCE_TYPE { NO_SOURCE , INVALID_FIXED_SOURCE_TYPE } ;

/// SOLVER block

enum WG_SOLVE_TYPE{ FP_SWEEPS, FP_DSA, KRYLOV_SWEEPS, KRYLOV_DSA, INVALID_WG_SOLVE_TYPE};

#endif