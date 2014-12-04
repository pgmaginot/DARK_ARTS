#ifndef Inputs_Allowed_h
#define Inputs_Allowed_h

/// Restart_type simulations allowed
enum RESTART_TYPE {MESH_REFINEMENT , RESTART , INVALID_RESTART_TYPE};

/// SPATIAL_DISCRETIZATION block
enum QUADRATURE_TYPE { GAUSS , LOBATTO , EQUAL_SPACED , INVALID_QUADRATURE_TYPE};

enum MATRIX_INTEGRATION { SELF_LUMPING , TRAD_LUMPING , EXACT , INVALID_MATRIX_INTEGRATION};

enum OPACITY_TREATMENT { MOMENT_PRESERVING , SLXS, INTERPOLATING, INVALID_OPACITY_TREATMENT};

/// ANGULAR_DISCRETIZATION block

enum ANGULAR_QUADRATURE_TYPE { GAUSS_ANGLE, LOBATTO_ANGLE, INVALID_ANGULAR_QUADRATURE_TYPE};

/// Introduced in the time TIME block

enum TIME_SOLVER {IMPLICIT_EULER , INVALID_TIME_SOLVER};

enum STARTING_METHOD { RAMP, EXPONENTIAL , VECTOR , INVALID_STARTING_METHOD };

/// Introduced in the REGIONS block

enum GRID_SPACING {EQUAL , LOG, INVALID_GRID_SPACING};

/// Introduced by the MATERIALS block
/**
  CONSTANT_XS - constant in space and temperature
  RATIONAL - function of temperature only
  TABLE_LOOKUP - function of temperature and group
  POLYNOMIAL_SPACE - function of space
*/
enum OPACITY_TYPE { CONSTANT_XS, RATIONAL, TABLE_LOOKUP , POLYNOMIAL_SPACE , INVALID_OPACITY_TYPE};

enum CV_TYPE { CONSTANT_CV , RATIONAL_CV , INVALID_CV_TYPE} ;

enum FIXED_SOURCE_TYPE { NO_SOURCE , MMS_SOURCE , INVALID_FIXED_SOURCE_TYPE } ;

/// For building Source if MMS_SOURCE_TYPE == MMS_SOURCE

enum TIME_MMS_TYPE { POLY_TIME, COS_TIME, INVALID_TIME_MMS_TYPE};

enum RADIATION_SPACE_MMS { RAD_POLY_SPACE , RAD_COS_SPACE , INVALID_RADIATION_SPACE_MMS };

enum RADIATION_ANGLE_MMS { MMS_ISOTROPIC , MMS_ANGLE_POLY , INVALID_RADIATION_ANGLE_MMS };

enum TEMPERATURE_SPACE_MMS { TEMP_POLY_SPACE , TEMP_COS_SPACE , INVALID_TEMPERATURE_SPACE_MMS }; 

/// Physical units, or unity for MMS / benchmark problems

enum UNITS_TYPE { CM_SH_KEV , UNITY, INVALID_UNITS_TYPE};

/// SOLVER block

enum WG_SOLVE_TYPE{ FP_SWEEPS, FP_DSA, KRYLOV_SWEEPS, KRYLOV_DSA, INVALID_WG_SOLVE_TYPE};

enum ARD_SOLVE_TYPE{ FP_NO_ACCEL, FP_LMFGA, KRYLOV_LMFGA, INVALID_ARD_SOLVE_TYPE};

/// IC_BC block
enum RADIATION_BC_TYPE{ VACUUM_BC, REFLECTIVE_BC, INCIDENT_BC, MMS_BC, INVALID_RADIATION_BC_TYPE};

enum INCIDENT_BC_VALUE_TYPE{ INCIDENT_CURRENT , INCIDENT_TEMPERATURE, INVALID_INCIDENT_BC_VALUE_TYPE};

enum BC_ENERGY_DEPENDENCE{ PLANCKIAN , INVALID_BC_ENERGY_DEPENDENCE};

enum BC_ANGLE_DEPENDENCE{ BC_ISOTROPIC , BC_GLANCE , BC_NORMAL, INVALID_BC_ANGLE_DEPENDENCE };

enum BC_TIME_DEPENDENCE{ BC_BURST , BC_CONSTANT , INVALID_BC_TIME_DEPENDENCE};



enum TEMPERATURE_IC_TYPE{ CONSTANT_TEMPERATURE_IC,  MMS_TEMPERATURE_IC, INVALID_TEMPERATURE_IC_TYPE};

enum RADIATION_IC_TYPE{ PLANCKIAN_IC, MMS_RADIATION_IC, INVALID_RADIATION_IC_TYPE};



#endif