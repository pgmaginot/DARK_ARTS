/** @file   Absorption_Opacity_Adams_Novak_Model.cc
  *   @author pmaginot
  *   @brief Implement the concrete Absorption_Opacity_Adams_Novak_Model class
  *   \f$ \sigma_a = \frac{c_1}{c_2 + T^P} \f$
*/
#include "Absorption_Opacity_Adams_Novak_Model.h"

Absorption_Opacity_Adams_Novak_Model::Absorption_Opacity_Adams_Novak_Model( const Input_Reader& input_reader) 
{

}

Absorption_Opacity_Adams_Novak_Model::~Absorption_Opacity_Adams_Novak_Model()
{

}

double Absorption_Opacity_Adams_Novak_Model::get_absorption_opacity(const int group, const double temperature, const double position)
{
  return 0.;
}


// /*!
   // The model opacity  given in "Asymptotic Analysis of a Computational Method
   // for Time- and Frequency- Dependent Radiative Transfer", Adams and Nowak,
   // Journal of Computational Physics 146, 366-403 (1998) is:

    // \f$ \sigma(E,T) = \sigma_0 frac(1-exp(-E/kT),E^3) \f$.

    // To find multigroup (Planck weighted) model opacities, we need

    // \f$ \int_E1^E2 {\sigma(E,T)*B(E,T)}dE / \int_E1^E2 {B(E,T)}dE \f$

    // which yields

    // \f$ 2*\sigma_0/(h^3 c^2)*(exp(-E1/kT)-exp(-E2/kT) /
    // \int_E1^E2 {B(E,T)}dE \f$.
// */

// inline double Opacities::get_model_opacity(double temperature, double E_min,
  // double E_max, opac_type type)
// {
  // // Note that our units for B are NOT per steradian

  // double mod;

  // if( type == OT_COUNT )
    // type = OT_default;

  // assert( temperature >= 0.0 );
  // assert( E_min >= 0.0 );
  // assert( E_max >= E_min );

  // if( type == OT_model_ross )
    // throw Dark_Arts_Exception(SUPPOSRT_OBJECT, "Rosseland-weighted model opacity has not been implemented");

  // double K  = BOLTZMANN_CONSTANT; // eV K^-1
  // double H  = PLANCK_CONSTANT;    // eV s
  // double C  = SPEED_OF_LIGHT;     // cm s^-1
  // double Pi = PI;

  // Planck planck1;
  // double Bg = planck1.integrate_B(temperature, E_min, E_max);

  // if (Bg == 0)
  // {
    // mod = OPACITY_MAX;
    // // I need to do this integral properly; but 10^6 is better than:
    // // double E = .5*(E_min + E_max);
    // // mod = mod_const*(1-exp(-E/(k*temperature))*std::pow(E,-3));
  // }
  // else
  // {
    // double opacity = 8 * Pi * mod_const * std::pow(H,-3.0) * std::pow(C,-2.0) *
     // (exp(-E_min/(K*temperature)) - exp(-E_max/(K*temperature)));
    // mod = opacity / Bg ;
  // }
  // return mod;
// }
