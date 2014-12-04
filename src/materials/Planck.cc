// Planck.cc
//
// 11/4/2014 PGM- Changed to return answers in per steradians 
//
// 10/13/2010 WDH
//
// Modifed integrate_B and integrate_dBdT to return values that are NOT per
// steradian (just like the integrate_B_grey and integrate_dBdT_grey).
//
//
//
//Fri, Dec-17, 2010
//Using this class obtained from PDT.
//Has done some modifications to be compatible to my code.
//==============================================================================

#include "Planck.h"


using std::vector;

Planck::Planck(double accuracy_parameter, const Input_Reader& input_reader, const double sn_weights)
:
  m_sn_weight{sn_weights}
{
  accuracy = accuracy_parameter;
  pi = 3.14159265358979323846;
  if(input_reader.use_weird_units() )
  {
    if( input_reader.get_units_type() == UNITY )
    {
      h = 1.;    
      k = 1.;      
      c = 1.;   
      a = 1.;    
    }
    else
    {
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Other Unit Types Not coded in Planck.cc");
    }
  }
  else
  {
    h = 4.1356668e-15;     // ev s
    //h=PLANCK_CONSTANT_KeV_SH;
    k = 8.6173423e-5;        // ev K^-1
    //k=BOLTZMANN_CONSTANT_JK_KEV;
    c = 299792458.0*100;    // cm s^-1
    //c=SPEED_OF_LIGHT_SH;
    
    a = (8.0*pow(pi,5)*pow(k,4))/(15.0*pow(h,3)*pow(c,3)); // eV/cm^3-K^4
    //a=0.00472218;
  }

	gauss_quad();

}

double Planck::get_c(void) const
{
  return c;
}

double Planck::getTempForEnergy(double E) const
{
	if( E < 0 )
	{
		stringstream str;
		str << "Planckian integration requires temperature and energy to be positive.";
		//throw CommonException(str, CET_INTERNAL_ERROR);
	}
	return pow(E/(a*c),0.25); // Note that this is NOT a per steradian quantity
}

double Planck::integrate_B_grey(double T) const
{
	if( T < 0 )
	{
		stringstream str;
		str << "Planckian integration requires temperature to be positive.";
		//throw CommonException(str, CET_INTERNAL_ERROR);
	}
	return a*c*pow(T,4)/m_sn_weight; /// Note this is a per steradian quantity
}

double Planck::integrate_B(double T, double E_min, double E_max) const
{
	if( T == 0 || E_min == E_max )
		return 0;
	if( T < 0 || E_min < 0 || E_max < E_min )
	{
		stringstream str;
		str << "Planckian integration requires temperature and energy to be "
			<< "positive.";
		//throw CommonException(str, CET_INTERNAL_ERROR);
	}

	// unitless representation of energy and temperature
	double z1 = E_min / (k * T);
	double z2 = E_max / (k * T);

	// Bg is the group averaged planck function
	double Bg;


	// integrate by gaussian quadrature
	// ============
	if ( z2 <= .7 )
	{

		// accumulator and mapping values
		double gauss = 0;
		double g_mid = .5 * (E_max + E_min);
		double g_map = .5 * (E_max - E_min);

		for (unsigned int i=0; i<points.size(); i++)
			gauss += g_map * weights[i] * get_B(T, g_mid + g_map * points[i]);

		Bg = gauss;

	}
	// integrate by infinite sum
	// =================
	else if ( z1 >= .5 )
	{

		int N = 32;
		double SUM_one = 0;
		double SUM_two = 0;

		// determine number of terms necessary
		// the expression in the conditional expresses a bound on error
		SUM_one = exp(-z1) * (z1*z1*z1 + 3*z1*z1 + 6*z1 + 6);

		while ( exp(-(N+1)*z1)/(1-exp(-z1)) * pow((double)N+1,-4) * (pow((N+1)*z1,3) +
			3*pow((N+1)*z1,2) + 6*(N+1)*z1 + 6) / SUM_one  >  accuracy)
			N++;

		SUM_one=0;
		// compute integral from E_min to infinity
		for (int n=N; n!=0; n--)
		{
			SUM_one += exp(-n*z1)/pow((double)n,4) * (pow(n*z1,3)+3*pow(n*z1,2)+6*n*z1+6);
			SUM_two += exp(-n*z2)/pow((double)n,4) * (pow(n*z2,3)+3*pow(n*z2,2)+6*n*z2+6);
		}

		Bg = 2 * pow(k*T,4) * (SUM_one-SUM_two) / (pow(h,3) * pow(c,2));
	}
	// split interval
	// ==============
	else
	{

		// split point, .6 for now
		z1 = .6;
		double gauss = 0;
		double g_mid = .5 * (z1*k*T + E_min);
		double g_map = .5 * (z1*k*T - E_min);

		for (unsigned int i=0; i<points.size(); i++)
		{
			gauss += g_map * weights[i] * get_B(T, g_mid + g_map * points[i]);
		}

		int N = 32;
		double SUM_one = 0;
		double SUM_two = 0;

		SUM_one = exp(-z1) * (z1*z1*z1 + 3*z1*z1 + 6*z1 + 6);
		while ( exp(-(N+1)*z1)/(1-exp(-z1)) * pow((double)N+1,-4) * (pow((N+1)*z1,3) +
			3*pow((N+1)*z1,2) + 6*(N+1)*z1 + 6) / SUM_one  >  accuracy)
			N++;

		SUM_one=0;
		for (int n=N; n!=0; n--)
		{
			SUM_one += exp(-n*z1)/pow((double)n,4) * (pow(n*z1,3)+3*pow(n*z1,2)+6*n*z1+6);
			SUM_two += exp(-n*z2)/pow((double)n,4) * (pow(n*z2,3)+3*pow(n*z2,2)+6*n*z2+6);
		}

		Bg = gauss + 2 * pow(k*T,4)*(SUM_one-SUM_two) / (pow(h,3)*pow(c,2));

	}

	return Bg*4.0*pi/m_sn_weight; // Note that this is a per steradian quantity

} // Planck::integrate_B()


double Planck::integrate_dBdT_grey(double T) const
{
	return 4.0*a*c*pow(T,3)/m_sn_weight; // Note that this is a per steradian quantity
}


double Planck::integrate_dBdT(double T, double E_min, double E_max) const
{
	if( T == 0 || E_min == E_max )
		return 0;
	if( T < 0 || E_min < 0 || E_max < E_min )
	{
		stringstream str;
		str << "Integration of the temperature derivative of the Planckian "
			<< "requires temperature and energy to be positive.";
		//throw CommonException(str, CET_INTERNAL_ERROR);
	}

	// unitless representation of energy and temperature
	double z1 = E_min / (k * T);
	double z2 = E_max / (k * T);

	// dBgdT is the integral-dE of the partial-dT
	double dBgdT;


	// integrate by gaussian quadrature
	// ============
	if ( z2 <= .7 )
	{

		double gauss = 0;
		double g_mid = .5 * (E_max + E_min);
		double g_map = .5 * (E_max - E_min);


		for (unsigned int i=0; i<points.size(); i++)
			gauss += g_map * weights[i] * get_dBdT(T, g_mid + g_map * points[i]);

		dBgdT = gauss;

	}
	// integrate by infinite sum
	// =================
	else if ( z1 >= .5 )
	{

		int N = 32;
		double SUM_one = 0;
		double SUM_two = 0;

		SUM_one = exp(-z1) * (pow(z1,4) + 4*pow(z1,3) + 12*z1*z1 + 24*z1 + 24);
		while ( exp(-(N+1)*z1)/(1-exp(-z1)) * pow((double)N+1,-4) * (pow((N+1)*z1,4) +
			4*pow((N+1)*z1,3) + 12*pow((N+1)*z1,2) + 24*(N+1)*z1 + 24)
			/ SUM_one  >  accuracy)
			N++;
		SUM_one=0;

		for (int n=N; n!=0; n--)
		{
			SUM_one += exp(-n*z1)/pow((double)n,4) *(pow(n*z1,4) + 4*pow(n*z1,3) + 12*pow(n*z1,2) + 24*n*z1 + 24);
			SUM_two += exp(-n*z2)/pow((double)n,4) *
				(pow(n*z2,4) + 4*pow(n*z2,3) + 12*pow(n*z2,2) + 24*n*z2 + 24);
		}

		dBgdT = 2 * pow(k,4) * pow(T,3) * (SUM_one-SUM_two) / (pow(h,3) * pow(c,2));

	}
	// split interval
	// ==============
	else
	{

		z1 = .6;
		double gauss = 0;
		double g_mid = .5 * (z1*k*T + E_min);
		double g_map = .5 * (z1*k*T - E_min);

		for (unsigned int i=0; i<points.size(); i++)
		{
			gauss += g_map * weights[i] * get_dBdT(T, g_mid + g_map * points[i]);
		}

		int N = 32;
		double SUM_one = 0;
		double SUM_two = 0;

		SUM_one = exp(-z1) * (pow(z1,4) + 4*pow(z1,3) + 12*z1*z1 + 24*z1 + 24);
		while ( exp(-(N+1)*z1)/(1-exp(-z1)) * pow((double)N+1,-4) * (pow((N+1)*z1,4) +
			4*pow((N+1)*z1,3) + 12*pow((N+1)*z1,2) + 24*(N+1)*z1 + 24)
			/ SUM_one  >  accuracy)
			N++;

		SUM_one=0;
		for (int n=N; n!=0; n--)
		{
			SUM_one += exp(-n*z1)/pow((double)n,4) *
				(pow(n*z1,4) + 4*pow(n*z1,3) + 12*pow(n*z1,2) + 24*n*z1+24);
			SUM_two += exp(-n*z2)/pow((double)n,4) *
				(pow(n*z2,4) + 4*pow(n*z2,3) + 12*pow(n*z2,2) + 24*n*z2+24);
		}

		dBgdT = gauss + 2*pow(k,4)*pow(T,3)*(SUM_one-SUM_two) / (pow(h,3)*pow(c,2));

	}

	return dBgdT*4.0*pi/m_sn_weight; // Note that this is a per steradian quantity

} // Planck::integrate_dBdT() remake
void Planck::gauss_quad()
{
	// This function sets the points and weights for integrating the planck
	// function from low_bound to high_bound (in eV) by gaussian quadrature.
	// n specifies the order of the quadrature.

	// description of method
	// ====================================================
	//
	// A gaussian quadrature of order n is given by n points and n associated
	// weights.  The points are a direct mapping of the roots of the nth order
	// legendre polynomial, and the weights, which sum to the interval of
	// integration, are a function of the root values and the first derivative
	// of the legendre polynomial
	//
	// ====================================================


	// declaration of variables
	// ====================================================
	//

	// a quadrature order of 12 is sufficient to integrate B for z < 1
	// for more on this, see B_poly_z.cc
	const short unsigned int order = 12;

	// resize the vectors
	points.resize(order);
	weights.resize(order);

	// integer division sets midpoint given even or odd n
	short unsigned int midpoint = (order + 1) / 2;

	// sum of weights to normalize at end of function
	long double weight_sum = 0;

	// j, j-1, and j-2 order legendre polynomials at mu
	long double p_j, p_jminus1, p_jminus2;

	// first derivative of nth order legendre
	long double p_deriv;

	// mu represents the independent for the legendre polynomial, which is
	// energy mapped to (-1, 1).
	// we will iterate on mu to find the roots of the nth legendre polynomial
	long double mu, old_mu;

	// loop indices.  i represents quad points, j the order of leg polynomials
	short unsigned int i, j;

	// convergence parameter and flag
	// *** does this map straight? ***
	const double tolerance = accuracy;
	bool converged = false;

	// ====================================================


	// computations
	// ====================================================
	//

	for (i=0; i < midpoint; i++)
	{
		// set guess for mu; you'll have to ask kt about the logic!
		mu = cos(3.14159265358979323846264338 * (i + 0.75) / (order + 0.5));

		// loop until p_j at mu is sufficiently small to represent a root
		converged = false;
		while (!converged)
		{
			p_jminus1 = 0;
			p_j = 1;

			// set the nth order legendre polynomial at mu by recursion
			for (j=1; j<=order; j++)
			{
				p_jminus2 = p_jminus1;
				p_jminus1 = p_j;

				// this is the recursion relation for legendre polynomials
				p_j = ( (2*j - 1) * mu * p_jminus1 - (j - 1) * p_jminus2) / (j);
			}

			p_deriv = j * (mu * p_j - p_jminus1) / (mu*mu - 1);

			// iterating over mu to find the ith root
			old_mu = mu;
			mu = old_mu - p_j / p_deriv;

			if ( fabsl(mu - old_mu) < tolerance)
				converged = true;

		} // while !converged

		// set points utilizing symmetry
		points[i] = - mu;
		points[order-1 - i] = mu;

		// set weights, also symmetric
		weights[i] = 1 / ( (1 - mu*mu) * p_deriv*p_deriv);
		weights[order-1 - i] = weights[i];

		weight_sum += weights[i] + weights[order - 1 - i];
		if (i == order - 1 - i)
			weight_sum -= weights[i];


	} // for i to midpoint

	// normailize weights so that they sum to the length of the interval
	for (i=0; i<order; i++)
		weights[i] *= 2 / weight_sum;

	// ====================================================

} // Planck::gauss_quad()


// class Planck

