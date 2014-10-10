/** @file   Angular_Quadrature.cc
  *   @author pmaginot
  *   @brief Implement the Angular_Quadrature class, stores discrete ordinates, weights, Legendre polynomials
 */

#include "Angular_Quadrature.h"

Angular_Quadrature::Angular_Quadrature(const Input_Reader& input_reader, const Quadrule_New& quad_fun)
:   
  m_n_dir{input_reader.get_number_of_angles() },
  m_n_groups{input_reader.get_number_of_groups() },   
  m_n_legendre_moments{input_reader.get_number_of_legendre_moments() }
{
  ANGULAR_QUADRATURE_TYPE quad_type = input_reader.get_angular_quadrature_type();
 
  m_mu.resize(m_n_dir,0.);
  m_w.resize(m_n_dir,0.);
  
  const int n_leg_evals = m_n_legendre_moments * m_n_dir;
  m_legendre_poly.resize(n_leg_evals,0.);
  
  if(quad_type == GAUSS_ANGLE)
  {
    quad_fun.legendre_dr_compute( m_n_dir, m_mu, m_w );
  }
  else if(quad_type == LOBATTO_ANGLE)
  {
    quad_fun.lobatto_compute( m_n_dir, m_mu, m_w );
  }
  
  m_sum_w = 0.;
  for(int d = 0;d<m_n_dir;d++)
    m_sum_w += m_w[d];
  
  /// calculate legendre polynomials of discrete ordinates
  Legendre_Poly_Evaluation leg_poly;
  for(int d=0; d<m_n_dir ; d++)
  {
    leg_poly.get_evaluated_legendre_polynomials( m_mu[d] , m_n_legendre_moments - 1 , d*m_n_legendre_moments, m_legendre_poly );
  }

  /// add in \f$ \frac{2l + 1}{2} \f$
  int cnt = 0;
  for(int d=0; d<m_n_dir ; d++)
  {
    for(int l=0;l<m_n_legendre_moments;l++)
    {
      m_legendre_poly[cnt] *= (2.*double(l) + 1.)/m_sum_w;
      cnt++;
    }
  }
}
    
int Angular_Quadrature::get_number_of_dir(void) const
{
  return m_n_dir;
}

int Angular_Quadrature::get_number_of_groups(void) const{
  return m_n_groups;
}

int Angular_Quadrature::get_number_of_leg_moments(void) const
{
  return m_n_legendre_moments;
}

double Angular_Quadrature::get_leg_poly(const int dir, const int mom) const
{
  return m_legendre_poly[m_n_legendre_moments*dir + mom];
}

double Angular_Quadrature::get_mu(const int dir) const
{
  return m_mu[dir];
}

double Angular_Quadrature::get_w(const int dir) const
{
  return m_w[dir];
}

double Angular_Quadrature::get_sum_w(void) const
{
  return m_sum_w;
}

double Angular_Quadrature::calculate_boundary_conditions(const int dir, const int grp, const double time) const
{  
  return 1.;
}