#include "Materials.h"

/** @file   Materials.cc
  *   @author pmaginot
  *   @brief Create material property objects for each region, and evaluate material properties
  *   material properties are considered to be \f$ \sigma_a,~\sigma_s,~C_v \f$ and the Planckian emission term
  *   Materials class will also handle the case when the number of integration points is not equal to the number
  *   of cross section evaluation points (e.g. cell-wise constant xs vs. self-lumping vs. moment preserving vs. 
  *   cross section interpolation not equal to temperature interpolation
*/

/* *****************************************************************
*   Public Member functions
*  ****************************************************************/

Materials::Materials( const Input_Reader& input_reader, 
  const Fem_Quadrature& fem_quadrature, 
  const Cell_Data& cell_data,
  const Angular_Quadrature& angular_quadrature)
:
  m_planck( 1.0E-15 , input_reader, angular_quadrature.get_sum_w() ),
  m_num_materials(input_reader.get_number_of_materials()),
  m_n_xs_quad( fem_quadrature.get_number_of_xs_point() ),
  m_n_el_cell( fem_quadrature.get_number_of_interpolation_points() ),
  m_dx(-1.),
  m_cell_data( cell_data),
  m_n_source_pts( fem_quadrature.get_number_of_source_points() ),
  m_left_edge_set(false),
  m_right_edge_set(false),
  m_current_material(-1)
{ 
  try{
    /// resize material property objects
    m_abs_opacities.resize(m_num_materials);
    m_scat_opacities.resize(m_num_materials);
    m_cv_obj.resize(m_num_materials);
    m_source_t.resize(m_num_materials);
    m_source_i.resize(m_num_materials);  
    
    /// load material property objects for each different material
    load_materials(input_reader , angular_quadrature);
      
    /// create smart pointers that look at input once, then select cross section treatment type
    switch( input_reader.get_opacity_treatment() )
    {
      case SLXS:
      {
        m_xs_treatment = std::make_shared<XS_Treatment_SLXS>(fem_quadrature) ;
        break;
      }
      case INTERPOLATING:
      {
        m_xs_treatment =  std::make_shared<XS_Treatment_Interpolating>(fem_quadrature );
        break;
      }
      case MOMENT_PRESERVING:
      {
        m_xs_treatment =  std::make_shared<XS_Treatment_Moment_Preserving>(fem_quadrature , input_reader) ;
        break;
      }
      case INVALID_OPACITY_TREATMENT:
      {
        throw Dark_Arts_Exception( SUPPORT_OBJECT ,  "Invalid OPacity treament still exists in Materials initializations" );
        break;
      }    
    }
    
    /// get cross section evaluation quadrature points, dfem evaluations at integration points and cell edges
    fem_quadrature.get_xs_eval_points(m_xs_eval_quad);  
    fem_quadrature.get_dfem_at_xs_eval_points(m_dfem_at_xs);    
    fem_quadrature.get_dfem_at_edges(m_dfem_at_left_edge,m_dfem_at_right_edge);  
    fem_quadrature.get_dfem_at_source_points(m_dfem_at_source);
        
    /// allocate space for physical position and temperature at material property evaluation points
    m_xs_position.resize(m_n_xs_quad,0.);  
    m_t_at_xs_eval_points.resize(m_n_xs_quad,0.);
    /// allocate scratch vector space for material property evaluations
    m_mat_property_evals.resize(m_n_xs_quad,0.);
    
    /// qet source moment quadrature points and allocate space for source moment evaluations
    fem_quadrature.get_source_points(m_source_quad);
    m_position_at_source_quad.resize(m_n_source_pts,0.);
    m_temperature_at_source_quad.resize(m_n_source_pts,0.);
    
    /// get energy bounds 
    if(angular_quadrature.get_number_of_groups() > 1)
    {
      /// having energy group bounds only makes sense if this is a grey problem
      input_reader.get_lower_energy_bounds(m_grp_e_min);
      input_reader.get_upper_energy_bounds(m_grp_e_max);
    }
  }
  catch(const Dark_Arts_Exception& da_exception)
  {
    da_exception.message();
  }
}

void Materials::calculate_left_edge_temp_and_position(const int cell_num, const Eigen::VectorXd& temperature)
{
  m_current_material = m_cell_data.get_cell_material_number(cell_num);
  
  /// calculate the local position at the xs evalaution points
  m_x_left = m_cell_data.get_cell_left_edge(cell_num);;
      
  m_t_left_bound = 0.; 
  
  for(int dfem_t = 0; dfem_t < m_n_el_cell ; dfem_t++)
    m_t_left_bound += temperature(dfem_t)*m_dfem_at_left_edge[dfem_t];
   
  m_left_edge_set = true;
  
  return;
}

void Materials::calculate_right_edge_temp_and_position(const int cell_num, const Eigen::VectorXd& temperature)
{
  m_current_material = m_cell_data.get_cell_material_number(cell_num);
  
  /// calculate the local position at the xs evalaution points
  m_x_right = m_cell_data.get_cell_left_edge(cell_num) + m_cell_data.get_cell_width(cell_num);

  m_t_right_bound = 0.;  
  for(int dfem_t = 0; dfem_t < m_n_el_cell ; dfem_t++)
    m_t_right_bound += temperature(dfem_t)*m_dfem_at_right_edge[dfem_t];
    
  m_right_edge_set = true;
   
  return;
}

void Materials::calculate_local_temp_and_position(const int cell_num, const Eigen::VectorXd& temperature)
{
  /// determine what material we are in
  m_current_material = m_cell_data.get_cell_material_number(cell_num);
  
  /// calculate the local position at the xs evalaution points
  m_dx = m_cell_data.get_cell_width(cell_num);
  double xL = m_cell_data.get_cell_left_edge(cell_num);
  
  /// for evaluating material properties
  for(int i=0; i< m_n_xs_quad; i++)
    m_xs_position[i] = xL + m_dx/2.*(1. + m_xs_eval_quad[i]);
  
  /// for evaluating driving source moments
  for(int i=0; i< m_n_source_pts; i++)
    m_position_at_source_quad[i] = xL + m_dx/2.*(1. + m_source_quad[i]);
   
  /// calculate the local temperature at the xs evaluation points
  /// basis function evaluation is laid out in a vector: B_{1,1} ... B_{1,N_eval} B_{2,1} ...
  for(int i=0; i < m_n_xs_quad; i++)
    m_t_at_xs_eval_points[i] = 0.;    
  
  for(int i=0; i < m_n_source_pts; i++)
    m_temperature_at_source_quad[i] = 0.;
    
  int p=0;
  int src_p = 0;
  for(int dfem_t = 0; dfem_t < m_n_el_cell ; dfem_t++)
  {
    /// temperature is an Eigen::VectorXd
    double t_el = temperature(dfem_t);
     for(int i=0; i< m_n_xs_quad; i++)
    {
      m_t_at_xs_eval_points[i] += t_el*m_dfem_at_xs[p];
      p++;
    }   
    for(int j=0; j< m_n_source_pts; j++)
    {
      m_temperature_at_source_quad[j] += t_el*m_dfem_at_source[src_p];
      src_p++;
    }      
  }
  return;
}
     
void Materials::get_sigma_a(const int grp, std::vector<double>& sig_a)
{
  for(int i=0; i < m_n_xs_quad ; i++)
  {
    m_mat_property_evals[i] = m_abs_opacities[m_current_material]->get_absorption_opacity(
      grp, m_t_at_xs_eval_points[i], m_xs_position[i]);      
  }
  
  m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, sig_a);
    
  return;
}

void Materials::clear_left_edge_set(void){ m_left_edge_set = false; return; }

void Materials::clear_right_edge_set(void){ m_right_edge_set = false; return; }
  
double Materials::get_left_sigma_a(const int grp)
{
  if(!m_left_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate sigma_a on left edge without left edge position and temperature being set");
    
  return m_abs_opacities[m_current_material]->get_absorption_opacity(grp, m_t_left_bound, m_x_left);
}

double Materials::get_right_sigma_a(const int grp)
{
  if(!m_right_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate sigma_a on right edge without right edge position and temperature being set");
    
  return m_abs_opacities[m_current_material]->get_absorption_opacity(grp, m_t_right_bound, m_x_right);
}

double Materials::get_left_sigma_s(const int grp,const int l_mom)
{
  if(!m_left_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate sigma_s on left edge without left edge position and temperature being set");
    
  return m_scat_opacities[m_current_material]->get_scattering_opacity( l_mom, grp, m_t_left_bound, m_x_left);
}

double Materials::get_right_sigma_s(const int grp, const int l_mom)
{
  if(!m_right_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate sigma_s on right edge without right edge position and temperature being set");
    
  return m_scat_opacities[m_current_material]->get_scattering_opacity( l_mom, grp, m_t_right_bound, m_x_right);
}

double Materials::get_left_cv(void)
{
  if(!m_left_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate CV on left edge without left edge position and temperature being set");
    
  return m_cv_obj[m_current_material]->get_cv(m_x_left, m_t_left_bound);
}

double Materials::get_right_cv(void)
{
  if(!m_right_edge_set)
    throw Dark_Arts_Exception(SUPPORT_OBJECT, "Attempting to evaluate CV on right edge without right edge position and temperature being set");
    
  return m_cv_obj[m_current_material]->get_cv(m_x_right, m_t_right_bound);
}

  /// calculate \f$ \sigma_s \f$ for all DFEM integration points, groups, and scattering moments for cell cell_num
void Materials::get_sigma_s(const int grp, const int l_mom, std::vector<double>& sig_s)
{  
  for(int i=0; i < m_n_xs_quad; i++)
  {  
    m_mat_property_evals[i] = m_scat_opacities[m_current_material]->get_scattering_opacity(
      l_mom, grp, m_t_at_xs_eval_points[i], m_xs_position[i]);
  }
      
  m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, sig_s);

  return;
}

/// calculate \f$ C_v \f$ for all DFEM integration points and groups for cell cell_num
void Materials::get_cv(std::vector<double>& cv)
{
  for(int i=0; i < m_n_xs_quad; i++)
  {
    m_mat_property_evals[i] = m_cv_obj[m_current_material]->get_cv(m_xs_position[i], m_t_at_xs_eval_points[i]);
  }  
  /// no special mapping necessary, just calculate material cv (at dfem integration points) directly (no group or legendre moment dependence)
  m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, cv);

  return;
}

/** evaluate temperature  driving source moments at a fine mesh number of points.
  Do not apply any type of opacity_treatment.  Exact integration works best.
  V_Matrix_Construction has called this function and will assemble the source moments
*/
void Materials::get_temperature_source(const double time, std::vector<double>& t_source_evals)
{
  for(int i=0; i < m_n_source_pts; i++)
  {  
    // std::cout << "x eval in materials: " << m_position_at_source_quad[i] << std::endl;
    t_source_evals[i] = m_source_t[m_current_material]->get_temperature_source(m_position_at_source_quad[i], time);
  }
  return;
}

/** evaluate intensity  driving source moments at a fine mesh number of points.
  Do not apply any type of opacity_treatment.  Exact integration works best.
  V_Matrix_Construction has called this function and will assemble the source moments
*/
void Materials::get_intensity_source(const double time, const int grp, const int dir, std::vector<double>& i_source_evals)
{
  for(int i=0; i < m_n_source_pts; i++)
    i_source_evals[i] = m_source_i[m_current_material]->get_intensity_source(m_position_at_source_quad[i], grp, dir, time);
  
  return;
}

void Materials::get_mf_planck(const Eigen::VectorXd& t_eval_vec, const int grp, Eigen::VectorXd& planck)
{
  for(int i=0 ; i < m_n_el_cell ; i++)
    planck(i) = m_planck.integrate_B(t_eval_vec(i), m_grp_e_min[grp], m_grp_e_max[grp] );

  return;
}

void Materials::get_grey_planck(const Eigen::VectorXd& t_eval_vec, Eigen::VectorXd& planck)
{
  for(int i=0 ; i < m_n_el_cell ; i++)
    planck(i) = m_planck.integrate_B_grey( t_eval_vec(i) );
    
  return;
}

double Materials::get_mf_planck(const double t_eval, const int grp)
{
  return m_planck.integrate_B(t_eval, m_grp_e_min[grp], m_grp_e_max[grp] );
}

double Materials::get_grey_planck(const double t_eval)
{
  return m_planck.integrate_B_grey(t_eval);
}

void Materials::get_mf_planck_derivative(const Eigen::VectorXd& t_eval_vec, const int grp, Eigen::MatrixXd& d_planck)
{
  for(int i=0 ; i < m_n_el_cell ; i++)
    d_planck(i,i) = m_planck.integrate_dBdT(t_eval_vec(i), m_grp_e_min[grp], m_grp_e_max[grp] );

  return;
}

void Materials::get_grey_planck_derivative(const Eigen::VectorXd& t_eval_vec, Eigen::MatrixXd& d_planck)
{
  for(int i=0 ; i < m_n_el_cell ; i++)
    d_planck(i,i) = m_planck.integrate_dBdT_grey( t_eval_vec(i) );

  return;
}

double Materials::get_c(void)
{
  return m_planck.get_c() ;
}

/* *****************************************************************
*   Private Member functions
*  ****************************************************************/

void Materials::load_materials(const Input_Reader& input_reader, const Angular_Quadrature& angular_quadrature)
{
  for(int mat_num=0; mat_num< m_num_materials ; mat_num++)
  {
    /// absorption opacity
    OPACITY_TYPE abs_op_type = input_reader.get_abs_opacity_type(mat_num);
    if(abs_op_type == CONSTANT_XS)
    {
      m_abs_opacities[mat_num] = std::make_shared<Absorption_Opacity_Constant>( input_reader, mat_num) ;
    }
    else if(abs_op_type == RATIONAL)
    {
      m_abs_opacities[mat_num] = std::make_shared<Absorption_Opacity_Rational>( input_reader, mat_num ) ;
    }
    else if(abs_op_type == TABLE_LOOKUP)
    {
      m_abs_opacities[mat_num] = std::make_shared<Absorption_Opacity_Table>( input_reader, mat_num );
    }
    else if(abs_op_type == POLYNOMIAL_SPACE)
    {
      m_abs_opacities[mat_num] = std::make_shared<Absorption_Opacity_Polynomial_Space>( input_reader, mat_num );
    }
    else
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Requesting to create an undefined VAbsorption_Opactiy object");
    
    /// scattering opacity
    OPACITY_TYPE scat_op_type = input_reader.get_scat_opacity_type(mat_num);
    if(scat_op_type == CONSTANT_XS)
    {
      m_scat_opacities[mat_num] = std::make_shared<Scattering_Opacity_Constant>( input_reader, mat_num) ;
    }
    else if(scat_op_type == RATIONAL)
    {
      m_scat_opacities[mat_num] = std::make_shared<Scattering_Opacity_Rational>( input_reader, mat_num) ;
    }
    else if(scat_op_type == TABLE_LOOKUP)
    {
      throw Dark_Arts_Exception(SUPPORT_OBJECT,"No table look-up scattering opacities.  In Matierals constructor");
    }
    else if(scat_op_type == POLYNOMIAL_SPACE)
    {
      m_scat_opacities[mat_num] = std::make_shared<Scattering_Opacity_Polynomial_Space>( input_reader, mat_num );
    }
    else
      throw Dark_Arts_Exception( SUPPORT_OBJECT , "Requesting to create an undefined VScattering_Opactiy object");
    
    /// cv
    CV_TYPE cv_type = input_reader.get_cv_type(mat_num);
    if(cv_type == CONSTANT_CV)
    {
      m_cv_obj[mat_num] = std::make_shared<Cv_Constant>( input_reader, mat_num)   ;
    }
    else if(cv_type == RATIONAL_CV)
    {
      m_cv_obj[mat_num] = std::make_shared<Cv_Rational>( input_reader, mat_num)   ;
    }
    else if(cv_type == POLYNOMIAL_CV)
    {
      m_cv_obj[mat_num] = std::make_shared<Cv_Polynomial_Temperature>( input_reader, mat_num );
    }
    else
      throw Dark_Arts_Exception(SUPPORT_OBJECT , "Requesting to initialize an undefined VCv object");
    
    /// radiation source
    FIXED_SOURCE_TYPE rad_src_type = input_reader.get_radiation_source_type(mat_num);
    switch(rad_src_type)
    {
      case NO_SOURCE:
      {
        m_source_i[mat_num] = std::make_shared<Source_I_None>( input_reader, mat_num)   ;
        break;
      }
      case MMS_SOURCE:
      {
        m_source_i[mat_num] = std::make_shared<Source_I_MMS>( input_reader , 
          angular_quadrature , m_abs_opacities , m_scat_opacities, mat_num , m_planck)   ;
        break;
      }
      case CONSTANT_SOURCE:
      {
        m_source_i[mat_num] = std::make_shared<Source_I_Constant>( input_reader, mat_num)   ;
        break;
      }
      case INVALID_FIXED_SOURCE_TYPE:
      {
        throw Dark_Arts_Exception(INPUT, "Trying to create an invalid radiation source for Materials object");
        break;      
      }
    }
    /// temperature source
    FIXED_SOURCE_TYPE temp_src_type = input_reader.get_temperature_source_type(mat_num);
    switch(temp_src_type)
    {
      case NO_SOURCE:
      {
      /// default is to set source to 0
        m_source_t[mat_num] = std::make_shared<Source_T_None>( input_reader, mat_num)   ;
        break;
      }
      case CONSTANT_SOURCE:
      {
      /// default is to set source to 0
        m_source_t[mat_num] = std::make_shared<Source_T_Constant>( input_reader, mat_num)   ;
        break;
      }
      case MMS_SOURCE:
      {
        m_source_t[mat_num] = std::make_shared<Source_T_MMS>( input_reader , angular_quadrature,
          m_abs_opacities , m_cv_obj , mat_num , m_planck) ;
        break;
      }
      case INVALID_FIXED_SOURCE_TYPE:
      {
        throw Dark_Arts_Exception (SUPPORT_OBJECT,  "Attempting to create an invalid temperature source in Materials");
        break;
      }
    }
  }
}

double Materials::get_cell_width(void)
{
  return m_dx;
}

