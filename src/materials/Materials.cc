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

Materials::Materials( const Input_Reader& input_reader, const Fem_Quadrature& fem_quadrature, Cell_Data* cell_ptr)
{
  /// Loop over the number of materials, load each object type as a appropriate
  m_num_materials = input_reader.get_number_of_materials();
  
  m_abs_opacities.resize(m_num_materials);
  m_scat_opacities.resize(m_num_materials);

  m_cv_obj.resize(m_num_materials);
  m_source_t.resize(m_num_materials);
  m_source_i.resize(m_num_materials);  
  
  load_materials(input_reader);
  
  /// create a smart pointers that looks at input once, then select cross section generating functions
  switch( input_reader.get_opacity_treatment() )
  {
    case SLXS:
    {
      m_xs_treatment = std::shared_ptr<V_XS_Treatment>( new XS_Treatment_SLXS(fem_quadrature) );
      break;
    }
    case INTERPOLATING:
    {
      m_xs_treatment =  std::shared_ptr<V_XS_Treatment>( new XS_Treatment_Interpolating(fem_quadrature) );
      break;
    }
    case MOMENT_PRESERVING:
    {
      m_xs_treatment =  std::shared_ptr<V_XS_Treatment>( new XS_Treatment_Moment_Preserving(fem_quadrature , input_reader) );
      break;
    }
    case INVALID_OPACITY_TREATMENT:
    {
      std::cerr << "Invalid OPacity treament still exists in Materials initializations\n" ;
      exit(EXIT_FAILURE);
      break;
    }    
  }
  
  fem_quadrature.get_xs_eval_points(m_xs_eval_quad);
  m_n_xs_quad = m_xs_eval_quad.size();
  fem_quadrature.get_dfem_at_xs_eval_points(m_dfem_at_xs);
  
  m_n_el_cell = fem_quadrature.get_number_of_interpolation_points();
  
  fem_quadrature.get_dfem_at_edges(m_dfem_at_left_edge,m_dfem_at_right_edge);
  
  m_xs_position.resize(m_n_xs_quad,0.);
  
  /// data used to access arrays
  m_n_groups = input_reader.get_number_of_groups();
  m_n_l_mom = input_reader.get_number_of_legendre_moments();
  m_n_dfem_integration_points = fem_quadrature.get_number_of_integration_points();
  
  m_sig_a.resize(m_n_groups*m_n_dfem_integration_points,0.);
  m_sig_s.resize(m_n_groups*m_n_dfem_integration_points*m_n_l_mom,0.);;
  m_cv.resize(m_n_dfem_integration_points,0.);
  m_sig_a_boundary.resize(2*m_n_groups,0.);
  m_sig_s_boundary.resize(2*m_n_groups*m_n_l_mom,0.);
  m_cv_boundary.resize(2,0.);
  m_mat_property_evals.resize(m_n_xs_quad,0.);
  m_t_at_xs_eval_points.resize(m_n_xs_quad,0.);
  
  cell_data_ptr = cell_ptr;
}

void Materials::calculate_local_temp_and_position(const int cell_num, const Temperature_Data& temperature)
{
  /// determine what material we are in
  m_current_material = cell_data_ptr->get_cell_material_number(cell_num);
  
  /// calculate the local position at the xs evalaution points
  double xL = cell_data_ptr->get_cell_width(cell_num);
  double dx = cell_data_ptr->get_cell_left_edge(cell_num);
  
  for(int i=0; i< m_n_xs_quad; i++)
  {
    m_xs_position[i] = xL + dx/2.*(1. + m_xs_eval_quad[i]);
  }
  
  m_x_left = xL;
  m_x_right = xL + dx;
  
  /// calculate the local temperature at the xs evaluation points
  /// basis function evaluation is laid out in a vector: B_{1,1} ... B_{1,N_eval} B_{2,1} ...
  for(int i=0; i < m_n_xs_quad; i++)
    m_t_at_xs_eval_points[i] = 0.;
    
  m_t_left_bound = 0.; 
  m_t_right_bound = 0.;
  
  
  int p=0;
  for(int dfem_t = 0; dfem_t < m_n_el_cell ; dfem_t++)
  {
    double t_el = temperature.get_temperature(dfem_t, cell_num);
    m_t_left_bound += t_el*m_dfem_at_left_edge[dfem_t];
    m_t_right_bound += t_el*m_dfem_at_right_edge[dfem_t];
    for(int i=0; i< m_n_xs_quad; i++)
    {
      m_t_at_xs_eval_points[i] += t_el*m_dfem_at_xs[p];
      p++;
    }
  }
  
  return;
}
   
void Materials::update_sigma_a(void)
{
  int cnt = 0;
  for(int g=0; g<m_n_groups ; g++)
  {
    for(int i=0; i < m_n_xs_quad ; i++)
    {
      m_mat_property_evals[cnt] = m_abs_opacities[m_current_material]->get_absorption_opacity(
        g, m_t_at_xs_eval_points[i], m_xs_position[i]);
      cnt++;
    }
    
    m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, m_mat_mapped);
    /** group g opacities have now been calculated at dfem integration points, load into vector that
      is laid out at g=0 to g=m_n_groups { j=0 ... N_P DFEM integration points }
    */
    for(int j=0; j<m_n_dfem_integration_points;j++)
      m_sig_a[g*m_n_dfem_integration_points + j] = m_mat_mapped[j];
    
    /// calculate edge values
    m_sig_a_boundary[2*g] = m_abs_opacities[m_current_material]->get_absorption_opacity(
        g, m_t_left_bound, m_x_left);
    m_sig_a_boundary[2*g+1] = m_abs_opacities[m_current_material]->get_absorption_opacity(
        g, m_t_right_bound, m_x_right);
  }

  return;
}
  /// calculate \f$ \sigma_s \f$ for all DFEM integration points, groups, and scattering moments for cell cell_num
void Materials::update_sigma_s(void)
{
  // for(int l=0; l< m_n_leg_mom; l++)
  // {
    // for(int g=0; g<m_n_groups ; g++)
    // {
      // for(int i=0; i < m_n_xs_quad i++)
      // {  
        // m_mat_property_evals[cnt] = m_scat_opacities[m_current_material]->get_scat_opacity(
          // l, g, m_t_at_xs_eval_points[i], m_xs_position[i]);
      // }
      
      // m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, m_mat_mapped);
      // /** moment l, group g opacities have now been calculated at dfem integration points, load into vector that
        // is laid out at l=0 to l = m_n_leg_mom {g=0 to g=m_n_groups { j=0 ... N_P DFEM integration points }}
      // */
      // for(int j=0; j<m_n_integration_points;j++)
        // m_sig_s[l* m_n_groups*m_n_integration_points + g*m_n_integration_points + j] = m_mat_mapped[j];
    // }
  // }

  return;
}

/// calculate \f$ C_v \f$ for all DFEM integration points and groups for cell cell_num
void Materials::update_cv(void)
{
  // for(int i=0; i < m_n_xs_quad i++)
  // {
    // m_mat_property_evals[cnt] = m_cv_obj[m_current_material]->get_cv(m_xs_position[i], m_t_at_xs_eval_points[i]);
  // }  
  // m_xs_treatment->calculate_xs_at_integration_points(m_mat_property_evals, m_cv);

  return;
}

// void Materials::get_rx_sigma_a(const double dx, const double xL, const int mat_num, const int grp, 
  // const std::vector<double>& temperature, std::vector<double>& sig_a)
// {
  // m_sig_a_wiz->query_sig_a(dx, xL, mat_num, grp, temperature, sig_a);
  // return;
// }  
// void Materials::get_rx_sigma_s(const double dx, const double xL, const int mat_num, const int grp, 
  // const int l_mom, const std::vector<double>& temperature, std::vector<double>& sig_s)
// {
  // return;
// }
// void Materials::get_rx_cv(const double dx, const double xL, const int mat_num, 
  // const std::vector<double>& temperature, std::vector<double>& cv)
// {
  // return;
// }

/* *****************************************************************
*   Private Member functions
*  ****************************************************************/

void Materials::load_materials(const Input_Reader& input_reader)
{
  for(int mat_num=0; mat_num< m_num_materials ; mat_num++)
  {
    /// absorption opacity
    OPACITY_TYPE abs_op_type = input_reader.get_abs_opacity_type(mat_num);
    if(abs_op_type == CONSTANT_XS)
    {
      m_abs_opacities[mat_num] = std::shared_ptr<VAbsorption_Opacity> 
        (new Absorption_Opacity_Constant( input_reader, mat_num)  ) ;
    }
    else if(abs_op_type == RATIONAL)
    {
      m_abs_opacities[mat_num] = std::shared_ptr<VAbsorption_Opacity> 
        (new Absorption_Opacity_Rational( input_reader, mat_num ) );
    }
    else if(abs_op_type == TABLE_LOOKUP)
    {
      m_abs_opacities[mat_num] = std::shared_ptr<VAbsorption_Opacity> 
        (new Absorption_Opacity_Table( input_reader, mat_num ));
    }
    
    /// scattering opacity
    OPACITY_TYPE scat_op_type = input_reader.get_scat_opacity_type(mat_num);
    if(scat_op_type == CONSTANT_XS)
    {
      m_scat_opacities[mat_num] = std::shared_ptr<VScattering_Opacity> 
        (new Scattering_Opacity_Constant( input_reader, mat_num)  ) ;
    }
    else if(scat_op_type == RATIONAL)
    {
      m_scat_opacities[mat_num] = std::shared_ptr<VScattering_Opacity>
        (new Scattering_Opacity_Rational( input_reader, mat_num)  ) ;
    }
    else if(scat_op_type == TABLE_LOOKUP)
    {
      m_scat_opacities[mat_num] = std::shared_ptr<VScattering_Opacity>
        (new Scattering_Opacity_Table( input_reader, mat_num)  ) ;
    }
    
    /// cv
    CV_TYPE cv_type = input_reader.get_cv_type(mat_num);
    if(cv_type == CONSTANT_CV)
    {
      m_cv_obj[mat_num] = std::shared_ptr<VCv> (new Cv_Constant( input_reader, mat_num)  ) ;
    }
    
    /// radiation source
    FIXED_SOURCE_TYPE rad_src_type = input_reader.get_radiation_source_type(mat_num);
    if(rad_src_type == NO_SOURCE)
    {
      m_source_i[mat_num] = std::shared_ptr<VSource_I> (new Source_I_Constant( input_reader, mat_num)  ) ;
    }
    
    /// temperature source
    FIXED_SOURCE_TYPE temp_src_type = input_reader.get_temperature_source_type(mat_num);
    if(temp_src_type == NO_SOURCE)
    {
      /// default is to set source to 0
      m_source_t[mat_num] = std::shared_ptr<VSource_T> (new Source_T_Constant( input_reader, mat_num)  ) ;
    }
  }
}

