#include "Materials.h"

Materials::Materials( const Input_Reader& input_reader)
{
  /// Loop over the number of materials, load each object type as a appropriate
  m_num_materials = input_reader.get_number_of_materials();
  
  m_abs_opacities.resize(m_num_materials);
  m_scat_opacities.resize(m_num_materials);

  m_cv.resize(m_num_materials);
  m_source_t.resize(m_num_materials);
  m_source_i.resize(m_num_materials);  
}

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
      m_cv[mat_num] = std::shared_ptr<VCv> (new Cv_Constant( input_reader, mat_num)  ) ;
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