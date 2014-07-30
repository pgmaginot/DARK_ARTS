#include "Materials.h"

Materials::Materials( const Input_Reader& input_reader)
{
  /// Loop over the number of materials, load each object type as a appropriate
  m_num_materials = input_reader.get_number_of_materials();
  for(int mat_num=0; mat_num< m_num_materials ; mat_num++)
  {
    /// absorption opacity
    OPACITY_TYPE abs_op_type = input_reader.get_abs_opacity_type(mat_num);
    if(abs_op_type == CONSTANT_XS)
    {
    
    }
    else if(abs_op_type == RATIONAL)
    {
    
    }
    else if(abs_op_type == TABLE_LOOKUP)
    {
    
    }
    
    /// scattering opacity
    OPACITY_TYPE scat_op_type = input_reader.get_scat_opacity_type(mat_num);
    if(scat_op_type == CONSTANT_XS)
    {
    
    }
    else if(scat_op_type == RATIONAL)
    {
    
    }
    else if(scat_op_type == TABLE_LOOKUP)
    {
    
    }
    
    /// cv
    CV_TYPE cv_type = input_reader.get_cv_type(mat_num);
    if(cv_type == CONSTANT_CV)
    {
    
    }
    
    /// radiation source
    FIXED_SOURCE_TYPE rad_src_type = input_reader.get_radiation_source_type(mat_num);
    if(rad_src_type == NO_SOURCE)
    {
    
    }
    
    /// temperature source
    FIXED_SOURCE_TYPE temp_src_type = input_reader.get_temperature_source_type(mat_num);
    if(temp_src_type == NO_SOURCE)
    {
    
    }
  }
}