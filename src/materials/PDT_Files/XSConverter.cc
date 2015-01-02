// Copyright (c) 2000-2008, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.

// If written permission is obtained for redistribution or further use, the
// following conditions must be met:

// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.

// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.

// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


// =============================================================================
//
// 01/18/2013 hayes
//
// Added converter from ISOTXS text files to PDT format
//
  /*!
    \file    Converter.cc
    \author  Michael Adams, Aaron Roney, Bruno Turcksin
    \date    June 2010

    \brief
     This program converts data files of several formats to the standard PDT
   data format.  Tags are "format", "file", and "pdtfile".  "format" is the
   format of the source data; "file" is the source file; and "pdtfile" is the
   name of the file to be written.
     For neutron and gamme data, we also offer the capability to split and/or
   coalesce energy groups.  For example, all thermal groups can be lumped
   together, or a sliver group can be created to capture spectral lines.
  */
// =============================================================================


#include "XSConverter.h"
#include <cstdlib>
#include <iostream>
#include <string>

using std::string;
using std::fstream;
using std::stringstream;
using std::endl;
using std::cin;
using std::cout;
using std::ofstream;
using std::map;
using std::exception;
using std::make_pair;
using std::noskipws;
using std::skipws;
using std::vector;
using std::ifstream;
using std::setw;
using std::setprecision;
using std::pair;

// =============================================================================
//    Converter Main
// =============================================================================
int main(int argc, char** argv)
{
  string format;
  string file;
  string pdtfile;
  bool b_arg4 = false;

  // Converter input file fields.
  bool pdtConverterInput = false;
  vector<float> clydes;
  vector<command> manipulationCommands;

  if(argc == 2) // Using converter input file.
  {
    pdtConverterInput = true;
    ifstream fin;

    try
    {
      char dummy[1000];
      int streamsize;

      fin.open(argv[1]);
      streamsize = 1000;

      // Get line 1 arguments (format + ',' + file + ',' + pdtfile + '\n')
      fin.getline(dummy, streamsize, ',');
      format = dummy;
      fin.getline(dummy, streamsize, ',');
      file = dummy;
      fin.getline(dummy, streamsize);
      pdtfile = dummy;

      // Get line 2 arguments (clyde1,clyde2,clyde3,etc.)
      string clydeLine;
      fin.getline(dummy, streamsize);
      clydeLine = dummy;
      stringstream clydeLineStream;
      clydeLineStream << clydeLine;

      while(clydeLineStream.getline(dummy, streamsize, ','))
      {
        stringstream clydeStream;
        float clyde;
        clydeStream << dummy;
        clydeStream >> clyde;
        clydes.push_back(clyde);
      }

      // Get line 2 arguments (command1,command2,command3,etc.)
      string commandLine;
      fin.getline(dummy, streamsize);
      commandLine = dummy;
      stringstream commandLineStream;
      commandLineStream << commandLine;

      while(commandLineStream.getline(dummy, streamsize, ','))
      {
        stringstream commandStream;
        commandStream << dummy;
        command cmd;
        commandStream >> cmd.commandType;
        commandStream >> cmd.groupType;
        commandStream >> cmd.operator1;
        commandStream >> cmd.operator2;
        manipulationCommands.push_back(cmd);
      }
    }
    catch(exception& e)
    {
      std::cout << std::endl << std::endl
      << "***ERROR***" << std::endl
      << "Standard exception caught; reported message follows:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Terminating program." << std::endl;
      return 1;
    }
  }
  else if(argc == 4 || argc == 5) // Using normal operation.
  {
    format = argv[1];
    file = argv[2];
    pdtfile = argv[3];
    b_arg4 = false;

    if( argv[4] )
      b_arg4 = true;
  }
  else // Incorrect number.
  {
    stringstream str;
    str << "ERROR:  XSConverter requires three arguments:  source format, "
      << "source file, and output file." << endl;
    str << "ERROR:  XSConverter requires one argument: converter file." << endl;
    throw CommonException( str, CET_INPUT_ERROR );
  }

  try
  {
    if ( format == "pdt" )
    {
      if( b_arg4 )
        std::cerr << "*** WARNING ***  Fourth argument not needed for PDT "
          << "format. " << std::endl;
      Converter converter;
      converter.info.sourcefile = file;
      converter.info.format = MF_PDT;
      converter.read_pdt(file);
      converter.write_pdt(pdtfile);
    }
    else if( format == "crash" )
    {
      double opac_max;
      std::cout << "Enter maximum opacity value in cm^2/g:  ";
      std::cin >> opac_max;
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_opac();
      converter.info.type = MDT_THERMAL;
      converter.info.units = MDU_cmspg;
      converter.info.format = MF_CRASH;
      converter.info.microscopic = true;
      converter.read_crash(file, opac_max);
      converter.write_pdt(pdtfile);
    }
    else if( format == "tops" )
    {
      double e_max, opac_max;
      if( ! b_arg4 ) {
        std::cerr << "*** WARNING ***  Fourth argument needed for TOPS format."
          << endl << "Enter highest energy group boundary in eV:  ";
        std::cin >> e_max;
      }
      else
        e_max = atof( argv[4] );
      std::cout << "Enter maximum opacity value in cm^2/g:  ";
      std::cin >> opac_max;
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_opac();
      converter.info.type = MDT_THERMAL;
      converter.info.units = MDU_cmspg;
      converter.info.format = MF_TOPS;
      converter.info.microscopic = true;
      converter.read_tops(file, opac_max);
      converter.fix_energy(e_max);
      converter.reverse_opacities();
      converter.write_pdt(pdtfile);
    }
    else if( format == "sesame" )
    {
      if( b_arg4 )
        std::cerr << "*** WARNING ***  Fourth argument not needed for SESAME "
          << "format." << std::endl;
      double opac_max;
      std::cout << "Enter maximum opacity value in cm^2/g:  ";
      std::cin >> opac_max;
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_opac();
      converter.info.type = MDT_THERMAL;
      converter.info.units = MDU_cmspg;
      converter.info.format = MF_SESAME;
      converter.info.multigroup = false;
      converter.info.microscopic = true;
      converter.read_sesame(file, opac_max);
      converter.reverse_opacities();
      converter.write_pdt(pdtfile);
    }
    else if( format == "csv" )
    {
      if( b_arg4 )
        std::cerr << "*** WARNING ***  Fourth argument not needed for CSV "
          << "format." << std::endl;
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_opac();
      converter.info.type = MDT_THERMAL;
      converter.info.units = MDU_cmspg;
      converter.info.format = MF_CSV;
      converter.info.multigroup = false;
      converter.info.microscopic = true;
      converter.read_csv(file);
      converter.write_pdt(pdtfile);
    }
    else if( format == "cx_test" )
    {
      if( b_arg4 )
        std::cerr << "*** WARNING ***  Fourth argument not needed for cx_test "
          << "format." << std::endl;
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_cx();
      converter.info.type = MDT_NEUTRON;
      converter.info.units = MDU_cms;
      converter.info.format = MF_CX_TEST;
      converter.info.microscopic = true;
      converter.read_test(file);
      converter.write_pdt(pdtfile);
    }
    else if( format == "ampx_neutron" || format == "ampx_coupled")
    {
      if(pdtConverterInput)
      {
        vector<float>::iterator clyde = clydes.begin();
        for( ; clyde != clydes.end(); ++clyde)
        {
          Converter converter;
          converter.info.sourcefile = file;

          // Clyde does not need assignment due to iterator.
          std::cout << "Auto selection of clyde " << *clyde << "." << std::endl;
          cout << "clyde = " << *clyde << endl;
          if( format == "ampx_neutron" )
          {
            converter.set_ptrs_cx();
            converter.info.type = MDT_NEUTRON;
            converter.info.units = MDU_barns;
            converter.info.format = MF_AMPX_COUPLED;
            converter.info.microscopic = true;
            converter.info.format = MF_AMPX_NEUTRON;
            converter.read_ampx_neutron(file, *clyde);
          }
          if( format == "ampx_coupled" )
          {
            converter.set_ptrs_cx();
            converter.info.type = MDT_COUPLED_NG;
            converter.info.units = MDU_barns;
            converter.info.format = MF_AMPX_COUPLED;
            converter.info.microscopic = true;
            converter.read_ampx_coupled(file, *clyde);
          }

          // Group splitting and coalescing option.
          vector<command>::iterator cmd = manipulationCommands.begin();
          for( ; cmd != manipulationCommands.end(); ++cmd)
          {
            if((*cmd).commandType == "SPLIT")
            {
              stringstream ss;
              ss << (*cmd).commandType << " " << (*cmd).groupType << " "
                 << (*cmd).operator1 << " " << (*cmd).operator2;
              string commandString;
              std::getline(ss, commandString);
              converter.group_manipulation(commandString);

              vector<command>::iterator adjusterator = cmd + 1;
              for( ; adjusterator != manipulationCommands.end(); ++adjusterator)
              {
                if((*adjusterator).operator1 > (*cmd).operator1)
                {
                  if((*adjusterator).commandType == "SPLIT")
                    (*adjusterator).operator1++;
                  else if((*adjusterator).commandType == "COMB")
                  {
                    (*adjusterator).operator1++;
                    (*adjusterator).operator2++;
                  }
                }
              }
            }
            else if((*cmd).commandType == "COMB")
            {
              int start = (*cmd).operator1;
              int end = (*cmd).operator2;
              int count = 0;

              while(start < end)
              {
                stringstream ss;
                ss << (*cmd).commandType << " " << (*cmd).groupType << " "
                   << start << " " << (start + 1);
                string commandString;
                std::getline(ss, commandString);
                converter.group_manipulation(commandString);

                end--;
                count++;
              }

              vector<command>::iterator adjusterator = cmd + 1;
              for( ; adjusterator != manipulationCommands.end(); ++adjusterator)
              {
                if((*adjusterator).operator1 > (*cmd).operator2)
                {
                  if((*adjusterator).commandType == "SPLIT")
                    (*adjusterator).operator1 -= count;
                  else if((*adjusterator).commandType == "COMB")
                  {
                    (*adjusterator).operator1 -= count;
                    (*adjusterator).operator2 -= count;
                  }
                }
              }
            }
          }

          stringstream ss;
          ss << pdtfile << "_" << *clyde;
          string filename;
          ss >> filename;

          converter.write_pdt(ss.str());
        }
      }
      else
      {
        Converter converter;
        converter.info.sourcefile = file;

        float clyde;
        if( ! b_arg4 ) {
          std::cerr << "*** WARNING ***  Fourth argument needed for AMPX "
            << "format." << endl << "Enter Clyde number of desired material: ";
         std::cin >> clyde;
        }
        else
          clyde = (float) atof( argv[4] );

        if(format == "ampx_neutron")
        {
          converter.set_ptrs_cx();
          converter.info.type = MDT_NEUTRON;
          converter.info.units = MDU_barns;
          converter.info.format = MF_AMPX_NEUTRON;
          converter.info.microscopic = true;
          converter.read_ampx_neutron(file, clyde);
        }
        else if(format == "ampx_coupled")
        {
          converter.set_ptrs_cx();
          converter.info.type = MDT_COUPLED_NG;
          converter.info.units = MDU_barns;
          converter.info.format = MF_AMPX_COUPLED;
          converter.info.microscopic = true;
          converter.read_ampx_coupled(file, clyde);
        }

        // Group splitting and coalescing option.
        char response = 'y';
        while(response == 'y')
        {
          std::cout << std::endl << "Do you want to run the group splitting "
          << "and coalescing module (y/n)?: ";
          cin >> response;
          if(response == 'y')
            converter.group_manipulation();
        }

        converter.write_pdt(pdtfile);
      }
    }
    else if( format == "matxs" )
    {
      //Whether to combine the xs into one
      // condense == 0: no condensation
      // condense >  0: condense and only output that matrix
      // condense <  0: condense and output all transfer matrices
      int opt_combine;
      if( ! b_arg4 )
  opt_combine = 1;
      else
  opt_combine = (int) atof( argv[4] );

      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_cx();
      converter.info.units = MDU_barns;
      converter.info.format = MF_MATXS;
      converter.info.microscopic = true;
      converter.read_matxs(file, pdtfile, opt_combine);
    }
    else if( format == "isotxs" )
    {
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_cx();
      converter.info.units = MDU_barns;
      converter.info.format = MF_ISOTXS;
      converter.info.type = MDT_NEUTRON;
      converter.info.microscopic = true;
      converter.read_isotxs(file, pdtfile);
    }
    else if( format == "dragon" )
    {
      Converter converter;
      converter.info.sourcefile = file;
      converter.set_ptrs_cx();
      converter.info.type = MDT_NEUTRON;
      converter.info.units = MDU_icm;
      converter.info.format = MF_DRAGON;
      converter.info.microscopic = false;
      converter.read_dragon(file, pdtfile);
    }
    else if (format == "cepxs")
    {
        if (b_arg4)
            std::cerr << "*** WARNING ***  Fourth argument not needed for "
              << "cx_test format." << std::endl;
        Converter converter;
        converter.info.sourcefile = file;
        converter.set_ptrs_cx();
        converter.info.type = MDT_ELECTRON;
        converter.info.units = MDU_icm;
        converter.info.format = MF_CEPXS;
        converter.info.microscopic = false;
        converter.read_cepxs(file);
        converter.write_pdt(pdtfile);
    }
    else if (format == "cepxs-bfp")
    {
        if (b_arg4)
            std::cerr << "*** WARNING *** Fourth argument not needed for "
              << "cx_test format." << std::endl;
        Converter converter;
        converter.info.sourcefile = file;
        converter.set_ptrs_cx();
        converter.info.type = MDT_ELECTRON;
        converter.info.units = MDU_icm;
        converter.info.format = MF_CEPXS_BFP;
        converter.info.microscopic = false;
        converter.read_cepxs_bfp(file);
        converter.write_pdt(pdtfile);
    }
    else
    {
      std::cout << "Invalid file format specified.  Arguments are "
                << "format, file to read, and file to write." << endl
                << "Format options currently are tops, sesame, csv, cx_test, "
                << "ampx_neutron, ampx_coupled, cepxs, cepxs-bfp, and isotxs" << endl
                << "You specified " << format << ", " << file << ", and "
                << pdtfile << '.' << endl;
      return 1;
    }
  }
  catch( CommonException& ce )
  {
    ce.print();
    std::cout << "Terminating program." << std::endl;
    exit(-1);
  }
  catch( std::exception& e )
  {
    std::cout << std::endl << std::endl
      << "ERROR:" << std::endl
      << "Standard exception caught; reported message follows:" << std::endl;
    std::cout << e.what() << std::endl;
    std::cout << "Terminating program." << std::endl;
    exit(-1);
  }

  return 0;

}
