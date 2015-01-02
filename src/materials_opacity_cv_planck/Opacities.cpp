#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <cmath>

#include <map>
#include <set>
#ifndef Data_Transfer_Object_h
#define Data_Transfer_Object_h
#include "Data_Transfer_Object.h"
#endif
#ifndef Constants_h
#define Constants_h
#include "Constants.h"
#endif
#ifndef Math_Utils_h
#define Math_Utils_h
#include "Math_Utils.h"
#endif
#ifndef Common_Utils_h
#define Common_Utils_h
#include "Common_Utils.h"
#endif
#ifndef Opacities_h
#define Opacities_h
#include "Opacities.h"
#endif
#ifndef Opacity_File_Details_h
#define Opacity_File_Details_h
#include "Opacity_File_Details.h"
#endif
#include <stdlib.h>

using std::cout;
using std::stringstream;
using std::vector;


void Opacities::initialize(string filename){
	//read_file(filename);
}

Opacity_File_Details Opacities::getEnergyGrpDetails(string filename,Opacity_File_Details ret_Opacity_File_Details)
{

	set <double> grp_set;

	int num_temps;
	int num_dens;
	int num_groups;

	ifstream fin( filename.c_str() );
	if( !fin.is_open() ) {
				
	}
	double conversion_factor;
	char line[250];
	string word;
	string sType="Macroscopic";
	double dDensity;
	fin.getline(line, 250);
	fin.getline(line, 250);

	// "This file is a " multigroup / single-group
	fin >> word >> word >> word >> word;
	string multi;
	fin >> multi;
	
	fin.getline(line, 250);
	fin >> num_temps >> word >> num_dens >> word >> word >> num_groups;
		
	ret_Opacity_File_Details.iNoOfDensities=num_dens;
	ret_Opacity_File_Details.iNoOfEnergyGrps=num_groups;
	ret_Opacity_File_Details.iNoOfTempEntries=num_temps;

	if(num_temps==1)
	{
		//It is a constant opacity problem.
		//Opacities dont depend on temperature.
		ret_Opacity_File_Details.isConstantOpacity='y';		
	}
	else
		ret_Opacity_File_Details.isConstantOpacity='n';
	

	vector<double> temps( num_temps,0.0 );
	set<double> temps_set;	
		
	vector<double> densities( num_dens,0.0 );
	vector<double> group_structure( num_groups+1,0.0 );
	vector<string> groups_mat(num_groups,"");
	
	//vector<vector<double> > sii_centers_curr(num_temps,num_groups);
        vector<vector<double> > sii_centers_curr(num_temps, vector<double>(num_groups, 0.0));
	
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin >> sType;
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int t=0; t!=num_temps; ++t )
	{
		fin >> temps[t];
		temps_set.insert(temps[t]);
	}
	//get Temperatures
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int d=0; d!=num_dens; ++d )
		fin >> densities[d];
	
	if(sType.compare("Macroscopic")==0){
		//Macroscopic

	}else{
		//Microscopic

	}
	
	stringstream strStrm;
	string tempStr="";
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int e=0; e!=num_groups+1; ++e )
	{
		fin >> group_structure[e];
		grp_set.insert(group_structure[e]);
	}

	int iLoop=0;
	set<double>::iterator it;
	vector<double> grp_vector( num_groups+1,0.0 );
	for (it=grp_set.begin(); it!=grp_set.end(); it++)
	{
		grp_vector[iLoop]=*it;
		iLoop++;
	}

	ret_Opacity_File_Details.grp_vector=grp_vector;
	ret_Opacity_File_Details.density_mat=densities;

	return ret_Opacity_File_Details;
	
}
//Opacity_Object Opacities::getOpacityData1(map<string,Opacity_Object> opacitiesMap,vector<double> temperature_mat,string strEnergyGrp,double dConvertedTemp,double dDensity)
//{
//	
//	int num_temps=temperature_mat.size();
//	stringstream strStmOpacKey;	
//	strStmOpacKey << strEnergyGrp<<":"<<dConvertedTemp<<":"<<dDensity;
//	string strOpacitykey=strStmOpacKey.str();
//	Opacity_Object opObj;
//
//	if(temperature_mat[0]==-1)
//	{
//		stringstream strStmOpacKey1;	
//		strStmOpacKey1 << strEnergyGrp<<":"<<"-1"<<":"<<dDensity;
//		string strOpacitykey1=strStmOpacKey1.str();
//		return (opacitiesMap.find(strOpacitykey1))->second;
//	}
//
//	if((opacitiesMap.find(strOpacitykey)) != opacitiesMap.end())
//		return (opacitiesMap.find(strOpacitykey))->second;
//	else 
//	{
//				
//		if(dConvertedTemp > temperature_mat[num_temps-1])
//		{
//			stringstream strStmOpacKeyRight;	
//			strStmOpacKeyRight << strEnergyGrp<<":"<<temperature_mat[num_temps-1]<<":"<<dDensity;
//			string strOpacitykeyRight=strStmOpacKeyRight.str();
//			Opacity_Object opDefaultRight;
//			opDefaultRight=(opacitiesMap.find(strOpacitykeyRight))->second;
//
//			
//			/*double dFracChange=((dConvertedTemp)/temps1[num_temps-1]);
//			double dSigA=dFracChange*opDefaultRight.dSigmaA;
//			double dSigT=dFracChange*opDefaultRight.dSigmaT;
//			double dSigS=dFracChange*opDefaultRight.dSigmaS;
//			opObj.dDensity=dDensity;
//			opObj.dSigmaA=dSigA;
//			opObj.dSigmaS=dSigS;
//			opObj.dSigmaT=dSigT;
//			opObj.dTemperature=dConvertedTemp;*/
//
//
//			opObj.dDensity=opDefaultRight.dDensity;
//			opObj.dSigmaA=opDefaultRight.dSigmaA;
//			opObj.dSigmaS=opDefaultRight.dSigmaS;
//			opObj.dSigmaT=opDefaultRight.dSigmaT;
//			opObj.dTemperature=dConvertedTemp;
//
//
//		}
//		else if(dConvertedTemp < temperature_mat[0])
//		{
//			stringstream strStmOpacKeyLeft;	
//			strStmOpacKeyLeft << strEnergyGrp<<":"<<temperature_mat[0]<<":"<<dDensity;
//			string strOpacitykeyLeft=strStmOpacKeyLeft.str();
//			Opacity_Object opDefaultLeft;
//			opDefaultLeft=(opacitiesMap.find(strOpacitykeyLeft))->second;
//
//			/*double dFracChange=((dConvertedTemp)/temps1[0]);
//			double dSigA=dFracChange*opDefaultLeft.dSigmaA;
//			double dSigT=dFracChange*opDefaultLeft.dSigmaT;
//			double dSigS=dFracChange*opDefaultLeft.dSigmaS;
//			opObj.dDensity=dDensity;
//			opObj.dSigmaA=dSigA;
//			opObj.dSigmaS=dSigS;
//			opObj.dSigmaT=dSigT;
//			opObj.dTemperature=dConvertedTemp;*/
//
//			opObj.dDensity=opDefaultLeft.dDensity;
//			opObj.dSigmaA=opDefaultLeft.dSigmaA;
//			opObj.dSigmaS=opDefaultLeft.dSigmaS;
//			opObj.dSigmaT=opDefaultLeft.dSigmaT;
//			opObj.dTemperature=dConvertedTemp;
//		}
//		else
//		{	
//			
//
//			for(int iCntLoop=1;iCntLoop < num_temps;iCntLoop++)
//			{
//				if((temperature_mat[iCntLoop-1] < dConvertedTemp)&&(temperature_mat[iCntLoop] > dConvertedTemp))
//				{
//					stringstream strStmOpacKeyLeft;	
//					strStmOpacKeyLeft << strEnergyGrp<<":"<<temperature_mat[iCntLoop-1]<<":"<<dDensity;
//					string strOpacitykeyLeft=strStmOpacKeyLeft.str();
//					stringstream strStmOpacKeyRight;	
//					strStmOpacKeyRight << strEnergyGrp<<":"<<temperature_mat[iCntLoop]<<":"<<dDensity;
//					string strOpacitykeyRight=strStmOpacKeyRight.str();
//					Opacity_Object opDefaultLeft;
//					
//					opDefaultLeft=(opacitiesMap.find(strOpacitykeyLeft))->second;
//					Opacity_Object opDefaultRight;
//					opDefaultRight=(opacitiesMap.find(strOpacitykeyRight))->second;
//
//					double dFracChange=((dConvertedTemp-temperature_mat[iCntLoop-1])/(temperature_mat[iCntLoop]-temperature_mat[iCntLoop-1]));
//					double dSigA=opDefaultLeft.dSigmaA+(dFracChange*(opDefaultRight.dSigmaA-opDefaultLeft.dSigmaA));						
//					double dSigT=opDefaultLeft.dSigmaT+(dFracChange*(opDefaultRight.dSigmaT-opDefaultLeft.dSigmaT));
//					double dSigS=opDefaultLeft.dSigmaS+(dFracChange*(opDefaultRight.dSigmaS-opDefaultLeft.dSigmaS));
//					opObj.dDensity=dDensity;
//					opObj.dSigmaA=dSigA;
//					opObj.dSigmaS=dSigS;
//					opObj.dSigmaT=dSigT;
//					opObj.dTemperature=dConvertedTemp;
//								
//					
//
//					break;
//				}
//			}
//		}
//		
//				
//
//		/*opDefault.dSigmaA=3000;
//		opDefault.dSigmaS=2000;
//		opDefault.dSigmaT=5000;*/
//		
//		return opObj;
//	}
//}
//
Opacity_Object Opacities::getOpacityData(Opacity_File_Details opacityFileDetails,vector<double> temperature_mat,string strEnergyGrp,double dConvertedTemp,double dDensity)
{

	map<string,Opacity_Object> opacitiesMap;

	Opacity_Object opObject1;
	vector<Opacity_Object> opacities_vector_1d(opacityFileDetails.iNoOfDensities,opObject1);
	vector<vector<Opacity_Object> > opacities_vector_2d(opacityFileDetails.iNoOfEnergyGrps,opacities_vector_1d);
	vector<vector<vector<Opacity_Object> > > opacities_vector(opacityFileDetails.iNoOfTempEntries,opacities_vector_2d);
	opacities_vector=opacityFileDetails.opacities_vector;
	
	int num_temps=temperature_mat.size();
	/*stringstream strStmOpacKey;	
	strStmOpacKey << strEnergyGrp<<":"<<dConvertedTemp<<":"<<dDensity;
	string strOpacitykey=strStmOpacKey.str();*/
	Opacity_Object opObj;
	Opacity_Object opObjFromVector;

	/*if(temperature_mat[0]==-1)
	{
		stringstream strStmOpacKey1;	
		strStmOpacKey1 << strEnergyGrp<<":"<<"-1"<<":"<<dDensity;
		string strOpacitykey1=strStmOpacKey1.str();
		return (opacitiesMap.find(strOpacitykey1))->second;
	}*/
	
	int iTempIndex=-1;
	int iDensityIndex=-1;
	int iEnergGrpIndex=-1;

	stringstream strStmTemp;
	strStmTemp << dConvertedTemp;
	if(((opacityFileDetails.tempMap).find(strStmTemp.str())) != (opacityFileDetails.tempMap).end())
		iTempIndex=(opacityFileDetails.tempMap).find(strStmTemp.str())->second;
	strStmTemp.str("");

	stringstream strStmDensity;

	//Call get interpolated density index based on  input density
	iDensityIndex=getDensityIndex(dDensity,opacityFileDetails);
	//strStmDensity << dDensity;

	//if(((opacityFileDetails.densityMap).find(strStmDensity.str())) != (opacityFileDetails.densityMap).end())
		//iDensityIndex=(opacityFileDetails.densityMap).find(strStmDensity.str())->second;
	//strStmDensity.str("");

	if(((opacityFileDetails.energyGrpMap).find(strEnergyGrp)) != (opacityFileDetails.energyGrpMap).end())
		iEnergGrpIndex=opacityFileDetails.energyGrpMap.find(strEnergyGrp)->second;

	if(iTempIndex != -1 && iDensityIndex != -1 && iEnergGrpIndex != -1)
	{
		//Opacities for the specific temperature, density and energy group combination is present.
		Opacity_Object opObjFromVector=opacities_vector[iTempIndex][iEnergGrpIndex][iDensityIndex];
		opObj.dDensity=opObjFromVector.dDensity;
		opObj.dSigmaA=opObjFromVector.dSigmaA;
		opObj.dSigmaS=opObjFromVector.dSigmaS;
		opObj.dSigmaT=opObjFromVector.dSigmaT;
		opObj.dTemperature=dConvertedTemp;
	}
	else
	{
		if(iTempIndex == -1 && iDensityIndex != -1 && iEnergGrpIndex != -1)
		{
			//Linear Interpolation of Temperature values.
			if(dConvertedTemp > temperature_mat[num_temps-1])
			{
				strStmTemp << temperature_mat[num_temps-1];
				iTempIndex=(opacityFileDetails.tempMap).find(strStmTemp.str())->second;
				opObjFromVector=opacities_vector[iTempIndex][iEnergGrpIndex][iDensityIndex];
				opObj.dDensity=opObjFromVector.dDensity;
				opObj.dSigmaA=opObjFromVector.dSigmaA;
				opObj.dSigmaS=opObjFromVector.dSigmaS;
				opObj.dSigmaT=opObjFromVector.dSigmaT;
				opObj.dTemperature=dConvertedTemp;
			}
			else if(dConvertedTemp < temperature_mat[0])
			{
				strStmTemp << temperature_mat[0];
				iTempIndex=(opacityFileDetails.tempMap).find(strStmTemp.str())->second;
				opObjFromVector=opacities_vector[iTempIndex][iEnergGrpIndex][iDensityIndex];
				opObj.dDensity=opObjFromVector.dDensity;
				opObj.dSigmaA=opObjFromVector.dSigmaA;
				opObj.dSigmaS=opObjFromVector.dSigmaS;
				opObj.dSigmaT=opObjFromVector.dSigmaT;
				opObj.dTemperature=dConvertedTemp;
			}
			else
			{
				for(int iCntLoop=1;iCntLoop < num_temps;iCntLoop++)
				{
					if((temperature_mat[iCntLoop-1] < dConvertedTemp)&&(temperature_mat[iCntLoop] > dConvertedTemp))
					{
						int iTempLeftIndex;
						stringstream strStmOpacKeyLeft;	
						strStmOpacKeyLeft << temperature_mat[iCntLoop-1];
						iTempLeftIndex=(opacityFileDetails.tempMap).find(strStmOpacKeyLeft.str())->second;

						int iTempRightIndex;
						stringstream strStmOpacKeyRight;	
						strStmOpacKeyRight << temperature_mat[iCntLoop];
						iTempRightIndex=(opacityFileDetails.tempMap).find(strStmOpacKeyRight.str())->second;

						Opacity_Object opDefaultLeft;					
						opDefaultLeft=opacities_vector[iTempLeftIndex][iEnergGrpIndex][iDensityIndex];
						Opacity_Object opDefaultRight;
						opDefaultRight=opacities_vector[iTempRightIndex][iEnergGrpIndex][iDensityIndex];

						double dFracChange=((dConvertedTemp-temperature_mat[iCntLoop-1])/(temperature_mat[iCntLoop]-temperature_mat[iCntLoop-1]));
						double dSigA=opDefaultLeft.dSigmaA+(dFracChange*(opDefaultRight.dSigmaA-opDefaultLeft.dSigmaA));						
						double dSigT=opDefaultLeft.dSigmaT+(dFracChange*(opDefaultRight.dSigmaT-opDefaultLeft.dSigmaT));
						double dSigS=opDefaultLeft.dSigmaS+(dFracChange*(opDefaultRight.dSigmaS-opDefaultLeft.dSigmaS));
						opObj.dDensity=dDensity;
						opObj.dSigmaA=dSigA;
						opObj.dSigmaS=dSigS;
						opObj.dSigmaT=dSigT;
						opObj.dTemperature=dConvertedTemp;
					

						break;
					}
				}
			}
		}
		else
		{
		}
	}

	return opObj;	
}

int Opacities::getDensityIndex(double dDensity,Opacity_File_Details opacityFileDetails)
{
	map<string,int> densityMap=opacityFileDetails.densityMap;
	map<string,int>::iterator it;
	double dCurrDensity;
	double dPrevDensity=0;
	double dIndex=0;
	Common_Utils common_Utils;
	vector<double> density_mat(densityMap.size(),0.0);

	for( map<string,int>::iterator it=densityMap.begin(); it!=densityMap.end(); ++it){
		dCurrDensity=atof(((*it).first).c_str());
		dIndex=(*it).second;
		density_mat[dIndex]=dCurrDensity;
	}
	
	for(int iLoop=0 ; iLoop < density_mat.size(); iLoop++){
		dCurrDensity=density_mat[iLoop];
		if(dDensity >= dPrevDensity && dDensity <= dCurrDensity)
		{
			if(((dDensity-dPrevDensity) < (dCurrDensity-dDensity)) && (dPrevDensity !=0 ))
			{
				return iLoop-1;
			}
			else
			{
				int retIndex=densityMap.find(common_Utils.convert_double_to_string(dCurrDensity))->second;
				return iLoop;
			}
		}
		else
		{
			dPrevDensity=dCurrDensity;
		}
	}

	return (density_mat.size()-1);
}


Opacity_File_Details Opacities::getopacitiesMapFromFile(string strFilename)
{

	Opacity_File_Details opacityFileDetails;

	string filename=strFilename;
	set<double>::iterator it;
	Opacity_Object opObj1;
	Opacity_Object opObj2;
	Opacity_Object opObj3;

	int num_temps;
	int num_dens;
	int num_groups;
	
	ifstream fin( filename.c_str() );
	if( !fin.is_open() ) {				
	}
	double conversion_factor;
	char line[250];
	string word;
	string sType="Macroscopic";
	fin.getline(line, 250);
	fin.getline(line, 250);

	// "This file is a " multigroup / single-group
	fin >> word >> word >> word >> word;
	string multi;
	fin >> multi;
	
	fin.getline(line, 250);
	fin >> num_temps >> word >> num_dens >> word >> word >> num_groups;

	map<string,int> densityMap;
	map<string,int> energyGrpMap;
	map<string,int> tempMap;
	//Create a 3-dimensional matrix with density in innermost vector followed
	// by energy _groups and then temperature.
	Opacity_Object opObject;
	vector<Opacity_Object> opacities_vector_1d(num_dens,opObject);
	vector<vector<Opacity_Object> > opacities_vector_2d(num_groups,opacities_vector_1d);
	vector<vector<vector<Opacity_Object> > > opacities_vector(num_temps,opacities_vector_2d);
		
	vector<double> temps( num_temps,0.0 );
	set<double> temps_set;
	double dDensity;
	vector<double> group_structure( num_groups+1,0.0 );
	vector<string> groups_mat(num_groups,"");	

	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin >> sType;
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int t=0; t!=num_temps; ++t )
	{
		fin >> temps[t];
		temps_set.insert(temps[t]);		
	}
	
	int iTempLoop=0;
	stringstream strStrmTemp1;
	for (it=temps_set.begin(); it!=temps_set.end(); it++)
	{
		strStrmTemp1 << *it;
		tempMap[strStrmTemp1.str()]=iTempLoop;
		iTempLoop++;
		strStrmTemp1.str("");
	}

	//get Temperatures
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int d=0; d!=num_dens; ++d )
	{
		fin >> dDensity;
		strStrmTemp1 << dDensity;
		densityMap[strStrmTemp1.str()]=d;
		strStrmTemp1.str("");
	}
	if(sType.compare("Macroscopic")==0){
		//Macroscopic
	}else{
		//Microscopic
	}
	
	stringstream strStrm;
	string tempStr="";
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int e=0; e!=num_groups+1; ++e )
	{
		fin >> group_structure[e];
	}
	for( int e=0; e < num_groups; e++ )
 
	{
		strStrm << (group_structure[e+1]) << "-" << group_structure[e];
		tempStr=strStrm.str();
		Opacity_Object opObj;
		groups_mat[e]=tempStr;
		energyGrpMap[tempStr]=(e);
		strStrm.str("");
	}
			
	fin.getline(line, 250);
	fin.getline(line, 250);
	double dTemp_temp;
	double dDensity_temp;

	Opacity_Object op_temp1;
	Opacity_Object op_temp2;
	map<string,Opacity_Object>::iterator it_temp1;
	stringstream strStrmTemp;
	string strTemp="";

	map<string,int>::iterator itTempIndex;
	map<string,int>::iterator itDensityIndex;
	stringstream strTempStrm;
	stringstream strDensityStrm;
	int iTempIndex;
	int iDensityIndex;
	
	int iIterTermination=1;
	while(fin.eof()==false)
	{
		fin >> word >> word >> dTemp_temp >> word>> word >>dDensity_temp;
		fin.getline(line, 250);
		fin.getline(line, 250);
		int iCrossSectionType;
		int iTermination=0;

		strTempStrm << dTemp_temp;
		itTempIndex=tempMap.find(strTempStrm.str());
		strTempStrm.str("");
		iTempIndex=itTempIndex->second;

		strDensityStrm << dDensity_temp;
		itDensityIndex=densityMap.find(strDensityStrm.str());
		strDensityStrm.str("");
		iDensityIndex=itDensityIndex->second;
				
		while(iTermination < 2)
		{
			fin >> word >> iCrossSectionType;		
			fin.getline(line, 250);
			if(iCrossSectionType==3011 && word.compare("MT")==0)
			{				
				for( int e=0; e!=num_groups; ++e )					
				{
					op_temp1=opacities_vector[iTempIndex][e][iDensityIndex];
					double dSigmA;
					fin >> dSigmA;		
					op_temp1.dSigmaA=dSigmA;
					op_temp1.dSigmaT=dSigmA;
					op_temp1.dSigmaS=op_temp1.dSigmaT-dSigmA;
					/*op_temp1.dSigmaA=100000;
					op_temp1.dSigmaT=7800000;
					op_temp1.dSigmaS=7700000;*/
					op_temp1.dTemperature=dTemp_temp;
					op_temp1.dDensity=dDensity_temp;
					opacities_vector[iTempIndex][e][iDensityIndex]=op_temp1;					
				}
				iTermination++;
			}else if(iCrossSectionType==3012 && word.compare("MT")==0)				
			{
				for( int e=0; e!=num_groups; ++e )
				{
					double dSigmT;
					fin >> dSigmT;					
				}
				iTermination++;			
			}else{
				double num;
				if((iCrossSectionType==3001 || iCrossSectionType==3002) && word.compare("MT")==0)
				{
					fin >> num;
				}
				else
				{
					for( int e=0; e!=num_groups; ++e )
						fin >> num;
				}
			}			
		}
	}

	fin.close();

	opacityFileDetails.densityMap=densityMap;
	opacityFileDetails.energyGrpMap=energyGrpMap;
	opacityFileDetails.tempMap=tempMap;
	opacityFileDetails.opacities_vector=opacities_vector;
	
	return opacityFileDetails;
}

vector<double> Opacities::getTemperatureVector(string strFileName)
{
	
	int num_temps;
	int num_dens;
	int num_groups;

	string filename=strFileName;

	ifstream fin( filename.c_str() );
	if( !fin.is_open() ) {
				
	}
	double conversion_factor;
	char line[250];
	string word;
	string sType="Macroscopic";
	fin.getline(line, 250);
	fin.getline(line, 250);

	// "This file is a " multigroup / single-group
	fin >> word >> word >> word >> word;
	string multi;
	fin >> multi;
	
	fin.getline(line, 250);
	fin >> num_temps >> word >> num_dens >> word >> word >> num_groups;
		
	vector<double> temps( num_temps,0.0 );
	set<double> temps_set;
	vector<double> densities( num_dens,0.0 );
	vector<double> group_structure( num_groups+1,0.0 );
	vector<string> groups_mat(num_groups,"");
	vector<double> sigmaA_mat( num_groups,0.0 );
	vector<double> sigmaS_mat( num_groups,0.0 );
	vector<double> sigmaT_mat( num_groups,0.0 );

	vector<vector<double> > sii_centers_curr(num_temps,vector<double>(num_groups,0.0));
	

	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin >> sType;
	fin.getline(line, 250);
	fin.getline(line, 250);
	fin.getline(line, 250);
	for( int t=0; t!=num_temps; ++t )
	{
		fin >> temps[t];
		temps_set.insert(temps[t]);
	}
	set<double>::iterator it;
	vector<double> temps1( num_temps,0.0 );
	int iLoop=0;
	for (it=temps_set.begin(); it!=temps_set.end(); it++)
	{
		temps1[iLoop]=*it;
		iLoop++;
	}
	fin.close();
	return temps1;
}




Opacities::Opacities(void)
{

}


Opacities::~Opacities(void)
{
}

	




