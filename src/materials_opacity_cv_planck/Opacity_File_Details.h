#ifndef Rad_Helper_h
#define Rad_Helper_h
#include "Rad_Helper.h"
#endif

#ifndef Opacity_Object_h
#define Opacity_Object_h
#include "Opacity_Object.h"
#endif

#include <vector>
#include <map>

using std::vector;
using std::map;

class Opacity_File_Details
{
public:
	int iNoOfTempEntries;
	int iNoOfEnergyGrps;
	int iNoOfDensities;
	//For a constant density
	//double dDensity;
	vector<double> grp_vector;
	vector<double> density_mat;
	char isConstantOpacity;

	map<string,int> densityMap;
	map<string,int> energyGrpMap;
	map<string,int> tempMap;
	vector<vector<vector<Opacity_Object> > > opacities_vector;
	
	Opacity_File_Details(void);
	~Opacity_File_Details(void);
};

