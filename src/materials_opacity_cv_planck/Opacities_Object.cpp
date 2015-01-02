#include "Opacities_Object.h"
#include "Opacity_Object.h"
#include <set>
#include <map>

int no_of_temperatures;
map<string,Opacity_Object> opacitiesMap;
set<double> temps_set;

Opacities_Object::Opacities_Object(void)
{
}


Opacities_Object::~Opacities_Object(void)
{
}
