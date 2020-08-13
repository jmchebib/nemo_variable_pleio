/** $Id: MPStatHandler.cc,v 1.13 2015-07-13 08:52:57 fred Exp $
 *
 *  @file MPStatHandler.cc
 *  Nemo2
 *
 *   Copyright (C) 2006-2015 Frederic Guillaume
 *   frederic.guillaume@ieu.uzh.ch
 *
 *   This file is part of Nemo
 *
 *   Nemo is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   Nemo is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  
 *  Created on @date 11.07.05
 *  @author fred
 */

#include <sstream>
#include "MPStatHandler.h"
#include "metapop.h"
#include "output.h"
#include "tstring.h"

// ------------------------------------------------------------------------------

//                             StatHandler

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool MPStatHandler::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-MPStatHandler::setStatRecorders ");
#endif
  
  if(token.compare("extrate") == 0) {
    
    add("Extinction rate","extrate",ALL,0,0,&MPStatHandler::getObsrvdExtinctionRate,0,0,
        &MPStatHandler::setObsrvdExtinctionRate);
		
  } else if(token.compare("pop") == 0) {
    
    setStatsForPop();
    
  } else if(token.compare("pop.patch") == 0) {
    
    setStatsForPopPerPatch();
    
  } else if(token.compare("off.fem.patch") == 0) {

    addIndNumPerPatch(FEM, OFFSPRG);
    
  } else if(token.compare("off.mal.patch") == 0) {
  
    addIndNumPerPatch(MAL, OFFSPRG);

  } else if(token.compare("adlt.fem.patch") == 0) {
  
    addIndNumPerPatch(FEM, ADULTS);

  } else if(token.compare("adlt.mal.patch") == 0) {

    addIndNumPerPatch(MAL, ADULTS);

  } else if(token.compare("adlt.sexratio") == 0) {
    
    add("Adults sex ratio","adlt.sexratio",ADULTS,0,0,&MPStatHandler::getAdultSexRatio,0,0,0);
    
  } else if(token.compare("off.sexratio") == 0) {
    
    add("Offspring Sex Ratio","off.sexratio",OFFSPRG,0,0,&MPStatHandler::getOffsprgSexRatio,0,0,0);
    
  } else if(token.compare("migrants") == 0) {
    
    setStatsForMigrants();
    
  } else if(token.compare("migrants.patch") == 0) {
    
    setStatsForMigrantsPerPatch();
    
  } else if(token.compare("demography") == 0) {
    
    setStatsForDemography(OFFSPRG);
    setStatsForDemography(ADULTS);
    
  } else if(token.compare("off.demography") == 0) {
    
    setStatsForDemography(OFFSPRG);

  } else if(token.compare("adlt.demography") == 0) {

    setStatsForDemography(ADULTS);
    
  } else if(token.compare("fecundity") == 0) {
    
    add("Females assigned fecundity","adlt.femfec",ADULTS,1,0,0,&MPStatHandler::getMeanMatings,0,0);
    add("Females realized fecundity","adlt.femrealfec",ADULTS,1,0,0,&MPStatHandler::setReproductiveStats,0,0);
    add("Females reproductive var  ","adlt.femvarfec",ADULTS,0,0,&MPStatHandler::getReproductiveVar,0,0,0);
    add("Males realized fecundity  ","adlt.malrealfec",ADULTS,0,0,0,&MPStatHandler::setReproductiveStats,0,0);   
    add("Males reproductive var    ","adlt.malvarfec",ADULTS,0,0,&MPStatHandler::getReproductiveVar,0,0,0); 
    
  } else if(token.compare("kinship") == 0) {
    
    add("Proportion of full-sib offspring","off.fsib",OFFSPRG,3,0,0,&MPStatHandler::getSibProportion,0,&MPStatHandler::setKinship);
    add("Proportion of paternal half-sib ","off.phsib",OFFSPRG,2,0,0,&MPStatHandler::getSibProportion,0,0);
    add("Proportion of maternal half-sib ","off.mhsib",OFFSPRG,1,0,0, &MPStatHandler::getSibProportion,0,0);
    add("Proportion of non-sib offspring ","off.nsib",OFFSPRG,0,0,0,&MPStatHandler::getSibProportion,0,0);
    add("Proportion of selfed offspring  ","off.self",OFFSPRG,4,0,0,&MPStatHandler::getSibProportion,0,0);
    
  } else if(token == "pedigree") {
    
    add("Proportion of other deme mated","ped.outb",OFFSPRG,0,0,0,&MPStatHandler::getPedProportion,0,&MPStatHandler::setPedegreeCount);
    add("Proportion of same deme mated ","ped.outw",OFFSPRG,1,0,0,&MPStatHandler::getPedProportion,0,0);
    add("Proportion of half-sib mated ","ped.hsib",OFFSPRG,2,0,0, &MPStatHandler::getPedProportion,0,0);
    add("Proportion of full-sib mated ","ped.fsib",OFFSPRG,3,0,0,&MPStatHandler::getPedProportion,0,0);
    add("Proportion of selfed mated ","ped.self",OFFSPRG,4,0,0,&MPStatHandler::getPedProportion,0,0);    
    
  } else
    return false;
  
  return true;
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::setStatsForPop()
{
  setStatsForDemography(OFFSPRG);
  add("Offspring Sex Ratio","off.sexratio",OFFSPRG,0,0,&MPStatHandler::getOffsprgSexRatio,0,0,0);

  setStatsForDemography(ADULTS);
  add("Adults sex ratio","adlt.sexratio",ADULTS,0,0,&MPStatHandler::getAdultSexRatio,0,0,0);
  
  add("Extinction rate","extrate",ALL,0,0,&MPStatHandler::getObsrvdExtinctionRate,0,0,
      &MPStatHandler::setObsrvdExtinctionRate);
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::setStatsForDemography(age_t AGE)
{
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name = (AGE == ADULTS ? "Adults ":"Offspring ");
  
  add(name + "number", suffix + "nbr", AGE, AGE, 0, 0, &MPStatHandler::getPopulationSize,0,0);
  add(name + "fem nbr", suffix + "nbfem", AGE, FEM, AGE, 0, 0, &MPStatHandler::getMeanPatchSizePerSex,0);
  add(name + "mal nbr", suffix + "nbmal", AGE, MAL, AGE, 0, 0, &MPStatHandler::getMeanPatchSizePerSex,0);
  add(name + "density", suffix + "density", AGE, AGE, 0, 0, &MPStatHandler::getMeanPatchDensity,0,0);
  add(name + "density var", suffix + "dvar", AGE, AGE, 0, 0, &MPStatHandler::getMeanPatchDensityVariance,0,0);
  
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::setStatsForPopPerPatch()
{
  addIndNumPerPatch(FEM, OFFSPRG);
  addIndNumPerPatch(MAL, OFFSPRG);
  addIndNumPerPatch(FEM, ADULTS);
  addIndNumPerPatch(MAL, ADULTS);
  addPatchAge();
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::addIndNumPerPatch(sex_t SEX, age_t AGE)
{
  double (MPStatHandler::* setter) (unsigned int) = 0;
  
  if (SEX) {
    
    if (AGE == ADULTS) {
      setter = &MPStatHandler::getFemNumber;
    } else if (AGE == OFFSPRG) {
      setter =  &MPStatHandler::getOffFemNumber;
    }
  
  } else {
    
    if (AGE == ADULTS) {
      setter = &MPStatHandler::getMalNumber;
    } else if (AGE == OFFSPRG) {
      setter =  &MPStatHandler::getOffMalNumber;
    }
    
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name = (SEX == FEM ? suffix + "fem." : suffix + "mal.");
  string patch = "p1";
    
  add("", name + patch, AGE, 0, 0, 0, setter, 0, 0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    patch = "p" + tstring::int2str(i+1);
    add("", name + patch, AGE, i, 0, 0, setter, 0, 0);
  }
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::addPatchAge()
{
  
  std::ostringstream name, sub_name;
  
  add("Patch 1 age","age.patch1",ALL,0,0,0,&MPStatHandler::getPatchAge,0,0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    name<<"Patch "<<i+1<<" age";
    sub_name<<"age.patch"<<i+1;
    add(name.str(),sub_name.str(),ALL,i,0,0,&MPStatHandler::getPatchAge,0,0);
    name.str("");
    sub_name.str("");
  }
  
  add("Mean patch age","patch.avg.age",ALL,0,0,&MPStatHandler::getMeanPatchAge,0,0,0);
  add("Extinction rate","extrate",ALL,0,0,&MPStatHandler::getObsrvdExtinctionRate,0,0,
      &MPStatHandler::setObsrvdExtinctionRate);
}

// ----------------------------------------------------------------------------------------
void MPStatHandler::setStatsForMigrants()
{
  add("Emigrant nbr","emigrants",ALL,0,0,&MPStatHandler::getMeanEmigrantPerPatch,0,0,0);
  add("Imigrant nbr","imigrants",ALL,0,0,&MPStatHandler::getMeanImigrantPerPatch,0,0,0);
  add("Residant nbr","residents",ALL,0,0,&MPStatHandler::getMeanResidantPerPatch,0,0,0);
  add("Imigration rate","imigrate",ALL,0,0,&MPStatHandler::getMeanMigrantRatio,0,0,0);
  add("Coloniser nbr","colonisers",ALL,0,0,&MPStatHandler::getMeanKolonisersPerPatch,0,0,0);
  add("Colonisation rate","colonrate",ALL,0,0,&MPStatHandler::getMeanKolonisersProportion,0,0,0);
}
// ----------------------------------------------------------------------------------------
void MPStatHandler::setStatsForMigrantsPerPatch()
{  std::ostringstream name, sub_name;
  
  add("Emigrants patch 1","emigr.p1", ALL,0,0,0,&MPStatHandler::getEmigrantInPatch,0,0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    name<<"Emigrants patch "<<i+1;
    sub_name<<"emigr.p"<<i+1;
    add(name.str(),sub_name.str(),ALL,i,0,0,&MPStatHandler::getEmigrantInPatch,0,0);
    name.str("");
    sub_name.str("");
  }
  
  add("Residants patch 1","resid.p1", ALL,0,0,0,&MPStatHandler::getResidantInPatch,0,0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    name<<"Residants patch "<<i+1;
    sub_name<<"resid.p"<<i+1;
    add(name.str(),sub_name.str(),ALL,i,0,0,&MPStatHandler::getResidantInPatch,0,0);
    name.str("");
    sub_name.str("");
  }
  
  add("Imig. rate patch 1","imrate.p1", ALL,0,0,0,&MPStatHandler::getImigrateInPatch,0,0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    name<<"Imig. rate patch "<<i+1;
    sub_name<<"imrate.p"<<i+1;
    add(name.str(),sub_name.str(),ALL,i,0,0,&MPStatHandler::getImigrateInPatch,0,0);
    name.str("");
    sub_name.str("");
  }
  
  add("Colons patch 1","colo.p1", ALL,0,0,0,&MPStatHandler::getKolonisersInPatch,0,0);
  for(unsigned int i = 1; i < _pop->getPatchNbr(); i++) {
    name<<"Colons patch "<<i+1;
    sub_name<<"colo.p"<<i+1;
    add(name.str(),sub_name.str(),ALL,i,0,0,&MPStatHandler::getKolonisersInPatch,0,0);
    name.str("");
    sub_name.str("");
  }
}
