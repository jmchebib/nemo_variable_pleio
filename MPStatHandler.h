/** $Id: MPStatHandler.h,v 1.10 2015-07-13 08:52:58 fred Exp $
*
*  @file MPStatHandler.h
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
*  created on @date 11.07.2005
* 
*  @author fred
*/

#ifndef MPSTATHANDLER_H
#define MPSTATHANDLER_H

#include <list>
#include "stathandler.h"

class Individual;
class Metapop;
/**@todo Add pedigree stats*/
/**A StatHandler for the Metapop SimComponent.*/
class MPStatHandler : public StatHandler<MPStatHandler> {

  double meanEmigrant,meanImigrant,meanResidant,meanKolonisers,meanDeadDisp;
  double ObservedExtinctionRate;
  double _sib_prop[5], _ped_prop[5];
  double _var_reprod_success;
  
public:
  
  MPStatHandler( ) { }
  
  virtual ~MPStatHandler() {}

  virtual bool setStatRecorders(std::string& token);
  
  void   setStatsForDemography(age_t AGE);
  void   setStatsForPop();
  void   setStatsForPopPerPatch();
  void   setStatsForMigrants();
  void   setStatsForMigrantsPerPatch();
  void   addIndNumPerPatch(sex_t SEX, age_t AGE);
  void   addPatchAge();
  ///@name Migration
  ///@{
  double getMeanEmigrantPerPatch           ();
  double getMeanImigrantPerPatch           ();
  double getMeanMigrantRatio               ();
  double getMeanResidantPerPatch           ();
  double getMeanKolonisersProportion       ();
  double getMeanKolonisersPerPatch         ();
  double getEmigrantInPatch  (unsigned int i);
  double getResidantInPatch  (unsigned int i);
  double getImigrateInPatch  (unsigned int i);
  double getKolonisersInPatch(unsigned int i);
  ///@}
  ///@name Patch extinction
  ///@{
  void   setObsrvdExtinctionRate         ();
  double getObsrvdExtinctionRate         ()  {return ObservedExtinctionRate;}
  double get_isAlive                     ();
  double getPatchAge                     (unsigned int i);
  double getMeanPatchAge                 ();
  ///@}
  ///@name Demography
  ///@{
  double getAdultSexRatio           ();
  double getOffsprgSexRatio         ();
  double getPatchSize               (unsigned int age, unsigned int patch);
  double getPopulationSize          (unsigned int age);
  double getMeanPatchSize           (unsigned int age);
  double getMeanPatchSizePerSex     (unsigned int sex, unsigned int age);
  double getMeanPatchDensity        (age_t AGE);
  double getMeanPatchDensityVariance(age_t AGE);
  
  double getFemNumber               (unsigned int i);
  double getMalNumber               (unsigned int i);
  double getOffFemNumber            (unsigned int i);
  double getOffMalNumber            (unsigned int i);
  double getMeanAssignedFecundity   (unsigned int sex);
  double getMeanMatings             (unsigned int sex);
  double setReproductiveStats       (unsigned int sex);
  double getReproductiveVar         () {return _var_reprod_success;}
  ///@}
  ///@name Kinship
  ///@{
  void   setKinship                 ();
  void   setKinClassCounter         (Individual *I1, Individual *I2);
  double getSibProportion           (unsigned int i) {return _sib_prop[i];}
  ///@}
  ///@name Pedegree
  ///@{
  void   setPedegreeCount           ();
  double getPedProportion           (unsigned int i) {return _ped_prop[i];}
  ///@}
};

#endif
