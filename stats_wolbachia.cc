/**  $Id: stats_wolbachia.cc,v 1.7 2015-07-13 08:52:56 fred Exp $
*
*  @file stats_wolbachia.cc
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
*  Created on @date 05.08.2004
*
*  @author fred

*/
#include <cmath>
#include "ttwolbachia.h"
#include "metapop.h"

// ----------------------------------------------------------------------------------------
void TTWolbachiaSH::setInfectionStats ( )
{
  double Fsize = 0, Msize = 0, val, local_inf;
  double *all;
  Patch* crnt_patch;
    
  Fsize = _pop->size(FEM, ADULTS);
  Msize = _pop->size(MAL, ADULTS);

  if( Fsize != 0 )
    all = new double [ _pop->getPatchNbr() ];
  else {
    _Fmean = nanf("NULL");
    _Mmean = nanf("NULL");
    _extrate = nanf("NULL");
    _var = nanf("NULL");
    return;
  }
    
  _Fmean = _Mmean = 0;
  _extrate = 0;

  //females:
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
    crnt_patch = _pop->getPatch(i);
    local_inf = 0;
    all[i] = 0;
    unsigned int fsize = crnt_patch->size(FEM, ADLTx);
    for(unsigned int j = 0; j < fsize; ++j) {
      val = (double)*(bool*)(crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_TTidx));
      local_inf += val;
    }
    
    if(local_inf != 0) {
      all[i] = local_inf / fsize;
      _Fmean += all[i];
    } else  _extrate++;
      
  }
  
  _Fmean = (Fsize != 0 ? _Fmean / (_pop->getPatchNbr() - _extrate) : nanf("NULL"));
  
  _extrate /= _pop->getPatchNbr();
 
  _var = 0;
  double local_mean = 0;
  
  for (unsigned int i = 0; i < _pop->getPatchNbr(); i++)
    local_mean += all[i];
  
  local_mean /= _pop->getPatchNbr();
  
  //hacked from the GSL:
  for (unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    const long double delta = (all[i] - local_mean);
    _var += (delta * delta - _var) / (i + 1);
  }

 //males
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
    crnt_patch = _pop->getPatch(i);
    for(unsigned int j = 0, msize = crnt_patch->size(MAL, ADLTx); j < msize; ++j)
      _Mmean += (double)*(bool*)(crnt_patch->get(MAL, ADLTx, j)->getTraitValue(_TTidx));
  }
  
  _Mmean = (Msize != 0 ? _Mmean / Msize : nanf("NULL"));
  
  delete [] all;
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getMeanOffsprgInfection(unsigned int sex)
{
  unsigned int i,j;
  double indNbr = 0, mean = 0;
  Patch* crnt_patch;
  
  
  if(sex) {//females:
    for(i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      for(j = 0; j < crnt_patch->size(FEM, OFFSx); ++j)
        mean += (double)*(bool*)crnt_patch->get(FEM, OFFSx, j)->getTraitValue(_TTidx);
	}
    indNbr = _pop->size(FEM, OFFSPRG);
  } else { //males
    for(i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      for(j = 0; j < crnt_patch->size(MAL, OFFSx); ++j)
        mean += (double)*(bool*)crnt_patch->get(MAL, OFFSx, j)->getTraitValue(_TTidx);
	}
    indNbr = _pop->size(MAL, OFFSPRG);
  }
  
  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));  
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getMeanFemaleInfection_perPatch( unsigned int i)
{
  double indNbr = 0, mean = 0;
  Patch* crnt_patch = _pop->getPatch(i);
  
  indNbr = crnt_patch->size(FEM, ADLTx);
  for(unsigned int j = 0; j < indNbr; ++j)
    mean += (double)*(bool*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_TTidx);
  
  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getMeanMaleInfection_perPatch( unsigned int i)
{
  double indNbr = 0, mean = 0;
  Patch* crnt_patch = _pop->getPatch(i);
  
  indNbr = crnt_patch->size(MAL, ADLTx);
  for(unsigned int j = 0; j < indNbr; ++j)
    mean += (double)*(bool*)crnt_patch->get(MAL, ADLTx, j)->getTraitValue(_TTidx);
  
  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getMeanOffsprgFemaleInfection_perPatch( unsigned int i)
{
  double indNbr = 0, mean = 0;
  Patch* crnt_patch = _pop->getPatch(i);
  
  indNbr = crnt_patch->size(FEM, OFFSx);
  
  for(unsigned int j = 0; j < indNbr; ++j)
      mean += (double)*(bool*)crnt_patch->get(FEM, OFFSx, j)->getTraitValue(_TTidx);

  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getMeanOffsprgMaleInfection_perPatch( unsigned int i)
{
  double indNbr = 0, mean = 0;
  Patch* crnt_patch = _pop->getPatch(i);
  
  indNbr = crnt_patch->size(MAL, OFFSx);
  
  for(unsigned int j = 0; j < indNbr; ++j)
    mean += (double)*(bool*)crnt_patch->get(MAL, OFFSx, j)->getTraitValue(_TTidx);

  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
double TTWolbachiaSH::getIcompatibleMatingFreq()
{
//  double sum = 0;
//  Patch* crnt_patch;
//  
//  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
//	crnt_patch = _pop->getPatch(i);
//	for(unsigned int j = 0; j < crnt_patch->size(FEM, ADLTx);++j)
//	  sum += crnt_patch->get(FEM, ADLTx, j)->getFecundity();
//  }
//  return (double)Patch::DeadOffsprgs/sum;
  return nanf("NULL");
}

