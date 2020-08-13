/** $Id: stats_disp.cc,v 1.7 2015-07-13 08:52:58 fred Exp $
 *
*  @file stats_disp.cc
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
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*  
*  Created on @date 23.01.2004.
*  @author fred
**/

#include <cmath>
#include "ttdispersal.h"


// ------------------------------------------------------------------------------

//                          dispersal analysis

// ----------------------------------------------------------------------------------------
// getMeanDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getMeanDispRate()
{
  double meanDisp = 0;
  unsigned int fem = 0, mal = 0, tot = 0,i,j;
  Patch* crnt_patch;
  
  _meanFemDisp = 0;
  _meanMalDisp = 0;
  
  for( i = 0; i < _pop->getPatchNbr(); ++i) {
	crnt_patch = _pop->getPatch(i);
    for( j = 0; j < crnt_patch->size(FEM, ADLTx); ++j) {
      _meanFemDisp += *(double*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_fdispIdx);
      fem++;
    }
    for( j = 0; j < crnt_patch->size(MAL, ADLTx);++j) {
      _meanMalDisp += *(double*)crnt_patch->get(MAL, ADLTx, j)->getTraitValue(_mdispIdx);
      mal++;
    }
  }
  meanDisp = _meanFemDisp + _meanMalDisp;
  tot = fem + mal;
  _meanFemDisp = (fem != 0 ? _meanFemDisp/fem : 0.0);
  _meanMalDisp = (mal != 0 ? _meanMalDisp/mal : 0.0);
  
  return (tot != 0 ? meanDisp/tot : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getMeanDispRate(bool sex)
{
  double meanDispRate = 0;
  unsigned int i,j,indnbr = 0;
  Patch* crnt_patch;
  
  if(sex) {
	for( i=0;i<_pop->getPatchNbr();++i) {
	  crnt_patch = _pop->getPatch(i);
	  for( j=0;j<crnt_patch->size(FEM, ADLTx);++j) {
		meanDispRate += *(double*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue( _fdispIdx);
		indnbr++;
	  }
	}
  } else {
	for( i=0;i<_pop->getPatchNbr();++i) {
	  crnt_patch = _pop->getPatch(i);
	  for( j=0;j<crnt_patch->size(MAL, ADLTx);++j) {
		meanDispRate += *(double*)crnt_patch->get(MAL, ADLTx, j)->getTraitValue(_mdispIdx);
		indnbr++;
	  }
	}
  }
  
  return (indnbr != 0 ? meanDispRate/indnbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getOffsprgMeanDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getOffsprgMeanDispRate()
{
  double meanDisp = 0;
  unsigned int fem = 0, mal = 0, tot = 0;
  Patch* crnt_patch;
  
  _meanOffFemDisp = 0;
  _meanOffMalDisp = 0;
  
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
    
	crnt_patch = _pop->getPatch(i);
    
    for(unsigned int j = 0, size = crnt_patch->size(FEM, OFFSx); j < size; ++j)
      _meanOffFemDisp += *(double*)crnt_patch->get( FEM, OFFSx, j )->getTraitValue( _fdispIdx );
    
    for(unsigned int j = 0, size = crnt_patch->size(MAL, OFFSx); j < size; ++j)
      _meanOffMalDisp += *(double*)crnt_patch->get( MAL, OFFSx, j )->getTraitValue( _mdispIdx );

  }

  meanDisp = _meanOffFemDisp + _meanOffMalDisp;
  _meanOffFemDisp = ( (fem = _pop->size(FEM, OFFSPRG)) != 0 ? _meanOffFemDisp/fem : 0.0);
  _meanOffMalDisp = ( (mal = _pop->size(MAL, OFFSPRG)) != 0 ? _meanOffMalDisp/mal : 0.0);
  
  return ( (tot = _pop->size(OFFSPRG)) != 0 ? meanDisp/tot : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getOffsprgMeanDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getOffsprgMeanDispRate(bool sex)
{
  double meanDispRate = 0;
  unsigned int indnbr = 0;
  int TType = (sex ? _fdispIdx : _mdispIdx);
  Patch* crnt_patch;
  
  if(sex) {
    
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i){
      
      crnt_patch = _pop->getPatch(i);
      
      for(unsigned int j = 0, size = crnt_patch->size(FEM, OFFSx); j < size; ++j)
        meanDispRate += *(double*)crnt_patch->get( FEM, OFFSx, j )->getTraitValue(TType);
    }    
    indnbr = _pop->size(FEM, OFFSPRG);
    
  } else {
    
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i){
      
      crnt_patch = _pop->getPatch(i);
      
      for(unsigned int j = 0, size = crnt_patch->size(MAL, OFFSx); j < size; ++j)
        meanDispRate += *(double*)crnt_patch->get( MAL, OFFSx, j )->getTraitValue(TType);
    }    
    indnbr = _pop->size(MAL, OFFSPRG);
  }
  
  return (indnbr != 0 ? meanDispRate/indnbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanFemDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getMeanFemDispRate()
{
  double meanDispRate = 0;
  unsigned int indnbr = 0;
  Patch* crnt_patch;
  
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
    
	crnt_patch = _pop->getPatch(i);
    
    for(unsigned int j = 0, size = crnt_patch->size(FEM, ADLTx); j < size; ++j) 
	  meanDispRate += *(double*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_fdispIdx);
  }
  indnbr = _pop->size(FEM, ADULTS);
  return (indnbr != 0 ? meanDispRate/indnbr : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanMalDispRate
// ----------------------------------------------------------------------------------------
double TTDispersalSH::getMeanMalDispRate()
{
  double meanDispRate = 0;
  unsigned int indnbr = 0;
  Patch* crnt_patch;
  
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	
    crnt_patch = _pop->getPatch(i);
	
    for(unsigned int j = 0, size = crnt_patch->size(MAL, ADLTx); j < size; ++j) 
	  meanDispRate += *(double*)crnt_patch->get(MAL, ADLTx, j)->getTraitValue(_mdispIdx);
  }
  indnbr = _pop->size(MAL, ADULTS);
  return (indnbr != 0 ? meanDispRate/indnbr : nanf("NULL"));
}

