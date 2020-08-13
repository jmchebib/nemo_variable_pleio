/** $Id: stats_coa.cc,v 1.9 2015-07-13 08:52:58 fred Exp $
 *
 *  @file stats_coa.cc
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
 *  Created on @date 23.01.2004
 *
 *  @author fred
 */

#include <iostream>
#include <cmath>
#include "ttneutralgenes.h"
#include "individual.h"

using namespace std;

// ------------------------------------------------------------------------------

//                          coancestries analysis
// ----------------------------------------------------------------------------------------
// coancestry
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::Coancestry(void** ind1, void** ind2, unsigned int nb_locus)
{
  unsigned int p = 0;
  unsigned char **seq1, **seq2;
  
  seq1 = (unsigned char**)ind1;
  seq2 = (unsigned char**)ind2;
  
  for (unsigned int k = 0; k < nb_locus; ++k)
    p += !(seq1[0][k]^seq2[0][k]) + !(seq1[0][k]^seq2[1][k]) + !(seq1[1][k]^seq2[0][k]) + !(seq1[1][k]^seq2[1][k]);
  
  return (double)p/(4.0*nb_locus);
}
// ----------------------------------------------------------------------------------------
// setCoaMatrix
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setCoaMatrix (age_idx age_pos, unsigned char dim)
{
  unsigned int Fsize,Msize,tot_size,size_i,size_l, wt;
  unsigned int nb_coeff = 0;
  unsigned int i, j, k, l;
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  Patch *P1, *P2;
  
  if(_coa_matrix == NULL)
    
    _coa_matrix = new TMatrix(patchNbr, patchNbr);
  
  else if( _coa_matrix->length() != patchNbr * patchNbr)
    
    _coa_matrix->reset(patchNbr, patchNbr);
  
  _coa_matrix->assign(0);
  
  if(dim & 1) {
    
    _mean_theta = 0;
    
    wt = 0;
    
    //first fill the diagonale: within deme coancestry (theta)
    for(i = 0; i < patchNbr; ++i) {
      
      P1 = _pop->getPatch(i);
      Fsize = P1->size(FEM, age_pos);
      Msize = P1->size(MAL, age_pos);
      
      tot_size = Fsize + Msize;
      
      if(tot_size != 0) {
        
        //fem-fem coa
          for(j = 0; j < Fsize; ++j)
            for(k = j+1; k < Fsize; ++k)
              _coa_matrix->plus(i, i, Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 P1->get(FEM, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 nb_locus));
        //mal-mal coa
        for(j = 0; j < Msize; ++j)
          for(k = j+1; k < Msize; ++k)
            _coa_matrix->plus(i, i, Coancestry(P1->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                               P1->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                               nb_locus));
        //fem-mal coa
        for(j = 0; j < Fsize; ++j)
          for(k = 0; k < Msize; ++k)
            _coa_matrix->plus(i, i, Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                               P1->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                               nb_locus));
        
        _coa_matrix->divide(i, i, tot_size*(tot_size -1)/2.0); 
      }//end if
      
      _mean_theta += tot_size * _coa_matrix->get(i, i);
      
      wt += tot_size;
      
    }//end for patchNbr
    
    //weighted average:
    _mean_theta /= wt;
    
  } //end if diag
  
  if(dim & 2) {
    //fill the first upper half of the matrix: between deme coancestry (alpha)
    
    _mean_alpha = 0;
    
    wt = 0;
    
    for(i = 0; i < patchNbr-1; ++i) {
      
      P1 = _pop->getPatch(i);
      
      for(l = i+1; l < patchNbr; ++l) {
        
        P2 = _pop->getPatch(l);
        
        tot_size = P1->size(age_pos) * P2->size(age_pos);
        
        if(tot_size != 0) {
          //females i vs. females l
          size_i = P1->size(FEM, age_pos);
          size_l = P2->size(FEM, age_pos);
          nb_coeff = size_i * size_l;
          
          for(j = 0; j < size_i; ++j)
            for(k = 0; k < size_l; ++k)
              _coa_matrix->plus(i, l, Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 P2->get(FEM, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 nb_locus));
          //females i vs. males l
          size_l = P2->size(MAL, age_pos);
          nb_coeff += size_i * size_l;
          
          for(j = 0; j < size_i; ++j)
            for(k = 0; k < size_l; ++k)
              _coa_matrix->plus(i, l, Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 P2->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 nb_locus));
          //males i vs. males l
          size_i = P1->size(MAL, age_pos);
          size_l = P2->size(MAL, age_pos);
          nb_coeff += size_i * size_l;
          
          for(j = 0; j < size_i; ++j)
            for(k = 0; k < size_l; ++k)
              _coa_matrix->plus(i, l, Coancestry(P1->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 P2->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 nb_locus));
          //males i vs. females l
          size_l = P2->size(FEM, age_pos);
          nb_coeff += size_i * size_l;
          
          for(j = 0; j < size_i; ++j)
            for(k = 0; k < size_l; ++k)
              _coa_matrix->plus(i, l, Coancestry(P1->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 P2->get(FEM, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                                 nb_locus));
          _coa_matrix->divide(i, l, nb_coeff);
        }//endif
        
        _mean_alpha += tot_size * _coa_matrix->get(i, l);
        
        wt += tot_size;
        
      }//end for P2
    }//end for P1 
    
    //weighted average:
    _mean_alpha /= wt;
    
  }//end if upper half
  
}
// ----------------------------------------------------------------------------------------
// setAdults_Theta
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setAdults_Theta()
{
  unsigned int Fsize,Msize,i,j,k, FFsize, MMsize, FMsize;
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  Patch *P1;
  double mean = 0, grand_mean = 0;
  
  Theta_FF = 0;
  Theta_MM = 0;
  Theta_FM = 0;
  
  _mean_theta = 0;
  
  for(i = 0; i < patchNbr; ++i) {
    
    P1 = _pop->getPatch(i);
    Fsize = P1->size(FEM, ADLTx);
    Msize = P1->size(MAL, ADLTx);
    
    FFsize = Fsize*(Fsize-1)/2;
    MMsize = Msize*(Msize-1)/2;
    FMsize = Fsize*Msize;
    
    grand_mean = 0;
    
    if(Fsize != 0 || Msize != 0) {
      
      if(Fsize != 0) {
        mean = 0;
        for(j = 0; j < Fsize-1;++j)
          for(k = j+1; k < Fsize;++k)
            mean += Coancestry(P1->get(FEM, ADLTx, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               P1->get(FEM, ADLTx, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               nb_locus);
        Theta_FF += mean/FFsize;
        grand_mean += mean;
      }
            
      if(Msize != 0) {
        mean = 0;
        for(j = 0; j < Msize-1;++j)
          for(k = j+1; k < Msize;++k)
            mean += Coancestry(P1->get(MAL, ADLTx, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               P1->get(MAL, ADLTx, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               nb_locus);
        Theta_MM += mean/MMsize;
        grand_mean += mean;
      }
      
      
      if(Fsize != 0 && Msize != 0) {
        mean = 0;
        for(j = 0; j < Fsize;++j)
          for(k = 0; k < Msize;++k)
            mean += Coancestry(P1->get(FEM, ADLTx, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               P1->get(MAL, ADLTx, k)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               nb_locus);
        Theta_FM += mean/FMsize;
        grand_mean += mean;
      }
      _mean_theta += grand_mean / (FFsize + MMsize + FMsize);
    }
  }
  
  _mean_theta /= patchNbr;
  Theta_FF /= patchNbr;
  Theta_MM /= patchNbr;
  Theta_FM /= patchNbr;
}
// ------------------------------------------------------------------------------

//                          kinship analysis
// ----------------------------------------------------------------------------------------
// setSibsStats
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setSibStats()
{
  unsigned int i,j,k;
  unsigned int patchNbr = this->_pop->getPatchNbr(), Fsize, Msize;
  Patch* current_patch;
  Individual *I1, *I2;
  
  for(i = 0;i < 4; ++i) {
    _sib_prop[i] = 0.0;
    _sib_coa[i] = 0.0;
  }
  
  for( i = 0; i < patchNbr; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    if ( (Fsize = current_patch->size(FEM, OFFSx)) != 0) {
      
      for(j = 0; j < Fsize -1; ++j) {
        
        I1 = current_patch->get(FEM, OFFSx, j);
        
        for(k = j+1; k < Fsize; ++k) {
          
          I2 = current_patch->get(FEM, OFFSx, k);
          
          setSibCoa(I1, I2);
        }
      }
    }
    
    if ( (Msize = current_patch->size(MAL, OFFSx)) != 0) {
      
      for(j = 0; j < Msize -1; ++j) {
        
        I1 = current_patch->get(MAL, OFFSx, j);
        
        for(k = j+1; k < Msize; ++k) {
          
          I2 = current_patch->get(MAL, OFFSx, k);
          
          setSibCoa(I1, I2);
        }
      }
    }
    
    //male-female
    for(j = 0; j < Msize; ++j) {
      
      I1 = current_patch->get(MAL, OFFSx, j);
      
      for(k = 0; k < Fsize; ++k) {
        
        I2 = current_patch->get(FEM, OFFSx, k);
        
        setSibCoa(I1, I2);
      }
    }
    
  }
  
  double tot = _sib_prop[0] + _sib_prop[1] + _sib_prop[2] + _sib_prop[3];
  
  for(i = 0 ; i < 4; ++i) {
    _sib_coa[i] = ( (_sib_prop[i] != 0) ? _sib_coa[i]/_sib_prop[i] : nanf("NULL"));
    _sib_prop[i] /= tot;
  }
}
// ----------------------------------------------------------------------------------------
// setSibCoa
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setSibCoa(Individual *I1, Individual *I2)
{
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  //non sibs [3]
  if((I1->getMotherID() != I2->getMotherID()) && (I1->getFatherID() != I2->getFatherID())) {
    _sib_prop[3]++;
    _sib_coa[3] += Coancestry(I1->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              I2->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              nb_locus);
  }
  //maternal half sibs [2]
  else if((I1->getMotherID() == I2->getMotherID()) && (I1->getFatherID() != I2->getFatherID())) {
    _sib_prop[2]++;
    _sib_coa[2] += Coancestry(I1->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              I2->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              nb_locus);
  }
  //paternal half sibs [1]
  else if((I1->getMotherID() != I2->getMotherID()) && (I1->getFatherID() == I2->getFatherID())) {
    _sib_prop[1]++;
    _sib_coa[1] += Coancestry(I1->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              I2->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              nb_locus);
  }
  //full sibs [0]
  else if((I1->getMotherID() == I2->getMotherID()) && (I1->getFatherID() == I2->getFatherID())) {
    _sib_prop[0]++;
    _sib_coa[0] += Coancestry(I1->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              I2->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                              nb_locus);
  }
}  
