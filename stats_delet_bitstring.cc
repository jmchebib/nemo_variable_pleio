/** $Id: stats_delet_bitstring.cc,v 1.7 2015-07-13 08:52:58 fred Exp $
*
*  @file stats_delet_bitstring.cc
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
*  @author fred
 **/

#include <iostream>
#include <cmath>
#include "ttdeletmutations_bitstring.h"
#include "output.h"

// ------------------------------------------------------------------------------

//                        Inbreeding analysis

// ----------------------------------------------------------------------------------------
// getMeanFecWithPatchMate
// ----------------------------------------------------------------------------------------
double TTDeletMutBitstrSH::getMeanFecWithPatchMate(bool HOME)
{
  double m_fec = 0;
  unsigned int nbfem = 0;
  Patch* patch;
  Individual* fem;
  
  if(HOME) {
    for(unsigned int i = 0, pnb = _pop->getPatchNbr(); i < pnb; ++i) {
      patch = _pop->getPatch(i);
      for(unsigned int j = 0, size = patch->size(FEM, ADLTx); j < size; ++j) {
        fem = patch->get(FEM, ADLTx, j);
        m_fec += fem->getFecWithHomePatchMate();
        nbfem += (fem->getLocalMatings() != 0);
      }
    }
  } else {
    for(unsigned int i = 0, pnb = _pop->getPatchNbr(); i < pnb; ++i) {
      patch = _pop->getPatch(i);
      for(unsigned int j = 0, size = _pop->getPatch(i)->size(FEM, ADLTx); j < size; ++j) {
        fem = patch->get(FEM, ADLTx, j);
        m_fec += fem->getFecWithOtherPatchMate();
        nbfem += (fem->getMatings(0) != 0);
      }
    }
  }
  return (nbfem != 0 ? (double) m_fec/nbfem : 0.0);
}
// ----------------------------------------------------------------------------------------
// getHeterosis
// ----------------------------------------------------------------------------------------
double TTDeletMutBitstrSH::getHeterosis()
{
  double Btheta = getMeanFecWithPatchMate(true);
  double Balpha = getMeanFecWithPatchMate(false);

  return (Balpha != 0 ? (1 - (Btheta/Balpha)): nanf("NULL"));
}
// ------------------------------------------------------------------------------

//                        Viability analysis

// ----------------------------------------------------------------------------------------
// setViability
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setViability(age_idx agex)
{
  Individual *ind;
  unsigned int ped_class;
  Patch *patch;
  
  for(unsigned i = 0; i < 5; i++) {
    _viability[i] = 0;
    _SibProps[i] = 0;
  }
  
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
    
    patch = _pop->getPatch(i);
    
    for(unsigned int j = 0, size = patch->size(MAL, agex); j < size; ++j) {
      
      ind = patch->get(MAL, agex, j);
      
      ped_class = ind->getPedigreeClass();
      
      _viability[ped_class] += *(double*)ind->getTrait(_SHLinkedTraitIndex)->getValue();
      
      _SibProps[ped_class]++;
    }
    
    for(unsigned int j = 0, size = patch->size(FEM, agex); j < size; ++j) {
      
      ind = patch->get(FEM, agex, j);
      
      ped_class = ind->getPedigreeClass();
      
      _viability[ped_class] += *(double*)ind->getTrait(_SHLinkedTraitIndex)->getValue();
      
      _SibProps[ped_class]++;
    }
    
  }
  
  double tot_ind = _SibProps[0] + _SibProps[1] + _SibProps[2] + _SibProps[3] + _SibProps[4];
  _meanViab = _viability[0] + _viability[1] + _viability[2] + _viability[3] + _viability[4];
  _meanViab /= tot_ind;
  
  for(unsigned i = 0; i < 5; i++)
    _viability[i] = (_SibProps[i] != 0 ? _viability[i] / _SibProps[i] : nanf("NULL"));
  
  for(unsigned i = 0; i < 5; i++)
    _SibProps[i] /= tot_ind;

}
// ----------------------------------------------------------------------------------------
// getMeanViability
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setMeanViability(age_idx agex)
{
  unsigned int ind_cnt = 0;
  Patch *patch;
  _meanViab = 0;
  for(unsigned int i = 0; i<_pop->getPatchNbr(); ++i) {
    patch = _pop->getPatch(i);
    for(unsigned int j = 0, size = patch->size(FEM, agex); j < size; ++j) {
      _meanViab += *(double*)patch->get(FEM, ADLTx, j)->getTraitValue(_SHLinkedTraitIndex);
      ind_cnt++;
    }
    for(unsigned int j = 0, size = patch->size(MAL, agex); j < size; ++j) {
      _meanViab += *(double*)patch->get(MAL, ADLTx, j)->getTraitValue(_SHLinkedTraitIndex);
      ind_cnt++;
    }
  }
  
  _meanViab = ( ind_cnt != 0 ? _meanViab/ind_cnt : nanf("NULL"));
}
// ------------------------------------------------------------------------------

//                      Mutation load analysis

// ----------------------------------------------------------------------------------------
// getLoad
// ----------------------------------------------------------------------------------------
double TTDeletMutBitstrSH::getLoad()
{
  unsigned int nb_patch=0;
  double mean_load=0;

  for(unsigned int i = 0; i<_pop->getPatchNbr(); i++){
    
    if(_pop->getPatch(i)->size(FEM, ADLTx) != 0) {

      mean_load += getPatchLoad(i);

      nb_patch++;
    }
  }

 return (nb_patch != 0 ? mean_load/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getPatchLoad
// ----------------------------------------------------------------------------------------
double TTDeletMutBitstrSH::getPatchLoad(unsigned int i)
{
  Individual *ind;
  double fec, matings;
  double rel_fec, max_fec = 0, mean_patch_fec = 0;
  
  for(unsigned int j = 0, size = _pop->getPatch(i)->size(FEM, ADLTx); j < size; j++) {
    ind = _pop->getPatch(i)->get(FEM, ADLTx, j);
    fec = (double)ind->getTotRealizedFecundity();
    matings = (double)ind->getTotMatings();
    rel_fec = (matings != 0 ? fec/matings : 0.0);
    mean_patch_fec += rel_fec;
    if(rel_fec > max_fec) max_fec = rel_fec;
  }
  
  mean_patch_fec /= _pop->getPatch(i)->size(FEM, ADLTx);
  
  return ( max_fec != 0 ? (max_fec - mean_patch_fec)/max_fec : nanf("NULL"));
}
// ------------------------------------------------------------------------------

//                    Deleterious mutations analysis

// ----------------------------------------------------------------------------------------
// TTDeletMutBitstrSH::setOffsprgDeletStats
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setDeletStats(age_t AGE)
{
  unsigned int i, j, k, nb_ind = 0, nb_patch = 0;
  unsigned int nb_locus = _SHLinkedTrait->get_nb_locus();
  unsigned int psize=0;
  TTDeletMutations_bitstring* trait;
  Patch* current_patch;
  age_idx agex = (AGE == ADULTS ? ADLTx : OFFSx);

  _deletHtzLoci = 0;
  _deletHmzLoci = 0;
  _deletAllCount = 0;
  _Hs = 0;
  _Ht = 0;
  _fixLocPerPatch = 0;
  _segrLocPerPatch = 0;

  if(_deletFreqTable != NULL)
    delete [] _deletFreqTable;

  _deletFreqTable = new double [nb_locus];

  double patch_freqbylocus[nb_locus];

  for(i = 0; i < nb_locus; ++i)
    _deletFreqTable[i] = 0;


  for(i = 0; i < _pop->getPatchNbr(); ++i) {

    current_patch = _pop->getPatch(i);

    if( (psize = current_patch->size(AGE)) == 0) continue;

    nb_patch++;

    nb_ind += psize;

    for(j = 0; j < nb_locus; ++j)
      patch_freqbylocus[j] = 0;
    
    //Females
    for(j = 0; j < current_patch->size(FEM, agex); ++j) {

      trait = dynamic_cast< TTDeletMutations_bitstring* > ( current_patch->get( FEM, agex, j )->getTrait( _SHLinkedTraitIndex ) );

      _deletHmzLoci += trait->get_nb_hmz_mutations();

      _deletHtzLoci += trait->get_nb_htz_mutations();

      for(k = 0; k < nb_locus; ++k)
        patch_freqbylocus[k] += trait->get_nb_mut_atLocus(k);

    }//end for females
    
    //Males
    for(j = 0; j < current_patch->size(MAL, agex); ++j) {
      
      trait = dynamic_cast< TTDeletMutations_bitstring* > ( current_patch->get( MAL, agex, j )->getTrait( _SHLinkedTraitIndex ) );
      
      _deletHmzLoci += trait->get_nb_hmz_mutations();
      
      _deletHtzLoci += trait->get_nb_htz_mutations();
      
      for(k = 0; k < nb_locus; ++k)
        patch_freqbylocus[k] += trait->get_nb_mut_atLocus(k);
      
    }//end for males
    
    double dipl_size = psize*2.0;
    for(j = 0; j < nb_locus; ++j) {
      //aggregate for the pop freq by locus:
      _deletFreqTable[j] += patch_freqbylocus[j];
      //aggregate mut nbr for the population tot mut nbr:
      _deletAllCount += patch_freqbylocus[j];
      //get the local Patch mut freq by locus:
      patch_freqbylocus[j] /= dipl_size;
      //compute the expected local htz = 2pq:
      _Hs += patch_freqbylocus[j] * (1.0 - patch_freqbylocus[j]);
    }

    for(j = 0; j < nb_locus; ++j) {
      //segregating and fixed loci in the Patch:
      if(patch_freqbylocus[j] == 1.0)
        _fixLocPerPatch++;
      else if(patch_freqbylocus[j] != 0)
        _segrLocPerPatch++;
    }

  }//end for Patch

  double dipl_size = nb_ind*2.0;
  
  _freq = 0;
  _fixloc = _segrloc = 0;
  
  for(i = 0; i < nb_locus; ++i){
    
    _deletFreqTable[i] /= dipl_size;
    
    _Ht += _deletFreqTable[i] * (1.0 - _deletFreqTable[i]);//2pq: expected htz
      
    _freq += _deletFreqTable[i];
    
    _fixloc += (_deletFreqTable[i] == 1);
    
    _segrloc += (_deletFreqTable[i] < 1 && _deletFreqTable[i] != 0);
  }
  
  //factorization not done in previous loops:
  _Hs *= 2;
  _Ht *= 2;
  
  _deletAllCount /= nb_ind;
  _freq /= nb_locus;
  _Hs /= nb_locus * nb_patch;
  _Ht /= nb_locus;
  _Ho = _deletHtzLoci / (nb_ind * nb_locus);
  _Hmz = _deletHmzLoci/ (nb_ind * nb_locus);
  _fixLocPerPatch /= nb_patch;
  _segrLocPerPatch /= nb_patch;

  setLethalEquivalents(AGE);
  
  setFst(AGE);
  
  delete [] _deletFreqTable;
  _deletFreqTable = 0;
}
// ----------------------------------------------------------------------------------------
// setLetahlEquivalents()
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setLethalEquivalents(age_t AGE)
{
  _letheq = 0;
  
  if(_deletFreqTable == 0) 
    fatal("allele frequency table not set when computing lethal equivalents\n");
  
  if(_isContinuousEffect) {
    float* s = _SHLinkedTrait->get_s_continous();
    for(int i = 0, nloc = _SHLinkedTrait->get_nb_locus(); i < nloc; i++)
      _letheq += _deletFreqTable[i] * s[i];
  } else
    _letheq = ((_deletHtzLoci/2 + _deletHmzLoci) * _SHLinkedTrait->get_strength())
              / _pop->size(AGE);
}
// ----------------------------------------------------------------------------------------
// TTDeletMutBitstrSH::getFst
// ----------------------------------------------------------------------------------------
void TTDeletMutBitstrSH::setFst(age_t AGE)
{
  unsigned i, nbpatch = 0;
  double Hsnei, Htnei, harmonic = 0, size;
  Patch* current_patch;

  for (i = 0; i < _pop->getPatchNbr(); ++i){

    current_patch = _pop->getPatch(i);

    if( (size = (double)current_patch->size(AGE)) != 0){
      nbpatch++;
      harmonic += 1.0/size;
    }
  }

  harmonic = (double)nbpatch/harmonic;

  Hsnei = harmonic/(harmonic-1.0)*(_Hs - ( _Ho/(2.0*harmonic) ) );
  Htnei = _Ht + ( Hsnei/(harmonic*nbpatch) ) - ( _Ho/(2.0*harmonic*nbpatch) );

  _Fst = 1.0 - (Hsnei/Htnei);
}

