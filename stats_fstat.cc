/** $Id: stats_fstat.cc,v 1.14 2016-02-03 14:12:18 fred Exp $
*  
*  @file stats_fstat.cc
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
#include "metapop.h"
#include "ttneutralgenes.h"
#include <assert.h>

// ------------------------------------------------------------------------------

//                          fstat's analysis

// ----------------------------------------------------------------------------------------
// allocateTable
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::allocateTables (unsigned int loci, unsigned int all)
{
  unsigned int nb_patch = _pop->getPatchNbr();
  unsigned int **sizes;
  
  sizes = new unsigned int * [nb_patch];
  
  for(unsigned int i = 0; i < nb_patch; ++i) {
    sizes[i] = new unsigned int [loci];
    for(unsigned int j = 0; j < loci; ++j)
      sizes[i][j] = all;
  }
  
  _alleleCountTable.allocate(nb_patch, loci, sizes);
  
  _alleleFreqTable.allocate(nb_patch, loci, sizes);
  
  _globalAlleleFreq.reset(loci, all);
  
  for(unsigned int i = 0; i < nb_patch; ++i)
    delete [] sizes[i];
  delete [] sizes;
  
  //reset the time info, to force recalculation after allocate
  _table_set_age = NONE;
  _table_set_gen = 0;
  _table_set_repl = 0;
}
// ----------------------------------------------------------------------------------------
// setAlleleTables
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setAlleleTables(age_t AGE)
{  
  
  if(_table_set_age == AGE 
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate())
    return;
  
  unsigned char** genes;
  unsigned int nb_locus = _SHLinkedTrait->get_locus_num(), nb_allele = _SHLinkedTrait->get_allele_num();
  unsigned int patch_size, patchNbr = _pop->getPatchNbr();
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  Patch* crnt_patch;
  
  //check if the population size has changed:
  if(_alleleCountTable.getNumGroups() != patchNbr)
    allocateTables(_SHLinkedTrait->get_locus_num(), _SHLinkedTrait->get_allele_num());
  
  _alleleCountTable.init(0);

  //counting the copies of each allele present in each patch:
  for(unsigned int k = 0; k < patchNbr; ++k) {

    crnt_patch = _pop->getPatch(k);

    for(unsigned int i = 0, size = crnt_patch->size(FEM, age_pos); i < size; ++i) {

      genes = (unsigned char**)crnt_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int j = 0; j < nb_locus; ++j)
        _alleleCountTable.increment(k, j, genes[0][j]);
      
      for (unsigned int j = 0; j < nb_locus; ++j)
        _alleleCountTable.increment(k, j, genes[1][j]);
      
    } 
    
    for(unsigned int i = 0, size = crnt_patch->size(MAL, age_pos); i < size; ++i) {
      
      genes = (unsigned char**)crnt_patch->get(MAL, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int j = 0; j < nb_locus; ++j)
        _alleleCountTable.increment(k, j, genes[0][j]);
      
      for (unsigned int j = 0; j < nb_locus; ++j)
        _alleleCountTable.increment(k, j, genes[1][j]);
      
    }
  }
  
  //allelic frequencies:
  for (unsigned int i = 0; i < patchNbr; ++i) {
    
    patch_size = _pop->getPatch(i)->size(age_pos) * 2;
    
    if (patch_size) {
      
    for (unsigned int l = 0; l < nb_locus; ++l)
      for (unsigned int u = 0; u < nb_allele; ++u)
        _alleleFreqTable.set(i, l, u, (double)_alleleCountTable.get(i, l, u) / (double)patch_size);
    
    } else {
      for (unsigned int l = 0; l < nb_locus; ++l)
        for (unsigned int u = 0; u < nb_allele; ++u)
          _alleleFreqTable.set(i, l, u, 0 );
    }

  }

  //population-wide allelic frequencies:  
  unsigned int tot_size = _pop->size(AGE) * 2;

  _globalAlleleFreq.assign(0);
  
  for(unsigned int i = 0; i < patchNbr; i++)
    for (unsigned int l = 0; l < nb_locus; ++l)
      for (unsigned int u = 0; u < nb_allele; ++u)
        _globalAlleleFreq.plus(l, u, _alleleCountTable.get(i, l, u));
  
  for (unsigned int l = 0; l < nb_locus; ++l)
    for (unsigned int u = 0; u < nb_allele; ++u)
      _globalAlleleFreq.divide(l, u, tot_size);
  
  _table_set_age = AGE;
  
  _table_set_gen = _pop->getCurrentGeneration();
  
  _table_set_repl = _pop->getCurrentReplicate();
}
// ----------------------------------------------------------------------------------------
// setHeteroTable
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setHeteroTable (age_t AGE)
{
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
  Patch* patch;
  unsigned int isHetero;
  unsigned char** genes;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
    
  unsigned int **sizes;
  
  sizes = new unsigned int * [patchNbr];
  
  for(unsigned int i = 0; i < patchNbr; ++i) {
    sizes[i] = new unsigned int [nb_locus];
    for(unsigned int j = 0; j < nb_locus; ++j)
      sizes[i][j] = nb_allele;
  }
  
  _heteroTable.update(patchNbr, nb_locus, sizes);
  _heteroTable.init(0);
  
  
  for (unsigned int i = 0; i < patchNbr; ++i) {
    
    patch = _pop->getPatch(i);
    
    if(!patch->size(age_pos)) continue;
    
    for(unsigned int j = 0, size = patch->size(FEM, age_pos); j < size; ++j) {
      
      genes = (unsigned char**)patch->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < nb_locus; ++l) {
        isHetero = (genes[0][l] != genes[1][l]);
        _heteroTable.plus(i, l, genes[0][l], isHetero);
        _heteroTable.plus(i, l, genes[1][l], isHetero);
      }
    }
    
    for(unsigned int j = 0, size = patch->size(MAL, age_pos); j < size; ++j) {
      
      genes = (unsigned char**)patch->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < nb_locus; ++l) {
        isHetero = (genes[0][l] != genes[1][l]);
        _heteroTable.plus(i, l, genes[0][l], isHetero);
        _heteroTable.plus(i, l, genes[1][l], isHetero);
      }
    }
  }  
  for(unsigned int i = 0; i < patchNbr; ++i) delete [] sizes[i];
  delete [] sizes;

}
// ----------------------------------------------------------------------------------------
// setHeteroTable
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setHeterozygosity (age_t AGE)
{
  unsigned int patchNbr  = this->_pop->getPatchNbr();
  unsigned int nb_locus  = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
  Patch* patch;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  setHeteroTable(AGE);
  
  for (unsigned int i = 0, size; i < patchNbr; ++i) {
    
    patch = _pop->getPatch(i);
    size = patch->size(age_pos);
    
    if(!size) continue;

    for (unsigned int l = 0; l < nb_locus; ++l) {
      for (unsigned int u = 0; u < nb_allele; ++u) {
        _heteroTable.divide(i, l , u, size);
      }
    }
  }
  
  
}
// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstat(age_t AGE)
{
  double harmonic = 0, nbpatch = 0, nbind;
  unsigned int patchNbr = _pop->getPatchNbr();
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  
  //check if the population size has changed:
  if(_alleleCountTable.getNumGroups() != patchNbr)
    allocateTables(_SHLinkedTrait->get_locus_num(), _SHLinkedTrait->get_allele_num());
  
  setAlleleTables(AGE);
  
  setLociDivCounter(AGE);

  //harmonic mean of patch sizes:
  for (unsigned int i = 0; i < patchNbr; ++i){
    
    nbind = _pop->size(AGE, i);
    if(nbind != 0){
      nbpatch++;
      harmonic += 1.0/nbind;
    }
  }
  
  harmonic = nbpatch / harmonic;
  
  _ho = setHo(age_pos);
  _hs = setHs(age_pos);
  _ht = setHt(age_pos);
  //Nei's corrections:
  _hsnei = (nbpatch != 0 ? harmonic/(harmonic-1.0)*(_hs-(_ho/(2.0*harmonic))) : nanf("NULL") );
  _htnei = (nbpatch != 0 ? _ht + (_hsnei/(harmonic*nbpatch))
            -(_ho/(2.0*harmonic*nbpatch)) : nanf("NULL") );
  
  _fis = ( _hsnei ? 1.0-(_ho/_hsnei) : nanf("NULL") );
  _fit = ( _htnei ? 1.0-(_ho/_htnei) : nanf("NULL") );
  _fst = ( _htnei ? 1.0-(_hsnei/_htnei) : nanf("NULL") );
}
// ----------------------------------------------------------------------------------------
// setLociDivCounter
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setLociDivCounter (age_t AGE)
{
  unsigned int i, j, k;
  unsigned int nbLoc = _SHLinkedTrait->get_locus_num(), nbAll = _SHLinkedTrait->get_allele_num();
  unsigned int nbpatch = 0, patchNbr = _pop->getPatchNbr();
  double patch_mean, pop_mean = 0;
  bool **pop_div;
  
  //number of alleles per locus, Patch and pop counters:
  pop_div = new bool * [nbLoc];
  
  for(i = 0; i < nbLoc; ++i) {
    pop_div[i] = new bool [nbAll];
    for(j = 0; j < nbAll;++j)
      pop_div[i][j] = 0;
  }
  
  for(k = 0; k < patchNbr; ++k) 
    nbpatch += (_pop->size(AGE, k) != 0);
  
  for(k = 0; k < patchNbr; ++k) {
	
    patch_mean = 0;
	
    for(i = 0; i < nbLoc; ++i)	  
	  for(j = 0; j < nbAll;++j) {
		patch_mean += (_alleleCountTable.get(k,i,j) != 0);
        pop_div[i][j] |= (_alleleCountTable.get(k,i,j) != 0);
      }
        //add mean nb of alleles per locus for Patch k to the pop mean
        pop_mean += patch_mean/nbLoc;
  }
  
  _nb_all_local = (nbpatch ? pop_mean/nbpatch : nanf("NULL"));

  _nb_all_global = 0;
  
  for(i = 0; i < nbLoc; ++i)	  
    for(j = 0; j < nbAll;++j)
      _nb_all_global += pop_div[i][j];
  
  _nb_all_global /= nbLoc;
  
  for(i = 0; i < nbLoc; ++i)
    delete [] pop_div[i];
  delete [] pop_div;
      
  //number of fixed loci, local and global counters:
  _fix_loc_local = 0;
  
  for (k = 0; k < patchNbr; ++k)
    for (i = 0; i < nbLoc; ++i)
      for (j = 0; j < nbAll; ++j)        
        _fix_loc_local += ( _alleleFreqTable.get(k, i, j) == 1 );
  
  _fix_loc_local /= nbpatch;
  
  _fix_loc_global = 0;
      
  //globally:  
  for (i = 0; i < nbLoc; ++i)
    for (j = 0; j < nbAll; ++j)  
      _fix_loc_global += ( _globalAlleleFreq.get(i, j) == 1 );

}
// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::setHo (age_idx age_pos)
{
  unsigned int nloc = _SHLinkedTrait->get_locus_num(), nbpatch = _pop->getPatchNbr(), psize;
  double indloc = 0, hetero = 0;
  unsigned char** genes;
  Patch* current_patch;
  
  for (unsigned int i = 0; i < nbpatch; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    psize = current_patch->size(FEM, age_pos);
    
    indloc += psize;
    
    for(unsigned int j = 0; j < psize; ++j) {
      
      genes = (unsigned char**)current_patch->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int k = 0; k < nloc; ++k)  hetero += (genes[0][k] != genes[1][k]);
    }
  }
  //doing some unrolling, loops for the males here:
  for (unsigned int i = 0; i < nbpatch; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    psize = current_patch->size(MAL, age_pos);
    
    indloc += psize;
    
    for(unsigned int j = 0; j < psize; ++j) {
      
      genes = (unsigned char**)current_patch->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int k = 0; k < nloc; ++k)  hetero += (genes[0][k] != genes[1][k]);
    }
  }
  
  indloc *= nloc;
  
  return (indloc != 0 ? hetero/indloc : 0.0);
}

// ----------------------------------------------------------------------------------------
// setHs
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::setHs (age_idx age_pos)
{
  unsigned int i, k, x, nb_patch=0;
  double hs = 0;
  unsigned int nbLoc = _SHLinkedTrait->get_locus_num(), nbAll = _SHLinkedTrait->get_allele_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  deque<double>genehet(nbLoc,1);
  deque<double>hspop(patchNbr,0);
  double freq;
  Patch* current_patch;
  
  //compute the expected heterozygosity for each patch and return average:
  for (i = 0; i < patchNbr; ++i) {
    
    current_patch = _pop->getPatch(i);
		
    genehet.assign(nbLoc, 1);
    
    if(current_patch->size(age_pos) != 0) {
      
      nb_patch++;
      
      for (k = 0; k < nbLoc; ++k) {
        
        for (x = 0; x < nbAll; ++x) {
          
          freq = _alleleFreqTable.get(i, k, x);
          
          freq *= freq; //squared frequencies (expected _homozygosity)
          
          genehet[k] -= freq; //1 - sum of p2 = expected heterozygosity
        }
        //expected heterozygosity:
        hspop[i] += genehet[k];
      }
    }//end_if size!=0
    
    hs += hspop[i];
  }//end patch loop
  
  return (nb_patch !=0 ? hs/(nbLoc*nb_patch) : 0.0);
}
// ----------------------------------------------------------------------------------------
//  setHt
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::setHt (age_idx age_pos)
{
  double ht = 0;
  unsigned int nbLoc = _SHLinkedTrait->get_locus_num(), nbAll = _SHLinkedTrait->get_allele_num();
  deque<double> genehet(nbLoc, 1);
  double freq;
  
  //get global allele frequencies per locus
  for (unsigned int l = 0; l < nbLoc; ++l) {
    for (unsigned int u = 0; u < nbAll; ++u) {
      
      freq = _globalAlleleFreq.get(l, u);
      
      freq *= freq; //squared frequencies
      
	  genehet[l] -= freq;  //1 - sum of p2 = expected heterozygosity
    }
    
    ht += genehet[l];
  }
  
  return (ht/nbLoc);
}

// ----------------------------------------------------------------------------------------
// setFstat2
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstat2(age_t AGE)
{
  double harmonic = 0, nbpatch = 0, nbind;
  unsigned int patchNbr = _pop->getPatchNbr();
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int nloc = _SHLinkedTrait->get_locus_num();
  
  //check if the population size has changed:
  if(_alleleCountTable.getNumGroups() != patchNbr)
    allocateTables(_SHLinkedTrait->get_locus_num(), _SHLinkedTrait->get_allele_num());
  
  setAlleleTables(AGE);
  
  setLociDivCounter(AGE);
  
  //harmonic mean of patch sizes:
  for (unsigned int i = 0; i < patchNbr; ++i){
    
    nbind = _pop->size(AGE, i);
    if(nbind != 0){
      nbpatch++;
      harmonic += 1.0/nbind;
    }
  }
  
  harmonic = nbpatch / harmonic;
  
  deque<double> ho = setHo2(age_pos);
  deque<double> hs = setHs2(age_pos);
  deque<double> ht = setHt2(age_pos);
  
  _fis = _fit = _fst = 0;
  
  double hsnei, htnei;
  
  for (unsigned int l = 0; l < nloc; ++l) {
  
  //Nei's corrections:
  hsnei = (nbpatch != 0 ? harmonic/(harmonic-1.0)*(hs[l]-(ho[l]/(2.0*harmonic))) : nanf("NULL") );
  htnei = (nbpatch != 0 ? ht[l] + (hsnei/(harmonic*nbpatch))
            -(ho[l]/(2.0*harmonic*nbpatch)) : 0 );
  _hsnei += hsnei;
  _htnei += htnei;
//    _hsnei += hs[l];
//    _htnei += ht[l];
    _ho += ho[l];
//  _fis = ( _hsnei ? 1.0-(_ho/_hsnei) : nanf("NULL") );
//  _fit = ( _htnei ? 1.0-(_ho/_htnei) : nanf("NULL") );
    if(ht[l] == 0) nloc--; //would mean locus is fixed
    _fst += ( ht[l] != 0 ? 1.0-(hsnei/htnei) : 0);

//    if(ht[l] == 0) nloc--;
//    
//    _fst += (ht[l] != 0 ? 1 - (hs[l]/ht[l]) : 0);
  }
  _hsnei /= nloc;
  _htnei /= nloc;
  _ho /= nloc;
  _fst /= nloc;
}
// ----------------------------------------------------------------------------------------
// setHo2
// ----------------------------------------------------------------------------------------
deque<double> TTNeutralGenesSH::setHo2 (age_idx age_pos)
{
  unsigned int nloc = _SHLinkedTrait->get_locus_num(), nbpatch = _pop->getPatchNbr();
  unsigned int msize, fsize, psize;
  deque<double> hetero(nloc,0);
  unsigned char** genes;
  Patch* current_patch;
    
  for (unsigned int i = 0; i < nbpatch; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    fsize = current_patch->size(FEM, age_pos);
    
    for(unsigned int j = 0; j < fsize; ++j) {
      
      genes = (unsigned char**)current_patch->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int k = 0; k < nloc; ++k)  hetero[k] += (genes[0][k] != genes[1][k]);
    }
        
    msize = current_patch->size(MAL, age_pos);
    
    for(unsigned int j = 0; j < msize; ++j) {
      
      genes = (unsigned char**)current_patch->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int k = 0; k < nloc; ++k)  hetero[k] += (genes[0][k] != genes[1][k]);
    }
    
  }

  psize = _pop->size(age_pos);
  
  if(psize != 0) 
    for (unsigned int k = 0; k < nloc; ++k)  
      hetero[k] /= psize;
  
  return hetero;
}

// ----------------------------------------------------------------------------------------
// setHs
// ----------------------------------------------------------------------------------------
deque<double> TTNeutralGenesSH::setHs2 (age_idx age_pos)
{
  unsigned int i, k, x, nb_patch = 0;
  unsigned int nbLoc = _SHLinkedTrait->get_locus_num(), nbAll = _SHLinkedTrait->get_allele_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  deque<double>genehet(nbLoc,1);
  deque<double>hsloc(nbLoc,0);
  double freq;
  Patch* current_patch;
  
  for (i = 0; i < patchNbr; ++i) {
    if(_pop->size(age_pos, i) != 0)
      nb_patch++;
  }
  
  //compute the expected heterozygosity for each patch and return average:
  for (k = 0; k < nbLoc; ++k) {    
    
    for (i = 0; i < patchNbr; ++i) {
      
      current_patch = _pop->getPatch(i);
      
      genehet[k] = 1;
      
      if(current_patch->size(age_pos) != 0) {
        
        for (x = 0; x < nbAll; ++x) {
          
          freq = _alleleFreqTable.get(i, k, x);
          
          freq *= freq; //squared frequencies (expected _homozygosity)
          
          genehet[k] -= freq; //1 - sum of p2 = expected heterozygosity
        }
        
      }//end_if size!=0
      
      hsloc[k] += genehet[k];
      
    }//end patch loop
    
    if(nb_patch !=0) hsloc[k] /= nb_patch;
    
  }//end loci loop
  
  return  hsloc;
}
// ----------------------------------------------------------------------------------------
//  setHt
// ----------------------------------------------------------------------------------------
deque<double> TTNeutralGenesSH::setHt2 (age_idx age_pos)
{
  unsigned int nbLoc = _SHLinkedTrait->get_locus_num(), nbAll = _SHLinkedTrait->get_allele_num();
  deque<double> genehet(nbLoc, 1);
  double freq;
  
  //get global allele frequencies per locus
  for (unsigned int l = 0; l < nbLoc; ++l) {
  
    for (unsigned int u = 0; u < nbAll; ++u) {
      
      freq = _globalAlleleFreq.get(l, u);
      
      freq *= freq; //squared frequencies
      
      genehet[l] -= freq;  //1 - sum of p2 = expected heterozygosity
    }
  }
  
  return genehet;
}


// ----------------------------------------------------------------------------------------
// setFstMatrix
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstMatrix(age_t AGE, unsigned char dim)
{
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();

  setAlleleTables(AGE);
  
  if(_fst_matrix == NULL)
    
    _fst_matrix = new TMatrix(patchNbr, patchNbr);
  
  else if( _fst_matrix->length() != patchNbr * patchNbr)
    
    _fst_matrix->reset(patchNbr, patchNbr);
  
  _fst_matrix->assign(nanf("NULL"));
  
  //init
  double *pop_weights = new double[patchNbr];
  double *pop_sizes = new double[patchNbr];
  double **numerator = new double*[patchNbr];
  for(unsigned int i = 0; i < patchNbr; i++) numerator[i] = new double [patchNbr];
  double tot_size;
  double numerator_W = 0;
  double denominator = 0;
  double sum_weights = 0;
  
  tot_size = _pop->size(AGE) * 2; 
  
  for(unsigned int i = 0; i < patchNbr; ++i) {
    pop_sizes[i] = _pop->size(AGE, i) * 2; 
    pop_weights[i] = pop_sizes[i] - (pop_sizes[i] * pop_sizes[i] / tot_size); //n_ic in Weir & Hill 2002
    sum_weights += pop_weights[i];
    for(unsigned int j = 0; j < patchNbr; j++)
      numerator[i][j] = 0;
  }
  
   double p, pq, var, num;

   for (unsigned int i = 0; i < patchNbr; ++i) {
     
     if( !pop_sizes[i] ) continue;
     
     for (unsigned int l = 0; l < nb_locus; ++l) {
       
       for (unsigned int u = 0; u < nb_allele; ++u) {
         
         p = _alleleFreqTable.get(i, l, u); //p_liu
         
         pq = p * (1 - p);
         
         var = p - _globalAlleleFreq.get(l, u); //(p_liu - pbar_u)^2 
         
         var *= var;
         
         num = pq * pop_sizes[i] / (pop_sizes[i] -1);

         numerator[i][i] += num;
         
         numerator_W += num * pop_sizes[i]; //see equ. 9, Weir & Hill 2002
         
         denominator += pop_sizes[i] * var + pop_weights[i] * pq; //common denominator
         
       } // end for allele
     }// end for locus
   }//end for pop

  for (unsigned int i = 0; i < patchNbr; ++i) {
    if( !pop_sizes[i] ) continue;
    _fst_matrix->set(i, i, 1 - (numerator[i][i] * sum_weights / denominator) );
  }   
   _fst_WH = 1 - ((numerator_W * sum_weights) / (denominator * tot_size)); //equ. 9 Weir & Hill 2002
   
   //pairwise Fst:
   if(dim & 2) {
     double pi, pj;
     for (unsigned int l = 0; l < nb_locus; ++l)
       for (unsigned int u = 0; u < nb_allele; ++u)
         for (unsigned int i = 0; i < patchNbr - 1; ++i) {
           if( !pop_sizes[i] ) continue;
           for (unsigned int j = i + 1; j < patchNbr; ++j) {
             if( !pop_sizes[j] ) continue;
             pi = _alleleFreqTable.get(i, l, u);
             pj = _alleleFreqTable.get(j, l, u);
             numerator[i][j] += pi * (1 - pj) + pj * (1 - pi); //equ. 7 of Weir & Hill 2002
           }
         }
     for (unsigned int i = 0; i < patchNbr - 1; ++i){
       if( !pop_sizes[i] ) continue;
       for (unsigned int j = i + 1; j < patchNbr; ++j){
         if( !pop_sizes[j] ) continue;
         _fst_matrix->set(i, j, 1 - ( (numerator[i][j] * sum_weights) / (2 * denominator)) );
       }
     }
   }
 
  
  delete [] pop_weights;
  delete [] pop_sizes;
  for(unsigned int i = 0; i < patchNbr; i++) delete [] numerator[i];
  delete [] numerator;
}
// ----------------------------------------------------------------------------------------
// setFstatWeirCockerham
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstatWeirCockerham(age_t AGE)
{
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
  //age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  
  setAlleleTables(AGE);
  setHeteroTable(AGE);
  
  //init
  double *pop_sizes = new double [patchNbr];
  double tot_size, inv_ntot;
  double sum_weights = 0;
  double nbar, nc, inv_nbar;
  unsigned int extantPs = 0;
  
  tot_size = _pop->size(AGE);
  
  for(unsigned int i = 0; i < patchNbr; i++) {
    pop_sizes[i] = _pop->size(AGE, i); 
    if(pop_sizes[i]) {
      extantPs++;
      sum_weights += (pop_sizes[i] * pop_sizes[i] / tot_size);
    }
  }
  nbar = tot_size/extantPs;
  nc = (tot_size - sum_weights)/(extantPs-1);
  inv_nbar = 1/(nbar - 1);
  inv_ntot = 1/tot_size;
  
  double var;
  double s2, pbar, hbar;
  double s2_denom = 1.0/((extantPs-1)*nbar),
		  r = (double)(extantPs-1)/extantPs,
		  hbar_factor=(2*nbar-1)/(4*nbar);
  double a = 0, b = 0, c = 0, x;
  
  for (unsigned int l = 0; l < nb_locus; ++l) {
    
    for (unsigned int u = 0; u < nb_allele; ++u) {
      
      s2 = pbar = hbar = 0;
      
      for (unsigned int i = 0; i < patchNbr; ++i) {
        
        var = _alleleFreqTable.get(i, l, u) - _globalAlleleFreq.get(l, u); //(p_liu - pbar_u)^2 
        
        var *= var;
        
        s2 += var * pop_sizes[i];
        
        hbar += _heteroTable.get(i, l, u);
        
      }//end for pop
      
      s2 *= s2_denom;
      pbar = _globalAlleleFreq.get(l, u);
      hbar *= inv_ntot;
      
      x = pbar * (1 - pbar) - r * s2;
      a += s2 - inv_nbar*( x - 0.25 * hbar);
      b += x - hbar_factor * hbar;
      c += hbar;
      
      
    } // end for allele
    
  }// end for locus
    
  a *= nbar/nc;
  b *= nbar/(nbar - 1);
  c *= 0.5;
  
  _fst_WC = a / (a + b + c);
  _fit_WC = (a + b) / (a + b + c);
  _fis_WC = b / (b + c);
  
  delete [] pop_sizes;   
}
// ----------------------------------------------------------------------------------------
// setFstatWeirCockerham_MS
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstatWeirCockerham_MS(age_t AGE)
{ /**Computes W&C F-stats using the Mean Squares approach, similar to the implementation in Hierfstat.
     This code gives the exact same results as the other version.
  */
	unsigned int patchNbr = this->_pop->getPatchNbr();
	unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
	unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
	//age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

	setAlleleTables(AGE);
	setHeteroTable(AGE);

	//init
	double *pop_sizes = new double [patchNbr];
	double tot_size;
	double sum_weights = 0;
	double nc;
	unsigned int extantPs = 0;

	tot_size = _pop->size(AGE);

	for(unsigned int i = 0; i < patchNbr; i++) {
		pop_sizes[i] = _pop->size(AGE, i);
		if(pop_sizes[i]) {
			extantPs++;
			sum_weights += (pop_sizes[i] * pop_sizes[i] / tot_size);
		}
	}

	nc = (tot_size - sum_weights)/(extantPs-1);

//	unsigned int np = extantPs;
	unsigned int npl = extantPs; //all loci typed in all patches

	//p = _alleleFreqTable
	//pb = _globalAlleleFreq

	unsigned int *alploc = new unsigned int [nb_locus];

	unsigned int **alploc_table = new unsigned int* [nb_locus];

	for(unsigned int i = 0; i < nb_locus; ++i)
		alploc_table[i] = new unsigned int[nb_allele];

	unsigned int tot_num_allele = 0;

	for(unsigned int l = 0; l < nb_locus; ++l){

		alploc[l] = 0;

		for(unsigned int cnt,  a = 0; a < nb_allele; ++a) {

			cnt=0;

			for(unsigned int i = 0; i < patchNbr; i++) {

				 cnt += _alleleCountTable.get(i,l,a);

			}
			alploc_table[l][a] = (cnt != 0);
			alploc[l] += (cnt != 0);
		}

		tot_num_allele += alploc[l];
	}

	//n, and nal are given by pop_sizes, same num ind typed at all loci in each patch
	//nc is the same for each locus
	//nt is given by tot_size, same tot num of ind typed for all loci

	//SSG: het/2 for each allele
	double *SSG = new double[tot_num_allele];
	double *SSP = new double[tot_num_allele];
	double *SSi = new double[tot_num_allele];

	unsigned int all_cntr = 0;

	double het, freq, var;

	for(unsigned int l = 0; l < nb_locus; ++l) {

		for(unsigned int a = 0; a < nb_allele & all_cntr < tot_num_allele; ++a) {

			if(alploc_table[l][a] == 0) continue; //do not consider alleles not present in the pop

			SSG[all_cntr] = 0;
			SSi[all_cntr] = 0;
			SSP[all_cntr] = 0;

			for(unsigned int p = 0; p < patchNbr; ++p){

				if(!_pop->size(AGE, p)) continue; //skip empty patches

				het = _heteroTable.get(p, l, a);

				freq = _alleleFreqTable.get(p, l, a);

				var = freq - _globalAlleleFreq.get(l, a); //(p_liu - pbar_u)^2

				var *= var;

				SSG[all_cntr] += het;

				SSi[all_cntr] += 2*pop_sizes[p]*freq*(1-freq) - het/2;

				SSP[all_cntr] += 2*pop_sizes[p]*var;
			}

			all_cntr++;
		}

	}


	assert(all_cntr == tot_num_allele);

	double *MSG = new double[tot_num_allele];
	double *MSP = new double[tot_num_allele];
	double *MSI = new double[tot_num_allele];
	double *sigw = new double[tot_num_allele];
	double *siga = new double[tot_num_allele];
	double *sigb = new double[tot_num_allele];

//	double *FST_pal = new double[tot_num_allele];
//	double *FIS_pal = new double[tot_num_allele];

	double SIGA = 0, SIGB = 0, SIGW = 0;

	for(unsigned int i = 0; i < tot_num_allele; ++i){

		MSG[i] = SSG[i] / (2 * tot_size);
		sigw[i] = MSG[i]; //wasted!

		MSP[i] = SSP[i] / (npl-1);

		MSI[i] = SSi[i]/ (tot_size - npl);

		sigb[i] = 0.5*(MSI[i] - MSG[i]);

		siga[i] = (MSP[i] - MSI[i])/(2*nc);

//		FST_pal[i] = siga[i]/(siga[i]+sigb[i]+sigw[i]);
//		FIS_pal[i] = sigb[i]/(sigb[i]+sigw[i]);

		SIGA += siga[i];
		SIGB += sigb[i];
		SIGW += sigw[i];
	}

	//per locus stats:
	if(_fst_WC_loc) delete [] _fst_WC_loc; _fst_WC_loc = new double [nb_locus];
	if(_fis_WC_loc) delete [] _fis_WC_loc; _fis_WC_loc = new double [nb_locus];
	if(_fit_WC_loc) delete [] _fit_WC_loc; _fit_WC_loc = new double [nb_locus];

	double lsiga, lsigb, lsigw;

//	cout<<"  computing sigma per locus\n";

	for(unsigned int allcntr = 0, i = 0; i < nb_locus; ++i) {

		lsiga = lsigb = lsigw = 0;

		for(unsigned int l = 0; l < alploc[i]; ++l) {

			lsiga += siga[allcntr];
			lsigb += sigb[allcntr];
			lsigw += sigw[allcntr];

			allcntr++;

		}

		_fst_WC_loc[i] = lsiga /(lsiga + lsigb + lsigw);
		_fis_WC_loc[i] = lsigb /(lsigb + lsigw);
		_fit_WC_loc[i] = (lsiga +lsigb) /(lsiga + lsigb + lsigw);

	}


	// Total F-stats
	_fst_WC = SIGA / (SIGA + SIGB + SIGW);
	_fit_WC = (SIGA + SIGB) / (SIGA + SIGB + SIGW);
	_fis_WC = SIGB / (SIGB + SIGW);

	delete[]pop_sizes;
	delete[]alploc;
	for(unsigned int i = 0; i < nb_locus; ++i)
		delete[]alploc_table[i];
	delete[]alploc_table;
	delete[]SSG;
	delete[]SSi;
	delete[]SSP;
	delete[]MSG;
	delete[]MSI;
	delete[]MSP;
	delete[]sigw;
	delete[]siga;
	delete[]sigb;

}
// ----------------------------------------------------------------------------------------
// setFstMatrix
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFst_li(unsigned int N, unsigned int L, double **array)
{
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
  
  assert(N == patchNbr);
  assert(L == nb_locus);
  
  setAlleleTables(ADULTS);
  
  double *pop_weights = new double[patchNbr];
  double *pop_sizes = new double[patchNbr];
  double *denominator = new double[nb_locus];
  double sum_weights = 0;
  
  double tot_size = _pop->size(ADULTS) * 2; //diploids!
  
  for(unsigned int i = 0; i < patchNbr; ++i) {
    pop_sizes[i] = _pop->size(ADULTS, i) * 2; 
    pop_weights[i] = pop_sizes[i] - (pop_sizes[i] * pop_sizes[i] / tot_size); //n_ic in Weir & Hill 2002
    sum_weights += pop_weights[i];
    for(unsigned int j = 0; j < nb_locus; ++j)
      array[i][j] = nanf("NULL");
  }
  
  for(unsigned int l = 0; l < nb_locus; ++l)
    denominator[l] = 0;
  
  double p, pq, var;
  
  for (unsigned int i = 0; i < patchNbr; ++i) {
    
    if( !pop_sizes[i] ) continue;
    
    for (unsigned int l = 0; l < nb_locus; ++l) {
      
      array[i][l] = 0;
      
      for (unsigned int u = 0; u < nb_allele; ++u) {
        
        p = _alleleFreqTable.get(i, l, u); //p_liu
        
        pq = p * (1 - p);
        
        var = p - _globalAlleleFreq.get(l, u); //(p_liu - pbar_u)
        
        var *= var;
        
        array[i][l] += pq * pop_sizes[i] / (pop_sizes[i] -1);
        
        denominator[l] += pop_sizes[i] * var + pop_weights[i] * pq;
        
      } // end for allele
    }// end for locus
  }//end for pop
  
  
  for (unsigned int i = 0; i < patchNbr; ++i) {
    if( !pop_sizes[i] ) continue;
    for (unsigned int l = 0; l < nb_locus; ++l) 
      array[i][l] = 1 - ( array[i][l] * sum_weights / denominator[l]);
  }
 
  delete [] pop_weights;
  delete [] pop_sizes;
  delete [] denominator;
}

// ----------------------------------------------------------------------------------------
// setNeiGeneticDistance
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setNeiGeneticDistance(age_t AGE)
{
  //see: Nei, M. 1978. Genetics 89:583-590
  unsigned int patchNbr = this->_pop->getPatchNbr();
  unsigned int nb_locus = this->_SHLinkedTrait->get_locus_num();
  unsigned int nb_allele = this->_SHLinkedTrait->get_allele_num();
  
  setAlleleTables(AGE);
  
  if(_D == NULL) 
    
    _D = new TMatrix(patchNbr, patchNbr);
  
  else if( _D->length() != patchNbr*patchNbr )
    
    _D->reset(patchNbr, patchNbr);
  
  _D->assign(nanf("NULL"));
  
  double num, denom, p1, p2, sum;
  double *sum_square = new double [patchNbr];
  double *sample_size = new double [patchNbr];
  
  for(unsigned int i = 0; i < patchNbr; ++i) {
    
    sample_size[i] = 2 * _pop->size(AGE, i);
    
    sum_square[i] = 0;
    
    if( !sample_size[i] ) continue;
    
    for(unsigned int l = 0; l < nb_locus; ++l) {
      
      sum = 0;
      
      for(unsigned int u = 0; u < nb_allele ; ++u) {
        p1 = _alleleFreqTable.get(i, l, u);
        sum += p1*p1;
      }
      
      sum_square[i] += (sample_size[i] * sum -1) / (sample_size[i] -1); //unbiased estimate, equ. 6 in Nei 1978 (Genetics)
    }
  }
  
  unsigned int pairs = 0;
  double Dpair;
  _meanD = 0;
  
  for(unsigned int i = 0; i < patchNbr-1; ++i) {
    
    if( !sample_size[i] ) continue;
    
    for(unsigned int j = i + 1; j < patchNbr; ++j) {

      if( !sample_size[j] ) continue;
      
      num = 0;
      pairs++;
      
      for(unsigned int l = 0; l < nb_locus; ++l) {
        for(unsigned int u = 0; u < nb_allele ; ++u) {
          p1 = _alleleFreqTable.get(i, l, u);
          p2 = _alleleFreqTable.get(j, l, u);
          num += p1*p2;
        }
      }      
      
      denom = sqrt(sum_square[i] * sum_square[j]);
      Dpair = -log(num/denom);
      _meanD += Dpair;
      _D->set(i, j, Dpair);
    }
  }
  _meanD /= pairs;
  
  delete [] sum_square;
  delete [] sample_size;
}
// ----------------------------------------------------------------------------------------
// getDxy
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::getDxy(unsigned int age_class)
{
  double D = 0;
  unsigned int p = 0;
  
  age_t AGE = (static_cast<age_idx>(age_class) == OFFSx ? OFFSPRG : ADULTS);
  
  for (unsigned int p1 = 0; p1 < _pop->getPatchNbr(); ++p1) {
    
    if( !_pop->size( AGE , p1) ) continue;
      
    for (unsigned int p2 = p1 + 1; p2 < _pop->getPatchNbr(); ++p2) {
      
      if( !_pop->size( AGE , p2) ) continue;

      D += getDxyPerPatch(static_cast<age_idx>(age_class), p1, p2);
    
      p++;
    }
  }
  
  return (p != 0 ? D/p : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getDxy
// ----------------------------------------------------------------------------------------
double TTNeutralGenesSH::getDxyPerPatch (age_idx age, unsigned int patch1, unsigned patch2)
{
  double D = 0;
  Patch* patch_1 = _pop->getPatchPtr(patch1);
  Patch* patch_2 = _pop->getPatchPtr(patch2);
  unsigned int size_1 = patch_1->size(FEM, age);
  unsigned int size_2 = patch_2->size(FEM, age);
  unsigned int N = 0;//4*size_1 * size_2;// (2N)^2, tot num of haplotype comparisons
  
  
  unsigned int num_loc = _SHLinkedTrait->get_locus_num();
  
  unsigned char **seq_1, **seq_2;
  
  //females:
  for (unsigned int i = 0; i < size_1; ++i) {
    
    seq_1 = (unsigned char**)patch_1->get(FEM, age, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
    
    for (unsigned int j = 0; j < size_2; ++j) {
      
      seq_2 = (unsigned char**)patch_2->get(FEM, age, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < num_loc; ++l) {
        D += (seq_1[0][l] != seq_2[0][l]);
        D += (seq_1[1][l] != seq_2[1][l]);
        D += (seq_1[0][l]!=seq_2[1][l]);
        D += (seq_1[1][l]!=seq_2[0][l]);
      } 
      N += 4;
    }
    
    for (unsigned int j = 0; j < patch_2->size(MAL, age); ++j) {
      
      seq_2 = (unsigned char**)patch_2->get(MAL, age, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < num_loc; ++l) {
        D += (seq_1[0][l]!=seq_2[0][l]);
        D += (seq_1[1][l]!=seq_2[1][l]);
        D += (seq_1[0][l]!=seq_2[1][l]);
        D += (seq_1[1][l]!=seq_2[0][l]);
      } 
      N += 4;
      
    }
  }
  
  //same for the males:
  size_1 = patch_1->size(MAL, age);
  size_2 = patch_2->size(MAL, age);
  
//  N += 4 * size_1 * size_2;
  
//  if( N == 0) return 0; //needless to go further
  
  for (unsigned int i = 0; i < size_1; ++i) {
    
    seq_1 = (unsigned char**)patch_1->get(MAL, age, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
    
    for (unsigned int j = 0; j < size_2; ++j) {
      
      seq_2 = (unsigned char**)patch_2->get(MAL, age, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < num_loc; ++l) {
        D += seq_1[0][l]!=seq_2[0][l];
        D += seq_1[1][l]!=seq_2[1][l];
        D += seq_1[0][l]!=seq_2[1][l];
        D += seq_1[1][l]!=seq_2[0][l];
      } 
      N += 4;
    }
    
    for (unsigned int j = 0; j < patch_2->size(FEM, age); ++j) {
      
      seq_2 = (unsigned char**)patch_2->get(FEM, age, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
      
      for (unsigned int l = 0; l < num_loc; ++l) {
        D += seq_1[0][l]!=seq_2[0][l];
        D += seq_1[1][l]!=seq_2[1][l];
        D += seq_1[0][l]!=seq_2[1][l];
        D += seq_1[1][l]!=seq_2[0][l];
      } 
      N += 4;
      
    }
  }
  
  if( N == 0) return 0;

  return D/N;
}
