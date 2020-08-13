/**  $Id: patch.cc,v 1.6 2015-07-13 08:52:59 fred Exp $
*
*  @file patch.cc
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
*  @author fred
*/

#include <iostream>
#include "metapop.h"
#include "Uniform.h"
#include "output.h"

using namespace std;

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Patch* Patch::init(unsigned int nbfem, unsigned int nbmal, unsigned int id)
{
  _ID = id;
  _KFem = nbfem;
  _KMal = nbmal;
  _K = _KFem + _KMal;
  _isExtinct = false;
  _age = 0;
  
  reset_counters();
 
  reset_containers();
  
  return this;
}
// ----------------------------------------------------------------------------------------
// reset_counters
// ----------------------------------------------------------------------------------------
void Patch::reset_counters()
{
  nbEmigrant = 0;
  nbImigrant = 0;
  nbPhilopat = 0;
  nbKolonisers = 0;
}
// ----------------------------------------------------------------------------------------
// reset_containers
// ----------------------------------------------------------------------------------------
void Patch::reset_containers()
{  
  for(unsigned int i=0; i < _nb_age_class; i++) {
    _containers[MAL][i].assign( _KMal, 0 );
    _containers[FEM][i].assign( _KFem, 0 );
    _sizes[MAL][i] = 0;
    _sizes[FEM][i] = 0;
    _capacities[MAL][i] = _KMal;
    _capacities[FEM][i] = _KFem;
  }
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_t AGE, Metapop* pop)
{
  unsigned int mask = 1;

  for(unsigned int i = 0; i < _nb_age_class; i++) {
    if( (mask & AGE) != 0) setNewGeneration(static_cast<age_idx>(i), pop);
    mask<<=1;
  }
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_idx AGE, Metapop* pop)
{  
  Individual *new_ind;

  //--------------------------------------------------------------------
  //if too much females in the Patch, flush them into the RecyclingPOOL
  if(size(FEM, AGE) > 0) flush(FEM, AGE, pop);
  
  for(unsigned int i = 0; i < _KFem; i++) {
    new_ind = pop->makeNewIndividual(0,0,FEM,_ID);
    new_ind->create_first_gen();
    add(FEM, AGE, new_ind);
  }
  
  //--------------------------------------------------------------------
  //males: same as for the females....
  if(size(MAL, AGE) > 0) flush(MAL, AGE, pop);
  
  for(unsigned int i = 0; i < _KMal; i++) {
    new_ind = pop->makeNewIndividual(0,0,MAL,_ID);
    new_ind->create_first_gen();
    add(MAL, AGE, new_ind);
  }
  
}
// ----------------------------------------------------------------------------------------
// ~Patch
// ----------------------------------------------------------------------------------------
Patch::~Patch()
{
//#ifdef _DEBUG_
//  message("Patch::~Patch\n");
//#endif
  
  for (unsigned int i = 0; i < 2; ++i)
    for(unsigned int j = 0; j < _nb_age_class; ++j)
      for(unsigned int k = 0; k < _sizes[i][j] ; ++k)
        delete _containers[i][j][k];
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Patch::show_up()
{
  message("Patch %i:\n  age: %i; K: %i, K_fem: %i; K_mal: %i\n",_ID, _age, _K, _KFem, _KMal);
  for(unsigned int j = 0; j < _nb_age_class; ++j)
    message("   age class %i: females: %i; males: %i\n", j, _sizes[FEM][j], _sizes[MAL][j]);
}

