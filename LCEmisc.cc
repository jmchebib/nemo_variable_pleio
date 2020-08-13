/** $Id: LCEmisc.cc,v 1.22 2016-09-28 13:59:08 fred Exp $
 *
 *  @file LCEmisc.cc
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
 *  created on @date 16.06.2005
 * 
 *  @author fred
 */

#include <cmath>
#include <deque>
#include "LCEmisc.h"
#include "metapop.h"
#include "Uniform.h"
#include "output.h"
#include "tstring.h"

// ------------------------------------------------------------------------------

//                           LCE_Aging

// ----------------------------------------------------------------------------------------
// LCE_Aging::execute
// ----------------------------------------------------------------------------------------
void LCE_Aging::execute ()
{  
#ifdef _DEBUG_
  message("LCE_Aging::execute (Patch nb: %i offsprg nb: %d adlt nb: %d; ",_popPtr->getPatchNbr()
          ,_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  unsigned int nbInd = 0;
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    patch->flush(ADLTx, _popPtr);
    
    nbInd = 0;
    
    while(nbInd++ < patch->get_KFem() && patch->size(FEM, OFFSx) != 0) 
      patch->move(FEM, OFFSx, ADLTx, RAND::Uniform( patch->size(FEM, OFFSx) ) );
    
    nbInd = 0;
    
    while(nbInd++ < patch->get_KMal() && patch->size(MAL, OFFSx) != 0)
      patch->move(MAL, OFFSx, ADLTx, RAND::Uniform( patch->size(MAL, OFFSx) ) );
    
    patch->flush(OFFSx, _popPtr);
    
    //set the Patch extinction and age tags:
    if(patch->size(ADULTS) == 0) {
      patch->set_isExtinct(true);
      patch->set_age(0);
    } else {
      patch->set_isExtinct(false);
      patch->set_age( patch->get_age() + 1 );
    }
  }
#ifdef _DEBUG_
  message("after: %i)\n",_popPtr->size( ));
#endif
}
// ------------------------------------------------------------------------------

//                             LCE_Regulation

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::execute ()
{
#ifdef _DEBUG_
  message("LCE_Regulation::execute (Patch nb: %i, offsprg: %i, adults: %i ",_popPtr->getPatchNbr()
          ,_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  Patch* patch;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i ++) {
    
    patch = _popPtr->getPatch(i);
    
    regulatePatch(patch, OFFSx, FEM);
    regulatePatch(patch, OFFSx, MAL);
    regulatePatch(patch, ADLTx, FEM);
    regulatePatch(patch, ADLTx, MAL);
    
  }
  
#ifdef _DEBUG_
  message("after: offsprg: %i, adults: %i)\n",_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  
}
// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulatePatch (Patch* patch, age_idx age, sex_t sex)
{
  unsigned int K = patch->get_K(sex);
  
  while( patch->size(sex, age) > K ) 
    _popPtr->recycle( patch->remove(sex, age, RAND::Uniform( patch->size(sex, age) ) ) );
}
// ------------------------------------------------------------------------------

//                         LCE_Patch_Extinction

//-----------------------------------------------------------------------------
// LCE_Patch_Extinction
//-----------------------------------------------------------------------------
LCE_Patch_Extinction::LCE_Patch_Extinction( ) : LifeCycleEvent("extinction",""), _Xtion_rate(0),
_harvest_size(0), _harvest_proportion(0), _harvest_size_varies(0), _by_size(0), 
_by_proportion(0), _harvest_dist_stdev(0), _extinction_threshold(0), _rand_size_fct(0)
{
  ParamUpdater< LCE_Patch_Extinction > * updater =
  new ParamUpdater< LCE_Patch_Extinction > (&LCE_Patch_Extinction::setParameters);
  add_parameter("extinction_rate", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_threshold", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_size", INT, 0, 0, 0, 0, updater);
  add_parameter("extinction_proportion", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_size_distribution", STR, 0, 0, 0, 0, updater);
  add_parameter("extinction_size_dist_stdev", DBL, 0, 0, 0, 0, updater);
  add_parameter("extinction_size_dist_shape", DBL, 0, 0, 0, 0, updater);
}

//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::init
//-----------------------------------------------------------------------------
bool LCE_Patch_Extinction::setParameters ()
{
  
  if(get_parameter("extinction_rate")->isSet()) {
    
    if(!_Xtion_rate) _Xtion_rate = new TMatrix();
    
    if(!set_matrix_param(_Xtion_rate, "extinction_rate")) return false;
    
  } else {
    if(_Xtion_rate) delete _Xtion_rate; _Xtion_rate = 0;
  }
  
  if(get_parameter("extinction_size")->isSet()) {
    
    if(!_harvest_size) _harvest_size = new TMatrix();
    
    if(!set_matrix_param(_harvest_size, "extinction_size")) return false;
    
    _by_size = true;
    
  } else {
    if(_harvest_size) delete _harvest_size; _harvest_size = 0;
    _by_size = false;
  }
  
  if(get_parameter("extinction_proportion")->isSet()) {
    if(!_harvest_proportion) _harvest_proportion = new TMatrix();
    
    if(!set_matrix_param(_harvest_proportion, "extinction_proportion")) return false;
    
    _by_proportion = true;
    
  } else {
    if(_harvest_proportion) delete _harvest_proportion; _harvest_proportion = 0;
    _by_proportion = false;
  }
  
  _extinction_threshold = get_parameter_value("extinction_threshold");
  
  if( !_Xtion_rate && !_by_size && !_by_proportion) {
    error("Please give one of the following parameter: \"extinction_rate\", \"extinction_size\", or \"extinction_proportion\".\n");
    return false;
  }
  else if(_by_size && _by_proportion) {
    warning("Both \"extinction_size\" and \"extinction_proportion\" are set, using sizes only.\n");
    _by_proportion = false;
  }
  
  if(get_parameter("extinction_size_distribution")->isSet()) {
    
    if(!_by_size) {
      error("\"extinction_size_distribution\" is set but the \"extinction_size\" parameter is not!\n");
      return false;
    }
    
    _harvest_distribution = _paramSet->getArg("extinction_size_distribution");
    _harvest_size_varies = true;
    _harvest_dist_stdev = get_parameter_value("extinction_size_dist_stdev");
    
    if(_harvest_distribution.compare("poisson") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_poisson;
    
    else if(_harvest_distribution.compare("uniform") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_uniform;
    
    else if(_harvest_distribution.compare("normal") == 0) {
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_gaussian;
      
      if(_harvest_dist_stdev == -1) {
        error("Standard deviation of the normal distribution for the harvesting size distribution is missing!\n");
        return false;
      }
      
    } else if(_harvest_distribution.compare("lognormal") == 0) {
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_lognormal;
      
      if(_harvest_dist_stdev == -1) {
        error("Standard deviation of the lognormal distribution for the harvesting size distribution is missing!\n");
        return false;
      }
      
    } else if(_harvest_distribution.compare("exponential") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_exp;
    
    else {
      error("Distribution \"%s\" is not a valid option for \"harvest_size_distribution\"\n",
            _harvest_distribution.c_str());
      return false;
    }
    
  }
  return true;
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::set_matrix_param
//-----------------------------------------------------------------------------
bool LCE_Patch_Extinction::set_matrix_param (TMatrix *mat, string name)
{
  double value;
  Param* param = get_parameter(name);
  
  if(param->isMatrix()) {
    
    param->getMatrix(mat);
    
    if(mat->getNbRows() > 1) {
      error("The \"%s\" matrix must be a one-dimensional array.\n", name.c_str());
      return false;
    }
    
    if(mat->getNbCols() != _popPtr->getPatchNbr()) {
      error("The length of the \"%s\" array must be equal to the number of patches.\n", name.c_str());
      return false;
    }
    
  } else {
    value = param->getValue();
    mat->reset(1, _popPtr->getPatchNbr());
    mat->assign(value);
  }
  return true; 
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::execute
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::execute ()
{
#ifdef _DEBUG_
  message("LCE_Patch_Extinction::execute ");
  unsigned int cnt = 0;
#endif
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    if(_by_size || _by_proportion) {
      do_remove(OFFSx, patch);
      do_remove(ADLTx, patch);
    } else if(_Xtion_rate)
      if( RAND::Uniform() < _Xtion_rate->get(0, i) )
        do_flush(patch);
    
    if(_extinction_threshold != -1) {
      if( _extinction_threshold < 1 && (double)patch->size(ALL)/patch->get_K() < _extinction_threshold ) 
        do_flush(patch);
      else if( patch->size(ALL) < _extinction_threshold )
        do_flush(patch);
    }
#ifdef _DEBUG_
    cnt += (patch->get_isExtinct());
#endif
  }
#ifdef _DEBUG_
  message("(%i extinct patches)\n",cnt);
#endif  
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::do_flush
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::do_flush (Patch *patch)
{
  patch->flush(_popPtr);
  patch->set_isExtinct(true);
  patch->set_age(0);
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::do_remove
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::do_remove(age_idx AGE, Patch* patch)
{
  unsigned int remove_size;
  sex_t sex;
  
  //check if probability of event is set, and if removal will happen
  if(_Xtion_rate) if( RAND::Uniform() > _Xtion_rate->get(0, patch->getID()) ) return;
  
  if(patch->size(AGE) != 0) {
    
    remove_size = get_harvest_size(AGE, patch);
    
    for(unsigned int i = 0; i < remove_size; ++i) {
      
      sex = (sex_t)RAND::RandBool(); 
      
      if(patch->size(sex, AGE) != 0)
        _popPtr->recycle( patch->remove( sex, AGE, (unsigned int)RAND::Uniform(patch->size(sex, AGE)) ) );
      else //we already know here that the patch is not empty
        _popPtr->recycle( patch->remove( (sex_t)!sex, AGE, (unsigned int)RAND::Uniform(patch->size( (sex_t)!sex, AGE) ) ) );
      
      if(patch->size(AGE) == 0) break;
    }
    //    cout<<"--removed "<<remove_size<<" individuals in age class "<<AGE<<", patch "<< patch->getID()<<" size = "<<patch->size(AGE)<<endl;
  }
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::get_harvest_size
//-----------------------------------------------------------------------------
unsigned int LCE_Patch_Extinction::get_harvest_size (age_idx AGE, Patch *patch)
{
  
  if( _by_size ) {
    
    if(_harvest_size_varies)  
      return (this->*_rand_size_fct) (_harvest_size->get(0, patch->getID()));
    
    else return (unsigned int)_harvest_size->get(0, patch->getID());
    
  } else if( _by_proportion ) {
    
    return (unsigned int)(_harvest_proportion->get(0, patch->getID()) * patch->size(AGE));
    
  }
  
  return 0;
}

// ------------------------------------------------------------------------------

//                             LCE_Cross
// ----------------------------------------------------------------------------------------
// LCE_Cross::execute
// ----------------------------------------------------------------------------------------
LCE_Cross::LCE_Cross() : LifeCycleEvent("cross",""), _nSire(0), _nDam(0), _nOffspring(0), 
_atGeneration(0), _reader(0)
{
  ParamUpdater< LCE_Cross > * updater = new ParamUpdater< LCE_Cross > (&LCE_Cross::setParameters);
  add_parameter("cross_num_sire", INT, 0, 0, 0, 0, updater);
  add_parameter("cross_num_dam", INT, 0, 0, 0, 0, updater);
  add_parameter("cross_num_offspring", INT, 0, 0, 0, 0, updater);
  add_parameter("cross_at_generation", INT, 1, 0, 0, 0, updater);
  add_parameter("cross_do_among_pop", BOOL, 0, 0, 0, 0, updater);
  add_parameter("cross_do_within_pop", BOOL, 0, 0, 0, 0, updater);
  add_parameter("cross_with_replacement", BOOL, 0, 0, 0, 0, updater);

  add_parameter("cross_with_pedigree_file", STR, 0, 0, 0, 0, updater);

}
// ----------------------------------------------------------------------------------------
// LCE_Cross::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Cross::setParameters ( )
{
  _nSire = (unsigned int)get_parameter_value("cross_num_sire");
  _nDam = (unsigned int)get_parameter_value("cross_num_dam");
  _nOffspring = (unsigned int)get_parameter_value("cross_num_offspring");
  _atGeneration = (unsigned int)get_parameter_value("cross_at_generation");
  
  if(get_parameter("cross_do_among_pop")->isSet())
    _doAmongPop = (unsigned int)get_parameter_value("cross_do_among_pop");
  else
    _doAmongPop = 0;
  
  string arg = get_parameter("cross_do_within_pop")->getArg();
  if (tstring::str2int(arg) == 1) {
    _doWithinPop = 1;
  } else if (tstring::str2int(arg) == 0) {
    _doWithinPop = 0;
  } else {
    _doWithinPop = 1;
  }
  
  //  if(get_parameter("cross_do_within_pop")->isSet())
  //    _doWithinPop = (unsigned int)get_parameter_value("cross_do_within_pop");
  //  else
  //    _doWithinPop = 1;
  
  if(get_parameter("cross_with_replacement")->isSet())
    _doReplace = (unsigned int)get_parameter_value("cross_with_replacement");
  else 
    _doReplace = 0;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::loadFileServices
// ----------------------------------------------------------------------------------------
void LCE_Cross::loadFileServices ( FileServices* loader )
{
	if( ! get_parameter("cross_with_pedigree_file")->isSet() ) return;

	if(_reader) delete _reader;
	_reader = new FHPedigreeReader(this);
	//set to read mode:
	_reader->set_isInputHandler(true);
	//attach to file manager:
	loader->attach_reader(_reader);

	_pedigree_file = get_parameter("cross_with_pedigree_file")->getArg();

	_reader->set(true, false, 1, _atGeneration, 0, "", this);

}
// ----------------------------------------------------------------------------------------
// LCE_Cross::execute
// ----------------------------------------------------------------------------------------
void LCE_Cross::execute ()
{
  Patch* patch;
  unsigned int nsire;
  deque<Individual*> males, females;
  
#ifdef _DEBUG_
  message("LCE_Cross::execute\n");
#endif
  
  if(_popPtr->getCurrentGeneration() != _atGeneration) return;
  
  if(get_parameter("cross_with_pedigree_file")->isSet())

	  generatePedigree();

  else {
  
	  if(_popPtr->size(OFFSPRG) != 0) {
		warning("Offspring are still present at time of crossing, flushing\n");
		_popPtr->flush(OFFSPRG);
	  }

	  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i ++) {

		patch = _popPtr->getPatch(i);

		if(patch->size(MAL, ADLTx) == 0) continue;

		if(_nSire > patch->size(MAL, ADLTx)) {
		  warning("LCE_Cross:: num_sire greater than actual number of males in patch %i, reseting.\n", patch->size(MAL, ADLTx));
		  nsire = patch->size(MAL, ADLTx);
		} else
		  nsire = _nSire;

		males.clear();
		patch->getCopy(MAL, ADLTx, males);

		if(_doAmongPop) sampleAmongPop(patch, males, nsire);

		males.clear();
		females.clear();

		if(patch->size(FEM, ADLTx) == 0) continue;


		if(_doWithinPop) {
		  patch->getCopy(MAL, ADLTx, males);
		  patch->getCopy(FEM, ADLTx, females);
		  sampleWithinPop(patch, males, females, nsire);
		}
	  }
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::generatePedigree
// ----------------------------------------------------------------------------------------
void LCE_Cross::generatePedigree ()
{
	Individual* ind;;

	//read the pedigree file only once per simulation:
	if(_popPtr->getCurrentReplicate() == 1) {

		_reader->FHread(_pedigree_file);

	}

	vector< unsigned long * > pedigree = _reader->getPedigree();

	_pedigree_pop.clear();

	//read in the individual ID of each record and create individuals
	for(unsigned int i = 0; i < pedigree.size(); ++i) {

//		cout<<"\nrecord #"<<i+1<<endl;

		if(_pedigree_pop.find(pedigree[i][0]) != _pedigree_pop.end())
			warning("found duplicated ID %i in pedigree, replacing it with latest record"
					, pedigree[i][0]);
		else
			_pedigree_pop[ pedigree[i][0] ] = _popPtr->getPrototypeClone();

//		cout<<"\nsetting new ind in pedigree: "<<pedigree[i][0]<<" "<<pedigree[i][1]
//											<<" "<<pedigree[i][2]<<"\n";
		ind = _pedigree_pop[ pedigree[i][0] ];

		ind->init(); //need to allocate traits' memory; warning: sets parents' ID to 0!!!
		ind->setID( pedigree[i][0] );
		ind->setFatherID( pedigree[i][1] );
		ind->setMotherID( pedigree[i][2] );

//		ind->show_up();
	}

	//create the genetics
	map< unsigned long, Individual*>::const_iterator REC = _pedigree_pop.begin();

	while(REC != _pedigree_pop.end()) {

		//recursively create ancestors of an individual, and the individual by inheritance
		// adding mutations as well
		if(!create_individual_ancestors( REC->second ))
			error("while creating individuals on a pedigree (ID=%i)\n",REC->first);

		REC++;
	}

	//deal with pedigreed population
	//for now, empty all patches and dump the ped pop in the first patch (offspring fems)
	_popPtr->flush();
	_popPtr->setCurrentAge(OFFSPRG);

	Patch* patch = _popPtr->getPatch(0);

	REC = _pedigree_pop.begin();
	while(REC != _pedigree_pop.end()) {

		patch->add(FEM, OFFSx, REC->second);

		REC++;
	}

}
// ----------------------------------------------------------------------------------------
// LCE_Cross::create_individual_ancestors
// ----------------------------------------------------------------------------------------
bool LCE_Cross::create_individual_ancestors(Individual* ind)
{
	Individual *mother = ind->getMother();
	Individual *father = ind->getFather();
	unsigned long sire = ind->getFatherID();
	unsigned long dam = ind->getMotherID();

	//the individual is set when the mum and dad ptrs are set, exit
	if(mother && father) return true;

	//check parents' ID for sampling in the whole pop; choose individuals in first patch
	if(dam == 0) {

		if(_popPtr->size(FEM, ADLTx, 0) != 0) {

			mother = _popPtr->get(FEM, ADLTx, RAND::Uniform( _popPtr->size(FEM, ADLTx, 0) ), 0);

		} else {
			mother = _popPtr->makeNewIndividual(0, 0, FEM, 0);
			mother->create_first_gen();
		}

		ind->setMother(mother);
		ind->setMotherID(mother->getID());

	} else if( !mother ) {

		if( create_individual_ancestors( _pedigree_pop[ dam ] ) ) {
			ind->setMother(_pedigree_pop[ dam ]);
		}
	}

	//father
	if(sire == 0) {

		if(_popPtr->size(MAL, ADLTx, 0) != 0 ) {

			father = _popPtr->get(MAL, ADLTx, RAND::Uniform( _popPtr->size(MAL, ADLTx, 0) ), 0);

		} else {
			father = _popPtr->makeNewIndividual(0, 0, MAL, 0);
			father->create_first_gen();
		}

		ind->setFather(father);
		ind->setFatherID(father->getID());

	} else if(!father) {

		if( create_individual_ancestors( _pedigree_pop[ sire ] ) )
			ind->setFather( _pedigree_pop[ sire ] );

	}

	mother = ind->getMother();
	father = ind->getFather();

	if(mother && father) {

//		cout << "create "<<ind->getID()<<" ("<<mother->getID()<<" + "<<father->getID()<<") "
//				<<"("<<dam<<" + "<<sire<<")\n";

		ind->create(true, true); //create with inheritance and mutation
	}

	return( mother && father ); //return false when parents are not created
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::sampleAmongPop
// ----------------------------------------------------------------------------------------
void LCE_Cross::sampleAmongPop(Patch* patch, deque<Individual*>& males, unsigned int nsire)
{
  unsigned int at, deme, npatch = _popPtr->getPatchNbr();
  Individual *sire, *dam;
  Patch* aimedPatch;
  
  for(unsigned int s = 0; s < nsire; ++s) {
    
    //sampling parents
    at = (unsigned int)RAND::Uniform( males.size() );
    
    sire = males[at];
    
    if( !_doReplace ) males.erase(males.begin() + at);
    
    for(unsigned int d = 0; d < _nDam; ++d) {
      //sample a deme
      do {
        deme = (unsigned int)RAND::Uniform(npatch);
      } while(deme == patch->getID() || patch->size(FEM, ADLTx) == 0);
      
      aimedPatch = _popPtr->getPatchPtr(deme);
      //sample a female
      dam = aimedPatch->get( FEM, ADLTx, (unsigned int)RAND::Uniform( aimedPatch->size(FEM, ADULTS) ) );
      
      if(dam == NULL) {d--; continue;}
      
      //mate and breed; ignore sex, put all offs in females' container
      for(unsigned int k = 0; k < _nOffspring; ++k) {
        
        patch->add(FEM, OFFSx, _popPtr->makeOffsprg(dam, sire, FEM, patch->getID()) );
        
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::sampleWithinPop
// ----------------------------------------------------------------------------------------
void LCE_Cross::sampleWithinPop(Patch* patch, deque<Individual*>& males, deque<Individual*>& females, unsigned int nsire)
{  
  unsigned int ndam, at;
  Individual *sire, *dam;
  
  if(!_doReplace && nsire * _nDam > patch->size(FEM,ADLTx)) {
    ndam = patch->size(FEM,ADLTx)/nsire;
    warning("LCE_Cross:: total num of dam greater than available females in patch %i, reseting num_dam to %i.\n",patch->getID(), ndam);
  } else
    ndam = _nDam;
  
  if(nsire < 2 || ndam == 0) return;
  
  for(unsigned int s = 0; s < nsire; ++s) {
    
    //sampling parents
    at = (unsigned int)RAND::Uniform( males.size() );
    
    sire = males[at];
    
    if( !_doReplace ) males.erase(males.begin() + at);
    
    for(unsigned int d = 0; d < ndam; ++d) {
      
      at = (unsigned int)RAND::Uniform( females.size() );
      
      dam = females[at];
      
      if( !_doReplace ) females.erase(females.begin() + at);
      
      //now mate and breed; ignore sex, put all offs in females' container
      for(unsigned int k = 0; k < _nOffspring; ++k) {
        
        patch->add(FEM, OFFSx, _popPtr->makeOffsprg(dam, sire, FEM, patch->getID()) );
        
      }
    }
  }
}
// ----------------------------------------------------------------------------------------

//                             FHPedigreeReader
// ----------------------------------------------------------------------------------------
FHPedigreeReader::FHPedigreeReader(LCE_Cross* event):
		EventFileHandler< LCE_Cross > (event, ".ped")
{

}
void FHPedigreeReader::FHread (string& filename)
{
	ifstream FILE(filename.c_str(),ios::in);

	if(!FILE) fatal("could not open pedigree input file \"%s\"\n",filename.c_str());

	string ID, sire, dam;

	unsigned int line = 0;

	unsigned long * record = 0;

	_pedigree.clear();

	message("\n++++ reading pedigree file: %s\n",filename.c_str());

	while(FILE){

		line++;

		FILE >> ID >> sire >> dam;

		if(!FILE || FILE.eof()) break;

		record = new unsigned long [3];

		if(tstring::isanumber(ID)) record[0] = tstring::str2ulong(ID);
		else if(tstring::isNA(ID)) record[0] = 0;
		else {
			error("%s: ID on line %i is not a number\n",filename.c_str(),line);
			delete [] record; record = 0;
			continue;
		}

		if(tstring::isanumber(sire)) record[1] = tstring::str2ulong(sire);
		else if(tstring::isNA(sire)) record[1] = 0;
		else {
			error("%s: 'sire' on line %i is not a number\n",filename.c_str(),line);
			delete [] record; record = 0;
			continue;
		}

		if(tstring::isanumber(dam)) record[2] = tstring::str2ulong(dam);
		else if(tstring::isNA(dam)) record[2] = 0;
		else {
			error("%s: 'dam' on line %i is not a number\n",filename.c_str(),line);
			delete [] record; record = 0;
			continue;
		}
//		cout<<line<<": "<<ID<<" "<<sire<<" "<<dam<<endl;

		_pedigree.push_back(record);
	}

	if(FILE) FILE.close();

	message("++++ pedigree has %i valid records\n", _pedigree.size());
}

// ----------------------------------------------------------------------------------------

//                             LCE_Resize
// ----------------------------------------------------------------------------------------
// LCE_Resize::execute
// ----------------------------------------------------------------------------------------
LCE_Resize::LCE_Resize() : LifeCycleEvent("resize",""), _genITER(0), _atGeneration(0),
_do_flush(0), _do_fill(0), _do_regulate(0), _setAge(0), _patchBackup(0)
{
  ParamUpdater< LCE_Resize > * updater = new ParamUpdater< LCE_Resize > (&LCE_Resize::setParameters);
  add_parameter("resize_at_generation", INT, 1, 0, 0, 0, updater);
  updater = new ParamUpdater< LCE_Resize > (&LCE_Resize::updateParameters);
  add_parameter("resize_patch_number", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_patch_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_female_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_male_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_age_class", STR, 0, 0, 0, 0, updater);
  add_parameter("resize_do_flush", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_do_regulate", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_do_fill", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_keep_patch", MAT, 0, 0, 0, 0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Resize::setParameters ( )
{
  
  _generations.clear();
  
  if(_paramSet->isMatrix("resize_at_generation") ) {
    
    TMatrix tmp;
    
    _paramSet->getMatrix("resize_at_generation", &tmp);
    
    for(unsigned int i = 0; i < tmp.getNbCols(); i++)
      _generations.push_back((int)tmp.get(0, i));
    
  } else
    _generations.push_back((int)get_parameter_value("resize_at_generation"));
  
  _genITER = _generations.begin();
  _atGeneration =  (*_genITER);
  
  return updateParameters();
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::updateParameters
// ----------------------------------------------------------------------------------------
bool LCE_Resize::updateParameters ( )
{
  string option = _paramSet->getArg("resize_age_class");
  
  if(option.length() == 0) _setAge = ALL;
  else if(option.compare("OFFSPRG") == 0 || option.compare("offspring") == 0 ||
          option.compare("0") == 0){
    _setAge = OFFSPRG;
  } else if(option.compare("ADULTS") == 0 || option.compare("adults") == 0 ||
            option.compare("1") == 0) {
    _setAge = ADULTS;
  } else if(option.compare("ALL") == 0 || option.compare("all") == 0) {
    _setAge = ALL;
  } else {
    warning("\"%s\" is not a valid option for parameter \"resize_age_class\".\n",option.c_str());
    _setAge = ALL;
  }
  
  if( !get_parameter("resize_patch_number")->isSet() &&
     !get_parameter("resize_patch_capacity")->isSet() &&
     !get_parameter("resize_female_capacity")->isSet() &&
     !get_parameter("resize_male_capacity")->isSet() &&
     !get_parameter("resize_keep_patch")->isSet()) {
    error("LCE_Resize:: at least one of the population size parameters must be specified.\n");
    return false;
  }
  
  _do_regulate = get_parameter("resize_do_regulate")->isSet();
  _do_flush = get_parameter("resize_do_flush")->isSet();
  _do_fill = get_parameter("resize_do_fill")->isSet();
  
  if(get_parameter("resize_keep_patch")->isSet()) {
    _paramSet->getMatrix("resize_keep_patch", &_patch2keep);
    if(_patch2keep.getNbCols() > _popPtr->getPatchNbr() && 
       _popPtr->getCurrentGeneration() == (unsigned)_atGeneration) {
      error("LCE_Resize:: more patches to keep than existing patches in the population!\n");
      return false;
    }
  } else
    _patch2keep.reset();
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::execute
// ----------------------------------------------------------------------------------------
void LCE_Resize::execute ()
{

  if(_atGeneration > 0)
    if(_popPtr->getCurrentGeneration() != (unsigned)_atGeneration) return;
  
  if(++_genITER != _generations.end())
    _atGeneration = (*_genITER);
  else {
    _genITER = _generations.begin();
    _atGeneration = (*_genITER);
  }
  
  // if(_everyGeneration > 0)
  //    if(_atGeneration > 0) {
  //      //in this case, the _atGeneration value is taken as the starting point
  //      if( _popPtr->getCurrentGeneration() < (unsigned)_atGeneration ||
  //         (_popPtr->getCurrentGeneration() - _atGeneration) % _everyGeneration ) return;
  //    } else if(_popPtr->getCurrentGeneration() % _everyGeneration) return;
  
#ifdef _DEBUG_
  message("\nLCE_Resize::execute at %i (Patch nb: %i offsprg nb: %i adlt nb: %i)"
          ,_popPtr->getCurrentGeneration(),_popPtr->getPatchArraySize(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
#endif
  
  if( !(_popPtr->getCurrentAge() & _setAge) || _popPtr->getCurrentAge() == NONE ) {
    warning("LCE_Resize::execute:: required aged class not present in population, exiting.\n");
    return;
  }
  //first reset the metapop size parameters:
  ParamSet *originalPSet = _popPtr->get_paramset();
  //make a backup copy of the initial parameters state to restore them later:
  ParamSet *pSet = new ParamSet(*originalPSet);
  
  if( get_parameter("resize_patch_number")->isSet() )
    pSet->set_param("patch_number", _paramSet->getArg("resize_patch_number"));
  
  if( get_parameter("resize_keep_patch")->isSet() ) {
    if(_patch2keep.getNbCols() <= _popPtr->getPatchNbr()) {
      pSet->set_param("patch_number", tstring::int2str(_patch2keep.getNbCols()));
    } else if(_patch2keep.getNbCols() > _popPtr->getPatchNbr() && 
              _popPtr->getCurrentGeneration() == (unsigned)_atGeneration)
      //this is a problem if one want to resize AND reorder... unless temporal arg used in pop params
      fatal("LCE_Resize::execute:: more patches to keep than present in the population!\n");
  }
  
  if( get_parameter("resize_patch_capacity")->isSet() )
    pSet->set_param("patch_capacity", _paramSet->getArg("resize_patch_capacity"));
  
  if( get_parameter("resize_female_capacity")->isSet() )
    pSet->set_param("patch_nbfem", _paramSet->getArg("resize_female_capacity"));
  
  if( get_parameter("resize_male_capacity")->isSet() )
    pSet->set_param("patch_nbmal", _paramSet->getArg("resize_male_capacity"));
  
  _popPtr->set_paramset(pSet);
  //just reset the patch number and patch capacities, doesn't touch the patch structure yet
  _popPtr->setPopulationParameters();
  
  //now restore the original param set from the backup, will be needed to start the next replicate:
  _popPtr->set_paramset(originalPSet);
  
  delete pSet;
  
  if(_do_flush) {
    
    //new populations will be empty and supernumerary patches are flushed and deleted
    buildNewPatchArrayNoBackup();
    
  } else {
    
    if(_patchBackup) {
      _patchBackup->clear();
      delete _patchBackup;
    }
    _patchBackup = new Patch();
    
    _patchBackup->init(_popPtr->getPatchKFem(), _popPtr->getPatchKMal(), 0);
    
    //should be used when existing patches are merged together (fusion) or split (fission)
    buildNewPatchArrayWithBackup();
    
  }
  
  // updatePatchCapacities();
  
  if(_do_regulate)
    if(_do_flush)  regulate( &LCE_Resize::regulateAgeClassNoBackup );
    else           regulate( &LCE_Resize::regulateAgeClassWithBackup );
  
  if(_do_fill)
    if(_do_flush) fillPop( &LCE_Resize::fillPatchNoBackup );
    else          fillPop( &LCE_Resize::fillPatchWithBackup );
  
#ifdef _DEBUG_
  message(" ---> (Patch nb: %i offsprg nb: %i adlt nb: %i, do_fill %i, do_flush %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ), _do_fill, _do_flush);
#endif
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::updatePatchCapacities
// ----------------------------------------------------------------------------------------
void LCE_Resize::updatePatchCapacities()
{
  TMatrix *cap = _popPtr->getPatchCapacities();
  
  if(cap->getNbCols() != _popPtr->getPatchNbr())
    fatal("LCE_Resize:: more patch capacities than patches in the population!\n");
  
  for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    _popPtr->getPatchPtr(i)->set_KFem((unsigned int)cap->get(FEM, i));
    _popPtr->getPatchPtr(i)->set_KMal((unsigned int)cap->get(MAL, i));
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::buildNewPatchArrayNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::buildNewPatchArrayNoBackup()
{
  if(_patch2keep.getNbCols() != 0)  removeDesignatedPatch(false);
  
  //
  //add or remove patches until the array size matches Metapop::_patchNbr
  //up to this point, the number of patches to keep could not be != from that nbre
  _popPtr->updatePatchArray();
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::buildNewPatchArrayWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::buildNewPatchArrayWithBackup()
{
  unsigned int patchNbr =  _popPtr->getPatchNbr();
  
  if(_patch2keep.getNbCols() != 0)  removeDesignatedPatch(true);
  
  //reset the right number of patches to match new parameters
  if(_popPtr->getPatchArraySize() > patchNbr) {
    
    while(_popPtr->getPatchArraySize() > patchNbr) {
      _popPtr->getPatchPtr(0)->copy2patch(_patchBackup);
      _popPtr->getPatchPtr(0)->clear();
      _popPtr->removePatch(0);
    }
  } else while(_popPtr->getPatchArraySize() < patchNbr) _popPtr->addPatch(new Patch());
  //set the patch ID and Ks:
  _popPtr->updatePatchState();
  //regulate the patches, backup the supernumerary individuals, will be used in empty patches
  if(_do_fill) regulate( &LCE_Resize::regulateAgeClassWithBackup );
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::removeDesignatedPatch
// ----------------------------------------------------------------------------------------
void LCE_Resize::removeDesignatedPatch(bool do_backup)
{
  deque< Patch* > to_keep;
  Patch* patch;
  unsigned int patchNbr = _popPtr->getPatchArraySize() ;
  unsigned int i, nKeep = _patch2keep.getNbCols();
  unsigned int id;
  
  if(patchNbr < nKeep) fatal("LCE_Resize::more patches to keep than available in the population\n");
  
  for(i = 0; i < nKeep; i++)
    to_keep.push_back( _popPtr->getPatchPtr((unsigned int)_patch2keep.get(0, i)-1 ) );
  
  //remove the patches we want to keep from the patch array, without touching the contents
  for(i = 0; i < nKeep; i++) {
    id = to_keep[i]->getID();
    for(unsigned int j = 0; j < _popPtr->getPatchArraySize(); ++j)
      if(_popPtr->getPatchPtr(j)->getID() == id){ _popPtr->removePatch(j); break;}
  }
  //delete the patches we don't want to keep, remove the individuals and backup them if needed
  for(i = 0; i < _popPtr->getPatchArraySize(); i++) {
    
    patch = _popPtr->getPatchPtr(i);
    
    if(do_backup) { patch->copy2patch(_patchBackup); patch->clear(); }
    
    _popPtr->deletePatch(i); //this deletes the patch and the individuals it contains
    
  }
  //rebuild the patch array with the patches we want, in the right order
  for(i = 0; i < nKeep; i++)
    _popPtr->addPatch(to_keep[i]);
  
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPop
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPop ( void (LCE_Resize:: *fillFuncPtr) (unsigned int p, age_idx age))
{
  age_t active_age = (_popPtr->getCurrentAge() & _setAge);
  
  if(active_age & OFFSPRG) { 
    for(unsigned int i=0; i < _popPtr->getPatchNbr(); ++i) 
      (this->*fillFuncPtr)(i, OFFSx); 
  }
  
  if(active_age & ADULTS) {
    for(unsigned int i=0; i < _popPtr->getPatchNbr(); ++i) 
      (this->*fillFuncPtr)(i, ADLTx); 
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPatchNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPatchNoBackup(unsigned int p, age_idx age)
{
  Patch *patch = _popPtr->getPatchPtr(p);
  while (patch->size(FEM, age) != patch->get_KFem()) {
    patch->add(FEM, age, _popPtr->makeNewIndividual(NULL, NULL, FEM, patch->getID()));
  }
  while (patch->size(MAL, age) != patch->get_KMal()) {
    patch->add(MAL, age, _popPtr->makeNewIndividual(NULL, NULL, MAL, patch->getID()));
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPatchWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPatchWithBackup(unsigned int p, age_idx age)
{
  Patch *patch = _popPtr->getPatchPtr(p);
  while (patch->size(FEM, age) != patch->get_KFem() && _patchBackup->size(FEM, age) != 0) {
    patch->add(FEM, age, _patchBackup->get(FEM, age, 0));
    _patchBackup->remove(FEM, age, 0);
  }
  while (patch->size(MAL, age) != patch->get_KMal() && _patchBackup->size(MAL, age) != 0) {
    patch->add(MAL, age, _patchBackup->get(MAL, age, 0));
    _patchBackup->remove(MAL, age, 0);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulate( void (LCE_Resize::* regFuncPtr) (Patch *patch, age_idx age))
{
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    patch = _popPtr->getPatchPtr(i);
    (this->*regFuncPtr)(patch, OFFSx);
    (this->*regFuncPtr)(patch, ADLTx);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateAgeClassWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulateAgeClassWithBackup(Patch *patch, age_idx age)
{
  unsigned int ind;
  
  while(patch->size(FEM, age) > patch->get_KFem()) {
    ind = RAND::Uniform( patch->size(FEM, age) ) ;
    _patchBackup->add(FEM, age, patch->get(FEM, age, ind));
    patch->remove(FEM, age, ind);
  }
  
  while(patch->size(MAL, age) > patch->get_KMal()) {
    ind = RAND::Uniform( patch->size(MAL, age) ) ;
    _patchBackup->add(MAL, age, patch->get(MAL, age, ind));
    patch->remove(MAL, age, ind);
  }
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateAgeClassNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulateAgeClassNoBackup(Patch *patch, age_idx age)
{
  unsigned int ind;
  
  while(patch->size(FEM, age) > patch->get_KFem()) {
    ind = RAND::Uniform( patch->size(FEM, age) ) ;
    _popPtr->recycle( patch->get(FEM, age, ind) );
    patch->remove(FEM, age, ind);
  }
  
  while(patch->size(MAL, age) > patch->get_KMal()) {
    ind = RAND::Uniform( patch->size(MAL, age) ) ;
    _popPtr->recycle( patch->get(MAL, age, ind) );
    patch->remove(MAL, age, ind);
  }
  
}
