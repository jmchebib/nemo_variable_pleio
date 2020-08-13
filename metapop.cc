/** $Id: metapop.cc,v 1.19 2016-09-28 14:57:03 fred Exp $
 *
 *  @file metapop.cc
 *  Nemo2
 *
 *  Copyright (C) 2006-2015 Frederic Guillaume
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
 *  Created on @date 07.08.2004
 *
 *  @author fred
 */

#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <dirent.h>
#include <map>
#include "Uniform.h"
#include "output.h"
#include "metapop.h"
#include "individual.h"
#include "lifecycleevent.h"
#include "simulation.h"
#include "simenv.h"

using namespace std;

extern MPIenv *_myenv;

// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
Metapop::Metapop() : _mpimgr(0), _statHandler(), _writer(0),_loader(), _source(0),
_source_preserve(0), _source_load(0), _source_replicates(0), _source_replicate_digits(0),
_source_start_at_replicate(1), _source_generation(0), _source_load_periodicity(0),
_patchNbr(0), _patchK(0), _patchKfem(0), _patchKmal(0), _generations(0), _replicates(0),
_currentGeneration(0), _currentReplicate(0), _currentAge(NONE), _requiredAge(0)
{
  set_paramset("population", true, this);
  ParamUpdater<Metapop>* upd = new ParamUpdater<Metapop>( &Metapop::updatePopulationParameters );
  add_parameter("patch_capacity",INT,false,false,0,0, upd);
  add_parameter("patch_nbfem",INT,false,false,0,0, upd);
  add_parameter("patch_nbmal",INT,false,false,0,0, upd);
  add_parameter("patch_number",INT,false,false,0,0, upd);
  
  add_parameter("pop_output", BOOL, false, false, 0, 0);
  add_parameter("pop_output_dir", STR, false, false, 0, 0);
  add_parameter("pop_output_logtime", INT, false, false, 0, 0);
  add_parameter("pop_output_patch_sample_size", INT, false, false, 0, 0);

  upd = new ParamUpdater<Metapop>( &Metapop::setSourceParameters );
  add_parameter("source_pop",STR,false,false,0,0, upd);
  add_parameter("source_file_type",STR,false,false,0,0, upd);
  add_parameter("source_preserve",BOOL,false,false,0,0, upd);
  add_parameter("source_replicates",INT,false,false,0,0, upd);
  add_parameter("source_replicate_digit",INT,false,false,0,0, upd);
  add_parameter("source_start_at_replicate",INT,false,false,0,0, upd);
  add_parameter("source_generation",INT,false,false,0,0, upd);
  add_parameter("source_fill_age_class",STR,false,false,0,0, upd);
}
// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
Metapop::~Metapop()
{
#ifdef _DEBUG_
  message("Metapop::~Metapop\n");
#endif
  clear();
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool Metapop::init( )
{
  if( !(_paramSet->isSet()) ) {
    error("parameters in \"population\" are not properly set!\n");
    return false;
  }
  
  if(!setParameters()) return false;
  
  buildPatchArray();
  
  //empty and clean the RecyclingPOOL, safer...
  purgeRecyclingPOOL();
  
  return true;
}
// ----------------------------------------------------------------------------------------
// setParameters
// ----------------------------------------------------------------------------------------
bool Metapop::setParameters()
{
  if(!setPopulationParameters()) return false;
  if(!setSourceParameters()) return false;
  return true;
}
// ----------------------------------------------------------------------------------------
// setPopulationParameters
// ----------------------------------------------------------------------------------------
bool Metapop::setSourceParameters()
{
  if(_paramSet->isSet("source_pop")) {
    _source_load = true;
    _source_name = _paramSet->getArg("source_pop");
    _source_preserve = _paramSet->isSet("source_preserve");
    _source_replicates = _paramSet->isSet("source_replicates") ? 
    (unsigned int)_paramSet->getValue("source_replicates") : 0;
    
    _source_replicate_digits = _paramSet->isSet("source_replicate_digit") ? 
    (unsigned int)_paramSet->getValue("source_replicate_digit") : 1;
    
    _source_start_at_replicate = _paramSet->isSet("source_start_at_replicate") ? 
    (unsigned int)_paramSet->getValue("source_start_at_replicate") : 1;
    
    _source_generation = _paramSet->isSet("source_generation") ? 
    (unsigned int)_paramSet->getValue("source_generation") : 0;
    
    _source_filetype = _paramSet->getArg("source_file_type");
    
    if(_source_filetype.length() == 0)  _source_filetype = ".bin";
    
    _source_required_age = _paramSet->getArg("source_fill_age_class");
    
    if(_source_required_age.length() == 0) _requiredAge = NONE;
    else if(_source_required_age.compare("OFFSPRG") == 0 || _source_required_age.compare("offspring") == 0 ||
            _source_required_age.compare("0") == 0){
      _requiredAge = OFFSPRG;
    } else if(_source_required_age.compare("ADULTS") == 0 || _source_required_age.compare("adults") == 0 ||
              _source_required_age.compare("1") == 0) {
      _requiredAge = ADULTS;
    } else
      _requiredAge = NONE; //for now...
    
    
  } else {
    _source = NULL;
    _source_preserve = false;
    _source_load = false;
    _source_replicates = 0;
    _source_generation = 0;
    _source_name = "";
    _source_filetype = ".bin";
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// setPopulationParameters
// ----------------------------------------------------------------------------------------
bool Metapop::setPopulationParameters ()
{

  if(_paramSet->isSet("patch_number")) 
    _patchNbr = (unsigned int)_paramSet->getValue("patch_number");
  else _patchNbr = 0;
  
  if(_paramSet->isSet("patch_capacity")) {
    
    if( _paramSet->isMatrix("patch_capacity") ) {
      
      setPatchCapacities("patch_capacity");
      
    } else if( !(_paramSet->isSet("patch_number")) ) {
      
      error("param \"patch_number\" is missing!\n");
      return false;
      
    } else {
      _patchK = (unsigned int)_paramSet->getValue("patch_capacity");
      _patchKfem = _patchKmal = _patchK/2;
      setPatchCapacities();
    }
    
  } else if(_paramSet->isSet("patch_nbfem") && _paramSet->isSet("patch_nbmal")) {
    
    if( !(_paramSet->isMatrix("patch_nbfem")) && !(_paramSet->isMatrix("patch_nbmal")) ) {
      
      if( !(_paramSet->isSet("patch_number")) ) {
        
        error("param \"patch_number\" is missing!\n");
        return false;
        
      } else {
        _patchKfem = (unsigned int)_paramSet->getValue("patch_nbfem");
        _patchKmal = (unsigned int)_paramSet->getValue("patch_nbmal");
        _patchK = _patchKfem + _patchKmal;
        setPatchCapacities();
      }
      
    } else {
      
      if( !(_paramSet->isMatrix("patch_nbfem")) && _paramSet->isMatrix("patch_nbmal") ) {
        _patchKfem = (unsigned int)_paramSet->getValue("patch_nbfem");
        setPatchCapacities(MAL,"patch_nbmal");
      } else if( _paramSet->isMatrix("patch_nbfem") && !(_paramSet->isMatrix("patch_nbmal")) ) {
        _patchKmal = (unsigned int)_paramSet->getValue("patch_nbmal");
        setPatchCapacities(FEM,"patch_nbfem");
      } else
        setPatchCapacities("patch_nbfem","patch_nbmal");
    }
    
  } else {
    error("population parameters are not properly set!\n");
    return false;
  }

  return true;
}
// ----------------------------------------------------------------------------------------
// updatePopulationParameters
// ----------------------------------------------------------------------------------------
bool Metapop::updatePopulationParameters ()
{
  if(!setPopulationParameters()) return false;
  updatePatchArray();
  //if we added patches, they will be empty
  //fusion of existing pop is not possible here, 
  //the individuals are flushed when deleting the patches
  return true;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void Metapop::loadFileServices ( FileServices* loader )
{

  if(get_parameter("pop_output")->isSet()){

	  if(_writer == 0) _writer = new MPFileHandler();

	  int size = 0;

	  if(get_parameter("pop_output_patch_sample_size")->isSet())
		  size = get_parameter_value("pop_output_patch_sample_size");

	  _writer->setOption(size);

	  Param* param = get_parameter("pop_output_logtime");

	  if(param->isMatrix()) {

		TMatrix temp;
		param->getMatrix(&temp);
		_writer->set_multi(true, true, 1, &temp, get_parameter("pop_output_dir")->getArg());

	  } else   //  rpl_per, gen_per, rpl_occ, gen_occ, rank (0), path, self-ref
		_writer->set(true, true, 1, (param->isSet() ? (int)param->getValue() : 0),
					 0, get_parameter("pop_output_dir")->getArg());

	  loader->attach(_writer);

  } else if(_writer) {
	  delete _writer;
	  _writer = NULL;
  }
}
// ----------------------------------------------------------------------------------------
// buildPatchArray
// ----------------------------------------------------------------------------------------
void Metapop::buildPatchArray()
{
  resizePatchArray();
  
  //set the population capacities:
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _vPatch[i]->flush(this);
    _vPatch[i]->init((unsigned int)_patchSizes.get(FEM,i), (unsigned int)_patchSizes.get(MAL,i), i);
  }
}
// ----------------------------------------------------------------------------------------
// resizePatchArray
// ----------------------------------------------------------------------------------------
void Metapop::resizePatchArray()
{
  //reset the right number of patches for the new simulation
  if(_vPatch.size() > _patchNbr) {
    while(_vPatch.size() > _patchNbr) {
      delete _vPatch[0];
      _vPatch.pop_front();
    }
  } else while(_vPatch.size() < _patchNbr) _vPatch.push_back(new Patch());
}
// ----------------------------------------------------------------------------------------
// updatePatchArray
// ----------------------------------------------------------------------------------------
void Metapop::updatePatchArray()
{
  //remove or add patches as needed:
  resizePatchArray();
  //reset the patch capacities and patch ID:
  updatePatchState();
}
// ----------------------------------------------------------------------------------------
// updatePatchState
// ----------------------------------------------------------------------------------------
void Metapop::updatePatchState()
{
  //reset the patch capacities and patch ID:
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _vPatch[i]->setID(i);
    _vPatch[i]->set_K((unsigned int)_patchSizes.get(FEM,i) + (unsigned int)_patchSizes.get(MAL,i));
    _vPatch[i]->set_KFem((unsigned int)_patchSizes.get(FEM,i));
    _vPatch[i]->set_KMal((unsigned int)_patchSizes.get(MAL,i));
  }
}
// ----------------------------------------------------------------------------------------
// setPatchCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setPatchCapacities()
{
  _patchSizes.reset(2, _patchNbr);
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _patchSizes.set(FEM, i, _patchKfem);
    _patchSizes.set(MAL, i, _patchKmal);
  }
}
// ----------------------------------------------------------------------------------------
// setPatchCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setPatchCapacities(string param)
{
  double* size_array;
  unsigned int Knum;
  TMatrix popK;
  
  _paramSet->getMatrix(param, &popK);
  
  Knum = popK.length();
  
  if(_patchNbr == 0 || _patchNbr < Knum)  _patchNbr = Knum;
  
  _patchSizes.reset(2, _patchNbr);
  
  size_array = popK.get();
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _patchSizes.set(FEM, i, (unsigned int)size_array[i % Knum]/2);
    _patchSizes.set(MAL, i, (unsigned int)size_array[i % Knum]/2);
  }
}
// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
void Metapop::setPatchCapacities(sex_t SEX, string param)
{
  double* size_array;
  unsigned int size_, Knum;
  TMatrix popK;
  
  _paramSet->getMatrix(param, &popK);
  
  Knum = popK.length();
  
  if(_patchNbr == 0 || _patchNbr < Knum)  _patchNbr = Knum;
  
  _patchSizes.reset(2, _patchNbr);
  
  size_array = popK.get();
  
  if(SEX == FEM)
    size_ = _patchKmal;
  else
    size_ = _patchKfem;
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _patchSizes.set(SEX, i, (unsigned int)size_array[i % Knum]);
    _patchSizes.set(!SEX, i, size_);
  }
}
// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
void Metapop::setPatchCapacities(string paramfem, string parammal)
{
  double *size_fem, *size_mal;
  unsigned int KFnum, KMnum;
  TMatrix popKfem, popKmal;
  
  _paramSet->getMatrix(paramfem,&popKfem);
  _paramSet->getMatrix(parammal,&popKmal);
  
  KFnum = popKfem.length();
  KMnum = popKmal.length();
  
  if(_patchNbr == 0 && KFnum != KMnum){
    warning("not same number of elements in females and males capacity matrices!\n");
    warning("setting the number of populations from size of longest capacity array (= %i).\n",max(KFnum, KMnum));
  }
  if(_patchNbr < max(KFnum, KMnum)) _patchNbr = max(KFnum, KMnum);
  
  _patchSizes.reset(2, _patchNbr);
  
  size_fem = popKfem.get();
  size_mal = popKmal.get();
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    _patchSizes.set(FEM, i, (unsigned int)size_fem[i % KFnum]);
    _patchSizes.set(MAL, i, (unsigned int)size_mal[i % KMnum]);
  }   
}
// ----------------------------------------------------------------------------------------
// store_data
// ----------------------------------------------------------------------------------------
void Metapop::store_data ( BinaryStorageBuffer* saver )
{
  unsigned int *sizes[3];
  unsigned char separator[2] = {'@','P'};
  
  sizes[0] = new unsigned int [_patchNbr]; 
  sizes[1] = new unsigned int [_patchNbr]; 
  sizes[2] = new unsigned int [_patchNbr];
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    sizes[0][i] = size(OFFSPRG,i);
    sizes[1][i] = size(FEM,ADULTS,i);
    sizes[2][i] = size(MAL,ADULTS,i);
  }
  //store the data, begin with pop separator and number of patches:
  saver->store(&separator, 2 * sizeof(unsigned char));
  
  saver->store(&_patchNbr, sizeof(unsigned int));
  //store the Patch sizes
  //offspring
  saver->store(sizes[0], _patchNbr * sizeof(unsigned int));
  //females adults
  saver->store(sizes[1], _patchNbr * sizeof(unsigned int));
  //males adults
  saver->store(sizes[2], _patchNbr * sizeof(unsigned int));
  
  int byte_count = saver->getByteLength();
  //record all individual informations: IDs, matings, etc.
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    //first offspring:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM,OFFSx); ++j)
      _vPatch[i]->get(FEM, OFFSx, j)->store_data(saver);	
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL,OFFSx); ++j)
      _vPatch[i]->get(MAL, OFFSx, j)->store_data(saver);	
    
    //then adults:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM, ADLTx); ++j)
      _vPatch[i]->get(FEM, ADLTx, j)->store_data(saver);
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL, ADLTx); ++j)
      _vPatch[i]->get(MAL, ADLTx, j)->store_data(saver);
  }
  
#ifdef _DEBUG_
  message("Metapop::store_data :stored %ikB of individual data (%i individuals)\n",
          (saver->getByteLength()-byte_count)/1024, size());
#endif
  
  //records the trait sequences:
  map<trait_t, TraitPrototype *> traits = this->getTraitPrototypes();
  map<trait_t, TraitPrototype *>::iterator tt = traits.begin();
  //  trait_t type;
  separator[1] = 'T'; //trait separator = '@T'
  
  byte_count = saver->getByteLength();
  
  while(tt != traits.end()) {
    saver->store(&separator, 2 * sizeof(unsigned char));
    //store the trait type:
    saver->store((void*)tt->first.c_str(), TRAIT_T_MAX); 
    
    //then ask the prototype to store its data:
    tt->second->store_data(saver);
    //store the traits data:
    store_trait(tt->second->get_index(), saver);
    
    tt++;
  }
#ifdef _DEBUG_
  message("Metapop::store_data :stored %ikB of traits data\n",
          (saver->getByteLength()-byte_count)/1024);
#endif
  
  for(unsigned int i = 0; i < 3; i++)
    delete [] sizes[i];
}
// ----------------------------------------------------------------------------------------
// store_trait
// ----------------------------------------------------------------------------------------
void Metapop::store_trait (int trait_idx, BinaryStorageBuffer* saver)
{  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    //first offspring:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM, OFFSx); ++j) {
      _vPatch[i]->get(FEM, OFFSx, j)->getTrait(trait_idx)->store_data(saver);	
    }
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL, OFFSx); ++j) {
      _vPatch[i]->get(MAL, OFFSx, j)->getTrait(trait_idx)->store_data(saver);	
    }
    //then adults:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM, ADLTx); ++j)
      _vPatch[i]->get(FEM, ADLTx, j)->getTrait(trait_idx)->store_data(saver);
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL, ADLTx); ++j)
      _vPatch[i]->get(MAL, ADLTx, j)->getTrait(trait_idx)->store_data(saver);
  }
}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
bool Metapop::retrieve_data ( BinaryStorageBuffer* loader )
{
#ifdef _DEBUG_ 
  message("Metapop::retrieve_data, %iB of data read so far\n",loader->getBytesOut()); 
#endif
  unsigned int *sizes[3], dummy_int;
  unsigned char separator[2];
  
  sizes[0] = new unsigned int [_patchNbr]; 
  sizes[1] = new unsigned int [_patchNbr]; 
  sizes[2] = new unsigned int [_patchNbr];
  
  loader->read(&separator, 2 * sizeof(unsigned char));
  
  if(separator[0] != '@' || separator[1] != 'P') {
    error("Binary file appears corrupted:\n >>>> Metapop::retrieve_data::wrong population seprarator\n");
    return false;
  }
  
  loader->read(&dummy_int, sizeof(unsigned int));
  
  if(dummy_int != _patchNbr) {
    error("Population in binary file differs from simulation settings:\n >>>> Metapop::retrieve_data:number of Patch differ from parameter value\n");
    return false;
  }
  //get the Patch sizes
  //offspring
  loader->read(sizes[0], _patchNbr * sizeof(unsigned int));
  //females adults
  loader->read(sizes[1], _patchNbr * sizeof(unsigned int));
  //males adults
  loader->read(sizes[2], _patchNbr * sizeof(unsigned int));
  
  Individual* ind;
  unsigned int bytes_cnt = loader->getBytesOut();
  //retrieve all individual informations: IDs, matings, etc.
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    //first, offspring:
    for(unsigned int j = 0; j < sizes[0][i]; ++j) {
      ind = this->getNewIndividual();
      ind->retrieve_data(loader);
      _vPatch[i]->add(ind->getSex(), OFFSx, ind);	
    }
    //then adult females:
    for(unsigned int j = 0; j < sizes[1][i]; ++j) {
      ind = this->getNewIndividual();
      ind->retrieve_data(loader);
      _vPatch[i]->add(FEM, ADLTx, ind);
    }
    //and adult males:
    for(unsigned int j = 0; j < sizes[2][i]; ++j) {
      ind = this->getNewIndividual();
      ind->retrieve_data(loader);
      _vPatch[i]->add(MAL, ADLTx, ind);
    }
    
  }
  
#ifdef _DEBUG_ 
  message("Metapop::retrieve_data::retrieved %ikB of ind data (%i individuals)\n",
          (loader->getBytesOut()-bytes_cnt)/1024, size()); 
#endif
  //retrieve traits sequences
  map<trait_t, TraitPrototype *>::iterator tt = _protoTraits.begin();
  
  char trait_name[6] = {'\0','\0','\0','\0','\0','\0'};
  
  unsigned int trait_cntr = 0;
  
  loader->read(&separator, 2 * sizeof(unsigned char));
  
  if(separator[0] != '@' || separator[1] != 'T') {
    error("Binary file appears corrupted:\n >>>> Metapop::retrieve_data::wrong trait seprarator\n");
    return false;
  }
  
  bytes_cnt = loader->getBytesOut();
  
  do {
    //get the trait type:
    loader->read(&trait_name[0], TRAIT_T_MAX);
    
    string dummy_trait(trait_name);
    
    //get the prototype:
    tt = _protoTraits.find(dummy_trait);
    
    if( tt == _protoTraits.end() ) {
	    error("Trait(s) in binary file differ from simulation settings:\n >>>> Metapop::retrieve_data::trait in file not present in protoype\n");
	    return false;
    }
    //then ask the prototype to retrieve its data:
    tt->second->retrieve_data(loader);
    
#ifdef _DEBUG_ 
    message("%iB of trait data read so far (trait %s)\n",loader->getBytesOut()-bytes_cnt,tt->first.c_str()); 
#endif
    
    //get the traits data:
    read_trait(tt->second->get_index(), loader);
    
    trait_cntr++;
    
    //get the next separator, should be present if right number of bytes have been read
    loader->read(&separator, 2 * sizeof(unsigned char));
    
    
  } while (separator[0] == '@' && separator[1] == 'T') ;
  
#ifdef _DEBUG_ 
  message("Metapop::retrieve_data::retrieved %ikB of trait data\n",(loader->getBytesOut()-bytes_cnt)/1024); 
#endif
  
  if(trait_cntr != _protoTraits.size()){
    error("Trait(s) in binary file differ from simulation settings:\n >>>> Metapop::retrieve_data::some traits are missing from binary file\n");
    return false;
  }
  
  if(separator[0] != '@') 
    error("Binary file appears corrupted:\n >>>> Metapop::retrieve_data::separator not found at end of pop record!\n");
  
  for(unsigned int i = 0; i < 3; i++)
    delete [] sizes[i];
  
  return true;
}
// ----------------------------------------------------------------------------------------
// read_trait
// ----------------------------------------------------------------------------------------
void Metapop::read_trait (int trait_idx, BinaryStorageBuffer* loader)
{ 
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    
    //first offspring:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM, OFFSx); ++j){
      _vPatch[i]->get(FEM, OFFSx, j)->getTrait(trait_idx)->retrieve_data(loader);	
    }
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL, OFFSx); ++j){
      _vPatch[i]->get(MAL, OFFSx, j)->getTrait(trait_idx)->retrieve_data(loader);	
    }
    
    //then adults:
    for(unsigned int j = 0; j < _vPatch[i]->size(FEM, ADLTx); ++j)
      _vPatch[i]->get(FEM, ADLTx, j)->getTrait(trait_idx)->retrieve_data(loader);
    
    for(unsigned int j = 0; j < _vPatch[i]->size(MAL, ADLTx); ++j)
      _vPatch[i]->get(MAL, ADLTx, j)->getTrait(trait_idx)->retrieve_data(loader);
  }
}
// ----------------------------------------------------------------------------------------
// setPopulation
// ----------------------------------------------------------------------------------------
void Metapop::setPopulation(unsigned int currentReplicate, unsigned int replicates)
{
#ifdef _DEBUG_
  message("Metapop::setPopulation: ");
  fflush(stdout);
#endif
  
  _currentReplicate = currentReplicate;
  
  //reset the population parameters, they may have changed during the life cycle
  setPopulationParameters();
  //reset the patch array, flush the patches and reset the patch capacities.
  buildPatchArray();
  
  _currentAge = NONE;
  
  //find the age class required to start the life cycle with:
  if(_requiredAge == NONE) 
    _requiredAge = SIMenv::MainSim->getFirstRequiredAgeInLifeCycle();
  
#ifdef _DEBUG_
  message("required age is: %i ", _requiredAge);
  fflush(stdout);
#endif
  //load first generation from a source population
  if(_source_load) {    
    
    //delete ind in recycler, needed when dealing with large source pop and lots of replicates:
    purgeRecyclingPOOL();
    
    _source_load_periodicity = (unsigned int)ceil((double)replicates/_source_replicates);
    
    if(currentReplicate == 1 ||
       _source_replicates >= replicates || //load from a different source every replicate
       !(_source_replicates != 0 && (currentReplicate - 1) % _source_load_periodicity) ) //change source every _source_replicates-nth replicates
      loadSourcePopulation();
    
    //load in preserve mode, i.e. make a copy of the saved pop
    if(_source_preserve) { 
      //basic checks:
      //not necessary anymore, after all...
//      if(_source->getPatchNbr() != _patchNbr)
//        fatal("loading population from source file failed: preserve mode forbids difference in number of patches\n");

      //have to remove check as K parameters might not reflect state of stored pop
      //this happens when temporal params are present in init files
//      if(_source->getPatchKFem() != _patchKfem || _source->getPatchKMal() != _patchKmal)
//        fatal("loading population from source file failed: preserve mode forbids difference in patch sizes\n");
      
      setPopulationFromSourceInPreserveMode();
      
    } else setPopulationFromSource();
    
  } else { //not source_load
    
    for(unsigned int i = 0; i < _patchNbr; ++i)
      _vPatch[i]->setNewGeneration( _requiredAge, this );
    
    _currentAge = _requiredAge;
    
  }
  
#ifdef _DEBUG_
  message("loaded age is: %i \n",_currentAge);
#endif
  
  if(_currentAge == NONE) warning("Metapop::setPopulation: generation 0 is empty!!\n");
}
// ----------------------------------------------------------------------------------------
// setPopulationFromSourceInPreserveMode
// ----------------------------------------------------------------------------------------
void Metapop::setPopulationFromSourceInPreserveMode()
{
  Patch* src_patch;
  age_t source_age;
  unsigned int num_patch = (unsigned int)min((double)_vPatch.size(), (double)_source->getPatchNbr());
   
  if(_requiredAge != NONE && _source->size(_requiredAge) == 0 && _source->size() != 0) {
    
    source_age = _source->getCurrentAge();
    
    warning("required age %i not present in source, using age class %i instead (preserve mode).\n"
            ,_requiredAge,source_age);
    
    for(unsigned int i = 0; i < num_patch; ++i) {
      src_patch = _source->getPatchPtr(i);
      fillPatchFromSource(FEM, src_patch, _vPatch[i], source_age);
      fillPatchFromSource(MAL, src_patch, _vPatch[i], source_age);
    }
    
  } else {
    
    if (_requiredAge == NONE) _requiredAge = ALL;
    
    for(unsigned int i = 0; i < num_patch; ++i) {
      src_patch = _source->getPatchPtr(i);
      if(_requiredAge & OFFSPRG) {
        fillPatchFromSource(FEM, src_patch, _vPatch[i], OFFSPRG);
        fillPatchFromSource(MAL, src_patch, _vPatch[i], OFFSPRG);
      }
      if(_requiredAge & ADULTS) {
        fillPatchFromSource(FEM, src_patch, _vPatch[i], ADULTS);
        fillPatchFromSource(MAL, src_patch, _vPatch[i], ADULTS);
      }
    }
  }
  
  if(size( ADULTS ) != 0)
    _currentAge |= ADULTS;
  
  if(size( OFFSPRG ) != 0)
    _currentAge |= OFFSPRG;
  
} 
// ----------------------------------------------------------------------------------------
// fillPatchFromSource
// ----------------------------------------------------------------------------------------
void Metapop::fillPatchFromSource(sex_t SEX, Patch* src, Patch* patch, age_t AGE)
{
  age_idx in_age = (AGE == OFFSPRG? OFFSx : ADLTx);
  age_idx to_age = (_requiredAge | AGE? in_age : ADLTx);//this to be sure we load some
  
  Individual* new_ind;
  
  for(unsigned int j = 0; j < src->size(SEX, in_age); ++j) {
    new_ind = getNewIndividual(); //this correctly sets pointers to TProto's, and params
    (*new_ind) = (*src->get(SEX, in_age, j)); //this should only copy genes
    patch->add(SEX , to_age, new_ind );
  }
}
// ----------------------------------------------------------------------------------------
// setPopulationFromSource
// ----------------------------------------------------------------------------------------
void Metapop::setPopulationFromSource()
{
  //load in non-preserve mode, i.e. use saved pop as a source of new individuals
  //individuals are drawn without replacement, sex and age are preserved.
  deque< Individual* > src_fem_pool;
  deque< Individual* > src_mal_pool;
  
  if( _requiredAge & OFFSPRG ) {
    
    if(_source->size(OFFSPRG) != 0)
      _source->getAllIndividuals(OFFSx, src_fem_pool, src_mal_pool);
    else {
      warning("source population does not contain offspring individuals, using adults instead.\n"); 
      //total size of source has been checked before
      _source->getAllIndividuals(ADLTx, src_fem_pool, src_mal_pool);
    }
    
    fillPopulationFromSource(OFFSx, FEM, src_fem_pool);
    fillPopulationFromSource(OFFSx, MAL, src_mal_pool);
    src_fem_pool.clear();
    src_mal_pool.clear();
    
    _currentAge |= OFFSPRG;
  }
  
  if( _requiredAge & ADULTS ) {
    
    if(_source->size(ADULTS) != 0)
      _source->getAllIndividuals(ADLTx, src_fem_pool, src_mal_pool);
    else {
      warning("source population does not contain adult individuals, using offspring instead.\n"); 
      //total size of source has been checked before
      _source->getAllIndividuals(OFFSx, src_fem_pool, src_mal_pool);
    }
    
    fillPopulationFromSource(ADLTx, FEM, src_fem_pool);
    fillPopulationFromSource(ADLTx, MAL, src_mal_pool);
    src_fem_pool.clear();
    src_mal_pool.clear();
    
    _currentAge |= ADULTS;
  }
  
} 
// ----------------------------------------------------------------------------------------
// fillPopulationFromSource
// ----------------------------------------------------------------------------------------
void Metapop::fillPopulationFromSource(age_idx AGE, sex_t SEX, deque<Individual*>& src_pool)
{
  unsigned int Ktot = 0;
  unsigned int ind_pos ;
  Individual* new_ind;
  
  for(unsigned int i = 0; i < _patchNbr; ++i) 
    Ktot += _vPatch[i]->get_K(SEX);
  
  if(src_pool.size() < Ktot)
    warning("Number of %s %s in source metapop (%i) not enough to fill current metapop (%i).\n",
            (SEX ? "female" : "male"), (AGE == OFFSx ? "offspring" : "adults"),
            src_pool.size(), Ktot);
  
  for(unsigned int i = 0; i < _patchNbr; ++i) {
    
    for(unsigned int j = 0, psize = _vPatch[i]->get_K(SEX); j < psize && src_pool.size() != 0; ++j) {
      
      new_ind = getNewIndividual();
      
      ind_pos = RAND::Uniform( src_pool.size() );
      
      _vPatch[i]->add( SEX, AGE, &( (*new_ind) = *(src_pool[ ind_pos ]) ) );
      
      src_pool.erase( src_pool.begin() + ind_pos );  //this is without replacement!!!!
      
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// loadSourcePopulation
// ----------------------------------------------------------------------------------------
void Metapop::loadSourcePopulation( )
{  
  std::string filename;
  
  if(_source_replicates != 0) {
    //find the replicate string of the source file
    ostringstream rpl;
    rpl.fill('0');
    rpl.width(_source_replicate_digits);
    
    if(_source_replicates >= _replicates)
      rpl << _currentReplicate + _source_start_at_replicate - 1;
    else if (_source_replicates != 1)
      //less replicates in source, reload the same repl every `source_replicates'
      rpl << ((_currentReplicate-1) / _source_load_periodicity) + 1;
    else
      rpl << 1; //a bit dumb but just in case...
    
    filename = _source_name + "_" + rpl.str() + _source_filetype;
    
  } else {
    
    filename = _source_name;
    
  }
  
  if(_source_filetype.compare(".bin") == 0) 
    loadPopFromBinarySource(filename);
  else
    loadPopFromTraitFile(filename);
  
  if(_source->size() == 0)  fatal("source population is empty!\n");
  else {
    if(_source->size(OFFSPRG) != 0) _source->setCurrentAge( (_source->getCurrentAge() | OFFSPRG) );
    if(_source->size(ADULTS)  != 0) _source->setCurrentAge( (_source->getCurrentAge() | ADULTS)  );
  }
  
#ifdef _DEBUG_
  message("source pop loaded: offspring: %i f, %i m; adults: %i f, %i m\n", _source->size(FEM,OFFSPRG),
          _source->size(MAL,OFFSPRG), _source->size(FEM,ADULTS), _source->size(MAL,ADULTS));
#endif
}
// ----------------------------------------------------------------------------------------
// loadPopFromBinarySource
// ----------------------------------------------------------------------------------------
void Metapop::loadPopFromBinarySource( string &filename )
{
  Individual  *new_ind, *src_tmp;
  
  if(!(_source = _loader.extractPop(filename,_source_generation,SIMenv::MainSim, this)))
    fatal("Metapop::loadPopFromBinarySource:could not extract pop from binary file \"%s\"\n",filename.c_str());
  
  //check traits parameters (mainly sizes)
  new_ind = getNewIndividual();
  src_tmp = _source->getNewIndividual();
  
  if( (*new_ind) != (*src_tmp) ) {
    fatal("loading population from source file failed: trait genetic architecture differs.\n");
  }
  
  delete new_ind;
  delete src_tmp;
  
//  if(_source_preserve) {
    
    //basic checks:
//    if(_source->getPatchNbr() != _patchNbr)
//      fatal("loading population from source file failed: preserve mode forbids difference in number of patches\n");
    
//    if(_source->getPatchKFem() != _patchKfem || _source->getPatchKMal() != _patchKmal)
//      fatal("loading population from source file failed: preserve mode forbids difference in patch sizes\n");
//  }
  
}

// ----------------------------------------------------------------------------------------
// loadPopFromTraitFile
// ----------------------------------------------------------------------------------------
void Metapop::loadPopFromTraitFile( string &filename )
{
  if(_source != NULL) delete _source;
  
  _source = new Metapop();
  _source->setMPImanager( _mpimgr );
  
  //the population's parameters from current params
  _source->set_paramsetFromCopy((*this->get_paramset()));
  
  //we need this to be able to create new individuals in the source:
  _source->init( );
  _source->makePrototype(this->getTraitPrototypes());
  
  //now read data from genotype file:
  FileServices *FS = SIMenv::MainSim->get_FileServices();
  
  FileHandler* file = FS->getReader(_source_filetype);
  
  if(!file) fatal("no file reader exists with source file extension \"%s\".\n",_source_filetype.c_str());
  
  file->set_pop_ptr(_source);
  
  file->FHread(filename);
  
  if(_source->size() == 0) fatal("source population in file \"%s\" is empty.\n",filename.c_str());
}
/***************************************************************
 * Metapop::reset
 * resets each Patch, all individuals are moved to the POOL
 */
void Metapop::reset()
{
  unsigned int i;
  
  for(i = 0; i < _patchNbr; ++i) {
    _vPatch[i]->flush(this);
    _vPatch[i]->reset_counters();
  }
}
// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
void Metapop::clear()
{
  unsigned int i;
  
  for(i = 0; i < _vPatch.size(); ++i) delete _vPatch[i];
  
  _vPatch.clear();
  
  _patchNbr = 0;
  
  //  purgeRecyclingPOOL();
  
  //  if(_source) delete _source;
  //  _source = NULL;
}
// ----------------------------------------------------------------------------------------
// setCurrentAge
// ----------------------------------------------------------------------------------------
void Metapop::setCurrentAge ( LifeCycleEvent* LCE )
{
  _currentAge ^= LCE->removeAgeClass();
  _currentAge |= LCE->addAgeClass();
}
// ----------------------------------------------------------------------------------------
// getAllIndividuals
// ----------------------------------------------------------------------------------------
void Metapop::getAllIndividuals(age_idx AGE, deque<Individual*>& fem_pool, deque<Individual*>& mal_pool)
{
  
  for(unsigned int i = 0; i < _vPatch.size(); ++i) {
    
    for(unsigned int j = 0, psize = _vPatch[i]->size(MAL, AGE); j < psize; ++j)
      mal_pool.push_back( _vPatch[i]->get(MAL, AGE, j) );
    
    for(unsigned int j = 0, psize = _vPatch[i]->size(FEM, AGE); j < psize; ++j)
      fem_pool.push_back( _vPatch[i]->get(FEM, AGE, j) );
    
  }
  
}

void Metapop::show_up()
{
  message("Metapop:\n");
  message("nbre. of patches: %i(%i)\n", _vPatch.size(),_patchNbr);
  message("population size : %i\n", size());
  message("K = %i, K_fem = %i, K_mal = %i\n",_patchK, _patchKfem, _patchKmal);
  message("patch capacities: \n");
  _patchSizes.show_up();
  message("Patches:\n");
  for(unsigned int i = 0; i < _vPatch.size(); i++)
    _vPatch[i]->show_up();
}

//------------------------------------------------------------------------------------------
//
//                                       MPFileHandler
//
//------------------------------------------------------------------------------------------
void MPFileHandler::FHwrite()
{
	if (!_pop->isAlive()) return; //don't write anything if pop is extinct


	std::string filename = get_filename();

	std::ofstream FILE (filename.c_str(), ios::out);

	if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());

	FILE << "ID dad mum sex age home pop\n";

	Patch* patch;

	if(_patch_sample_size == 0) {

		for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {

			patch = _pop->getPatch(i);

			for(unsigned int S = 0; S < 2; ++S){

				printNoSample((sex_t)S, OFFSx, patch, FILE);

				printNoSample((sex_t)S, ADLTx, patch, FILE);

			}
		}

	} else {

		for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {

			patch = _pop->getPatch(i);

			createAndPrintSample(OFFSx, patch, FILE);

			createAndPrintSample(ADLTx, patch, FILE);
		}
	}

	FILE.close();

}

void MPFileHandler::printNoSample (sex_t SEX, age_idx AGE, Patch* patch, ofstream& FH )
{
	Individual* ind;
	for(unsigned int j = 0; j < patch->size(SEX, AGE); ++j) {

		ind = patch->get(SEX, AGE, j);

		FH<< ind->getID() << " " << ind->getFatherID() << " " << ind->getMotherID()
			<< " " << SEX << " " << AGE << " " << ind->getHome()+1
			<< " " << patch->getID()+1 << endl;
	}
}

void MPFileHandler::createAndPrintSample (age_idx AGE, Patch* patch, ofstream& FH)
{
	Individual* ind;
	unsigned int size_fem = patch->size(FEM, AGE);
	unsigned int size_mal = patch->size(MAL, AGE);
	unsigned int sample_size;

	int* sample = NULL;

	//Females -----------------------------------------------------------------------------
	sample_size = (_patch_sample_size < size_fem ? _patch_sample_size : size_fem);

	if( sample_size != 0) {

		sample = new int[sample_size];

		RAND::Sample(0, size_fem, sample_size, sample, false);

		for(unsigned int j = 0; j < sample_size; ++j) {

			ind = patch->get(FEM, AGE, sample[j]);

			FH<< ind->getID() << " " << ind->getFatherID() << " " << ind->getMotherID()
				<< " " << FEM << " " << AGE << " " << ind->getHome()+1
				<< " " << patch->getID()+1 << endl;
		}
	}

	//Males -------------------------------------------------------------------------------
	sample_size = (_patch_sample_size < size_mal ? _patch_sample_size : size_mal);

	if( sample_size != 0) {

		if(sample) delete [] sample;

		sample = new int[sample_size];

		RAND::Sample(0, size_mal, sample_size, sample, false);

		for(unsigned int j = 0; j < sample_size; ++j) {

			ind = patch->get(MAL, AGE, sample[j]);

			FH<< ind->getID() << " " << ind->getFatherID() << " " << ind->getMotherID()
				<< " " << MAL << " " << AGE << " " << ind->getHome()+1
				<< " " << patch->getID()+1 << endl;
		}
	}

	if(sample) delete [] sample;
}
