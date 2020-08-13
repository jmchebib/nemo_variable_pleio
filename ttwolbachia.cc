/** $Id: ttwolbachia.cc,v 1.17 2016-10-31 14:37:18 fred Exp $
*
*
*  @file ttwolbachia.cc
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

#include <sstream>
#include "ttwolbachia.h"
#include "metapop.h"
#include "simenv.h"
#include "output.h"
#include "utils.h"

// ----------------------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------------------
TProtoWolbachia::TProtoWolbachia () 
: _transmit_rate(0), _stats(0)
{
  set_paramset("wolbachia", false, this);
  add_parameter("wolbachia_transmission_rate",DBL,true,false,0,1);
}
// ----------------------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------------------
TProtoWolbachia::TProtoWolbachia (const TProtoWolbachia& TP) 
: _transmit_rate(TP._transmit_rate), _stats(0)
{
  _paramSet = new ParamSet( *(TP._paramSet) ) ;
}
// ----------------------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------------------
TProtoWolbachia::~TProtoWolbachia() 
{
  if(_stats != NULL) delete _stats;
}
// ----------------------------------------------------------------------------------------
// TProtoWolbachia::loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoWolbachia::loadStatServices ( StatServices* loader )
{
  if(_stats != NULL)
    delete _stats;
  _stats = new TTWolbachiaSH(this);
  loader->attach(_stats);
}
// ----------------------------------------------------------------------------------------
// TTWolbachia::operator=
// ----------------------------------------------------------------------------------------
TTWolbachia& TTWolbachia::operator= (const TTrait& T)
{
  const TTWolbachia& TD = dynamic_cast<const TTWolbachia&> (T);
  if(this != &TD) {
    _is_infected = TD._is_infected;
//    _transmit_rate = TD._transmit_rate;
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// TTWolbachia::operator==
// ----------------------------------------------------------------------------------------
bool TTWolbachia::operator== (const TTrait& T)
{
  if(this->get_type().compare(T.get_type()) != 0) return false;
//  const TTWolbachia& TD = dynamic_cast<const TTWolbachia&> (T);
//  if(this != &TD) {
//    if(_transmit_rate != TD._transmit_rate) return false;
//  }
  return true;
}
// ----------------------------------------------------------------------------------------
// TTWolbachia::operator!=
// ----------------------------------------------------------------------------------------
bool TTWolbachia::operator!= (const TTrait& T)
{
  if(!((*this) == T))
    return true;
  else
    return false;
}

// ------------------------------------------------------------------------------

//                             TTWolbachiaSH

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTWolbachiaSH::setStatRecorders(string& token)
{
//  cout << "TTWolbachiaSH::setStatRecorders("<<token<<")\n";

  if(token == "wolbachia") {
    
	  add("Wolbachia female infection (offsprng)","off.fwoinf",OFFSPRG,1,0,0,&TTWolbachiaSH::getMeanOffsprgInfection,0,0);
	  add("Wolbachia male infection   (offsprng)","off.mwoinf",OFFSPRG,0,0,0,&TTWolbachiaSH::getMeanOffsprgInfection,0,0);
    //	  add("Wolbachia incompatible mating freq   ","off.incmating",OFFSPRG,0,&TTWolbachiaSH::getIcompatibleMatingFreq,0,0,0);
	  add("Wolbachia female infection (adults)  ","adlt.fwoinf",ADULTS,1,0,0,&TTWolbachiaSH::getMeanInfection,0,&TTWolbachiaSH::setInfectionStats);
	  add("Wolbachia male infection   (adults)  ","adlt.mwoinf",ADULTS,0,0,0,&TTWolbachiaSH::getMeanInfection,0,0);
    add("Wolbachia demic infection variance   ","wolb.infvar",ADULTS,0,0,&TTWolbachiaSH::getDemicInfectionVar,0,0,0);
    add("Wolbachia demic extinction rate      ","wolb.extrate",ADULTS,0,0,&TTWolbachiaSH::getDemicExtinctionRate,0,0,0);
    
  } else if(token == "wolbachia_perpatch") {
    
    ostringstream name, sub_name;
    
    //	for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    //	  _name<<"Wolbachia fem infctn: patch "<<i+1;
    //	  sub_name<<"off.p"<<i+1<<"fwoinf";
    //	  add(_name.str(),sub_name.str(),OFFSPRG,i,0,0,&TTWolbachiaSH::getMeanOffsprgFemaleInfection_perPatch,0);
    //	  _name.str("");
    //	  _name<<"Wolbachia mal infctn: patch "<<i+1;
    //	  sub_name.str("");
    //	  sub_name<<"off.p"<<i+1<<"mwoinf";
    //	  add(_name.str(),sub_name.str(),OFFSPRG,i,0,0,&TTWolbachiaSH::getMeanOffsprgMaleInfection_perPatch,0);
    //	  sub_name.str("");
    //	  _name.str("");
    //	}
    
    for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
      name<<"Patch "<<i+1;
      sub_name<<"adlt.p"<<i+1<<"fwoinf";
      add(name.str(),sub_name.str(),ADULTS,i,0,0,&TTWolbachiaSH::getMeanFemaleInfection_perPatch,0,0);
      //	  sub_name.str("");
      //	  sub_name<<"adlt.p"<<i+1<<"mwoinf";
      //	  add(_name.str(),sub_name.str(),ADULTS,i,0,0,&TTWolbachiaSH::getMeanMaleInfection_perPatch,0);
      sub_name.str("");
      name.str("");
    }
    
  }else
    return false;
  
  return true;
}


// ----------------------------------------------------------------------------------------

//                        LCE_Breed_Wolbachia

//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia
// ----------------------------------------------------------------------------------------
LCE_Breed_Wolbachia::LCE_Breed_Wolbachia () : LifeCycleEvent("breed_wolbachia",WOLB),
_incomp_cost(0), _fec_cost(0), _infected_fec(0), _inoculum_size(0), _inoculum_time(0), 
_writer(0)
{
  ParamUpdater< LCE_Breed_Wolbachia > * updater =
  new ParamUpdater< LCE_Breed_Wolbachia > (&LCE_Breed_Wolbachia::setParameters);
  add_parameter("wolbachia_fecundity_cost",DBL,true,true,-1,1, updater);
  add_parameter("wolbachia_incompatibility_cost",DBL,true,true,-1,1, updater);
  add_parameter("wolbachia_inoculum_size",MAT,true,false,0,0, updater);
  add_parameter("wolbachia_inoculum_time",INT,true,false,0,0, updater);
  add_parameter("wolbachia_model",INT,false,true,1,2, updater);
  add_parameter("wolbachia_output_dir", STR, false, false, 0, 0);
}

LCE_Breed_Wolbachia::~LCE_Breed_Wolbachia ( ) 
{
  if(_writer) delete _writer; 
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Wolbachia::setParameters ()
{
  if(!LCE_Breed_base::setParameters( )) return false;
  
  _fec_cost = this->get_parameter_value("wolbachia_fecundity_cost");
  _infected_fec = this->getMeanFecundity(0) * (1 - _fec_cost);
  _incomp_cost = this->get_parameter_value("wolbachia_incompatibility_cost");
  _inoculum_time = (unsigned int)this->get_parameter_value("wolbachia_inoculum_time");
  
  TMatrix tmp;

  get_parameter("wolbachia_inoculum_size")->getMatrix(&tmp);

  if( !_inoculum_size ) _inoculum_size = new TMatrix();

  if(tmp.ncols() == 2 && tmp.nrows() == _popPtr->getPatchNbr())
	  _inoculum_size->copy(tmp);
  else
	  setSpatialMatrix("wolbachia_inoculum_size", "2", &tmp, _inoculum_size, 2, _popPtr->getPatchNbr());

  if (get_parameter("wolbachia_model")->isSet()) {
    _model = (unsigned int)get_parameter_value("wolbachia_model");
  } else
    _model = 1;
  
  
  switch(_model) {
    case 1: 
      _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_1;
      break;
    case 2:
      _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_2;
      break;
    default:  _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_1;
  }
  
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::loadFileServices
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::loadFileServices ( FileServices* loader )
{
  if(_writer != NULL) delete _writer;
  
  _writer = new TTWolbachiaFH(this);
  
  _writer->set(false, false, SIMenv::getReplicates(), SIMenv::getGenerations(), 0,
               get_parameter("wolbachia_output_dir")->getArg(), this);
  
  loader->attach(_writer);
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::execute ()
{
  
#ifdef _DEBUG_
  message("LCE_Breed_Wolbachia::execute\n");
#endif
  double infection_status;
  
  if( _popPtr->getCurrentGeneration() == _inoculum_time) {
	   inoculate_wolbachia();}
  
  if( _popPtr->getCurrentGeneration() > _inoculum_time) {
    
    infection_status = hasInfectedFemale();
    
    if( infection_status == 0 || infection_status == 1) {
      _writer->record(_popPtr->getCurrentReplicate(), _popPtr->getCurrentGeneration(), infection_status);
      _popPtr->reset();
      return;
    }}
  
  
  if(_popPtr->size(OFFSPRG) != 0) {
    warning("offspring containers not empty at time of breeding, flushing.\n");
    _popPtr->flush(OFFSx);
  }
  
  (this->*_breed_func_ptr)();
  
  if (_popPtr->getCurrentGeneration() == SIMenv::getGenerations() ) {
    _writer->record(_popPtr->getCurrentReplicate(), _popPtr->getCurrentGeneration(), hasInfectedFemale());
  }
  
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::wolbachia_model_1
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::wolbachia_model_1 ()
{
  Patch* current_patch;
  unsigned int indexOfMother, nbBaby;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->size(FEM, ADLTx) == 0 || current_patch->size(MAL, ADLTx) == 0) continue;
    
    for(indexOfMother = 0; indexOfMother < current_patch->size(FEM, ADLTx); indexOfMother++) {
      
      MotherPtr = current_patch->get(FEM, ADLTx, indexOfMother);
      
      if(*(bool*)MotherPtr->getTraitValue(_LCELinkedTraitIndex))
        nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity( _infected_fec ) );
      else
        nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity(i) );
      //-----------------------------------------------------------------------
      while(nbBaby != 0) {
        
        FatherPtr = getFatherPtr(current_patch, MotherPtr, indexOfMother);
        
        NewOffsprg = _popPtr->makeOffsprg(MotherPtr,FatherPtr,(sex_t)RAND::RandBool(),i);
        
        if(!(*(bool*)NewOffsprg->getTraitValue(_LCELinkedTraitIndex)) &&
           *(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex)) {
          if(RAND::Uniform() < _incomp_cost) {
            _popPtr->recycle(NewOffsprg);
          } else
            current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        } else
          current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        
        nbBaby--;
      }//_END_WHILE nbBaby
    }//end_for indexOfMother
  }//end_for patch
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::wolbachia_model_2
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::wolbachia_model_2 ()
{
  Patch* current_patch;
  unsigned int indexOfMother, nbBaby;
  double fec;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->size(FEM, ADLTx) == 0 || current_patch->size(MAL, ADLTx) == 0) continue;
    
    for(indexOfMother = 0; indexOfMother < current_patch->size(FEM, ADLTx); indexOfMother++) {
      
      MotherPtr = current_patch->get(FEM, ADLTx, indexOfMother);
      
      FatherPtr = getFatherPtr(current_patch, MotherPtr, indexOfMother);
      
      if(*(bool*)MotherPtr->getTraitValue(_LCELinkedTraitIndex)) {
        if(*(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex))
          fec = _infected_fec * (1 - _incomp_cost);
        else
          fec = _infected_fec;
      } else if (*(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex)) {
        fec = getMeanFecundity(i) * (1 - _incomp_cost);
      } else {
        fec = getMeanFecundity(i);
      }
      
      nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity( fec ) );
      
      //-----------------------------------------------------------------------
      for(;nbBaby != 0;nbBaby--) {
        
        NewOffsprg = _popPtr->makeOffsprg(MotherPtr,FatherPtr,(sex_t)RAND::RandBool(),i);
        
        current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        
      }//end_for nbBaby
    }//end_for indexOfMother
  }//end_for patch
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::inoculate_wolbachia
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::inoculate_wolbachia ()
{
  Patch* current_patch;
  bool T = 1;
  
  for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->get_isExtinct()) continue;
    
      for (unsigned int j = 0;
    		  j < _inoculum_size->get(i, 0) && j < current_patch->size(FEM, ADLTx); j++) {
        current_patch->get(FEM, ADLTx, j)->setTrait(_LCELinkedTraitIndex, &T);
      }

    for (unsigned int j = 0;
    		j < _inoculum_size->get(i, 1) && j < current_patch->size(MAL, ADLTx); j++) {
        current_patch->get(MAL, ADLTx, j)->setTrait(_LCELinkedTraitIndex, &T);
    }
  }
//  if(!ok) fatal("could not inoculate wolbachia, check the inoculum size!\n");
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::getFemaleInfection
// ----------------------------------------------------------------------------------------
double LCE_Breed_Wolbachia::hasInfectedFemale ()
{
  double indNbr = 0, mean = 0, size;
  Patch* crnt_patch;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {
    crnt_patch = _popPtr->getPatch(i);
    size = crnt_patch->size(FEM, ADLTx);
    indNbr += size;
    for(unsigned int j = 0; j < size; ++j)
      mean += (double)*(bool*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_LCELinkedTraitIndex);
  }
  
  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}

// ----------------------------------------------------------------------------------------

//                             TTWolbachiaFH

// ----------------------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------------------
void TTWolbachiaFH::record (unsigned int repl, unsigned int gen, double infection)
{
  _times[repl] = gen;
  _rate.push_back(infection);
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTWolbachiaFH::FHwrite ()
{
  string filename = get_filename();
  
  ofstream FILE (filename.c_str(), ios::out); 
  
  if(!FILE) fatal("could not open wolbachia output file \"%s\"\n", filename.c_str());
  
  FILE << "replicate\t" << "generation\t" << "infection_rate\n";
  
  unsigned int i = 0;
  
  map<unsigned int, unsigned int>::iterator IT = _times.begin();
  
  while (IT != _times.end() && i < _rate.size() ) {
    FILE << IT->first <<"\t"<< IT->second << "\t" << _rate[i] << endl;
    i++;
    IT++;
  }
  
  FILE.close();
  
}



