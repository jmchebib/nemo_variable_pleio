/**  $Id: ttdispersal.cc,v 1.12 2016-03-03 15:18:47 fred Exp $
*
*  @file ttdispersal.cc
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
*  created on @date 05.08.2004
*
*  @author fred
*/
#include <string.h>
#include <cmath>
#include <iostream>
#include "ttdispersal.h"
#include "Uniform.h"
#include "output.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                               ****** TProtoDispersal ******

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoDispersal::TProtoDispersal(sex_t sex) 
: _mut_rate(0), _mut_mean(0.2), _init_rate_fem(-1.0), _init_rate_mal(-1.0), _gender(sex), _stats(0)
{
  _type = (_gender == FEM ? FDISP : MDISP);

  set_paramset("dispersal", false, this);
  
  add_parameter("disp_mutation_rate",DBL,true,true,0,1);
  add_parameter("disp_mutation_mean",DBL,true,true,0,1);
  add_parameter("disp_init_rate",DBL,false,true,0,1);
  add_parameter("disp_init_rate_fem",DBL,false,true,0,1);
  add_parameter("disp_init_rate_mal",DBL,false,true,0,1);
  add_parameter("disp_init_distribution", STR, false, false, 0, 0);
  add_parameter("disp_init_dist_params", MAT, false, false, 0, 0);
}
// ----------------------------------------------------------------------------------------
// copy cstor
// ----------------------------------------------------------------------------------------
TProtoDispersal::TProtoDispersal(const TProtoDispersal& TP)
: _mut_rate(TP._mut_rate), _mut_mean(TP._mut_mean), _init_rate_fem(TP._init_rate_fem), 
  _init_rate_mal(TP._init_rate_mal), _gender(TP._gender), _type(TP._type), _stats(0)
{
  _paramSet = new ParamSet( *(TP._paramSet) ) ;
}
// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TProtoDispersal::~TProtoDispersal()
{
  if(_stats != NULL) delete _stats; 
}
// ----------------------------------------------------------------------------------------
// setParameters
// ----------------------------------------------------------------------------------------
bool TProtoDispersal::setParameters ()
{
  _mut_rate = get_parameter_value("disp_mutation_rate");
  _mut_mean = get_parameter_value("disp_mutation_mean");
  
  if (( get_parameter("disp_init_rate")->isSet() || 
       get_parameter("disp_init_rate_fem")->isSet() || 
       get_parameter("disp_init_rate_mal")->isSet() ) )
  {
    
    if(get_parameter("disp_init_distribution")->isSet()) {
          warning("both \"disp_init_rate[...]\" and \"disp_init_distribution\" are set\n");
          warning("discarding \"disp_init_rate[...]\" parameters\n");
          return setRandom();
        }
    
    return setNonRandom();
    
  } else if (get_parameter("disp_init_distribution")->isSet()) {
  
    return setRandom();
    
  } else { //no init parameter specified
    
    Param* init_param = get_parameter("disp_init_distribution");
    
    init_param->setArg("uniform");
    init_param->setIsSet(true);
    
    init_param = get_parameter("disp_init_dist_params");
    
    init_param->setArg("{{0,0}}");
    init_param->setIsSet(true);
    
    return setRandom();
    
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// setNonRandom
// ----------------------------------------------------------------------------------------
bool  TProtoDispersal::setNonRandom()
{
  
  Param* init_param = get_parameter("disp_init_rate");
  
  if(init_param->isSet()) {
    
    _init_rate_fem = init_param->getValue();
    
    _init_rate_mal = init_param->getValue();
    
  } else if (get_parameter("disp_init_rate_fem")->isSet() &&
             get_parameter("disp_init_rate_mal")->isSet()) {
    
    _init_rate_fem = get_parameter_value("disp_init_rate_fem");
  
    _init_rate_mal = get_parameter_value("disp_init_rate_mal");
  
  } else {
    error("dispersal trait init parameters not correctly set\n");
    error("note that \"disp_init_rate_fem\" and \"disp_init_rate_mal\" must be set together\n");
    return false;
  }
  
  _init_random = false;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// setRandom
// ----------------------------------------------------------------------------------------
bool TProtoDispersal::setRandom()
{
  
  Param *init_param = get_parameter("disp_init_distribution");
  
  if (init_param->isSet()) {

    if (!get_parameter("disp_init_dist_params")->isSet()) {
      error("parameter \"disp_init_distribution\" requires \"disp_init_dist_params\" to be set as well.\n");
      return false;
    } else {
      get_parameter("disp_init_dist_params")->getMatrix(&_init_dist_param);
    }

    _init_dist = init_param->getArg();
    
    
  } else { //this should not happen...
    error("something got wrong with the dispersal trait init distribution parameters\n");
    return false;
  }
  
  _init_random = true;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoDispersal::loadStatServices ( StatServices* loader )
{
  //allocate the stat handler
  if(_stats == NULL)
    _stats = new TTDispersalSH(this);
  if(_stats != NULL) loader->attach(_stats);
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTDispersal* TProtoDispersal::hatch ( )
{
  TTDispersal* new_trait = new TTDispersal(_gender);
  
  new_trait->set_mut_rate(_mut_rate);
  new_trait->set_mut_mean(_mut_mean);
  new_trait->set_init_rate_fem(_init_rate_fem);
  new_trait->set_init_rate_mal(_init_rate_mal);
  new_trait->set_gender(_gender);
  new_trait->set_type(_type);
  new_trait->set_proto(this);
  
  return new_trait;
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                               ****** TTDispersal ******

// ----------------------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------------------
TTDispersal::TTDispersal(sex_t sex)
: _mut_rate(0), _mut_mean(0.2), _init_rate_fem(-1.0), _init_rate_mal(-1.0), _myProto(0), 
  _gender(sex)
{
  _type = (_gender == FEM ? FDISP : MDISP);
  _sequence[0] = _sequence[1] = 0;
}
// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTDispersal& TTDispersal::operator= (const TTrait& T)
{
  const TTDispersal& TD = dynamic_cast<const TTDispersal&> (T);
  
  if(this != &TD) {
    _gender = TD._gender;
    _type = TD._type;
    _sequence[0] = TD._sequence[0];
    _sequence[1] = TD._sequence[1];
    _phenotype = TD._phenotype;
  }
  
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTDispersal::operator== (const TTrait& T)
{
  if(_type.compare(T.get_type()) != 0) return false;

  const TTDispersal& TD = dynamic_cast<const TTDispersal&> (T);
  
  if(this != &TD) {
    if(_gender != TD._gender)  return false;
    if(_type != TD._type)      return false;
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTDispersal::operator!= (const TTrait& T)
{
  if( !((*this) == T) )
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
void TTDispersal::init_sequence ()
{
  TMatrix *params;
  unsigned int row;
  
  if ( _myProto->get_init_mode() ) { //this is random mode
    
    params = _myProto->get_init_dist_params();    
    
    if (params->getNbRows() == 2) // means different mean/sd for each sex
      row = _gender;
    else 
      row = 0;
    
    string option = _myProto->get_init_dist() ;
    
    if (option == "normal") {
      
      do { _sequence[0] = params->get(row, 0) + RAND::Gaussian(params->get(row, 1)); }
      while (_sequence[0] < 0.0 || _sequence[0] > 1.0);
      
      do { _sequence[1] = params->get(row, 0) + RAND::Gaussian(params->get(row, 1)); }
      while (_sequence[1] < 0.0 || _sequence[1] > 1.0);

    } else if (option == "uniform") {
      
      _sequence[0] = RAND::Uniform();
      _sequence[1] = RAND::Uniform();
      
    } else {
      fatal("dispersal init distribution \"%s\" is not implemented\n",_myProto->get_init_dist().c_str());
    }
    
    
  } else {
    
    _sequence[0] = (_gender == FEM ? _init_rate_fem : _init_rate_mal);
    _sequence[1] = (_gender == FEM ? _init_rate_fem : _init_rate_mal);
    
  }
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
void TTDispersal::inherit (TTrait* mother, TTrait* father)
{
  if(mother->get_type() != _type || father->get_type() != _type) 
	fatal("TTDispersal::inherit::wrong parent's trait type \n");
  
  _sequence[0] = mother->get_allele_value(0,RAND::RandBool());
  _sequence[1] = father->get_allele_value(0,RAND::RandBool());
}
// ----------------------------------------------------------------------------------------
// mutate
// ----------------------------------------------------------------------------------------
void TTDispersal::mutate ()
{
  double step;
  
  unsigned int nbMut = (unsigned int)RAND::Poisson(2*_mut_rate);
  
  unsigned int mut_allele;
  
  for(unsigned int i = 0; i < nbMut; i++) {
    
    step =  -1.0 * _mut_mean * log(1 - RAND::Uniform());
    
    mut_allele = RAND::RandBool();
    
    if(RAND::RandBool())
      _sequence[mut_allele] = ((_sequence[i] + step) >= 1.0 ? 1.0 : _sequence[mut_allele] + step);
    else
      _sequence[mut_allele] = ((_sequence[i] - step) <= 0.0 ? 0.0 : _sequence[mut_allele] - step);
    
  }

}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTDispersal::show_up ()
{
  set_value();
  message("\n  Trait type: %s\n\
       value: %f\n\
    _sequence: \n\
%f\n\
%f\n",get_type().c_str(), *(double*)getValue(), _sequence[0], _sequence[1]);
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTDispersalSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTDispersalSH::setStatRecorders ");
#endif
  if(token.compare("disp") == 0) {
	
	add("Dispersal rate (offspring)        ","off.disp",OFFSPRG,0,0,&TTDispersalSH::getOffsprgMeanDispRate,0,0,0);
	add("Dispersal rate (offspring females)","off.fdisp",OFFSPRG,0,0,&TTDispersalSH::getmeanOFD,0,0,0);
	add("Dispersal rate (offspring males)  ","off.mdisp",OFFSPRG,0,0,&TTDispersalSH::getmeanOMD,0,0,0);
	
	add("Dispersal rate (adults)       ","adlt.disp",ADULTS,0,0,&TTDispersalSH::getMeanDispRate,0,0,0);
	add("Dispersal rate (adult females)","adlt.fdisp",ADULTS,0,0,&TTDispersalSH::getmeanFD,0,0,0);
	add("Dispersal rate (adult males)  ","adlt.mdisp",ADULTS,0,0,&TTDispersalSH::getmeanMD,0,0,0);
	
  } else if(token.compare("adlt.disp") == 0) {
	
	add("Dispersal rate (adults)       ","adlt.disp",ADULTS,0,0,&TTDispersalSH::getMeanDispRate,0,0,0);
	add("Dispersal rate (adult females)","adlt.fdisp",ADULTS,0,0,&TTDispersalSH::getmeanFD,0,0,0);
	add("Dispersal rate (adult males)  ","adlt.mdisp",ADULTS,0,0,&TTDispersalSH::getmeanMD,0,0,0);
	
  } else if(token.compare("off.disp") == 0) {
	
	add("Dispersal rate (offspring)        ","off.disp",OFFSPRG,0,0,&TTDispersalSH::getOffsprgMeanDispRate,0,0,0);
	add("Dispersal rate (offspring females)","off.fdisp",OFFSPRG,0,0,&TTDispersalSH::getmeanOFD,0,0,0);
	add("Dispersal rate (offspring males)  ","off.mdisp",OFFSPRG,0,0,&TTDispersalSH::getmeanOMD,0,0,0);
	
  } else
	return false;
  
  return true;
}	
