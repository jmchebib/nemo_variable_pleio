/**  $Id: ttwolbachia.h,v 1.15 2016-10-31 14:37:18 fred Exp $
*
*  @file ttwolbachia.h
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
 

#ifndef TTWOLBACHIA_H
#define TTWOLBACHIA_H

#include "ttrait.h"
#include "lifecycleevent.h"
#include "filehandler.h"
#include "stathandler.h"
#include "binarystoragebuffer.h"
#include "Uniform.h"
#include "LCEbreed.h"

class TTWolbachiaSH;

// ------------------------------------------------------------------------------

//                                TTWolbachia 

// ------------------------------------------------------------------------------
/**Trait used to study the dynamics of spread of Wolbachia, an endosymbiotic parasite causing
 * cytoplasmic incompatibility. The trait state is given by a unique haploid "gene" with two
 * alleles, 0 = uninfected and 1 = infected.
 */
class TTWolbachia : public TTrait
{
private:
  double _transmit_rate;  
  bool _is_infected;
  
public:
	
    TTWolbachia () 
    : _transmit_rate(0), _is_infected(0) { }
  
  TTWolbachia (const TTWolbachia& T) 
    : _transmit_rate(T._transmit_rate), _is_infected(T._is_infected) { }
  
  
  virtual ~TTWolbachia				  () {}
  
  void set_transmit_rate              (double val) {_transmit_rate = val;}
  
  virtual void    init				    ()  {_is_infected = 0;}
  virtual void    init_sequence		()  {_is_infected = 0;}
  virtual void    reset				    ()  {_is_infected = 0;}
  virtual void    inherit         (TTrait* mother, TTrait* father)
  {  _is_infected = *(bool*)mother->getValue();  }
  virtual void    mutate          ()
  {  _is_infected &= (RAND::Uniform() < _transmit_rate);  }  
  virtual void*   set_trait			  (void* value)
  {  _is_infected = *(bool*)value; return &_is_infected; }
  virtual void    set_sequence    (void** seq)  { }
  virtual void    set_value       ()  { }
  virtual void*   getValue			  () const {return (void*)&_is_infected;}
  virtual trait_t get_type			  () const {return WOLB;}
  virtual void**  get_sequence    () const {return NULL;}  
  virtual double  get_allele_value(int loc, int all) {return _is_infected;}  
  virtual void    set_allele_value(unsigned int locus, unsigned int allele, double value)  {_is_infected=(bool)value;}  
  virtual void    show_up			    ()       {}
  virtual TTWolbachia*  clone		  ()       {return new TTWolbachia(*this);}
  virtual TTWolbachia& operator=  (const TTrait& T);
  virtual bool    operator==      (const TTrait& T);
  virtual bool    operator!=      (const TTrait& T);
    
  virtual void    store_data      (BinaryStorageBuffer* saver)  
  { 
    unsigned char dummy = static_cast<unsigned char>(_is_infected);
    saver->store(&dummy, 1);
  }
  virtual bool    retrieve_data   (BinaryStorageBuffer* reader)
  {
    unsigned char dummy;
    reader->read(&dummy, 1);
    _is_infected = dummy;
    return true;
  }
  
};

// ------------------------------------------------------------------------------

//                                 TProtoWolbachia  

// ------------------------------------------------------------------------------
/**Prototype of the Wolbachia trait.*/
class TProtoWolbachia : public TraitPrototype {
public:
  TProtoWolbachia();
  TProtoWolbachia(const TProtoWolbachia&);
  ~TProtoWolbachia();
  
  virtual void init(){}
  virtual void reset(){}
  virtual bool setParameters()
  { _transmit_rate =  get_parameter_value("wolbachia_transmission_rate"); return true;}
  
  virtual TTWolbachia* hatch() 
  {
    TTWolbachia* new_trait = new TTWolbachia();
    new_trait->set_transmit_rate(_transmit_rate);
    return new_trait;
  }
  
  virtual TProtoWolbachia* clone() {return new TProtoWolbachia(*this);}
  
  virtual trait_t get_type ( ) const {return WOLB;}
  
  virtual void store_data (BinaryStorageBuffer* saver)     {/*we have nothing to save...*/}
  virtual bool retrieve_data (BinaryStorageBuffer* reader) {return true;}
  
  virtual void    loadFileServices ( FileServices* loader ) {}
  virtual void    loadStatServices ( StatServices* loader );
  
private:
  
    double _transmit_rate;
  TTWolbachiaSH* _stats;
};

// ------------------------------------------------------------------------------

//                                TTWolbachiaSH 

// ------------------------------------------------------------------------------
/**StatHandler of the Wolbachia trait.*/
class TTWolbachiaSH : public TraitStatHandler<TProtoWolbachia, TTWolbachiaSH> {
  
  TProtoWolbachia* _trait;
  int _TTidx;
  double _Fmean, _Mmean, _var, _extrate;
public:
  
  TTWolbachiaSH (TProtoWolbachia* TT) 
  : TraitStatHandler<TProtoWolbachia, TTWolbachiaSH> (TT),
  _trait(TT), _TTidx(TT->get_index()), _Fmean(0), _Mmean(0), _var(0), _extrate(0) {}
  
  virtual ~TTWolbachiaSH ( ) { }
  
  virtual bool setStatRecorders (string& token);
  
  void   setInfectionStats                      ( );
  double getMeanInfection                       (unsigned int sex) {return ((bool)sex ? _Fmean : _Mmean);}
  double getMeanOffsprgInfection                (unsigned int sex);
  double getMeanFemaleInfection_perPatch        (unsigned int patch);
  double getMeanMaleInfection_perPatch          (unsigned int patch);
  double getMeanOffsprgFemaleInfection_perPatch (unsigned int patch);
  double getMeanOffsprgMaleInfection_perPatch   (unsigned int patch);
  double getIcompatibleMatingFreq				( );
  double getDemicInfectionVar                   ( ) {return _var;}
  double getDemicExtinctionRate                 ( ) {return _extrate;}
};



class TTWolbachiaFH;


/**Breeding LCE when individuals carry the Wolbachia endosymbiotic parasite. 
 * See the user manual for an explanation of the parameters defined here.
 */
class LCE_Breed_Wolbachia : public virtual LCE_Breed_base  {
  
  double _incomp_cost;
  double _fec_cost;
  double _infected_fec;
  TMatrix* _inoculum_size;
  unsigned int _inoculum_time;
  unsigned int _model;
  
  void (LCE_Breed_Wolbachia::* _breed_func_ptr) ();
  
  TTWolbachiaFH* _writer;
  
  void inoculate_wolbachia();
  double hasInfectedFemale ();
  
public:
	
	LCE_Breed_Wolbachia ( );
  
  virtual ~LCE_Breed_Wolbachia ( );
  
  void wolbachia_model_1 ();
  void wolbachia_model_2 ();
  
  virtual bool setParameters ();
  virtual void  execute ();
  
  virtual LifeCycleEvent* clone ( ) 
  { return new LCE_Breed_Wolbachia(); }
  
  
  virtual void loadFileServices ( FileServices* loader );
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
};

// ------------------------------------------------------------------------------

//                               ****** TTWolbachiaFH ******

// ------------------------------------------------------------------------------
/**FileHandler of the Wolbachia trait.
 Records the last generation of each replicate; i.e. the generation at which infection is 
 either lost, fixed, or incomplete due to generation time limit.
 */
class TTWolbachiaFH: public EventFileHandler<LCE_Breed_Wolbachia> {
  
  map< unsigned int, unsigned int > _times;
  vector< double > _rate;
  
public:
	
  TTWolbachiaFH ( LCE_Breed_Wolbachia* TP )
  : EventFileHandler<LCE_Breed_Wolbachia> (TP, ".wolb")
  { }
  
  //TTWolbachiaFH ( ) { }
  virtual ~TTWolbachiaFH ( ) { }
  
  void record (unsigned int repl, unsigned int gen, double infection);
  
  virtual void FHwrite  ();
  
  virtual void FHread (string& filename) {}
  
};

#endif //TTWOLBACHIA_H

