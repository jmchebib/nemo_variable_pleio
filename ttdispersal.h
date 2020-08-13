/**  $Id: ttdispersal.h,v 1.11 2016-03-03 10:45:47 fred Exp $
*
*  @file ttdispersal.h
*  Nemo2
*
*   Copyright (C) 2006-2014 Frederic Guillaume
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
*  @author fred
*/
 

#ifndef TTDISPERSALGENE_H
#define TTDISPERSALGENE_H

#include "ttrait.h"
#include "types.h"
#include "stathandler.h"
#include "metapop.h"
#include "binarystoragebuffer.h"

class TTDispersalSH;
class TProtoDispersal;

/**Evolving dispersal trait, codes for female (_type = FDISP) or male (_type = MDISP) sex-specific dispersal rates.**/ 
class TTDispersal : public TTrait
{
  /**The allelic mutation rate.
   * The mutation distribution is exponential, centered on the allelic value.*/
  double _mut_rate;
  /**The mean mutation step.**/
  double _mut_mean;
  /**Initial allele for female dispersal.**/
  double _init_rate_fem;
  /**Initial allele for male dispersal.**/
  double _init_rate_mal;
  
  TProtoDispersal* _myProto;
  
  /**the gender of the trait, will determine its type.**/
  sex_t   _gender;
  /** The trait's type.**/
  trait_t _type;
  /** One diploid locus coding for a sex-specific dispersal rate.**/
  double _sequence[2];
  double _phenotype;

public:
  /** @param sex determines the type of this trait (FDISP for female dispersal, MDISP for male dispersal) **/
  TTDispersal (sex_t sex);
  TTDispersal (const TTDispersal& TP)
    : _mut_rate(TP._mut_rate), _mut_mean(TP._mut_mean), _init_rate_fem(TP._init_rate_fem), 
    _init_rate_mal(TP._init_rate_mal), _myProto(TP._myProto), _gender(TP._gender), _type(TP._type)
  {}
  virtual ~TTDispersal () { }
  
  ///@name Setters
  ///@{
  void set_mut_rate                    (double val)        {_mut_rate = val;}
  void set_mut_mean                    (double val)        {_mut_mean = val;}
  void set_init_rate_fem               (double val)        {_init_rate_fem = val;}
  void set_init_rate_mal               (double val)        {_init_rate_mal = val;}
  void set_gender                      (sex_t val)         {_gender = val;}
  void set_type                        (trait_t val)       {_type = val;}
  void set_proto                       (TProtoDispersal* P){_myProto = P;}
  ///@}
  ///@name Implementations
  ///@{
  virtual void    init                 ( )                 { _sequence[0] = 0.0; _sequence[1] = 0.0; }
  virtual void    init_sequence        ( );
  virtual void    reset                ( )                 {init();}
  virtual void    inherit              (TTrait* mother, TTrait* father);
  virtual void    mutate               ( );
  virtual trait_t get_type             ( ) const           {return _type;}
  virtual void    set_value            ( )                 {_phenotype = (_sequence[0] + _sequence[1])/2.0;}
  /** @return the dispersal rate, mean of the 2 alleles beared at this dispersal locus. **/
  virtual void*   getValue             ( ) const           {return (void*)&_phenotype;}
  /** @return NULL, here the _sequence is not a dble ptr. **/
  virtual void**  get_sequence         ( ) const           {return 0;}
  virtual double  get_allele_value     (int loc, int all)  {return ( !(all<2) ? 0 : _sequence[all] );}
  virtual void    set_allele_value (unsigned int locus, unsigned int allele, double value)
    {assert(allele < 2); _sequence[allele] = value;}
  virtual void    set_sequence         (void** seq)        { }
  virtual void*   set_trait            (void* value)       {return value;}
  virtual void    show_up              ( );
  virtual TTDispersal*  clone          ( )                 {return new TTDispersal(*this);}
  virtual TTDispersal& operator= (const TTrait& TP);
  virtual bool operator== (const TTrait& TP);
  virtual bool operator!= (const TTrait& TP);
  
  //implements StorableComponent: 
  virtual void    store_data       (BinaryStorageBuffer* saver)  {saver->store(&_sequence, 2 * sizeof(double));}
  virtual bool    retrieve_data    (BinaryStorageBuffer* reader) {reader->read(&_sequence, 2 * sizeof(double));return true;}
  ///@}
};
// ------------------------------------------------------------------------------

//                                     TProtoDispersal 

// ------------------------------------------------------------------------------
/**Prototype of the evolving dispersal trait, defines the sex-specific trait type.*/
class TProtoDispersal : public TraitPrototype {
public:
  TProtoDispersal(sex_t sex);
  TProtoDispersal(const TProtoDispersal& TP);
  ~TProtoDispersal();
  //implements TraitPrototype:
  virtual   void              init (){}
  virtual   void             reset (){}
  virtual   TTDispersal*     hatch ();
  virtual   TProtoDispersal* clone () {return new TProtoDispersal(*this);}
  virtual   trait_t       get_type () const {return _type;}
  //implements StorableComponent: 
  virtual void    store_data       (BinaryStorageBuffer* saver)  {saver->store(&_type, sizeof(trait_t));saver->store(&_gender, sizeof(sex_t));}
  virtual bool    retrieve_data    (BinaryStorageBuffer* reader) {reader->read(&_type, sizeof(trait_t));reader->read(&_gender, sizeof(sex_t));return true;}
  //implements SimComponent:
  virtual bool setParameters();
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader );
  
  bool setNonRandom ();
  bool setRandom    ();
  
  bool get_init_mode () {return _init_random;}
  TMatrix* get_init_dist_params () {return &_init_dist_param;}
  string get_init_dist () {return _init_dist;}
  
private:
  /**The allelic mutation rate. **/
  double _mut_rate;
  /**Mean mutation step. **/
  double _mut_mean;
  /**Initial allele for female dispersal.**/
  double _init_rate_fem;
  /**Initial allele for male dispersal.**/
  double _init_rate_mal;
  
  string _init_dist;
  TMatrix _init_dist_param;
  bool _init_random;
  
  /**The gender of the trait, will determine its type.**/
  sex_t   _gender;
  /** The trait's type.**/
  trait_t _type;
  /** The trait's StatHandler. **/
  TTDispersalSH* _stats;
};

// ------------------------------------------------------------------------------

//                                      TTDispersalSH

// ------------------------------------------------------------------------------
/**The StatHandler for the evolving dispersal traits.**/
class TTDispersalSH: public StatHandler<TTDispersalSH> {
  
  TProtoDispersal* _trait;
  
  double _meanFemDisp, _meanMalDisp, _meanOffFemDisp, _meanOffMalDisp;
  
  int _fdispIdx, _mdispIdx;
    
public:
  
  TTDispersalSH (TProtoDispersal* TT) 
  : _trait(TT),_meanFemDisp(0), _meanMalDisp(0), _meanOffFemDisp(0), _meanOffMalDisp(0), _fdispIdx(-1), _mdispIdx(-1) { }
  
  virtual ~TTDispersalSH ( ) { }
  virtual void init ( ) 
  { 
    StatHandlerBase::init(); 
    _fdispIdx = _pop->getTraitIndex(FDISP); 
    _mdispIdx = _pop->getTraitIndex(MDISP);
  }

  virtual bool setStatRecorders (std::string& token);
  
  
  double getmeanOFD                        ()            {return _meanOffFemDisp;}
  double getmeanOMD                        ()            {return _meanOffMalDisp;}
  double getmeanFD                         ()            {return _meanFemDisp;}
  double getmeanMD                         ()            {return _meanMalDisp;}
  double getMeanDispRate                   ();
  double getOffsprgMeanDispRate            ();
  double getMeanDispRate                   (bool);
  double getOffsprgMeanDispRate            (bool);
  double getMeanFemDispRate                ();
  double getMeanMalDispRate                ();
  
};
#endif //TTDISPERSALGENE_H

