/**  $Id: ttbdmi.h,v 1.6 2016-09-28 15:00:27 fred Exp $
 *
 *  @file ttbdmi.h
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
 *  Created on @date 21.12.2010
 *  @author fred
 */

#ifndef _TTDBMI_H_
#define _TTDBMI_H_

#include "ttrait.h"
#include "types.h"
#include "stathandler.h"
#include "filehandler.h"
#include "lifecycleevent.h"
#include "metapop.h"
#include "bitstring.h"
#include "tmatrix.h"
#include "Uniform.h"

class TTBDMI_SH;
class TTBDMI_FH;
class TT_BDMI;
// ------------------------------------------------------------------------------

//   TProtoBDMI

// ------------------------------------------------------------------------------
class TProtoBDMI : public TTProtoWithMap {
  
  trait_t _type;
  
  bool _isHaploid;
  unsigned int _nb_locus;
  double _mut_rate;
  double _genomic_mut_rate;
  double _recomb_rate;
  double* _init_freq;
  bool    _isInitSet;
  static unsigned int _diploGenotTableCoding[3][3];
  
  void   (TT_BDMI::* _inherit_func_ptr)  (TTrait*, TTrait*);
  void   (TT_BDMI::* _mutation_func_ptr) (void);
  double (TT_BDMI::* _viability_func_ptr) (void);

  TMatrix *_genoTable;
  TTBDMI_SH* _stater;
  TTBDMI_FH* _writer;
  
public:
  
  TProtoBDMI();
  TProtoBDMI(const TProtoBDMI& TP);
  virtual ~TProtoBDMI();
  
  int             get_nb_locus     ( )  {return _nb_locus;}
  double          get_mut_rate     ( )  {return _mut_rate;}
  bool            isHaploid        ( )  {return _isHaploid;}
  void            set_init_freq    (double* val, unsigned int size);
  double          get_init_freq    (unsigned int i) {return _init_freq[i];}
  bool            isInitSet        ( )  {return _isInitSet;}

  double          getGenoFitnessHaplo (unsigned int row, unsigned int pos){
	  assert(pos < 4);
	  return _genoTable->get(row, pos);
  }
  
  double          getGenoFitnessDiplo (unsigned int row, unsigned int posA, unsigned int posB){
    return _genoTable->get(row, _diploGenotTableCoding[posA][posB] );
  }
  
  double          getGenoFitnessDiplo (unsigned int row, unsigned int pos){
	  assert(pos < 9);
      return _genoTable->get(row, pos);
  }

  void            setGenoFitnessValue(unsigned int row, unsigned int geno, double value)
  { _genoTable->set(row, geno, value); }

  void showGenoTable (unsigned int nrows);
  
  void            inherit          (sex_t SEX, bitstring* seq, bitstring** parent);
  
  ///@name Implementations
  ///@{  
  virtual   void            init     () {setParameters();}
  virtual   TTrait*         hatch    ();
  virtual   TraitPrototype* clone    () {return new TProtoBDMI(*this);}
  virtual   trait_t         get_type () const {return _type;}
  
  //implements StorageComponent
  virtual   void            store_data    (BinaryStorageBuffer* saver) {}
  virtual   bool            retrieve_data (BinaryStorageBuffer* reader) {return true;}
  //implements SimComponent
  virtual   bool            setParameters ( );
  virtual   void            loadFileServices ( FileServices* loader );
  virtual   void            loadStatServices ( StatServices* loader );
  ///@}
  
};
// ------------------------------------------------------------------------------

//  TT_BDMI

// ------------------------------------------------------------------------------
class TT_BDMI : public TTrait {

  TProtoBDMI* _myProto;
  
  bitstring* _sequence[2];
  
  double _phenotype;
  
  bool _isHaploid;
  
  unsigned int _nb_locus;
  double _mut_rate;
  double _genomic_mut_rate;
  double _recomb_rate;
  
  void   (TT_BDMI::* _inherit_func_ptr)   (TTrait*, TTrait*);
  void   (TT_BDMI::* _mutation_func_ptr)  (void);
  double (TT_BDMI::* _viability_func_ptr) (void);
  
  static unsigned int *_recomb_template, *_rSites;
  static unsigned char *_sites;
  static unsigned int _haploGenotCoding[2][2];
  static unsigned int _diploGenotCoding[2][2];

public:
  //cstor & dstor
  TT_BDMI()
  : _myProto(0), _phenotype(0), _isHaploid(0), _nb_locus(0), _mut_rate(0), _genomic_mut_rate(0),
  _recomb_rate(0), _inherit_func_ptr(0), _mutation_func_ptr(0), _viability_func_ptr(0)
  {_sequence[0] = _sequence[1] = NULL;}
  
  TT_BDMI(const TT_BDMI& T)
  : _myProto(T._myProto), _phenotype(0), _isHaploid(T._isHaploid), _nb_locus(T._nb_locus),
  _mut_rate(T._mut_rate), _genomic_mut_rate(T._genomic_mut_rate), _recomb_rate(T._recomb_rate), 
  _inherit_func_ptr(T._inherit_func_ptr), _mutation_func_ptr(T._mutation_func_ptr),
  _viability_func_ptr(T._viability_func_ptr)
  {_sequence[0] = _sequence[1] = NULL;}
  
  virtual ~TT_BDMI() { }
  
  ///@name Setters:
  ///@{
  void            set_nb_locus                  (int val)               {_nb_locus = val;}
  void            set_mut_rate                  (double val)            {_mut_rate = val;}
  void            set_geno_rate                 (double val)            {_genomic_mut_rate = val;}
  void            set_recomb_rate               (double val)            {_recomb_rate = val;}
  void            set_isHaploid                 (bool val)              {_isHaploid = val;}
  void            set_inherit_func_ptr          (void(TT_BDMI::* theFunc)(TTrait*, TTrait*)) 
                                                                        {_inherit_func_ptr = theFunc;}
  void            set_mutation_func_ptr         (void (TT_BDMI::* theFunc) (void))
                                                                        {_mutation_func_ptr = theFunc;}
  void            set_viability_func_ptr        (double (TT_BDMI::* theFunc) (void))
                                                                        {_viability_func_ptr = theFunc;}
  
  void            set_proto                     (TProtoBDMI* proto)    {_myProto = proto;}
  ///@}
  
  void            set_sequence                  (bitstring** seq);
  
  //inheritance routines:
  void            inherit_haplo                 (TTrait* mother, TTrait* father);
  void            inherit_diplo                 (TTrait* mother, TTrait* father);
  
  //mutation routines:
  void            mutate_haplo                  ( );
  void            mutate_diplo                  ( );
  
  double          viability_haplo               ( );
  double          viability_diplo               ( );
  
  unsigned int    get_num_mut_haplo             (unsigned int loc) {return (*_sequence[0])[loc];}
  unsigned int    get_num_mut_diplo             (unsigned int loc) {return (*_sequence[0])[loc]+(*_sequence[1])[loc];}
  
  
  ///@name Implementations
  ///@{
  virtual   void            init ();
  virtual   void            init_sequence ();
  virtual   void            reset ();
  virtual   void            inherit (TTrait* mother, TTrait* father){(this->* _inherit_func_ptr) (mother, father);}
  virtual   void            mutate       () {(this->*_mutation_func_ptr)();}
  virtual   void*           set_trait    (void* value) {return NULL;}
  virtual   void            set_sequence (void** seq) {}
  virtual   void            set_value    ()      {_phenotype = (this->*_viability_func_ptr)();}
  virtual   void*           getValue     () const {return (void*)&_phenotype;}
  virtual   trait_t         get_type     () const {return _myProto->get_type();}
  virtual   void**          get_sequence () const {return (void**)&_sequence[0];}
  virtual   double          get_allele_value   (int loc, int all);
  virtual   void            set_allele_value (unsigned int locus, unsigned int allele, double value);
  virtual   void            show_up  ();

  virtual   TT_BDMI*         clone () {return new TT_BDMI(*this);}
  virtual   TT_BDMI&         operator= (const TTrait&);
  virtual   bool             operator== (const TTrait&);
  virtual   bool             operator!= (const TTrait&);
  //implements StorableComponent:
  virtual void    store_data     (BinaryStorageBuffer* saver) {}
  virtual bool    retrieve_data  (BinaryStorageBuffer* reader) {return true;}
  ///@}
  
};

// ------------------------------------------------------------------------------

//  TTBDMI_FH

// ------------------------------------------------------------------------------
/** FileHandler for the DBMI trait.
    Records the complete genotype in a text file. */
class TTBDMI_FH : public TraitFileHandler< TProtoBDMI > {
  
public:
  
  TTBDMI_FH(TProtoBDMI* TP) : TraitFileHandler< TProtoBDMI > (TP,".dmi") {}
  virtual ~TTBDMI_FH ( ) { }
  
  void write_haplo (Patch* patch, sex_t SEX, age_idx AGE, ofstream& FH);
  void write_diplo (Patch* patch, sex_t SEX, age_idx AGE, ofstream& FH);
  
  virtual void FHwrite  ();
  virtual void FHread (string& filename) {}
  
};

// ------------------------------------------------------------------------------

//  TTBDMI_SH

// ------------------------------------------------------------------------------
/** StatHandler for the DBMI trait.
    Records the average allele frequencies at all loci, and the frequencies of
    incompatibilities (i.e. double heterozygotes).*/
class TTBDMI_SH : public TraitStatHandler< TProtoBDMI, TTBDMI_SH > {

  double _freq, _freqIcomp;
  double *_patchFreq, *_patchIcmp;
  
public:
  TTBDMI_SH(TProtoBDMI* TP) : TraitStatHandler< TProtoBDMI, TTBDMI_SH > (TP),
                               _patchFreq(0), _patchIcmp(0) {}
  
  virtual ~TTBDMI_SH(){if(_patchFreq)delete[]_patchFreq;if(_patchIcmp)delete[]_patchIcmp;}
  
  virtual bool setStatRecorders (std::string& token);
  
  void addStats          (age_t AGE);
  void setAdultStats     () ;
  void setOffsprgStats   ();
  void setStats          (age_idx agex, void(TTBDMI_SH::* cntFunc)(Patch*,sex_t,age_idx,double**,double**));
  void countAllele_haplo (Patch *patch, sex_t SEX, age_idx AGE, double **frqTab, double **icpTab);
  void countAllele_diplo (Patch *patch, sex_t SEX, age_idx AGE, double **frqTab, double **icpTab);
  double getFreq         () {return _freq;}
  double getFreqIcmp     () {return _freqIcomp;}
  double getPatchFreq    (unsigned int i) {return _patchFreq[i];}
  double getPatchIcmp    (unsigned int i) {return _patchIcmp[i];}
  
};

// ------------------------------------------------------------------------------

//  LCE_Init_BDMI

// ------------------------------------------------------------------------------
/** Allelic frequency initialiser for the DBMI trait.
    It executes at the first generation of each replicate only.*/
class LCE_Init_BDMI : public LifeCycleEvent {
  
  TMatrix _init_freq;
  unsigned int _nLocus;
  
public:
  
  LCE_Init_BDMI ( );
  
  virtual ~LCE_Init_BDMI ( ) { }
  
  bool setSpatialPattern(TMatrix& freq_mat, unsigned int patchNbr);
  bool setPatchFreq(TMatrix& freq_mat,TMatrix& pat_mat, unsigned int patchNbr);
  void init_value(sex_t SEX, age_idx age, unsigned int size, unsigned int deme);
  
  //LifeCycleEvent implementation:
  virtual void execute ();
  
  virtual LifeCycleEvent* clone ( ) {return new LCE_Init_BDMI();}
  
  virtual bool setParameters ();
  
  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
};

#endif
