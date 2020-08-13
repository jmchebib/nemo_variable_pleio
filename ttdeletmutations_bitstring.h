/**  $Id: ttdeletmutations_bitstring.h,v 1.16 2016-03-03 15:20:33 fred Exp $
*
*  @file ttdeletmutations_bitstring.h
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
*  Created on @date 22.03.2004
*  @author fred
*/

#ifndef TTDELETMUTATIONS_BITSTR_H
#define TTDELETMUTATIONS_BITSTR_H

#include "ttrait.h"
#include "types.h"
#include "stathandler.h"
#include "filehandler.h"
#include "datatable.h"
#include "metapop.h"
#include "bitstring.h"
#include "Uniform.h"

class TTDeletMutBitstrSH;
class TTDeletMutBitstrFH;
class TProtoDeletMutations_bitstring;


/**Bitstring implementation of TTDeletMutations with recombination.*/ 
class TTDeletMutations_bitstring : public TTrait {
private:
  
  TProtoDeletMutations_bitstring* _myProto;
  
  ///@name Parameters
  ///@{
  unsigned int _nb_locus;
  unsigned int _fitness_model;
  double _fitness_scaling_factor;
  double _init_freq;
  double _mut_rate;
  double _genomic_mut_rate;
  double _strength;
  double _dominance;
  bool   _continuous_effects;
  double (TTDeletMutations_bitstring::* _viability_func_ptr) (void);
  void   (TProtoDeletMutations_bitstring::* _inherit_func_ptr)(sex_t, bitstring*, bitstring**);
  void   (TTDeletMutations_bitstring::* _mutation_func_ptr)  (void);
  ///@}
  //globs:
  static  float** _effects;
  
  //counters:
  unsigned int _nb_mutations;
  unsigned int _nb_hmz_mutations;
  unsigned int _nb_htz_mutations;
  
  bitstring* sequence[2];
  bitstring* pleio_sequence[2];
  bitstring* mutcor_sequence[2];
  bitstring *_htz, *_hmz;
  double _phenotype;
  
  trait_t _type;
  
  void            set_nb_mutations              ( );
  void            set_nb_htz_mutations          ( );
  void            set_nb_hmz_mutations          ( );
  
  
public:
  
    //'tors:
    TTDeletMutations_bitstring () 
  : _myProto(0), _nb_locus(0), _fitness_model(0), _fitness_scaling_factor(1), _init_freq(0), 
  _mut_rate(0), _genomic_mut_rate(0), _strength(0), _dominance(0), _continuous_effects(0),
     _viability_func_ptr(0), _inherit_func_ptr(0),_mutation_func_ptr(0), _nb_mutations(0), 
     _nb_hmz_mutations(0),_nb_htz_mutations(0), _htz(0), _hmz(0), _phenotype(0), _type(DELE)
  {sequence[0] = sequence[1] = pleio_sequence[0] = pleio_sequence[1] = mutcor_sequence[0] = mutcor_sequence[1] = NULL;}
  
  TTDeletMutations_bitstring(const TTDeletMutations_bitstring& T) 
  : _myProto(T._myProto), _nb_locus(T._nb_locus), _fitness_model(T._fitness_model), 
  _fitness_scaling_factor(T._fitness_scaling_factor), 
  _init_freq(T._init_freq), _mut_rate(T._mut_rate), _genomic_mut_rate(T._genomic_mut_rate), 
  _strength(T._strength), _dominance(T._dominance), _continuous_effects(T._continuous_effects),
  _viability_func_ptr(T._viability_func_ptr),_inherit_func_ptr(T._inherit_func_ptr), 
  _mutation_func_ptr(T._mutation_func_ptr), _nb_mutations(0), _nb_hmz_mutations(0), 
  _nb_htz_mutations(0), _htz(0), _hmz(0), _phenotype(0), _type(DELE)
  {sequence[0] = sequence[1] = pleio_sequence[0] = pleio_sequence[1] = mutcor_sequence[0] = mutcor_sequence[1] = NULL;}
  
  virtual ~TTDeletMutations_bitstring() {reset();}
  
  
  ///@name Getters:
  ///@{
  unsigned int    get_nb_mutations              ( )                     {return _nb_mutations;}
  unsigned int    get_nb_mut_atLocus            (unsigned int loc)      {return (*sequence[0])[loc] + (*sequence[1])[loc];}
  unsigned int    get_nb_htz_mutations          ( )                     {return _nb_htz_mutations;}
  unsigned int    get_nb_hmz_mutations          ( )                     {return _nb_hmz_mutations;}
  bool            get_hmz_atLocus               (unsigned int loc)      {return (*sequence[0])[loc] & (*sequence[1])[loc];}
  bool            get_htz_atLocus               (unsigned int loc)      {return (*sequence[0])[loc] ^ (*sequence[1])[loc];}
  float**         get_effects                   ( )  const              {return _effects;}
  ///@}
  ///@name Setters:
  ///@{
  void            set_proto     (TProtoDeletMutations_bitstring* proto) {_myProto = proto;}
  void            set_nb_locus                  (int val)               {_nb_locus = val;}
  void            set_mut_rate                  (double val, int nloc)  {_mut_rate = val; _genomic_mut_rate = 2*nloc*_mut_rate;}
  void            set_strength                  (double val)            {_strength = val;}
  void            set_dominance                 (double val)            {_dominance = val;}
  void            set_viability_func_ptr        (unsigned int f_model, bool is_cont);
  void            set_inherit_func_ptr          (void(TProtoDeletMutations_bitstring::* theFunc)(sex_t, bitstring*, bitstring**)) 
                                                                        {_inherit_func_ptr = theFunc;}
  void            set_mutation_func_ptr         (unsigned int m_model);
  void            set_fitness_scaling_factor    (double val)            {_fitness_scaling_factor = val;}
  void            set_init_freq                 (double val)            {_init_freq = val;}
  void            set_fitness_model             (int val)               {_fitness_model = val;}
  void            set_continuous_effects        (bool val)              {_continuous_effects = val;}
  ///@}
  
  double          viability_multi               ( );
  double          viability_epist               ( );
  double          viability_multi_continuous    ( );
  double          viability_epist_continuous    ( );
  void            mutate_redraw                 ( );
  void            mutate_noredraw               ( );

  void            set_sequence                  (bitstring** seq);
  void            set_pleio_sequence            (bitstring** seq);
  void            set_mutcor_sequence           (bitstring** seq);
  //glob setters:
  static void     set_effects                   (float** fx);
  static void     set_recomb_template           (unsigned int size);
  
  ///@name Implementations
  ///@{
  virtual void    init                          ( );
  virtual void    init_sequence                 ( );
  virtual void	  reset                         ( );
  virtual void*   set_trait                     (void* value)           {return NULL;}
  virtual void**  get_sequence                  ( )        const        {return (void**)&sequence[0];}
  virtual void**  get_pleio_sequence            ( )        const        {return (void**)&pleio_sequence[0];}
  virtual void**  get_mutcor_sequence           ( )        const        {return (void**)&mutcor_sequence[0];}
  virtual double  get_allele_value              (int loc, int all);
  /*Be aware that the set_allele_value here changes the mutation effect for all individuals in the pop!!*/
  virtual void    set_allele_value              (unsigned int locus, unsigned int allele, double value);
  virtual void    set_sequence                  (void** seq)            {}
  virtual void    set_pleio_sequence            (void** seq)            {}
  virtual void    set_mutcor_sequence           (void** seq)            {}
  virtual trait_t get_type                      ( )        const        {return _type;}
  virtual void    inherit                       (TTrait* mother, TTrait* father);
  virtual void    mutate                        ( )                     {(this->*_mutation_func_ptr)();}
  virtual void    set_value                     ( );
  virtual void*   getValue                      ( )        const        {return (void*)&_phenotype;}
  virtual void    show_up                       ( );
  virtual TTDeletMutations_bitstring*   clone             ( )                     {return new TTDeletMutations_bitstring(*this);}
  virtual TTDeletMutations_bitstring& operator=(const TTrait& T);
  virtual bool operator==(const TTrait& T);
  virtual bool operator!=(const TTrait& T);
  
  //implements StorableComponent:
  virtual void    store_data     (BinaryStorageBuffer* saver); //   {saver->store(sequence[0], _nb_locus); saver->store(sequence[1], _nb_locus);}
  virtual bool    retrieve_data  (BinaryStorageBuffer* reader); //  {reader->read(sequence[0], _nb_locus); reader->read(sequence[1], _nb_locus);return true;}
  ///@}
};

// ------------------------------------------------------------------------------

//  TProtoDeletMutations

// ------------------------------------------------------------------------------
/**Prototype class of the bitstring-deleterious mutations trait class.*/
class TProtoDeletMutations_bitstring : public TTProtoWithMap {
public:  
  TProtoDeletMutations_bitstring ();
  TProtoDeletMutations_bitstring (const TProtoDeletMutations_bitstring& T);
  ~TProtoDeletMutations_bitstring ();
  
  int             get_nb_locus     ( )  {return _nb_locus;}
  double          get_mut_rate     ( )  {return _mut_rate;}
  double          get_strength     ( )  {return _strength;}
  double          get_dominance    ( )  {return _dominance;}
  int             get_dominance_model( )  {return _dominance_model;}
  bool            get_iscontinuous ( )  {return _continuous_effects;}
  
  void            set_effects         ( );
  double          set_effects_exp     ( )  {return RAND::Exponential(_strength);}
  double          set_effects_gamma   ( )  {return RAND::Gamma(_dist_p1, _dist_p2);}
  double          set_effects_lognorm ( )  {return RAND::LogNormal(_dist_p1, _dist_p2);}
  float*          get_s_continous     ( )  {return _effects[1];}
  float*          get_hs_continous    ( )  {return _effects[0];}
  
  void    inherit_low                  (sex_t SEX, bitstring* seq, bitstring** parent);
  void    inherit_free                 (sex_t SEX, bitstring* seq, bitstring** parent);

  bool   setSelectionParameters ();
  
  //implements TraitPrototype
  virtual   void                    init  (){setParameters();};
  
  virtual   void                    reset (){TTProtoWithMap::reset();}

  virtual   TTDeletMutations_bitstring*       hatch ();
  
  virtual   TProtoDeletMutations_bitstring*   clone ()      {return new TProtoDeletMutations_bitstring((*this));}
  
  virtual   trait_t                 get_type () const {return DELE;}
  //implements StorageComponent
  virtual   void                    store_data    (BinaryStorageBuffer* saver);
  
  virtual   bool                    retrieve_data (BinaryStorageBuffer* reader);
  //implements SimComponent
  virtual   bool                    setParameters ( );
  
  virtual   void                    loadFileServices ( FileServices* loader );
  
  virtual   void                    loadStatServices ( StatServices* loader );
  
private:
    
  unsigned int _nb_locus;
  unsigned int _fitness_model;
  unsigned int _mutation_model;
  int _dominance_model;
  double _fitness_scaling_factor;
  double _init_freq;
  double _mut_rate;
  double _strength;
  double _dominance;
  double _dist_p1; //parameter 1 of random effects distribution (mu or shape)
  double _dist_p2; //parameter 2 of random effects distribution (sigma or scale)
  bool   _continuous_effects;
  double (TTDeletMutations_bitstring::* _viability_func_ptr) (void);
  double (TProtoDeletMutations_bitstring::* _set_effects_func) (void);
  void (TProtoDeletMutations_bitstring::* _inherit_func_ptr) (sex_t, bitstring*, bitstring**);
  
  TTDeletMutBitstrSH* _stats;
  TTDeletMutBitstrFH* _writer;
  TTDeletMutBitstrFH* _reader;
  float** _effects;  
};

/**The StatHandler for TTDeletMutations_bitstring*/ 
class TTDeletMutBitstrSH : public TraitStatHandler<TProtoDeletMutations_bitstring, TTDeletMutBitstrSH> {
    
  double fecWithHomePatchMate, fecWithOtherPatchMate;
  //props and viability:
  double _SibProps[5], _viability[5], _meanViab;
  //0: outbred from migrants, 1: outbred from residants
  //2: inbred from half-sibs, 3: inbred from full-sibs
  //4: selfed offsprg
  
  double _deletHtzLoci, _deletHmzLoci, _fixLocPerPatch, _segrLocPerPatch;
  double _Ho, _Hs, _Ht, _Fst, _Hmz, _deletAllCount, _freq, _fixloc, _segrloc, _letheq;
  //table to store the deleterious mutations freq per locus:
  double *_deletFreqTable;
  
  bool _isContinuousEffect;
  
public:
	
  TTDeletMutBitstrSH ( TProtoDeletMutations_bitstring* TP ) 
  : TraitStatHandler<TProtoDeletMutations_bitstring, TTDeletMutBitstrSH> ( TP ), _deletFreqTable(0)
  { _isContinuousEffect = TP->get_iscontinuous(); }
  
  virtual ~TTDeletMutBitstrSH () {if(_deletFreqTable != 0) delete [] _deletFreqTable;}
  
  virtual bool setStatRecorders (std::string& token);
  
  void   setStatsForDeletMutations (age_t AGE);
  void   setViabStats           (age_t AGE);
  void   setDeletStats          (age_t AGE);
  void   setLethalEquivalents   (age_t AGE);
  void   setFst                 (age_t AGE);
  void   setAdultDeletStats              ()       {setDeletStats(ADULTS);}
  void   setOffsprgDeletStats            ()       {setDeletStats(OFFSPRG);}
  void   setViability        (age_idx agex);
  void   setAdultViab                    ()       {setViability(ADLTx);}
  void   setOffsprgViab                  ()       {setViability(OFFSx);}
  void   setMeanViability    (age_idx agex);
    
  
  double getMeanFecWithPatchMate(bool HOME);
  double getBtheta                       ()       {return fecWithHomePatchMate;}
  double getBalpha                       ()       {return fecWithOtherPatchMate;}
  double getHeterosis                    ();
  double getLoad                         ();
  double getPatchLoad      (unsigned int i);
  double getDeletAllFreq                 ()       {return _freq;}
  double getDeletAllHmz                  ()       {return _Hmz;}
  double getDeletAllHtz                  ()       {return _Ho;}
  double getFixedDeletLoci               ()       {return _fixloc;}
  double getFixedDeletLociPerPatch       ()       {return _fixLocPerPatch;}
  double getSegregatingDeletLoci         ()       {return _segrloc;}
  double getSegregatingDeletLociPerPatch ()       {return _segrLocPerPatch;}
  double getDeletAllPerGenome            ()       {return _deletAllCount;}
  double getHs                           ()       {return _Hs;}
  double getHt                           ()       {return _Ht;}
  double getFst                          ()       {return _Fst;}
  double getLethalEquivalents            ()       {return _letheq;}
  double getAdultsLetheq                 ()       {setDeletStats(ADULTS);  return _letheq;}
  double getOffsprgLetheq                ()       {setDeletStats(OFFSPRG); return _letheq;}
  
  double getViability      (unsigned int v)       {return _viability[v];}
  double getMeanViability                ()       {return _meanViab;}
  double getMeanViability  (unsigned int a)       {setMeanViability(static_cast<age_idx> (a)); return _meanViab;}
  
  double getSibProportions (unsigned int i)       {return _SibProps[i];}
  
};


/**The FileHandler associated with the TTDeletMutations_bitstring trait.
* Used to save genotypes in a text file.*/ 
class TTDeletMutBitstrFH : public TraitFileHandler<TProtoDeletMutations_bitstring> {
  
public:
    
  TTDeletMutBitstrFH (TProtoDeletMutations_bitstring* TP) 
  : TraitFileHandler<TProtoDeletMutations_bitstring> (TP,".del") { }
  
  virtual ~TTDeletMutBitstrFH ( ) { }
   
  virtual void FHwrite  ();
  virtual void FHread (string& filename);
};


#endif //TTDELETMUTATIONS_BITSTR_H
