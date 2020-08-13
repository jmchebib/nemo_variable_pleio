/**  $Id: ttneutralgenes.h,v 1.20 2017-06-20 08:35:36 fred Exp $
*
*  @file ttneutralgenes.h
*  Nemo2
*
*   Copyright (C) 2006-2012 Frederic Guillaume
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

#ifndef TTNEUTRALGENES_H
#define TTNEUTRALGENES_H

#include <cmath>
#include <vector>
#include "ttrait_with_map.h"
#include "types.h"
#include "filehandler.h"
#include "stathandler.h"
#include "datatable.h"
#include "metapop.h"
#include "binarystoragebuffer.h"

class TTNeutralGenesFH;
class TTNeutralGenesSH;
class TTNtrlPhenotyperFH;
class TProtoNeutralGenes;
// ----------------------------------------------------------------------------------------

//                                 N E U T R A L   T R A I T

// ----------------------------------------------------------------------------------------
/**Microsatellites genome.*/ 
class TTNeutralGenes : public TTrait 
{
private:
  TProtoNeutralGenes* _myProto;
  //parameters:
  unsigned int _allele_num;
  unsigned int _locus_num;
  unsigned int _ploidy;
  double _mut_rate;
  double _recomb_rate;
  unsigned int _mut_model;
  unsigned int _2L; //genome size (ploidy*num_locus)
  unsigned short _init_model;
  void (TTNeutralGenes::* _mutate_func_ptr) (void);
  void (TProtoNeutralGenes::* _inherit_func_ptr) (sex_t, unsigned char*, unsigned char**);

  //_sequence:
  unsigned char** _sequence; // ordered as [ploidy][locus]

  const trait_t _type;
   
public:
  
  TTNeutralGenes () 
    : _myProto(0),_allele_num(0),_locus_num(0),_ploidy(2),_mut_rate(0),_recomb_rate(0),_mut_model(0),
    _init_model(0), _mutate_func_ptr(0), _inherit_func_ptr(0), _sequence(0), _type(NTRL) { }
  
  TTNeutralGenes(const TTNeutralGenes& T) 
    : _myProto(T._myProto),
    _allele_num(T._allele_num), _locus_num(T._locus_num), _ploidy(2), _mut_rate(T._mut_rate),
    _recomb_rate(T._recomb_rate),_mut_model(T._mut_model), _2L(T._2L),
    _init_model(T._init_model),_mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr), 
    _sequence(0), _type(NTRL) { }
  
  virtual ~TTNeutralGenes ();
    
  ///@name Accessors
  ///@{
  unsigned int    get_ploidy           ( )					   {return _ploidy;}
  unsigned int    get_locus_num        ( )                     {return _locus_num;}
  unsigned int    get_allele_num       ( )                     {return _allele_num;}
  
  void            set_proto          (TProtoNeutralGenes* proto) {_myProto = proto;}
  void            set_locus_num        (int value)             {_locus_num = value;}
  void            set_allele_num       (int value)             {_allele_num = value;}
  void            set_mut_rate         (double value)          {_mut_rate = value;}
  void            set_2L               (unsigned int val)      {_2L = val;}
  void            set_recomb_rate      (double value)          {_recomb_rate = value;}
  void            set_mut_model        (int value)             {_mut_model = value;}
  void            set_init_model       (unsigned short val)    {_init_model = val;}
  void            set_mut_func_ptr     (void(TTNeutralGenes::* theFunc)(void))
  {_mutate_func_ptr = theFunc;}
  void            set_inherit_func_ptr (void(TProtoNeutralGenes::* theFunc)(sex_t, unsigned char*, unsigned char**))
  {_inherit_func_ptr = theFunc;}
  void            set_allele           (unsigned int loc, unsigned int al, unsigned char val) {_sequence[al][loc] = val;}
  ///@}
  
 // static void     set_recomb_template           (unsigned int size);
  
  ///@name Mutation models
  ///@{
  void            mutate_SSM           ( );
  void            mutate_KAM           ( );
  void            mutate_2all          ( );
  void            mutate_NULL          ( ) { }
  ///@}
  ///@name Implementations
  ///@{
  virtual TTNeutralGenes& operator= (const TTrait& T);
  virtual bool    operator== (const TTrait& T);
  virtual bool    operator!= (const TTrait& T);
  virtual void    init                 ( );
  virtual void    init_sequence        ( );
  virtual void    reset                ( );
  virtual void*   set_trait            (void* value)           {return NULL;}
  inline virtual void**  get_sequence  ( )  const              {return (void**)_sequence;}
  virtual double  get_allele_value     (int loc, int all)  
    {return ( !(loc<(int)_locus_num) || !(all<(int)_ploidy) ? 0 : (double)_sequence[all][loc]);}
  virtual void    set_allele_value (unsigned int locus, unsigned int allele, double value)
    {assert(locus < _locus_num && allele < 2); _sequence[allele][locus] = (unsigned char)value;}
  virtual void    set_sequence         (void** seq);
  inline virtual trait_t get_type      ( )  const              {return _type;}
  virtual void    set_value            ( )                     { }
  virtual void*   getValue             ( )	const              {return 0;}
  virtual void    inherit              (TTrait* mother, TTrait* father);
  virtual void    mutate               ( )                      {(this->*_mutate_func_ptr) ();}
  virtual void    show_up              ( );
  virtual TTNeutralGenes*  clone       ( )                      {return new TTNeutralGenes(*this);}
  
  //implements StorableComponent
  virtual void    store_data       (BinaryStorageBuffer* saver)  
  {
    for(unsigned int i = 0; i < _locus_num; ++i)
      for(unsigned int j = 0; j < _ploidy; ++j)
        saver->store(&_sequence[j][i], 1);
  }
  virtual bool    retrieve_data    (BinaryStorageBuffer* reader)
  { 
    for(unsigned int i = 0; i < _locus_num; ++i)
      for(unsigned int j = 0; j < _ploidy; ++j)
        reader->read(&_sequence[j][i], 1);
    return true;
  }
  ///@}
};  
// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   P R O T O T Y P E

// ----------------------------------------------------------------------------------------
/**Prototype class for the TTNeutralGenes trait class.**/
class TProtoNeutralGenes : public TTProtoWithMap {
private:
  unsigned int _allele_num;
  unsigned int _locus_num;
  unsigned int _ploidy;
  double _mut_rate;
  int _mut_model;
  unsigned short _init_model;
  double _recomb_rate;
  void (TTNeutralGenes::* _mutate_func_ptr) (void);
  void (TProtoNeutralGenes::* _inherit_func_ptr) (sex_t, unsigned char*, unsigned char**);
  
  vector< TTNeutralGenesFH* > _writers;
  TTNeutralGenesSH* _stats;
  const trait_t _type;
  
public:
  
  TProtoNeutralGenes ( );
  
  TProtoNeutralGenes(const TProtoNeutralGenes& T);
  
  virtual ~TProtoNeutralGenes ( );
  
  unsigned int     get_ploidy           ( )					   {return _ploidy;}
  unsigned int     get_locus_num        ( )            {return _locus_num;}
  unsigned int     get_allele_num       ( )            {return _allele_num;}
  
  TTNeutralGenesSH* get_stater () {return _stats;}

  void    inherit_low                  (sex_t SEX, unsigned char* seq, unsigned char** parent);
  void    inherit_free                 (sex_t SEX, unsigned char* seq, unsigned char** parent);

  //implementation of TraitPrototype:
  virtual void                     init  (){}
  
  virtual void                     reset (){TTProtoWithMap::reset();}
  
  virtual TTNeutralGenes*          hatch ();
  
  virtual TProtoNeutralGenes*      clone () {return new TProtoNeutralGenes(*this);}
  
  virtual   trait_t                 get_type () const {return _type;}
  //implementation of SimComponent:
  virtual bool setParameters();
  
  virtual void loadFileServices ( FileServices* loader );
  
  virtual void loadStatServices ( StatServices* loader );
  
  //implementation of StorableComponent:
  virtual void store_data    ( BinaryStorageBuffer* saver ) {saver->store(&_locus_num,sizeof(int));}
  
  virtual bool retrieve_data ( BinaryStorageBuffer* reader );
};

// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   F I L E   H A N D L E R

// ----------------------------------------------------------------------------------------

/**A file handler to save the neutral markers genotypes in the FSTAT format (extended). 
   By default, the file extension is ".dat" for the genotype file. It is changed to ".fsti"
   for the per-locus/per-patch Weir& Hill Fst's and ".freq" for the per-locus/per-patch allele
   frequencies. Also implements the FileHandler::FHread method to load a population's genotypes
   from an FSTAT file (using the 'source_pop' population parameter).
*/
class TTNeutralGenesFH: public TraitFileHandler<TProtoNeutralGenes> {
  
  string _output_option;

  void (TTNeutralGenesFH::* write_fct) ();
  
public:
	
  TTNeutralGenesFH ( TProtoNeutralGenes* TP )
    : TraitFileHandler<TProtoNeutralGenes> (TP, ".txt")
  { }
  
  virtual ~TTNeutralGenesFH ( ) { }

  virtual void FHwrite  (); // { (this->*write_fct)(); }
  
  virtual void FHread (string& filename);

  void write_TAB();
  void write_PLINK ();
  void write_PLINK_BED (ofstream &BED);
  void write_GENEPOP();
  void write_FSTAT();
  void write_Fst_i();
//  void write_freq();
  void write_varcompWC();
  
  void setOutputOption(string opt);
  void set_write_fct( void (TTNeutralGenesFH::* fct_ptr) () ) {write_fct = fct_ptr;}
  
};

// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   S T A T   H A N D L E R

// ----------------------------------------------------------------------------------------

/**The stat handler for neutral markers. */
class TTNeutralGenesSH: public TraitStatHandler<TProtoNeutralGenes, TTNeutralGenesSH> {
  
  DataTable< unsigned int > _alleleCountTable;
  DataTable< double > _alleleFreqTable;
  DataTable< double > _heteroTable;
  TMatrix _globalAlleleFreq;
  unsigned int _table_set_gen, _table_set_age, _table_set_repl;
  
  double Theta_FF, Theta_MM, Theta_FM;
  double _mean_theta, _mean_alpha;
  TMatrix *_coa_matrix;
  
  /**Kinship classes proportions*/
  double _sib_prop[4];
  double _sib_coa[4];
  
  /**F-statistics*/
  double _ho, _hs, _ht, _hsnei, _htnei, _nb_all_local, _nb_all_global,
    _fst, _fis, _fit, _fix_loc_local, _fix_loc_global;
  /**Weir & Hill (2002) F-stat estimates.*/
  double _fst_WH;
  /**Weir & Cockerham (1984) F-stat estimates.*/
  double _fst_WC, _fis_WC, _fit_WC;
  /**Per-locus F-stats (Weir&Cockerham).*/
  double *_fst_WC_loc, *_fis_WC_loc, *_fit_WC_loc;

  double _fst_W1, _fst_W2;
  /**Pairwise Fst matrix.*/
  TMatrix *_fst_matrix;
  
  //Nei's genetic distance:
  TMatrix *_D;
  double _meanD;
  
public:
  
	TTNeutralGenesSH (TProtoNeutralGenes* TP) 
    : TraitStatHandler<TProtoNeutralGenes, TTNeutralGenesSH>(TP), _table_set_gen(0), _table_set_age(0),
    _table_set_repl(0), _coa_matrix(0), _fst_WC_loc(0), _fis_WC_loc(0), _fit_WC_loc(0), _fst_matrix(0),
	_D(0)
	{ }
  
  virtual ~TTNeutralGenesSH ( ) 
  {
    if(_coa_matrix != NULL) delete _coa_matrix;
    if(_fst_matrix != NULL) delete _fst_matrix;
    if(_fst_WC_loc) delete[]_fst_WC_loc;
    if(_fis_WC_loc) delete[]_fis_WC_loc;
    if(_fit_WC_loc) delete[]_fit_WC_loc;
    if(_D != NULL) delete _D;
  }
  
  virtual void init ( ) ;
  
  virtual bool setStatRecorders (std::string& token);
  
  void    setFreqRecorders (age_t AGE);
  void    setFreqRecordersPerPatch (age_t AGE);
  void    setFstatRecorders (age_t AGE);
  void    setFstat2Recorders (age_t AGE);
  void    setFstatWCRecorders (age_t AGE);
  void    setCoaMatrixRecorders (age_t AGE, unsigned char dim);
  void    setFstMatrixRecorders (age_t AGE, unsigned char dim);
  void    setNeiGeneticDistanceRecorders (age_t AGE, bool pairwise);
  void    setDxyRecorders (age_t AGE, bool patchwise);
  
  ///@name Allele and genotype frequencies:
  ///@{
  void setAdultAlleleFreq     () {setAlleleTables(ADULTS);}
  void setOffspringAlleleFreq () {setAlleleTables(OFFSPRG);}\
  void setHeterozygosity      (age_t AGE);
  void setAdultHeterozygosity     () {setHeterozygosity(ADULTS);}
  void setOffspringHeterozygosity () {setHeterozygosity(OFFSPRG);}
  
  double getGlobalAlleleFreq (unsigned int loc, unsigned int all) {
    return _globalAlleleFreq.get(loc, all);
  }
  
  double getHeterozygosity (unsigned int loc) {
    double het = 0;
    for(unsigned int i = 0; i < _heteroTable.getNumGroups(); ++i )
      het += _heteroTable.get(i, loc, 0);
    return het/_heteroTable.getNumGroups(); //mean per patch heterozygosity
  }
  
  ///@}
  ///@name F-stats:
  ///@{
  void   setAlleleTables          (age_t AGE);
  void   setHeteroTable           (age_t AGE);
  void   allocateTables           (unsigned int loci, unsigned int all);
  
  /**Accessor to the table of allele frequencies, per patch.*/
  DataTable<double>* getAlleleFreqTable () {return &_alleleFreqTable;}
  
  DataTable< unsigned int >* getAlleleCountTable () {return &_alleleCountTable;}
  
  DataTable<double>*   getHeteroTable      () {return &_heteroTable;}
  
  /**Accessor to the table of allele frequencies in the whole population.*/
  TMatrix* getGlobalFreqs () {return &_globalAlleleFreq;}
  
  /**Computes the weighted within and between patch Fst's as well as the overall Fst (Theta).
    The method used here is that of Weir & Hill 2002, Ann. Rev. Genet. 36:721-750.
    The weighting is done for samples (patches) of unequal sizes.
    @param AGE the age class
    @param dim the dimension of the matrix to fill:
    - 1 = the diagonal (i.e. the wihtin patch Fst or theta_ii)
    - 2 = the upper half (i.e. the between patch Fst or theta_ii')
    - 3 = both
  */
  void   setFstMatrix             (age_t AGE, unsigned char dim);
  void   setAdultsFstMatrix       ()               {setFstMatrix(ADULTS, 3);}
  void   setAdultsFstWithin       ()               {setFstMatrix(ADULTS, 1);}
  void   setAdultsFstBetween      ()               {setFstMatrix(ADULTS, 2);}
  void   setOffsprgFstMatrix      ()               {setFstMatrix(OFFSPRG, 3);}
  void   setOffsprgFstWithin      ()               {setFstMatrix(OFFSPRG, 1);}
  void   setOffsprgFstBetween     ()               {setFstMatrix(OFFSPRG, 2);}
  /**Returns the weighted Fst using Weir & Hill (2002) method.
    This Fst is set by a previous call to setFstMatrix().*/
  double getWeightedFst           ()               {return _fst_WH;}
  /**Accessor to the Fst matrix as set by setFstMatrix().*/
  double getFst_ij                (unsigned int i)
  {
    unsigned int scale = (unsigned int)pow( 10.0, (int)log10((float)_fst_matrix->getNbCols()) + 1 );
    return _fst_matrix->get(i/scale, i%scale);
  }
  /**Computes the per-locus per-patch Fst values using Weir&Hill 2002 approach.*/
  void   setFst_li(unsigned int N, unsigned int L, double **array);
//  /**Computes raw Fst following the original definition (=var(p)/p_bar(1 - p_bar)).*/
//  double getFstWright(unsigned int i) {if(i == 1) return _fst_W1; else return _fst_W2;}
  
  /**Computes the F-statistics following Nei & Chesser (1983).*/
  void   setFstat                 (age_t AGE);
  void   setOffsprgFstat          ()               {setFstat(OFFSPRG);}
  void   setAdultsFstat           ()               {setFstat(ADULTS);}
  double setHo                    (age_idx age_pos);
  double setHs                    (age_idx age_pos);
  double setHt                    (age_idx age_pos);
  double getHsnei                 ()               {return _hsnei;}
  double getHtnei                 ()               {return _htnei;}
  double getHo                    ()               {return _ho;}
  double getHs                    ()               {return _hs;}
  double getHt                    ()               {return _ht;}
  double getFst                   ()               {return _fst;}
  double getFis                   ()               {return _fis;}
  double getFit                   ()               {return _fit;}
  /**New version of Nei & Chesser.*/
  void   setFstat2                 (age_t AGE);
  void   setOffsprgFstat2          ()               {setFstat2(OFFSPRG);}
  void   setAdultsFstat2           ()               {setFstat2(ADULTS);}
  deque<double> setHo2                    (age_idx age_pos);
  deque<double> setHs2                    (age_idx age_pos);
  deque<double> setHt2                    (age_idx age_pos);
  
  /**Computes the Weir & Cockerham (1984) Fstat values (Theta, F, and f).*/
  void   setFstatWeirCockerham    (age_t AGE);
  void   setFstatWeirCockerham_MS   (age_t AGE);
  void   setOffspringFstatWeirCockerham()          {setFstatWeirCockerham(OFFSPRG);}
  void   setAdultsFstatWeirCockerham()             {setFstatWeirCockerham(ADULTS);}
  double getFstWC                 ()               {return _fst_WC;}
  double getFisWC                 ()               {return _fis_WC;}
  double getFitWC                 ()               {return _fit_WC;}
  /**Sets the allelic diversity counters.*/
  void   setLociDivCounter        (age_t AGE);
  double getNbAllLocal            ()               {return _nb_all_local;}
  double getNbAllGlobal           ()               {return _nb_all_global;}
  double getFixLocLocal           ()               {return _fix_loc_local;}
  double getFixLocGlobal          ()               {return _fix_loc_global;}
  ///@}
  ///@name Coancestries
  ///@{
  /**Gives the coancestry (probability of identity by state) of two gene sequences.
     The probability returned is the average probability of having two identical alleles
     at a locus between the two sequences.
     @param ind1 first _sequence, treated as of type (unsigned char**)
     @param ind2 second _sequence, treated as of type (unsigned char**)
     @param nb_locus number of loci present in each _sequence
  */
  double Coancestry               (void** ind1, void** ind2, unsigned int nb_locus);
  /**Computes the within and between patches coancestry coefficients.
    @param age_pos the age class index
    @param dim the dimension of the matrix to fill:
               - 1 = the diagonal (i.e. the wihtin patch coancestries or theta's)
               - 2 = the upper half (i.e. the between patch coancestries or alpha's)
               - 3 = both
  */
  void   setCoaMatrix             (age_idx age_pos, unsigned char dim);
  void   setAdultsCoaMatrix       ()               {setCoaMatrix(ADLTx, 3);}
  void   setOffsprgCoaMatrix      ()               {setCoaMatrix(OFFSx, 3);}
  void   setAdultsCoaWithin       ()               {setCoaMatrix(ADLTx, 1);}
  void   setOffsprgCoaWithin      ()               {setCoaMatrix(OFFSx, 1);}
  void   setAdultsCoaBetween      ()               {setCoaMatrix(ADLTx, 2);}
  void   setOffsprgCoaBetween     ()               {setCoaMatrix(OFFSx, 2);}
  void   setAdults_Theta          ();
    
  /**Gets the given coancestry coefficient from the coancestry matrix.
    @param i combination of the row and column indexes (see setCoaMatrixRecorders()).
    \note the upper half and the diagonal of the matrix are filled, other positions are set to 0.
  */
  double getCoa                   (unsigned int i)
  {
    unsigned int scale = (unsigned int)pow( 10.0, (int)log10((float)_coa_matrix->getNbCols()) + 1 );
    return _coa_matrix->get(i/scale, i%scale);
  }
  double getMeanTheta             ()           {return _mean_theta;}
  double getMeanAlpha             ()           {return _mean_alpha;}
  /**Gives the mean within females coancestry coefficient.*/
  double getTheta_FF              ()		       {return Theta_FF;}
  /**Gives the mean within males coancestry coefficient.*/
  double getTheta_MM              ()		       {return Theta_MM;}
  /**Gives the mean between males and females coancestry coefficient.*/
  double getTheta_FM              ()		       {return Theta_FM;}
  void   setSibStats              ();
  void   setSibCoa                (Individual *I1, Individual *I2);
  double getSibProportions        (unsigned int i) {return _sib_prop[i];}
  double getSibCoaMeans           (unsigned int i) {return _sib_coa[i];}
  ///@}
  
  ///@name Nei's genetic distance:
  ///@{
  void    setAdltNeiGeneticDistance    ( ) {setNeiGeneticDistance(ADULTS);}
  void    setOffsprgNeiGeneticDistance ( ) {setNeiGeneticDistance(OFFSPRG);}
  void    setNeiGeneticDistance        (age_t AGE);
  double  getNeiGeneticDistance        (unsigned int i)  
  { 
    unsigned int scale = (unsigned int)pow( 10.0, (int)log10((float)_D->getNbCols()) + 1 );
    return _D->get(i/scale,i%scale);
  }
  double  getMeanNeiGeneticDistance    ( ) {return _meanD;}
  ///@}
  
  ///@name Sequence divergence: Dxy (Nei & Li 1979; Nei 1987)
  ///@{
  double getDxyOffspringPerPatch  (unsigned int patch1, unsigned patch2) {return getDxyPerPatch(OFFSx, patch1, patch2);}
  double getDxyAdultPerPatch      (unsigned int patch1, unsigned patch2) {return getDxyPerPatch(ADLTx, patch1, patch2);}
  double getDxyPerPatch           (age_idx age, unsigned int patch1, unsigned patch2);
  double getDxy                   (unsigned int age_class);
  ///@}
};


#endif //TTNEUTRALGENES_H

