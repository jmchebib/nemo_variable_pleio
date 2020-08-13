/**  $Id: ttneutralgenes.cc,v 1.30 2017-06-20 08:35:36 fred Exp $
 *
 *  @file ttneutralgenes.cc
 *  Nemo2
 *
 *   Copyright (C) 2006-2017 Frederic Guillaume
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
 *  Created on @date 08.05.2004
 *  @author fred
 */

#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "ttneutralgenes.h"
#include "Uniform.h"
#include "output.h"
#include "tstring.h"
#include <math.h>
#include "bitstring.h"
// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   P R O T O T Y P E

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoNeutralGenes::TProtoNeutralGenes( ) :_allele_num(0), _locus_num(0), _ploidy(2), _mut_rate(0),
_mut_model(0), _init_model(0), _mutate_func_ptr(0), _inherit_func_ptr(0),
_stats(0), _type(NTRL)
{
  set_paramset("neutralgenes", false, this);

  add_parameter("ntrl_loci",              INT, true, false, 0, 0);
  add_parameter("ntrl_init_model",        INT, false, true, 0, 3);
  add_parameter("ntrl_all",               INT, true, true, 2, 256);
  add_parameter("ntrl_mutation_rate",     DBL, true, true, 0, 1);
  add_parameter("ntrl_mutation_model",    INT, true, true, 0, 2);

  //genetic map parameters:
  TTProtoWithMap::addGeneticMapParameters("ntrl");

  //output
  add_parameter("ntrl_save_genotype",  STR, false, false, 0, 0);
  add_parameter("ntrl_save_fsti",      BOOL, false, false, 0, 0);
  add_parameter("ntrl_save_freq",      STR, false, false, 0, 0);
  add_parameter("ntrl_output_dir",     STR, false, false, 0, 0);
  add_parameter("ntrl_output_logtime", INT, false, false, 0, 0);
}
// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoNeutralGenes::TProtoNeutralGenes(const TProtoNeutralGenes& T) 
: _allele_num(T._allele_num), _locus_num(T._locus_num), _ploidy(2), _mut_rate(T._mut_rate),
  _mut_model(T._mut_model), _init_model(T._init_model),
  _mutate_func_ptr(T._mutate_func_ptr), _inherit_func_ptr(T._inherit_func_ptr), _stats(0), _type(NTRL)
{
  _paramSet = new ParamSet( *(T._paramSet) ) ;
}  
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
TProtoNeutralGenes::~TProtoNeutralGenes ()
{
  if(_stats != NULL)      {delete _stats; _stats = NULL;}
  for(unsigned int i = 0; i < _writers.size(); i++)
    delete _writers[i];
  _writers.clear();

}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TProtoNeutralGenes::setParameters ()
{
  _locus_num =  (unsigned int)get_parameter_value("ntrl_loci");
  _allele_num = (unsigned int)get_parameter_value("ntrl_all");
  _ploidy = 2;
  _mut_rate =  get_parameter_value("ntrl_mutation_rate");
  _mut_model = (unsigned int)get_parameter_value("ntrl_mutation_model");

  if(get_parameter("ntrl_init_model")->isSet())
    _init_model = (unsigned short)get_parameter_value("ntrl_init_model");
  else
    _init_model = 1; //maximum variance

  switch(_mut_model) {
    case 0:
      {
      _mutate_func_ptr = &TTNeutralGenes::mutate_NULL;
      break;
      }
    case 1:
      {
      _mutate_func_ptr = &TTNeutralGenes::mutate_SSM;
      break;
      }
    case 2:
      {
      _mutate_func_ptr = &TTNeutralGenes::mutate_KAM;
      break;
      }
    default:
      {
      error("wrong parameter value for parameter \"ntrl_mutation_model\", max is 2\n");
      break; //should return false
      }
  }
  //special case:
  if(_allele_num==2) _mutate_func_ptr = &TTNeutralGenes::mutate_2all;

  // recombination parameters, calling the genetic map
  _recombRate = get_parameter_value("ntrl_recombination_rate"); //member of TTProtoWithMap

  // we want to bypass the genetic map if loci are not linked
  if( _recombRate == 0.5 ) {

      _inherit_func_ptr = &TProtoNeutralGenes::inherit_free;

  } else {
      _inherit_func_ptr = &TProtoNeutralGenes::inherit_low;
      return TTProtoWithMap::setGeneticMapParameters ("ntrl");
  }

  return true;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TProtoNeutralGenes::loadFileServices  (FileServices* loader)
{ 
  TTNeutralGenesFH* writer;

  if(_writers.size() != 0) {
      for(unsigned int i = 0; i < _writers.size(); i++)
      delete _writers[i];
      _writers.clear();
  }
  // --- THE READER ---
  // note that this reader only reads in FSTAT format
  //always add the reader:
  writer = new TTNeutralGenesFH(this);
  //set to read mode:
  writer->set_isInputHandler(true);
  //attach to file manager:
  loader->attach_reader(writer);
  //add to list of FileHandlers, will be deleted upon destruction
  _writers.push_back( writer );

  // --- WRITERS ---
  Param* param = get_parameter("ntrl_output_logtime");
  TMatrix temp;
  bool isMatrix = 0;

  if(param->isMatrix()) {param->getMatrix(&temp); isMatrix = 1;}

  if( get_parameter("ntrl_save_genotype")->isSet() ) {

      writer = new TTNeutralGenesFH(this);

      if(isMatrix) writer->set_multi(true, true, 1, &temp, get_parameter("ntrl_output_dir")->getArg());
      //           rpl_per, gen_per, rpl_occ, gen_occ, rank, path, self-ref
      else writer->set(true, param->isSet(), 1, (param->isSet() ? (int)param->getValue() : 0), 0, get_parameter("ntrl_output_dir")->getArg(), this);

      if(get_parameter("ntrl_save_genotype")->getArg().length() == 0 || //default mode is tabular
        get_parameter("ntrl_save_genotype")->getArg() == "TAB" ||
        get_parameter("ntrl_save_genotype")->getArg() == "tab")
      writer->set_write_fct( &TTNeutralGenesFH::write_TAB );
      else if(get_parameter("ntrl_save_genotype")->getArg() == "FSTAT" || get_parameter("ntrl_save_genotype")->getArg() == "fstat") {
        writer->set_write_fct( &TTNeutralGenesFH::write_FSTAT );  writer->set_extension(".dat");
      } else if(get_parameter("ntrl_save_genotype")->getArg() == "GENEPOP" || get_parameter("ntrl_save_genotype")->getArg() == "genepop") {
        writer->set_write_fct( &TTNeutralGenesFH::write_GENEPOP );
      } else if(get_parameter("ntrl_save_genotype")->getArg() == "PLINK" || get_parameter("ntrl_save_genotype")->getArg() == "plink") {
        writer->set_write_fct( &TTNeutralGenesFH::write_PLINK );
      } else
      fatal("parameter \"ntrl_save_genotype\" options are \"TAB\", \"FSTAT\" or \"GENEPOP\".\n");

      loader->attach(writer);

      _writers.push_back( writer );
  }

  if( get_parameter("ntrl_save_fsti")->isSet() ) {

      writer = new TTNeutralGenesFH(this);

      if(isMatrix) writer->set_multi(true, true, 1, &temp, get_parameter("ntrl_output_dir")->getArg());
      //           rpl_per, gen_per, rpl_occ, gen_occ, rank, path, self-ref
      else writer->set(true, param->isSet(), 1, (param->isSet() ? (int)param->getValue() : 0), 0, get_parameter("ntrl_output_dir")->getArg(), this);

      writer->set_write_fct( &TTNeutralGenesFH::write_Fst_i );

      writer->set_extension(".fsti");

      loader->attach(writer);

      _writers.push_back( writer );
  }

  if( get_parameter("ntrl_save_freq")->isSet() ) {

      writer = new TTNeutralGenesFH(this);

      writer->set_extension(".freq");

      if (isMatrix) writer->set_multi(true, true, 1, &temp, get_parameter("ntrl_output_dir")->getArg());
      //           rpl_per, gen_per, rpl_occ, gen_occ, rank, path, self-ref
      else writer->set(true, param->isSet(), 1, (param->isSet() ? (int)param->getValue() : 0), 0, get_parameter("ntrl_output_dir")->getArg(), this);

      //check validity of output option:

      if(get_parameter("ntrl_save_freq")->getArg() == "allfreq") {

        warning("option \"allfreq\" to \"ntrl_save_freq\" is deprecated, using \"locus\" instead\n");

        writer->setOutputOption("locus");

      } else if(get_parameter("ntrl_save_freq")->getArg() == "vcomp") {

        warning("option \"vcomp\" to \"ntrl_save_freq\" is deprecated, using \"allele\" instead\n");

        writer->setOutputOption("allele");

      } else if( get_parameter("ntrl_save_freq")->getArg() != "1" &&
        get_parameter("ntrl_save_freq")->getArg() != "locus" &&
        get_parameter("ntrl_save_freq")->getArg() != "allele") {

        fatal("parameter \"ntrl_save_freq\" options are \"locus\" or \"allele\".\n");

      } else
      writer->setOutputOption(get_parameter("ntrl_save_freq")->getArg());

      writer->set_write_fct( &TTNeutralGenesFH::write_varcompWC );

      loader->attach(writer);

      _writers.push_back( writer );
  }

}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoNeutralGenes::loadStatServices  (StatServices* loader)
{
  //allocate the stat handler
  if(_stats == NULL)
    _stats = new TTNeutralGenesSH(this);

  if(_stats != NULL) {
      loader->attach(_stats);
  }
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTNeutralGenes* TProtoNeutralGenes::hatch ()
{
  TTNeutralGenes* kid = new TTNeutralGenes();

  kid->set_proto(this);
  kid->set_locus_num(_locus_num);
  kid->set_allele_num(_allele_num);
  kid->set_init_model(_init_model);
  kid->set_mut_model(_mut_model);
  kid->set_mut_rate(_mut_rate);
  kid->set_mut_func_ptr(_mutate_func_ptr);
  kid->set_2L(_ploidy*_locus_num);
  kid->set_inherit_func_ptr(_inherit_func_ptr);

  return kid;
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoNeutralGenes::inherit_free (sex_t SEX, unsigned char* seq, unsigned char** parent)
{
  for(unsigned int i = 0; i < _locus_num; ++i)
    seq[i] = parent[RAND::RandBool()][i];
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoNeutralGenes::inherit_low (sex_t SEX, unsigned char* seq, unsigned char** parent)
{
  register unsigned int prevLoc = 0, chrm_bloc;
  register bool flipper;

  vector< unsigned int >& recTable = _map.getRecLoci(SEX, _mapIndex);
  vector< bool > & firstRecPos = _map.getFirstRecPosition(SEX);

  unsigned int nbRec = recTable.size();

  //      cout << "TProtoNeutralGenes::inherit_low; sex="<<SEX<<" nb Rec = "<<nbRec-1<<" (recLoc[0]="<<recTable[0]<<")"<<endl;

  //  if (!nbRec) {
  //    memcpy(&seq[0], &parent[flipper][prevLoc], (recTable[rec] - prevLoc));
  //    return;
  //  }
  for(unsigned int c = 0, stride = 0, rec = 0; c < _numChromosome; ++c) {

      flipper = firstRecPos[c];

      chrm_bloc = stride + _numLociPerChrmsm[c]; //number of loci considered so far

      prevLoc = stride; //stride is the first locus of a chromosome

      //do recombination chromosome-wise
      //    cout<<"chrm "<<c<<" side="<<firstRecPos[c]<<endl;

      for(; rec < nbRec && recTable[rec] < chrm_bloc; rec++) {

        //      cout << " --copy from "<<prevLoc<<" to "<<recTable[rec]<<" (side="<<flipper<<")"<<endl;
        //           <<"("<<(recTable[rec] - prevLoc)<<" loc) on side "<<flipper<<endl;

        memcpy(&seq[prevLoc], &parent[flipper][prevLoc], (recTable[rec] - prevLoc));

        prevLoc = recTable[rec];

        flipper = !flipper;
      }
      //    cout << "copy end of chrmsm from "<<prevLoc<<" to "<<chrm_bloc
      //         <<" ("<<(chrm_bloc - prevLoc)<<" loc) on side "<<flipper<<endl;
      //    cout << " --copy from "<<prevLoc<<" to "<<chrm_bloc<<" (side="<<flipper<<")"<<endl;
      //copy what's left between the last x-over point and the end of the chrmsme
      memcpy(&seq[prevLoc], &parent[flipper][prevLoc], (chrm_bloc - prevLoc));

      stride += _numLociPerChrmsm[c];

  }
  //  cout<<" done"<<endl;
  //  cout << "parent chromosomes:\n --0:";
  //  for(unsigned int i = 0; i < _locus_num; ++i)
  //    cout << (int)parent[0][i];
  //  cout << "\n --1:";
  //  for(unsigned int i = 0; i < _locus_num; ++i)
  //    cout << (int)parent[1][i];
  //  cout << "\ngamete:\n --0:";
  //  for(unsigned int i = 0; i < _locus_num; ++i)
  //    cout << (int)seq[i];
  //  cout << "\n";

}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
bool TProtoNeutralGenes::retrieve_data ( BinaryStorageBuffer* reader )
{
  unsigned int dummy;
  reader->read(&dummy, sizeof(int));
  if(dummy != _locus_num ){
      error("TProtoDeletMutations::retrieve_data:nb locus in file differ from parameter value!\n");
      _locus_num = dummy;
  }
  return true;
}
// ----------------------------------------------------------------------------------------

//                                 N E U T R A L   T R A I T

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTNeutralGenes& TTNeutralGenes::operator= (const TTrait& T)
{
  const TTNeutralGenes& TN = dynamic_cast<const TTNeutralGenes&> (T);

  if(this != &TN) {

      _locus_num = TN._locus_num;
      _allele_num = TN._allele_num;
      _ploidy = 2; //by default
      _2L = TN._2L;

      reset();

      init();

      for(unsigned int i = 0; i < _locus_num; ++i) {
        //      for(int j = 0; j < _ploidy; ++j)
        _sequence[0][i] = TN._sequence[0][i];
        _sequence[1][i] = TN._sequence[1][i];
      }
  }

  return *this;
}
//-----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTNeutralGenes::operator== (const TTrait& T)
{
  if(_type.compare(T.get_type()) != 0) return false;
  const TTNeutralGenes& TN = dynamic_cast<const TTNeutralGenes&> (T);

  if(this != &TN) {
      if(_locus_num != TN._locus_num) return false;
      if(_allele_num != TN._allele_num) return false;
      //    if(_ploidy != TN._ploidy) return false;
  }
  return true;
}
//-----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTNeutralGenes::operator!= (const TTrait& T)
{
  if(!((*this) == T))
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTNeutralGenes::~TTNeutralGenes()
{
  if(_sequence != NULL){
      for(unsigned int i = 0; i < _ploidy; i++)
      delete [] _sequence[i];
      delete [] _sequence;
  }
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::init ()
{
  if(_sequence != NULL)
    fatal("TTNeutralGenes::init::_sequence is not NULL !\n");

  _sequence = new unsigned char* [_ploidy];

  for(unsigned int i = 0; i < _ploidy; i++) {
      _sequence[i] = new unsigned char[_locus_num];
      for(unsigned int j = 0; j < _locus_num; ++j)
      _sequence[i][j] = 1;
  }
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::init_sequence ()
{
  if(_sequence == NULL) { //that should never happen
    _sequence = new unsigned char* [_ploidy];
    for(unsigned int i = 0; i < _ploidy; i++)
      _sequence[i] = new unsigned char[_locus_num];
  }

  if(_init_model == 1) {

    for(unsigned int i = 0; i < _ploidy; ++i)
      for(unsigned int j = 0; j < _locus_num; j++) {
          if(_allele_num > 2)
            _sequence[i][j] = (unsigned char)RAND::Uniform(_allele_num);
          else
            _sequence[i][j] = (unsigned char)RAND::Bernoulli(0.5);
      }

  } else if (_init_model == 3 ) {

      for(unsigned int i = 0; i < _ploidy; ++i)
        for(unsigned int j = 0; j < _locus_num; j++)
          _sequence[i][j] = 0;

      mutate(); //do one round of mutation

  } else { //model 0; all monomorphic

    for(unsigned int i = 0; i < _ploidy; ++i)
      for(unsigned int j = 0; j < _locus_num; j++)
        _sequence[i][j] = 0;
  }
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::reset()
{
  if(_sequence != NULL) {

      for(unsigned int i = 0; i < _ploidy; i++)  delete [] _sequence[i];

      delete [] _sequence;

      _sequence = NULL;
  }
}
// ----------------------------------------------------------------------------------------
// set_sequence
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::set_sequence(void** seq)
{
  reset(); init();
  memcpy(_sequence[0],seq[0],_locus_num);
  memcpy(_sequence[1],seq[1],_locus_num);
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::show_up ()
{
  message("\n  Trait's type: ntrl\n\
       locus: %i\n\
     alleles: %i\n\
    _sequence:",_locus_num, _allele_num);

  for(unsigned int i = 0; (i < _locus_num && i < 10); i++)
    message("\n              %i %i",(int)_sequence[0][i],(int)_sequence[1][i]);
  if(_locus_num > 10) message("\n              ...");
  message("\n");
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::inherit (TTrait* mother, TTrait* father)
{
  unsigned char** mother_seq = (unsigned char**)mother->get_sequence();
  unsigned char** father_seq = (unsigned char**)father->get_sequence();

  (_myProto->* _inherit_func_ptr) (FEM, _sequence[FEM], mother_seq);

  (_myProto->* _inherit_func_ptr) (MAL, _sequence[MAL], father_seq);
}
// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::mutate_SSM ()
{
  unsigned int mutAll, mutLocus;
  bool direction;
  unsigned int NbMut = RAND::Binomial(_mut_rate, _2L);

  for(; NbMut != 0; NbMut--)
    {
      mutLocus = RAND::Uniform( _locus_num);
      mutAll = RAND::RandBool(); //RAND::Uniform( _ploidy); <--- only diploids for now
      direction = RAND::RandBool();
      //alleles values are from 0 to NtrlAll - 1 !!!
      if(direction && _sequence[mutAll][mutLocus] < _allele_num-1)
      _sequence[mutAll][mutLocus] += 1; //one step to the right
      else if(_sequence[mutAll][mutLocus] > 0) // !direction || all==_allele_num
      _sequence[mutAll][mutLocus] -= 1; //one step to the left
      else //!direction && all == 0
      _sequence[mutAll][mutLocus] += 1;
    }
}
// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void TTNeutralGenes::mutate_KAM ()
{
  unsigned int mutAll, mutLocus;
  unsigned char mut;
  unsigned int NbMut = RAND::Binomial(_mut_rate, _2L);

  for(; NbMut != 0; NbMut--) {
      mutLocus = RAND::Uniform( _locus_num);
      mutAll = RAND::RandBool(); //RAND::Uniform( _ploidy); <--- only diploids for now
      //assign an arbitrary allele value:
      do{
        mut = (unsigned char) (RAND::Uniform(_allele_num));
      } while (mut == _sequence[mutAll][mutLocus]);
      _sequence[mutAll][mutLocus] = mut;
  }
}
// ----------------------------------------------------------------------------------------
// mutate_2all
// ----------------------------------------------------------------------------------------
inline void TTNeutralGenes::mutate_2all ()
{
  register unsigned int mutAll, mutLocus;
  unsigned int NbMut = RAND::Binomial(_mut_rate, _2L);

  for( ; NbMut > 0; NbMut--) {
      mutLocus = RAND::Uniform( _locus_num);
      mutAll = RAND::RandBool(); //RAND::Uniform( _ploidy); <--- only diploids for now
      _sequence[mutAll][mutLocus] = !_sequence[mutAll][mutLocus];
  }
}

// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   F I L E   H A N D L E R

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::FHwrite  () 
{
  if(!_pop->isAlive()) return;

  (this->*write_fct)();
}
// ----------------------------------------------------------------------------------------
// write_TAB
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_TAB ()
{
  /**The file format is FSTAT-like, with age class and sex added after the pop id. The file
   * extension is ".dat".*/
  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned char** seq;
  Patch* current_patch;
  Individual *ind;

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

#ifdef _DEBUG_
  message("TTNeutralGenesFH::write_TAB (%s)\n",filename.c_str());
#endif

  ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open TAB output file!!\n");

  FILE<<"pop";

  for (unsigned int i = 0; i < nb_locus; ++i)
    for (unsigned int j = 0; j < ploidy; ++j)
      FILE<<" l"<<i+1<<"a"<<j+1;
  //add names for the three last fields:
  FILE<<" age sex home ped isMigrant father mother ID\n";

  for (unsigned int i = 0; i < patchNbr; ++i) {

      current_patch = _pop->getPatch(i);

      for (unsigned int j = 0; j < current_patch->size(FEM, OFFSx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(FEM, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE<<(unsigned int)(seq[l][k]+1)<<" ";
            }
        }
        FILE << OFFSx <<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
            << (ind->getFather() && ind->getMother() ?
              (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)
              <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;

      }

      for (unsigned int j = 0; j < current_patch->size(MAL, OFFSx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(MAL, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE<<(unsigned int)(seq[l][k]+1)<<" ";
            }
        }
        FILE << OFFSx <<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
            << (ind->getFather() && ind->getMother() ?
              (ind->getFather()->getHome() != i) + (ind->getMother()->getHome() != i) : 0)
              <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;

      }

      for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(FEM, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE<<(unsigned int)(seq[l][k]+1)<<" ";
            }
        }
        FILE << ADLTx <<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
            << (ind->getFather() && ind->getMother() ?
              (ind->getFather()->getHome() != i) + (ind->getMother()->getHome() != i) : 0)
              <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(MAL, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE<<(unsigned int)(seq[l][k]+1)<<" ";
            }
        }
        FILE << ADLTx <<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
            << (ind->getFather() && ind->getMother() ?
              (ind->getFather()->getHome() != i) + (ind->getMother()->getHome() != i) : 0)
              <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;
      }

  }

  FILE.close();

}
// ----------------------------------------------------------------------------------------
// write_PLINK
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_PLINK ()
{
  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned char** seq;
  Patch* current_patch;
  Individual *ind;

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + ".ped";

#ifdef _DEBUG_
  message("TTNeutralGenesFH::write_PLINK (%s)\n",filename.c_str());
#endif

  // the PED file -------------------------------------------------------------------------
  ofstream PED (filename.c_str(), ios::out);

  if(!PED) fatal("could not open plink .ped output file!!\n");

  // the FAM file -------------------------------------------------------------------------
  filename = get_path() + this->get_service()->getGenerationReplicateFileName() + ".fam";

  ofstream FAM (filename.c_str(), ios::out);

  if(!FAM) fatal("could not open plink .fam output file!!\n");

  char BASE[2] = {'A','G'};

  // BED files (binary output)
//  filename = get_path() + this->get_service()->getGenerationReplicateFileName() + ".bed";
//
//  ofstream BED (filename.c_str(), ios::out|ios::binary);
//
//  if(!BED) fatal("could not open plink .bed output file!!\n");
//
//  write_PLINK_BED(BED);


  unsigned long uniq_fam_id = 1;

  for (unsigned int i = 0; i < patchNbr; ++i) {

      current_patch = _pop->getPatch(i);

      for (unsigned int j = 0; j < current_patch->size(FEM, OFFSx); ++j) {

        ind = current_patch->get(FEM, OFFSx, j);

        PED<<"fam"<<uniq_fam_id<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
            <<" 2 -9 ";

        FAM<<"fam"<<uniq_fam_id++<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
           <<" 2 -9 "<<endl;

        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              PED<< BASE[ (unsigned int)(seq[l][k]) ]<<" ";
            }
        }

        PED <<std::endl;

      }

      for (unsigned int j = 0; j < current_patch->size(MAL, OFFSx); ++j) {

        ind = current_patch->get(MAL, OFFSx, j);

        PED<<"fam"<<uniq_fam_id<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
                <<" 1 -9 ";

        FAM<<"fam"<<uniq_fam_id++<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
               <<" 1 -9 "<<endl;

        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
          for (unsigned int l = 0; l < ploidy; ++l) {
            PED<< BASE[ (unsigned int)(seq[l][k]) ] <<" ";
          }
        }
        PED<<std::endl;

      }

      for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {

	        ind = current_patch->get(FEM, ADLTx, j);
        PED<<"fam"<<uniq_fam_id<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
            <<" 2 -9 ";

        FAM<<"fam"<<uniq_fam_id++<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
           <<" 2 -9 "<<endl;

        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              PED<< BASE[ (unsigned int)(seq[l][k]) ]<<" ";
            }
        }
        PED<<std::endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {

	        ind = current_patch->get(MAL, ADLTx, j);
        PED<<"fam"<<uniq_fam_id<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
            <<" 1 -9 ";

        FAM<<"fam"<<uniq_fam_id++<<" "<<ind->getID()<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()
                    <<" 1 -9 "<<endl;

        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              PED<< BASE[ (unsigned int)(seq[l][k]) ] <<" ";
            }
        }
        PED<<std::endl;
      }

  }

  // the MAP file -------------------------------------------------------------------------
  filename = get_path() + this->get_service()->getGenerationReplicateFileName() + ".map";

  ofstream MAP (filename.c_str(), ios::out);

  if(!MAP) fatal("could not open plink .map output file!!\n");

  // the BIM file -------------------------------------------------------------------------
  filename = get_path() + this->get_service()->getGenerationReplicateFileName() + ".bim";

  ofstream BIM (filename.c_str(), ios::out);

  if(!BIM) fatal("could not open plink .bim output file!!\n");

  double **map;
  map = new double* [nb_locus];
  for(unsigned int k = 0; k < nb_locus; ++k) map[k] = new double[2];

  if( _FHLinkedTrait->_map.getGeneticMap(_FHLinkedTrait->get_type(), map, nb_locus)) {

    for(unsigned int k = 0; k < nb_locus; ++k) {
      MAP<<map[k][0]+1<<" "<<"loc"<<k+1<<" "<<map[k][1]+1<<" "<<k+1<<endl;
      BIM<<map[k][0]+1<<" "<<"loc"<<k+1<<" "<<map[k][1]+1<<" "<<k+1<<" A G"<<endl; //A is the MAJ
    }
  }

  PED.close();
  FAM.close();
  MAP.close();
  BIM.close();
//  BED.close();

  for(unsigned int k = 0; k < nb_locus; ++k) delete [] map[k];
  delete [] map;
}
// ----------------------------------------------------------------------------------------
// write_PLINK_BED
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_PLINK_BED (ofstream &BED)
{
  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned int popsize = _pop->size();
  unsigned int n_blocs = (unsigned int)ceil(popsize / 4.0);
  Patch* current_patch;
  
  unsigned char BIN_GENO[4] = {0,1,2,3}; //00: homoz A/0; 01: missing; 10: heteroz; 11: homoz G/1

  unsigned char BIT_MASK = 0x3; //number 3 or '11'

  unsigned char MAGIC[4] = {0x6c,0x1b,0x01,'\0'};

  BED<<MAGIC;

  unsigned char BYTE = 0;
  
  bitstring DATA(n_blocs * 8); //one bloc is one byte
  
  for(unsigned int k = 0, pos; k < nb_locus; ++k) {

      DATA.reset();
      
      pos = 0;
      
      for (unsigned int i = 0; i < patchNbr; ++i) {

        current_patch = _pop->getPatch(i);

//        store_BED(FEM, OFFSx, BED, current_patch, DATA, &pos, k);
//        store_BED(MAL, OFFSx, BED, current_patch, DATA, &pos, k);
//        store_BED(FEM, ADLTx, BED, current_patch, DATA, &pos, k);
//        store_BED(MAL, ADLTx, BED, current_patch, DATA, &pos, k);
      }

  }
}

//void TTNeutralGenesFH::store_BED (sex_t SEX, age_idx age, ofstream &BED, Patch *patch,
//				  bitstring &DATA, unsigned int *pos, unsigned int locus)
//{
//  Individual *ind;
//  unsigned char** seq;
//
//  for (unsigned int j = 0; j < patch->size(SEX, age); ++j) {
//
//      ind = patch->get(SEX, age,j);
//      seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
//
//      if()
//  }
//
//}
// ----------------------------------------------------------------------------------------
// write_FSTAT
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_FSTAT ()
{
  /**The file format is FSTAT-like, with extra info about the individual added (age, sex, pedigree, origin). The file
   * extension is ".dat".*/
  unsigned int position;
  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_all   = _FHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned char** seq;
  Patch* current_patch;
  Individual *ind;

  position = nb_all > 99 ? 3 : 2; //assumes nb_all not sup. to 999

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

#ifdef _DEBUG_
  message("TTNeutralGenesFH::FHwrite (%s)\n",filename.c_str());
#endif

  ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open FSTAT output file!!\n");

  FILE<<patchNbr<<" "<< nb_locus + 4 <<" "<<nb_all<<" "<<position<<"\n";

  for (unsigned int i = 0; i < nb_locus; ++i)
    FILE<<"loc"<<i+1<<"\n";
  //add names for the three last fields:
  FILE<<"age\n"<<"sex\n"<<"ped\n"<<"origin\n";

  for (unsigned int i = 0; i < patchNbr; ++i) {

      current_patch = _pop->getPatch(i);

      for (unsigned int j = 0; j < current_patch->size(FEM, OFFSx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(FEM, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << OFFSx << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, OFFSx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(MAL, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << OFFSx << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(FEM, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << ADLTx << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {

        FILE<<i+1<<" ";
        ind = current_patch->get(MAL, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << ADLTx << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

  }

  FILE.close();

}
// ----------------------------------------------------------------------------------------
// write_GENEPOP
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_GENEPOP ()
{
  /**The file format is GENEPOP-like, with extra info about the individual added (age, sex, pedigree, origin). The file
   * extension is ".dat".*/
  unsigned int position;
  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_all   = _FHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned char** seq;
  Patch* current_patch;
  Individual *ind;

  position = nb_all > 99 ? 3 : 2; //assumes nb_all not sup. to 999

  if(ploidy == 1) position = 1;

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

#ifdef _DEBUG_
  message("TTNeutralGenesFH::FHwrite (%s)\n",filename.c_str());
#endif

  ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open FSTAT output file!!\n");

  FILE<<"Title line: "<<patchNbr<<" patches, "<<nb_locus<<" loci with "<<nb_all<<" alleles\n";

  for (unsigned int i = 0; i < nb_locus; ++i)
    FILE<<"loc"<<i+1<<", ";

  //add names for the three last fields:
  FILE<<"age, sex, ped, origin\n";

  for (unsigned int i = 0; i < patchNbr; ++i) {

      current_patch = _pop->getPatch(i);

      FILE<<"POP\n";

      for (unsigned int j = 0; j < current_patch->size(FEM, OFFSx); ++j) {

        FILE<<i+1<<", ";
        ind = current_patch->get(FEM, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << OFFSx << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, OFFSx); ++j) {

        FILE<<i+1<<", ";
        ind = current_patch->get(MAL, OFFSx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << OFFSx << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {

        FILE<<i+1<<", ";
        ind = current_patch->get(FEM, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << ADLTx << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

      for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {

        FILE<<i+1<<", ";
        ind = current_patch->get(MAL, ADLTx, j);
        seq = (unsigned char**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();

        for(unsigned int k = 0; k < nb_locus; ++k) {
            for (unsigned int l = 0; l < ploidy; ++l) {
              FILE.fill('0');
              FILE.width(position);
              FILE<<(unsigned int)(seq[l][k]+1);
            }
            FILE<<" ";
        }
        FILE << ADLTx << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
      }

  }

  FILE.close();

}
// ----------------------------------------------------------------------------------------
// setOutputOption
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::setOutputOption(string opt)
{

  if(opt != "1")
    _output_option = opt;
  else
    _output_option = "locus"; //default
}
// ----------------------------------------------------------------------------------------
// FHread
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::FHread (string& filename)
{ /** reads in only FSTAT format */
  unsigned int digit, nloci_infile, nall_infile, npatch_infile;
  //  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
  unsigned int nb_all   = _FHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();

  bool is_extended = false;

  ifstream FILE(filename.c_str(),ios::in);

  if(!FILE) fatal("could not open FSTAT input file \"%s\"\n", filename.c_str());

  FILE >> npatch_infile >> nloci_infile >> nall_infile >> digit;

  if(npatch_infile != patchNbr) fatal("# of patch in FSTAT file differs from simulation settings\n");

  if(nloci_infile != nb_locus && nloci_infile != nb_locus +4) {

      fatal("# of loci in FSTAT file (%i loci) differs from simulation settings (%i loci)\n",nloci_infile, nb_locus);

  } else if (nloci_infile == nb_locus + 4) {
      is_extended = true;
  } else
    is_extended = false;

  if(nall_infile != nb_all)  fatal("# of alleles in FSTAT file differs from simulation settings\n");

  digit = (int) pow(10.0,(double)digit);
  string loc_name;

  for(unsigned int i = 0; i < nloci_infile; ++i) {
      FILE >> loc_name; //we don't do anything with this here
  }


  unsigned int pop, genot, all0, all1, age, sex, ped, origin;
  age_idx agex;
  Individual *ind;
  unsigned char* seq[2];
  int lnbr = nloci_infile +2; //line counter
  seq[0] = new unsigned char[nb_locus];
  seq[1] = new unsigned char[nb_locus];

  while(FILE>>pop) {
      //check pop value
      if(pop > patchNbr) fatal("found an illegal pop identifier in FSTAT file (%i) at line %i",pop, lnbr);
      for(unsigned int i = 0; i < nb_locus; ++i) {
        FILE>>genot;

        all0 = genot/digit;
        all1 = genot%digit;

        if(all0 <= nb_all)
          seq[0][i] = all0 - 1;
        else {
            error("in FSTAT input file at line %i, locus %i : \
first allele value %d is greater than the max value specified (%i)!\n", lnbr, i+1, all0, nb_all);
            fatal("Please check the input file.\n");
        }

        if(all1 <= nb_all)
          seq[1][i] = all1 - 1;
        else {
            error("in FSTAT input file at line %i, locus %i : \
second allele value %i is greater than the max value specified (%i)!\n", lnbr, i+1, all1, nb_all);
            fatal("Please check the input file.\n");
        }
      }

      if(is_extended) {

        FILE >> age >> sex >> ped >> origin;
        //age index now saved in the FSTAT file
        agex = (static_cast<age_idx> (age) == ADLTx ? ADLTx : OFFSx);

      } else {
        // by default, individuals are adult females
        agex = ADLTx;
        sex = FEM;
        ped = 0;
        origin = 0;

      }


      ind = _pop->makeNewIndividual(0, 0, static_cast<sex_t> (sex), origin - 1);
      ind->setPedigreeClass((unsigned char)ped);
      ind->getTrait(_FHLinkedTraitIndex)->set_sequence((void**)seq);
      _pop->getPatch(pop-1)->add(static_cast<sex_t> (sex), agex, ind);

      lnbr++;

  }

  FILE.close();

  delete [] seq[0];
  delete [] seq[1];
}
// ----------------------------------------------------------------------------------------
// write_Fst_i
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_Fst_i ()
{  
  if(_pop->size(ADULTS) == 0) {
      warning("No adults in pop, not writing the Fst distribution to file.\n");
      return;
  }

  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  double **fst_i;
  bool added_stater = false;

  TTNeutralGenesSH* stater = _FHLinkedTrait->get_stater();

  if(stater == NULL) {
      stater = new TTNeutralGenesSH(_FHLinkedTrait);
      stater->allocateTables(_FHLinkedTrait->get_locus_num() , _FHLinkedTrait->get_allele_num());
      added_stater = true;
  }

  fst_i = new double* [patchNbr];
  for(unsigned int i = 0; i < patchNbr; i++)
    fst_i[i] = new double [nb_locus];

  stater->setFst_li(patchNbr, nb_locus, fst_i);

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

#ifdef _DEBUG_
  message("TTNeutralGenesFH::FHwrite (%s)\n",filename.c_str());
#endif

  ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open Fst_i output file!!\n");

  for(unsigned int i = 0; i < patchNbr; ++i)
    FILE<<"patch"<<i+1<<" ";

  FILE<<endl;

  for(unsigned int j = 0; j < nb_locus; ++j) {

      for(unsigned int i = 0; i < patchNbr; i++)
      FILE<<fst_i[i][j]<<" ";

      FILE<<endl;
  }

  FILE.close();

  if(added_stater) delete stater;

  for(unsigned int i = 0; i < patchNbr; i++)
    delete [] fst_i[i];

  delete [] fst_i;

}

// ----------------------------------------------------------------------------------------
// write_varcompWC
// ----------------------------------------------------------------------------------------
void TTNeutralGenesFH::write_varcompWC ()
{  
  unsigned int nb_allele = _FHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  bool added_stater = false;
  age_t AGE = _pop->getCurrentAge();

  if(AGE == (ADULTS | OFFSPRG) || AGE == ALL) {
      warning("saving only offspring neutral allele stats\n");
      AGE = OFFSPRG;
  }

  age_idx age = (AGE == OFFSPRG ? OFFSx : ADLTx);

  TTNeutralGenesSH* stater = _FHLinkedTrait->get_stater();

  if(stater == NULL) {
      stater = new TTNeutralGenesSH(_FHLinkedTrait);
      stater->allocateTables(_FHLinkedTrait->get_locus_num() , _FHLinkedTrait->get_allele_num());
      added_stater = true;
  }

  stater->setAlleleTables(AGE);
  stater->setHeteroTable(AGE);

  TMatrix *globalFreq = stater->getGlobalFreqs();
  DataTable< double > *freqTable = stater->getAlleleFreqTable();
  DataTable< double > *heteroTable = stater->getHeteroTable();
  DataTable< unsigned int > *alleleCountTable = stater->getAlleleCountTable();

  std::string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

#ifdef _DEBUG_
  message("TTNeutralGenesFH::write_varcompWC (%s)\n",filename.c_str());
#endif

  //init
  double* pop_sizes = new double [patchNbr];
  double tot_size;
  double sum_weights = 0;
  double nc;
  unsigned int extantPs = 0;

  tot_size = _pop->size(AGE);

  for(unsigned int i = 0; i < patchNbr; i++) {
      pop_sizes[i] = _pop->size(AGE, i);
      if(pop_sizes[i]) {
        extantPs++;
        sum_weights += (pop_sizes[i] * pop_sizes[i] / tot_size);
      }
  }

  nc = (tot_size - sum_weights)/(extantPs-1);

  //      unsigned int np = extantPs;
  unsigned int npl = extantPs; //all loci typed in all patches

  //p = _alleleFreqTable
  //pb = _globalAlleleFreq

  //counting num alleles per locus

  unsigned int *alploc = new unsigned int [nb_locus];

  unsigned int **alploc_table = new unsigned int* [nb_locus];

  for(unsigned int i = 0; i < nb_locus; ++i)
    alploc_table[i] = new unsigned int[nb_allele];

  unsigned int tot_num_allele = 0;

  for(unsigned int l = 0; l < nb_locus; ++l){

      alploc[l] = 0;

      for(unsigned int cnt,  a = 0; a < nb_allele; ++a) {

        cnt=0;

        for(unsigned int i = 0; i < patchNbr; i++) {

            cnt += alleleCountTable->get(i,l,a);

        }
        alploc_table[l][a] = (cnt != 0);
        alploc[l] += (cnt != 0);
      }
      tot_num_allele += alploc[l];
  }

  //correspondance with hierfstat implementation:
  //n, and nal are given by pop_sizes, same num ind typed at all loci in each patch
  //nc is the same for each locus
  //nt is given by tot_size, same tot num of ind typed at all loci

  //SSG = het/2 for each allele
  double *SSG = new double[tot_num_allele];
  double *SSP = new double[tot_num_allele];
  double *SSi = new double[tot_num_allele];
  double *loc_id = new double[tot_num_allele];
  double *al_id = new double [tot_num_allele];

  unsigned int all_cntr = 0;

  double het, freq, var;

  //computing sum of squares

  for(unsigned int l = 0; l < nb_locus; ++l) {

      for(unsigned int a = 0; a < nb_allele & all_cntr < tot_num_allele; ++a) {

        if(alploc_table[l][a] == 0) continue; //do not consider alleles not present in the whole pop

        //store locus and all identifiers for output
        loc_id[all_cntr] = l+1;
        al_id[all_cntr] = a+1;

        SSG[all_cntr] = 0;
        SSi[all_cntr] = 0;
        SSP[all_cntr] = 0;

        for(unsigned int p = 0; p < patchNbr; ++p){

            if(!_pop->size(AGE, p)) continue; //skip empty patches

            het = heteroTable->get(p, l, a);

            freq = freqTable->get(p, l, a);

            var = freq - globalFreq->get(l, a); //(p_liu - pbar_u)^2

            var *= var;

            SSG[all_cntr] += het;

            SSi[all_cntr] += 2*pop_sizes[p]*freq*(1-freq) - het/2;

            SSP[all_cntr] += 2*pop_sizes[p]*var;
        }

        all_cntr++;
      }

  }

  assert(all_cntr == tot_num_allele);

  // open the file

  ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open neutral vcomp output file!!\n");

  // print column names

  if(_output_option == "allele") {

      FILE<<"locus allele pbar het siga sigb sigw Fst Fis";

      for (unsigned int p = 1; p <= patchNbr; ++p)
      FILE<<" het.p"<<p;

      for (unsigned int p = 1; p <= patchNbr; ++p)
      FILE<<" freq.p"<<p;

      FILE<<endl;

  } else {

      FILE<<"locus maj.al pbar.maj.al het siga sigb sigw Fst Fis";

      for (unsigned int p = 1; p <= patchNbr; ++p)
      FILE<<" het.p"<<p;

      //allele frequencies in each patch
      for (unsigned int p = 1; p <= patchNbr; ++p)
      //                  for(unsigned int u = 0; u < nb_allele-1; ++u) //skip last allele, can be deduced...
      FILE<<" freq.maj."<<"p"<< p;

      FILE<<endl;
  }

  //-----------------------------------------------------------------------------------------
  // allele specific stats:

  double *MSG = new double[tot_num_allele];
  double *MSP = new double[tot_num_allele];
  double *MSI = new double[tot_num_allele];
  //      double *sigw = new double[tot_num_allele];
  double *siga = new double[tot_num_allele];
  double *sigb = new double[tot_num_allele];

  for(unsigned int i = 0; i < tot_num_allele; ++i){

      MSG[i] = SSG[i] / (2 * tot_size);
      //            sigw[i] = MSG[i]; //wasted!

      MSP[i] = SSP[i] / (npl-1);

      MSI[i] = SSi[i]/ (tot_size - npl);

      sigb[i] = 0.5*(MSI[i] - MSG[i]);

      siga[i] = (MSP[i] - MSI[i])/(2*nc);

      if(_output_option == "allele") {
        FILE<< loc_id[i] << " " << al_id[i] << " "
            //global allele frequency
            << globalFreq->get(loc_id[i] - 1, al_id[i] - 1) << " "
            //average heterozygosity:
            << SSG[i]/tot_size << " "
            //variance components:
            << siga[i] << " " << sigb[i] << " " << MSG[i] << " "
            //Fst
            << siga[i]/(siga[i]+sigb[i]+MSG[i]) << " "
            //Fis
            << sigb[i]/(sigb[i]+MSG[i]) << " ";
        //per patch heterozygosity:
        for (unsigned int p = 0; p < patchNbr; ++p) {
            FILE << heteroTable->get(p, loc_id[i]-1, al_id[i]-1)/pop_sizes[p] << " ";
        }
        //per patch allele frequency:
        for (unsigned int p = 0; p < patchNbr; ++p) {
            FILE << freqTable->get(p, loc_id[i]-1, al_id[i]-1) << " ";
        }

        FILE<< endl;
      }
  }

  //-----------------------------------------------------------------------------------------
  // locus specific stats:

  if(_output_option == "locus") {


      double lsiga, lsigb, lsigw, max_all_frq, het;
      unsigned int maj_al;

      all_cntr = 0;

      deque <double> loc_het = stater->setHo2(age);

      for(unsigned int i = 0; i < nb_locus; ++i) {

        lsiga = 0;
        lsigb = 0;
        lsigw = 0;

        max_all_frq = 0;

        for(unsigned int l = 0; l < alploc[i]; ++l) {

            lsiga += siga[all_cntr];
            lsigb += sigb[all_cntr];
            lsigw += MSG[all_cntr];

            if(max_all_frq < globalFreq->get(i, al_id[all_cntr]-1 ) ) {
              max_all_frq = globalFreq->get(i, al_id[all_cntr]-1 );
              maj_al = al_id[all_cntr];
            }

            all_cntr++;

        }

        FILE << i+1 <<" "<< maj_al <<" "<< max_all_frq <<" ";
        FILE<< loc_het[i] <<" "<< lsiga << " " << lsigb <<" ";
        FILE<< lsigw << " "<< lsiga /(lsiga + lsigb + lsigw) <<" ";
        FILE<< lsigb /(lsigb + lsigw);
        // patch heterozygosities:

        for(unsigned int p = 0; p < patchNbr; ++p) {
            het = 0;
            for(unsigned int a = 0; a < nb_allele; ++a){
              het += heteroTable->get(p, i, a);
            }
            FILE << " " << het/(2.0*pop_sizes[p]);
        }
        // patch allele freq?
        for(unsigned int p = 0; p < patchNbr; ++p) {
            FILE<< " " << freqTable->get(p, i, maj_al-1);
        }
        FILE << endl;
      }//END for locus

  } //END per locus output




  FILE.close();

  delete[]pop_sizes;
  delete[]alploc;
  for(unsigned int i = 0; i < nb_locus; ++i)
    delete[]alploc_table[i];
  delete[]alploc_table;
  delete[]loc_id;
  delete[]al_id;
  delete[]SSG;
  delete[]SSi;
  delete[]SSP;
  delete[]MSG;
  delete[]MSI;
  delete[]MSP;
  //      delete[]sigw;
  delete[]siga;
  delete[]sigb;
}

// ----------------------------------------------------------------------------------------

//                     N E U T R A L   T R A I T   S T A T   H A N D L E R

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::init()
{
  StatHandler<TTNeutralGenesSH>::init();

  allocateTables(_SHLinkedTrait->get_locus_num(),_SHLinkedTrait->get_allele_num());
}
// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTNeutralGenesSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTNeutralGenesSH::setStatRecorders ");
#endif
  if(token.compare("coa") == 0) {

      add("Wtn Patch Coancestry (offsprg)","off.theta",OFFSPRG,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setOffsprgCoaMatrix);
      add("Btn Patch Coancestry (offsprg)","off.alpha",OFFSPRG,0,0,&TTNeutralGenesSH::getMeanAlpha,0,0,0);

      add("Wtn Patch Coancestry   (adult)","adlt.theta",ADULTS,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setAdultsCoaMatrix);
      add("Btn Patch Coancestry   (adult)","adlt.alpha",ADULTS,0,0,&TTNeutralGenesSH::getMeanAlpha,0,0,0);

  } else if(token.compare("adlt.coa") == 0) {

      add("Wtn Patch Coancestry   (adult)","adlt.theta",ADULTS,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setAdultsCoaMatrix);
      add("Btn Patch Coancestry   (adult)","adlt.alpha",ADULTS,0,0,&TTNeutralGenesSH::getMeanAlpha,0,0,0);

  } else if(token.compare("off.coa") == 0) {

      add("Wtn Patch Coancestry (offsprg)","off.theta",OFFSPRG,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setOffsprgCoaMatrix);
      add("Btn Patch Coancestry (offsprg)","off.alpha",OFFSPRG,0,0,&TTNeutralGenesSH::getMeanAlpha,0,0,0);

  } else if(token.compare("adlt.coa.persex") == 0) {

      add("Female Theta      (adult)","adlt.theta",ADULTS,0,0,
        &TTNeutralGenesSH::getMeanTheta, 0, 0, &TTNeutralGenesSH::setAdults_Theta);
      add("Female Theta      (adult)","adlt.thetaFF",ADULTS,0,0,&TTNeutralGenesSH::getTheta_FF,0,0,0);
      add("Male Theta        (adult)","adlt.thetaMM",ADULTS,0,0,&TTNeutralGenesSH::getTheta_MM,0,0,0);
      add("Female-Male Theta (adult)","adlt.thetaFM",ADULTS,0,0,&TTNeutralGenesSH::getTheta_FM,0,0,0);

  } else if(token.compare("adlt.coa.within") == 0) {

      add("Wtn Patch Coancestry   (adult)","adlt.theta",ADULTS,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setAdultsCoaWithin);

  } else if(token.compare("off.coa.within") == 0) {

      add("Wtn Patch Coancestry (offsprg)","off.theta",OFFSPRG,0,0,
        &TTNeutralGenesSH::getMeanTheta,0,0,&TTNeutralGenesSH::setOffsprgCoaWithin);

  } else if(token.compare("adlt.coa.between") == 0) {

      add("Btn Patch Coancestry   (adult)","adlt.alpha",ADULTS,0,0,
        &TTNeutralGenesSH::getMeanAlpha,0,0,&TTNeutralGenesSH::setAdultsCoaBetween);

  } else if(token.compare("off.coa.between") == 0) {

      add("Btn Patch Coancestry (offsprg)","off.alpha",OFFSPRG,0,0,
        &TTNeutralGenesSH::getMeanAlpha,0,0,&TTNeutralGenesSH::setOffsprgCoaBetween);

  } else if(token.compare("coa.matrix") == 0) {

      setCoaMatrixRecorders(OFFSPRG, 3);
      setCoaMatrixRecorders(ADULTS, 3);

  } else if(token.compare("off.coa.matrix") == 0) {

      setCoaMatrixRecorders(OFFSPRG, 3);

  } else if(token.compare("adlt.coa.matrix") == 0) {

      setCoaMatrixRecorders(ADULTS, 3);

  } else if(token.compare("coa.matrix.within") == 0) {

      setCoaMatrixRecorders(OFFSPRG, 1);
      setCoaMatrixRecorders(ADULTS, 1);

  } else if(token.compare("off.coa.matrix.within") == 0) {

      setCoaMatrixRecorders(OFFSPRG, 1);

  } else if(token.compare("adlt.coa.matrix.within") == 0) {

      setCoaMatrixRecorders(ADULTS, 1);

  } else if(token.compare("sibcoa") == 0) {

      add("Proportion of full-sib offspring","prop.fsib",OFFSPRG,0,0,
        0,&TTNeutralGenesSH::getSibProportions,0,&TTNeutralGenesSH::setSibStats);
      add("Proportion of paternal half-sib ","prop.phsib",OFFSPRG,1,0,0,&TTNeutralGenesSH::getSibProportions,0,0);
      add("Proportion of maternal half-sib ","prop.mhsib",OFFSPRG,2,0,0, &TTNeutralGenesSH::getSibProportions,0,0);
      add("Proportion of non-sib offspring ","prop.nsib",OFFSPRG,3,0,0,&TTNeutralGenesSH::getSibProportions,0,0);
      add("Coancestry of full-sib offspring","coa.fsib",OFFSPRG,0,0,0,&TTNeutralGenesSH::getSibCoaMeans,0,0);
      add("Coancestry of paternal half-sib ","coa.phsib",OFFSPRG,1,0,0,&TTNeutralGenesSH::getSibCoaMeans,0,0);
      add("Coancestry of maternal half-sib ","coa.mhsib",OFFSPRG,2,0,0,&TTNeutralGenesSH::getSibCoaMeans,0,0);
      add("Coancestry of non-sib offspring ","coa.nsib",OFFSPRG,3,0,0,&TTNeutralGenesSH::getSibCoaMeans,0,0);

  } else if (token == "ntrl.freq") {

      setFreqRecorders(OFFSPRG);
      setFreqRecorders(ADULTS);

  } else if (token == "off.ntrl.freq") {

      setFreqRecorders(OFFSPRG);

  } else if (token == "adlt.ntrl.freq") {

      setFreqRecorders(ADULTS);

      //  } else if (token == "ntrl.freq.patch") {
      //
      //    setFreqRecordersPerPatch(ALL);
      //
      //  } else if (token == "off.ntrl.freq.patch") {
      //
      //    setFreqRecordersPerPatch(OFFSPRG);
      //
      //  } else if (token == "adlt.ntrl.freq.patch") {
      //
      //    setFreqRecordersPerPatch(ADULTS);

  } else if(token.compare("offsprgfstat") == 0 || token.compare("off.fstat") == 0) {

      setFstatRecorders(OFFSPRG);

  } else if(token.compare("adultfstat") == 0 || token.compare("adlt.fstat") == 0) {

      setFstatRecorders(ADULTS);

  } else if(token.compare("fstat") == 0) {

      setFstatRecorders(ALL);

  } else if(token == "off.fstat2") {

      setFstat2Recorders(OFFSPRG);

  } else if(token == "adlt.fstat2") {

      setFstat2Recorders(ADULTS);

  } else if(token == "fstat2") {

      setFstat2Recorders(ALL);

  } else if(token.compare("fstWC") == 0 || token.compare("fstatWC") == 0) {

      setFstatWCRecorders(ALL);

  } else if(token.compare("off.fstWC") == 0 || token.compare("off.fstatWC") == 0) {

      setFstatWCRecorders(OFFSPRG);

  } else if(token.compare("adlt.fstWC") == 0 || token.compare("adlt.fstatWC") == 0) {

      setFstatWCRecorders(ADULTS);

  } else if(token.compare("weighted.fst") == 0) {

      setFstMatrixRecorders(OFFSPRG, 0);
      setFstMatrixRecorders(ADULTS, 0);

  } else if(token.compare("off.weighted.fst") == 0) {

      setFstMatrixRecorders(OFFSPRG, 0);

  } else if(token.compare("adlt.weighted.fst") == 0) {

      setFstMatrixRecorders(ADULTS, 0);

  } else if(token.compare("weighted.fst.matrix") == 0) {

      setFstMatrixRecorders(OFFSPRG, 3);
      setFstMatrixRecorders(ADULTS, 3);

  } else if(token.compare("off.weighted.fst.matrix") == 0) {

      setFstMatrixRecorders(OFFSPRG, 3);

  } else if(token.compare("adlt.weighted.fst.matrix") == 0) {

      setFstMatrixRecorders(ADULTS, 3);

  } else if(token.compare("weighted.fst.within") == 0) {

      setFstMatrixRecorders(OFFSPRG, 1);
      setFstMatrixRecorders(ADULTS, 1);

  } else if(token.compare("off.weighted.fst.within") == 0) {

      setFstMatrixRecorders(OFFSPRG, 1);

  } else if(token.compare("adlt.weighted.fst.within") == 0) {

      setFstMatrixRecorders(ADULTS, 1);

  } else if(token.compare("adlt.NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(ADULTS, true);

  } else if(token.compare("off.NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(OFFSPRG, true);

  } else if(token.compare("NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(ADULTS, true);
      setNeiGeneticDistanceRecorders(OFFSPRG, true);

  } else if(token.compare("mean.NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(ADULTS, false);
      setNeiGeneticDistanceRecorders(OFFSPRG, false);

  } else if(token.compare("adlt.mean.NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(ADULTS, false);

  } else if(token.compare("off.mean.NeiDistance") == 0) {

      setNeiGeneticDistanceRecorders(OFFSPRG, false);

  } else if(token.compare("Dxy") == 0) {

      setDxyRecorders(ADULTS, false);
      setDxyRecorders(OFFSPRG, false);

  } else if(token.compare("off.Dxy") == 0) {

      setDxyRecorders(OFFSPRG, false);

  } else if(token.compare("adlt.Dxy") == 0) {

      setDxyRecorders(ADULTS, false);

  } else if(token.compare("Dxy.patch") == 0) {

      setDxyRecorders(ADULTS, true);
      setDxyRecorders(OFFSPRG, true);

  } else if(token.compare("off.Dxy.patch") == 0) {

      setDxyRecorders(OFFSPRG, true);

  } else if(token.compare("adlt.Dxy.patch") == 0) {

      setDxyRecorders(ADULTS, true);

  } else
    return false;

  return true;
}
// ----------------------------------------------------------------------------------------
// setCoaMatrixRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setCoaMatrixRecorders (age_t AGE, unsigned char dim)
{
  std::ostringstream name, sub_name;

  void (TTNeutralGenesSH::* setter) () = (AGE == ADULTS ?
      ( dim & 3 ? &TTNeutralGenesSH::setAdultsCoaMatrix :
        dim & 2 ? &TTNeutralGenesSH::setAdultsCoaBetween :
            &TTNeutralGenesSH::setAdultsCoaWithin ) :
            ( dim & 3 ? &TTNeutralGenesSH::setOffsprgCoaMatrix :
              dim & 2 ? &TTNeutralGenesSH::setOffsprgCoaBetween :
                  &TTNeutralGenesSH::setOffsprgCoaWithin) );

  const char *prefix = (AGE == ADULTS ? "adlt." : "off.");

  unsigned int nbpatch =  _pop->getPatchNbr();
  unsigned int scale = (unsigned int)pow(10.0, (int)log10((float)nbpatch) + 1);

  if(dim & 1) {
      name<<"Wtn Patch Coancestry ("<<prefix<<")";
      sub_name<< prefix << "theta";
      add(name.str(), sub_name.str(),  AGE,0, 0, &TTNeutralGenesSH::getMeanTheta, 0, 0, setter);
      name.str("");
      sub_name.str("");
  }
  if(dim & 2){
      name<<"Btn Patch Coancestry ("<<prefix<<")";
      sub_name<< prefix << "alpha";
      add(name.str(), sub_name.str(),  AGE, 0, 0, &TTNeutralGenesSH::getMeanAlpha, 0, 0, (!(dim&1) ? setter : 0) );
      name.str("");
      sub_name.str("");
  }
  if(dim & 1) {
      for(unsigned int i = 0; i < nbpatch; ++i) {
        name<<"coancestry "<<i+1<<"."<<i+1;
        sub_name<< prefix << "coa" << i+1 << "." << i+1;
        add(name.str(), sub_name.str(),  AGE, i*scale + i, 0, 0, &TTNeutralGenesSH::getCoa, 0, 0);
        name.str("");
        sub_name.str("");
      }
  }
  if(dim & 2){
      for(unsigned int i = 0; i < nbpatch; ++i) {
        for(unsigned int j = i+1; j < nbpatch; ++j) {
            name<<"coancestry "<<i+1<< "." <<j+1;
            sub_name<< prefix << "coa" << i+1 << "." << j+1;
            add(name.str(), sub_name.str(),  AGE, i*scale + j, 0, 0, &TTNeutralGenesSH::getCoa, 0, 0);
            name.str("");
            sub_name.str("");
        }
      }
  }

}
// ----------------------------------------------------------------------------------------
// setFreqRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFreqRecorders (age_t AGE)
{
  string prefix = (AGE == OFFSPRG ? "off." : "adlt.");

  void (TTNeutralGenesSH::* setter) () = (AGE == ADULTS ? &TTNeutralGenesSH::setAdultAlleleFreq :
      &TTNeutralGenesSH::setOffspringAlleleFreq );

  unsigned int nb_allele = _SHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _SHLinkedTrait->get_locus_num();

  for (unsigned int l = 0; l < nb_locus; ++l) {
      for (unsigned int u = 0; u < nb_allele-1; ++u) {
        add("", prefix + "ntrl.l" + tstring::int2str(l+1) + ".a" + tstring::int2str(u+1), AGE, l, u,
            0, 0, &TTNeutralGenesSH::getGlobalAlleleFreq, setter);
      }
  }

  setter = (AGE == ADULTS ? &TTNeutralGenesSH::setAdultHeterozygosity :
      &TTNeutralGenesSH::setOffspringHeterozygosity );

  for (unsigned int l = 0; l < nb_locus; ++l) {
      add("", prefix + "ntrl.l" + tstring::int2str(l+1) + ".Het", AGE, l, 0, 0,
        &TTNeutralGenesSH::getHeterozygosity, 0, setter);
      //    add("", prefix + "ntrl.l" + tstring::int2str(l) + ".Hom", AGE, l, 0, 0,
      //        &TTNeutralGenesSH::getHomozygosity, 0, 0);
  }

}

// ----------------------------------------------------------------------------------------
// setFstatRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstatRecorders (age_t AGE)
{
  if(AGE & OFFSPRG) {
      add("Nbr of Alleles per Locus - local  (offsprg)","off.allnbp",OFFSPRG,0,0,
        &TTNeutralGenesSH::getNbAllLocal,0,0,&TTNeutralGenesSH::setOffsprgFstat);
      add("Nbr of Alleles per Locus - global (offsprg)","off.allnb",OFFSPRG,0,0,&TTNeutralGenesSH::getNbAllGlobal,0,0,0);
      add("Nbr of Fixed Loci per Patch  (offsprg)","off.fixlocp",OFFSPRG,0,0,&TTNeutralGenesSH::getFixLocLocal,0,0,0);
      add("Nbr of Fixed Loci in the Pop (offsprg)","off.fixloc",OFFSPRG,0,0,&TTNeutralGenesSH::getFixLocGlobal,0,0,0);
      add("Ho     (offsprg)","off.ho",OFFSPRG,0,0,&TTNeutralGenesSH::getHo,0,0,0);
      add("Hs-Nei (offsprg)","off.hsnei",OFFSPRG,0,0,&TTNeutralGenesSH::getHsnei,0,0,0);
      add("Ht-Nei (offsprg)","off.htnei",OFFSPRG,0,0,&TTNeutralGenesSH::getHtnei,0,0,0);
      add("Fis    (offsprg)","off.fis",OFFSPRG,0,0,&TTNeutralGenesSH::getFis,0,0,0);
      add("Fst    (offsprg)","off.fst",OFFSPRG,0,0,&TTNeutralGenesSH::getFst,0,0,0);
      add("Fit    (offsprg)","off.fit",OFFSPRG,0,0,&TTNeutralGenesSH::getFit,0,0,0);
  }

  if(AGE & ADULTS) {
      add("Nbr of Alleles per Locus - local  (adult)","adlt.allnbp",ADULTS,0,0,
        &TTNeutralGenesSH::getNbAllLocal,0,0,&TTNeutralGenesSH::setAdultsFstat);
      add("Nbr of Alleles per Locus - global (adult)","adlt.allnb",ADULTS,0,0,&TTNeutralGenesSH::getNbAllGlobal,0,0,0);
      add("Nbr of Fixed Loci per Patch  (adult)","adlt.fixlocp",ADULTS,0,0,&TTNeutralGenesSH::getFixLocLocal,0,0,0);
      add("Nbr of Fixed Loci in the Pop (adult)","adlt.fixloc",ADULTS,0,0,&TTNeutralGenesSH::getFixLocGlobal,0,0,0);
      add("Ho       (adult)","adlt.ho",ADULTS,0,0,&TTNeutralGenesSH::getHo,0,0,0);
      add("Hs-Nei   (adult)","adlt.hsnei",ADULTS,0,0,&TTNeutralGenesSH::getHsnei,0,0,0);
      add("Ht-Nei   (adult)","adlt.htnei",ADULTS,0,0,&TTNeutralGenesSH::getHtnei,0,0,0);
      add("Fis      (adult)","adlt.fis",ADULTS,0,0,&TTNeutralGenesSH::getFis,0,0,0);
      add("Fst      (adult)","adlt.fst",ADULTS,0,0,&TTNeutralGenesSH::getFst,0,0,0);
      add("Fit      (adult)","adlt.fit",ADULTS,0,0,&TTNeutralGenesSH::getFit,0,0,0);
  }
}

// ----------------------------------------------------------------------------------------
// setFstatRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstat2Recorders (age_t AGE)
{
  if(AGE & OFFSPRG) {
      add("Ho     (offsprg)","off.ho2",OFFSPRG,0,0,&TTNeutralGenesSH::getHo,0,0,
        &TTNeutralGenesSH::setOffsprgFstat2);
      add("Hs-Nei (offsprg)","off.hsnei2",OFFSPRG,0,0,&TTNeutralGenesSH::getHsnei,0,0,0);
      add("Ht-Nei (offsprg)","off.htnei2",OFFSPRG,0,0,&TTNeutralGenesSH::getHtnei,0,0,0);
      add("Fst    (offsprg)","off.fst2",OFFSPRG,0,0,&TTNeutralGenesSH::getFst,0,0,0);
  }

  if(AGE & ADULTS) {
      add("Ho       (adult)","adlt.ho2",ADULTS,0,0,&TTNeutralGenesSH::getHo,0,0,
        &TTNeutralGenesSH::setAdultsFstat2);
      add("Hs-Nei   (adult)","adlt.hsnei2",ADULTS,0,0,&TTNeutralGenesSH::getHsnei,0,0,0);
      add("Ht-Nei   (adult)","adlt.htnei2",ADULTS,0,0,&TTNeutralGenesSH::getHtnei,0,0,0);
      add("Fst      (adult)","adlt.fst2",ADULTS,0,0,&TTNeutralGenesSH::getFst,0,0,0);
  }
}

// ----------------------------------------------------------------------------------------
// setFstatWCRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstatWCRecorders (age_t AGE)
{
  if(AGE & OFFSPRG) {

      add("Fis Weir & Cockerham","off.fis.WC",OFFSPRG,0,0,&TTNeutralGenesSH::getFisWC,0,0,&TTNeutralGenesSH::setOffspringFstatWeirCockerham);
      add("Fst Weir & Cockerham","off.fst.WC",OFFSPRG,0,0,&TTNeutralGenesSH::getFstWC,0,0,0);
      add("Fit Weir & Cockerham","off.fit.WC",OFFSPRG,0,0,&TTNeutralGenesSH::getFitWC,0,0,0);
  }

  if(AGE & ADULTS) {

      add("Fis Weir & Cockerham","adlt.fis.WC",ADULTS,0,0,&TTNeutralGenesSH::getFisWC,0,0,&TTNeutralGenesSH::setAdultsFstatWeirCockerham);
      add("Fst Weir & Cockerham","adlt.fst.WC",ADULTS,0,0,&TTNeutralGenesSH::getFstWC,0,0,0);
      add("Fit Weir & Cockerham","adlt.fit.WC",ADULTS,0,0,&TTNeutralGenesSH::getFitWC,0,0,0);

  }
}
// ----------------------------------------------------------------------------------------
// setFstMatrixRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setFstMatrixRecorders (age_t AGE, unsigned char dim)
{
  std::ostringstream name, sub_name;

  void (TTNeutralGenesSH::* setter) () = (AGE == ADULTS ?
      ( dim & 3 ? &TTNeutralGenesSH::setAdultsFstMatrix :
        dim & 2 ? &TTNeutralGenesSH::setAdultsFstBetween :
            &TTNeutralGenesSH::setAdultsFstWithin ) :
            ( dim & 3 ? &TTNeutralGenesSH::setOffsprgFstMatrix :
              dim & 2 ? &TTNeutralGenesSH::setOffsprgFstBetween :
                  &TTNeutralGenesSH::setOffsprgFstWithin) );

  const char *prefix = (AGE == ADULTS ? "adlt." : "off.");

  unsigned int nbpatch =  _pop->getPatchNbr();
  unsigned int scale = (unsigned int)pow(10.0, (int)log10((float)nbpatch) + 1);

  name<<"Weir&Hill weighted Fst ("<<prefix<<")";
  sub_name<< prefix << "fst.WH";
  add(name.str(), sub_name.str(),  AGE, 0, 0, &TTNeutralGenesSH::getWeightedFst, 0, 0, setter);
  name.str("");
  sub_name.str("");

  if(dim & 1) {
      for(unsigned int i = 0; i < nbpatch; ++i) {
        name<<"Weighted Fst "<<i+1<<"."<<i+1;
        sub_name<< prefix << "fst" << i+1 << "." << i+1;
        add(name.str(), sub_name.str(),  AGE, i*scale + i, 0, 0, &TTNeutralGenesSH::getFst_ij, 0,  0);
        name.str("");
        sub_name.str("");
      }
  }
  if(dim & 2){
      for(unsigned int i = 0; i < nbpatch; ++i) {
        for(unsigned int j = i+1; j < nbpatch; ++j) {
            name<<"Weighted Fst "<<i+1<< "." <<j+1;
            sub_name<< prefix << "fst" << i+1 << "." << j+1;
            add(name.str(), sub_name.str(),  AGE, i*scale + j, 0, 0, &TTNeutralGenesSH::getFst_ij, 0, 0);
            name.str("");
            sub_name.str("");
        }
      }
  }

}
// ----------------------------------------------------------------------------------------
// setNeiGeneticDistanceRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setNeiGeneticDistanceRecorders (age_t AGE, bool pairwise)
{
  string name, sub_name;
  unsigned int nbpatch =  _pop->getPatchNbr();
  string prefix = (AGE == ADULTS ? "adlt." : "off.");
  unsigned int scale = (unsigned int)pow(10.0, (int)log10((float)nbpatch) + 1);

  void (TTNeutralGenesSH::* setter) () = (AGE == ADULTS ?
      &TTNeutralGenesSH::setAdltNeiGeneticDistance
      : &TTNeutralGenesSH::setOffsprgNeiGeneticDistance );

  sub_name = prefix + "D";
  add("Average between pop Nei's D", sub_name,  AGE, 0, 0,
      &TTNeutralGenesSH::getMeanNeiGeneticDistance, 0, 0, setter);

  if(pairwise) {

      for(unsigned int i = 0; i < nbpatch -1; i++){
        for(unsigned int j = i+1; j < nbpatch; j++) {
            name = "Nei's D between pop" + tstring::int2str(i+1) + " and pop" + tstring::int2str(j+1);
            sub_name = prefix + "D.p" + tstring::int2str(i+1) + "p" + tstring::int2str(j+1);
            add(name, sub_name,  AGE, i*scale + j, 0, 0, &TTNeutralGenesSH::getNeiGeneticDistance, 0, 0);
        }
      }
  }
}

// ----------------------------------------------------------------------------------------
// setFreqRecorders
// ----------------------------------------------------------------------------------------
void TTNeutralGenesSH::setDxyRecorders (age_t AGE, bool patchwise)
{
  string prefix = (AGE == OFFSPRG ? "off." : "adlt.");
  age_idx age = (AGE == OFFSPRG ? OFFSx : ADLTx );
  string name = "Sequence divergence - Dxy";
  string sub_name = "Dxy";


  if (!patchwise) {

      add(name, prefix + sub_name, AGE, (unsigned int)age, 0, 0, &TTNeutralGenesSH::getDxy, 0, 0);

  } else {

      for (unsigned int p1 = 0; p1 < _pop->getPatchNbr(); ++p1) {
        for (unsigned int p2 = p1 + 1; p2 < _pop->getPatchNbr(); ++p2) {
            add(name, prefix + sub_name + ".p" + tstring::int2str(p1+1) + "p" + tstring::int2str(p2+1),
              AGE, p1, p2,
              0, 0, (AGE == OFFSPRG ? &TTNeutralGenesSH::getDxyOffspringPerPatch :
                  &TTNeutralGenesSH::getDxyAdultPerPatch),
                  0);
        }
      }

  }
}

