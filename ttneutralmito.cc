/**  $Id: ttneutralmito.cc,v 1.9 2016-03-03 15:22:12 fred Exp $
*
*  @file ttneutralgenes.cc
*  Nemo2
*
*   Copyright (C) 2012 Frederic Guillaume
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
*  Created on @date 04.04.2012
*  @author fred
*/

#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "ttneutralmito.h"
#include "Uniform.h"
#include "output.h"
#include "metapop.h"

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoNeutralMito::TProtoNeutralMito( ) :_allele_num(0), _locus_num(0), _mut_rate(0),
  _mut_model(0), _init_model(0), _mutate_func_ptr(0), _writer(0), _type("mito")
{
  set_paramset("neutralmito", false, this);
 
  add_parameter("mito_loci",           INT, true, false,0, 0);
  add_parameter("mito_all",            INT, true, true, 1, 256, 0);
  add_parameter("mito_init_model",     INT, false,true, 0, 2, 0);
  add_parameter("mito_mutation_rate",  DBL, true, true, 0, 1, 0);
  add_parameter("mito_mutation_model", INT, true, true, 0, 2, 0);
//  add_parameter("mito_recombination_rate",DBL, false,true, 0, 0.5, 0);
  
  add_parameter("mito_save_genotype",BOOL,false,false,0,0);
//  add_parameter("mito_save_fsti",BOOL,false,false,0,0);
//  add_parameter("mito_save_freq",STR,false,false,0,0);
  add_parameter("mito_output_dir",STR,false,false,0,0);
  add_parameter("mito_output_logtime",INT,false,false,0,0);
}
// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoNeutralMito::TProtoNeutralMito(const TProtoNeutralMito& T) 
: _allele_num(T._allele_num), _locus_num(T._locus_num), _mut_rate(T._mut_rate),
_mut_model(T._mut_model), _init_model(T._init_model), 
_mutate_func_ptr(T._mutate_func_ptr), _writer(0), _type("mito")
{
  _paramSet = new ParamSet( *(T._paramSet) ) ;
}  
// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TProtoNeutralMito::~TProtoNeutralMito ()
{
//  if(_stats != NULL)      {delete _stats; _stats = NULL;}
  if( _writer)    delete _writer;
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TProtoNeutralMito::setParameters ()
{
  _locus_num =  (unsigned int)get_parameter_value("mito_loci");
  _allele_num = (unsigned int)get_parameter_value("mito_all");
  _mut_rate =  get_parameter_value("mito_mutation_rate");
  _mut_model = (unsigned int)get_parameter_value("mito_mutation_model");
  
  if(get_parameter("mito_init_model")->isSet())
    _init_model = (unsigned short)get_parameter_value("mito_init_model");
  else
    _init_model = 1; //maximum variance
  
  switch(_mut_model) {
    case 0:
    {
      _mutate_func_ptr = &TTNeutralMito::mutate_NULL;
      break;
    }
    case 1:
    {
      _mutate_func_ptr = &TTNeutralMito::mutate_SSM;
      break;
    }
    case 2:
    {
      _mutate_func_ptr = &TTNeutralMito::mutate_KAM;
      break;
    }
    default:
    {
      error("wrong parameter value for parameter \"mito_mutation_model\", max is 2\n");
      break; //should return false
    }
  }

  //special case:
  if(_allele_num==2) _mutate_func_ptr = &TTNeutralMito::mutate_2all;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TProtoNeutralMito::loadFileServices  (FileServices* loader)
{ 
  Param* param = get_parameter("mito_output_logtime");
  int logtime = (param->isSet() ? (int)param->getValue() : 0);
    
  if(_writer != 0) {
    delete _writer;
    _writer = NULL;
  }
//  // --- THE READER ---
//  //always add the reader:
//  writer = new TTNeutralMitoFH(this);
//  //set to read mode:
//  writer->set_isInputHandler(true);
//  //attach to file manager:
//  loader->attach_reader(writer);
//  //add to list of FileHandlers, will be deleted upon destruction
//  _writers.push_back( writer );
//
//  // --- WRITERS ---
  if( get_parameter("mito_save_genotype")->isSet() ) {

    _writer = new TTNeutralMitoFH(this);

    //           rpl_per, gen_per, rpl_occ, gen_occ, rank, path, self-ref
    _writer->set(true, (logtime != 0), 1, logtime, 0, get_parameter("mito_output_dir")->getArg(), this);
    
    _writer->set_write_fct( &TTNeutralMitoFH::write_FSTAT );
    
    loader->attach(_writer);

//    _writers.push_back( writer );
  }
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoNeutralMito::loadStatServices  (StatServices* loader)
{
  //allocate the stat handler
//  if(_stats == NULL)
//    _stats = new TTNeutralMitoSH(this);
//  
//  if(_stats != NULL) {
//    loader->attach(_stats);
//  }
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTNeutralMito* TProtoNeutralMito::hatch ()
{
  TTNeutralMito* new_trait = new TTNeutralMito();
  
  new_trait->set_locus_num(_locus_num);
  new_trait->set_allele_num(_allele_num);
  new_trait->set_init_model(_init_model);
  new_trait->set_mut_model(_mut_model);
  new_trait->set_mut_rate(_mut_rate);
  new_trait->set_mut_func_ptr(_mutate_func_ptr);
  new_trait->set_mean_mut_num((double)_locus_num*_mut_rate);
  
  return new_trait;
}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
bool TProtoNeutralMito::retrieve_data ( BinaryStorageBuffer* reader )
{
  unsigned int dummy;
  reader->read(&dummy, sizeof(int));
  if(dummy != _locus_num ){
    error("TProtoNeutralMito::retrieve_data:nb locus in file differ from parameter value!\n");
    _locus_num = dummy;
  }
  return true;
}
// ----------------------------------------------------------------------------------------

//                            N E U T R A L mtD N A  T R A I T

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTNeutralMito& TTNeutralMito::operator= (const TTrait& T)
{
  const TTNeutralMito& TN = dynamic_cast<const TTNeutralMito&> (T);
  
  if(this != &TN) {
    _locus_num = TN._locus_num;
    _allele_num = TN._allele_num;
//    _init_model = TN._init_model;
//    _mut_rate = TN._mut_rate;
//    _mut_model = TN._mut_model;
//    _mutate_func_ptr = TN._mutate_func_ptr;
//    _mean_mut_num = TN._mean_mut_num;
   
    reset();
    
    init();
    
    for(unsigned int i = 0; i < _locus_num; ++i) 
      _sequence[i] = TN._sequence[i];
  }
  
  return *this;
}
//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTNeutralMito::operator== (const TTrait& T)
{
  if(_type.compare(T.get_type()) != 0) return false;
  const TTNeutralMito& TN = dynamic_cast<const TTNeutralMito&> (T);
  
  if(this != &TN) {
    if(_locus_num != TN._locus_num) return false;
  }
  return true;
}
//----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTNeutralMito::operator!= (const TTrait& T)
{
  if(!((*this) == T))
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTNeutralMito::~TTNeutralMito()
{
  if(_sequence != NULL) delete [] _sequence;
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTNeutralMito::init ()
{
  if(_sequence != NULL)
    fatal("TTNeutralMito::init::_sequence is not NULL !\n");

  _sequence = new unsigned char[_locus_num];
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
void TTNeutralMito::init_sequence ()
{
  if(_sequence == NULL) _sequence = new unsigned char[_locus_num];
  
  if(_init_model == 1)
    for(unsigned int j = 0; j < _locus_num; j++)
        _sequence[j] = (unsigned char)RAND::Uniform(_allele_num);
  else
    for(unsigned int j = 0; j < _locus_num; j++)
      _sequence[j] = 0;

}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void TTNeutralMito::reset()
{
  if(_sequence != NULL) {

    delete [] _sequence;

    _sequence = NULL;
  }
}
// ----------------------------------------------------------------------------------------
// set_sequence
// ----------------------------------------------------------------------------------------
void TTNeutralMito::set_sequence(void** seq)
{
//  reset(); init();
//  memcpy(_sequence[0],seq[0],_locus_num);
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
void TTNeutralMito::inherit (TTrait* mother, TTrait* father)
{
  unsigned char* mother_seq = dynamic_cast<TTNeutralMito*> (mother)->getHaploSequence();

  memcpy(_sequence, mother_seq, _locus_num);
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTNeutralMito::show_up ()
{
  message("\n  Trait's type: mito\n\
       locus: %i\n\
     alleles: %i\n\
    _sequence:",_locus_num, _allele_num);
  
  for(unsigned int i = 0; (i < _locus_num && i < 10); i++)
    message("\n              %i ",(int)_sequence[i]);
  if(_locus_num > 10) message("\n              ...");
  message("\n");
}
// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void TTNeutralMito::mutate_SSM ()
{
  unsigned int mutLocus;
  bool direction;
  unsigned int NbMut; // = gsl_ran_binomial(RAND::r, _mut_rate, _locus_num);

  for(NbMut = (unsigned int)RAND::Poisson(_mean_mut_num) ; NbMut != 0; NbMut--) {
//  for( ; NbMut > 0; NbMut--) {
    mutLocus = RAND::Uniform( _locus_num);
    direction = RAND::RandBool();
    //alleles values are from 0 to NtrlAll - 1 !!!
    if(direction && _sequence[mutLocus] < _allele_num-1)
      _sequence[mutLocus] += 1; //one step to the right
    else if(_sequence[mutLocus] > 0) // !direction || all==_allele_num
      _sequence[mutLocus] -= 1; //one step to the left
    else //!direction && all == 0
      _sequence[mutLocus] += 1;
  }
}
// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void TTNeutralMito::mutate_KAM ()
{
  unsigned int mutLocus;
  unsigned char mut;

  unsigned int NbMut; //= gsl_ran_binomial(RAND::r, _mut_rate, _locus_num);
  
  for(NbMut = RAND::Poisson(_mean_mut_num) ; NbMut != 0; NbMut--) {
//  for( ; NbMut > 0; NbMut--) {
    mutLocus = RAND::Uniform(_locus_num);
    //assign an arbytrary allele value:
    do{
      mut = (unsigned char) (RAND::Uniform(_allele_num));
    } while (mut == _sequence[mutLocus]);
    _sequence[mutLocus] = mut;
  }
}
// ----------------------------------------------------------------------------------------
// mutate_2all
// ----------------------------------------------------------------------------------------
inline void TTNeutralMito::mutate_2all ()
{
  register unsigned int mutLocus;

  unsigned int NbMut = RAND::Binomial(_mut_rate, _locus_num);
  
  for( ; NbMut > 0; NbMut--) {
    mutLocus = RAND::Uniform(_locus_num);
    _sequence[mutLocus] = !_sequence[mutLocus];
  }
}

// ------------------------------------------------------------------------------

//                             FSTAT file writer/

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTNeutralMitoFH::write_FSTAT ()
{
  if(!_pop->isAlive()) return;
  /**The file format is FSTAT-like, with age class and sex added after the pop id. The file
   * extension is ".dat".*/
  unsigned int position;
  unsigned int nb_all   = _FHLinkedTrait->get_allele_num();
  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
  unsigned int patchNbr = _pop->getPatchNbr();
  unsigned char* seq;
  Patch* current_patch;
  Individual *ind;
  
  position = nb_all > 100 ? 3 : 2; //assumes nb_all not sup. to 999
  
  std::string filename = get_path() + 
                          this->get_service()->getGenerationReplicateFileName() + get_extension();
  
#ifdef _DEBUG_
  message("TTNeutralMitoFH::FHwrite (%s)\n",filename.c_str());
#endif
  
  ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open FSTAT output file!!\n");
  
  FILE<<patchNbr<<" "<<nb_locus + 4<<" "<<nb_all<<" "<<position<<"\n";
  
  for (unsigned int i = 0; i < nb_locus; ++i)
    FILE<<"loc"<<i+1<<"\n";
  //add names for the three last fields:
  FILE<<"age\n"<<"sex\n"<<"ped\n"<<"origin\n";
  
  for (unsigned int i = 0; i < patchNbr; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    for (unsigned int j = 0; j < current_patch->size(FEM, OFFSx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(FEM, OFFSx, j);
      seq = dynamic_cast<TTNeutralMito*> ( ind->getTrait(_FHLinkedTraitIndex) )->getHaploSequence();
      
      for(unsigned int k = 0; k < nb_locus; ++k) {
        FILE.fill('0');
        FILE.width(position);
        FILE<<(unsigned int)(seq[k]+1)<<" ";
      }
      FILE << OFFSPRG << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
    }
    
    for (unsigned int j = 0; j < current_patch->size(MAL, OFFSx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(MAL, OFFSx, j);
      seq = dynamic_cast<TTNeutralMito*> ( ind->getTrait(_FHLinkedTraitIndex) )->getHaploSequence();
      
      for(unsigned int k = 0; k < nb_locus; ++k) {
        FILE.fill('0');
        FILE.width(position);
        FILE<<(unsigned int)(seq[k]+1)<<" ";
      }
      FILE << OFFSPRG << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
    }
    
    for (unsigned int j = 0; j < current_patch->size(FEM, ADLTx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(FEM, ADLTx, j);
      seq = dynamic_cast<TTNeutralMito*> ( ind->getTrait(_FHLinkedTraitIndex) )->getHaploSequence();
      
      for(unsigned int k = 0; k < nb_locus; ++k) {
        FILE.fill('0');
        FILE.width(position);
        FILE<<(unsigned int)(seq[k]+1)<<" ";
      }
      FILE << ADULTS << " "<< FEM << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
    }
    
    for (unsigned int j = 0; j < current_patch->size(MAL, ADLTx); ++j) {
      
      FILE<<i+1<<" ";
      ind = current_patch->get(MAL, ADLTx, j);
      seq = dynamic_cast<TTNeutralMito*> ( ind->getTrait(_FHLinkedTraitIndex) )->getHaploSequence();
      
      for(unsigned int k = 0; k < nb_locus; ++k) {
        FILE.fill('0');
        FILE.width(position);
        FILE<<(unsigned int)(seq[k]+1)<<" ";
      }
      FILE << ADULTS << " "<< MAL << " " << ind->getPedigreeClass() << " " << ind->getHome()+1<<endl;
    }
    
  }
  
  FILE.close();
  
}
// ----------------------------------------------------------------------------------------
// FHread
// ----------------------------------------------------------------------------------------
//void TTNeutralMitoFH::FHread (string& filename)
//{
//  unsigned int digit, nloci_infile, nall_infile, npatch_infile;
////  unsigned int ploidy   = _FHLinkedTrait->get_ploidy();
//  unsigned int nb_all   = _FHLinkedTrait->get_allele_num();
//  unsigned int nb_locus = _FHLinkedTrait->get_locus_num();
//  unsigned int patchNbr = _pop->getPatchNbr();
//  
//  bool is_extended = false;
//
//  ifstream FILE(filename.c_str(),ios::in);
//  
//  if(!FILE) fatal("could not open FSTAT input file \"%s\"\n", filename.c_str());
//  
//  FILE >> npatch_infile >> nloci_infile >> nall_infile >> digit;
//  
//  if(npatch_infile != patchNbr) fatal("# of patch in FSTAT file differs from simulation settings\n");
//  if(nloci_infile != nb_locus)
//    if(nloci_infile > nb_locus + 4 || nloci_infile < nb_locus)
//      fatal("# of loci in FSTAT file differs from simulation settings\n");
//    else if (nloci_infile == nb_locus + 4) is_extended = true;
//    else  is_extended = false;
//  
//  if(nall_infile != nb_all)  fatal("# of alleles in FSTAT file differs from simulation settings\n");
//  
//  digit = (int) pow(10.0,(double)digit);
//  string loc_name;
//  
//  for(unsigned int i = 0; i < nloci_infile; ++i) {
//    FILE >> loc_name;
//  }
//  
//  
//  unsigned int pop, genot, all0, all1, age, sex, ped, origin;
//  age_idx agex;
//  Individual *ind;
//  unsigned char* seq = new unsigned char[nb_locus];
//  int lnbr = nloci_infile +2;
//  
//  while(FILE>>pop) {
//
//    for(unsigned int i = 0; i < nb_locus; ++i) {
//      FILE>>genot;
//      
//      all0 = genot/digit;
//      all1 = genot%digit;
//      
//      if(all0 <= nb_all)
//        seq[i] = all0 - 1;
//      else {
//        error("in FSTAT input file at line %i, locus %i : \
//              first allele value %d is greater than the max value specified (%i)!\n", lnbr, i+1, all0, nb_all);
//        fatal("Please check the input file.\n");
//      }
//      
//    }
//    
//    if(is_extended) FILE >> age >> sex >> ped >> origin;
//    else {
//      age = ADULTS;
//      sex = FEM;
//      ped = 0;
//      origin = 0;
//    }
//    
//    agex = (age == ADULTS ? ADLTx : OFFSx);
//    
//    ind = _pop->makeNewIndividual(0, 0, static_cast<sex_t> (sex), origin - 1);
//    ind->setPedigreeClass((unsigned char)ped);
//    ind->getTrait(_FHLinkedTraitIndex)->set_sequence((void**)seq);
//    _pop->getPatch(pop-1)->add(static_cast<sex_t> (sex), agex, ind);
//    
//    lnbr++;
//    
//  }
//
//  FILE.close();  
//  
//  delete [] seq;
//}
//  
// ------------------------------------------------------------------------------

//                             StatHandler/

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::init()
//{
//  StatHandler<TTNeutralMitoSH>::init();
//
//  allocateTables(_SHLinkedTrait->get_locus_num(),_SHLinkedTrait->get_allele_num());
//}
// ----------------------------------------------------------------------------------------
// allocateTable
// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::allocateTables (unsigned int loci, unsigned int all)
//{
//  unsigned int nb_patch = _pop->getPatchNbr();
//  unsigned int **sizes;
//  
//  sizes = new unsigned int * [nb_patch];
//  
//  for(unsigned int i = 0; i < nb_patch; ++i) {
//    sizes[i] = new unsigned int [loci];
//    for(unsigned int j = 0; j < loci; ++j)
//      sizes[i][j] = all;
//  }
//
//  _alleleCountTable.allocate(nb_patch, loci, sizes);
//  
//  _alleleFreqTable.allocate(nb_patch, loci, sizes);
//  
//  _globalAlleleFreq.reset(loci, all);
//
//  for(unsigned int i = 0; i < nb_patch; ++i)
//    delete [] sizes[i];
//  delete [] sizes;
//}
// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
//bool TTNeutralMitoSH::setStatRecorders(std::string& token)
//{
//#ifdef _DEBUG_
//  message("-TTNeutralMitoSH::setStatRecorders ");
//#endif
//  if(token.compare("coa") == 0) {
//	
//	add("Wtn Patch Coancestry (offsprg)","off.theta",FLAT,OFFSPRG,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setOffsprgCoaMatrix);
//	add("Btn Patch Coancestry (offsprg)","off.alpha",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getMeanAlpha,0,0,0);
//    
//	add("Wtn Patch Coancestry   (adult)","adlt.theta",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setAdultsCoaMatrix);
//	add("Btn Patch Coancestry   (adult)","adlt.alpha",FLAT,ADULTS,0,&TTNeutralMitoSH::getMeanAlpha,0,0,0);
//    
//  } else if(token.compare("adlt.coa") == 0) {
//    
//	add("Wtn Patch Coancestry   (adult)","adlt.theta",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setAdultsCoaMatrix);
//	add("Btn Patch Coancestry   (adult)","adlt.alpha",FLAT,ADULTS,0,&TTNeutralMitoSH::getMeanAlpha,0,0,0);
//    
//  } else if(token.compare("off.coa") == 0) {
//    
//    add("Wtn Patch Coancestry (offsprg)","off.theta",FLAT,OFFSPRG,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setOffsprgCoaMatrix);
//	add("Btn Patch Coancestry (offsprg)","off.alpha",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getMeanAlpha,0,0,0);    
//   	
//  } else if(token.compare("adlt.coa.persex") == 0) {
//	
//	add("Female Theta      (adult)","adlt.theta",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getMeanTheta, 0, 0, &TTNeutralMitoSH::setAdults_Theta);
//	add("Female Theta      (adult)","adlt.thetaFF",FLAT,ADULTS,0,&TTNeutralMitoSH::getTheta_FF,0,0,0);
//	add("Male Theta        (adult)","adlt.thetaMM",FLAT,ADULTS,0,&TTNeutralMitoSH::getTheta_MM,0,0,0);
//	add("Female-Male Theta (adult)","adlt.thetaFM",FLAT,ADULTS,0,&TTNeutralMitoSH::getTheta_FM,0,0,0);
//    
//  } else if(token.compare("adlt.coa.within") == 0) {
//    
//	add("Wtn Patch Coancestry   (adult)","adlt.theta",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setAdultsCoaWithin);
//    
//  } else if(token.compare("off.coa.within") == 0) {
//    
//    add("Wtn Patch Coancestry (offsprg)","off.theta",FLAT,OFFSPRG,0,
//        &TTNeutralMitoSH::getMeanTheta,0,0,&TTNeutralMitoSH::setOffsprgCoaWithin);   
//    
//  } else if(token.compare("adlt.coa.between") == 0) {
//    
//	add("Btn Patch Coancestry   (adult)","adlt.alpha",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getMeanAlpha,0,0,&TTNeutralMitoSH::setAdultsCoaBetween);
//    
//  } else if(token.compare("off.coa.between") == 0) {
//    
//    add("Btn Patch Coancestry (offsprg)","off.alpha",FLAT,OFFSPRG,0,
//        &TTNeutralMitoSH::getMeanAlpha,0,0,&TTNeutralMitoSH::setOffsprgCoaBetween);   
//    
//  } else if(token.compare("coa.matrix") == 0) {
//    
//    setCoaMatrixRecorders(OFFSPRG, 3);
//    setCoaMatrixRecorders(ADULTS, 3);
//    
//  } else if(token.compare("off.coa.matrix") == 0) {
//        
//    setCoaMatrixRecorders(OFFSPRG, 3);
//    
//  } else if(token.compare("adlt.coa.matrix") == 0) {
//    
//    setCoaMatrixRecorders(ADULTS, 3);
//    
//  } else if(token.compare("coa.matrix.within") == 0) {
//    
//    setCoaMatrixRecorders(OFFSPRG, 1);
//    setCoaMatrixRecorders(ADULTS, 1);
//    
//  } else if(token.compare("off.coa.matrix.within") == 0) {
//    
//    setCoaMatrixRecorders(OFFSPRG, 1);
//    
//  } else if(token.compare("adlt.coa.matrix.within") == 0) {
//    
//    setCoaMatrixRecorders(ADULTS, 1);
//    
//  } else if(token.compare("sibcoa") == 0) {
//
//    add("Proportion of full-sib offspring","prop.fsib",FLAT,OFFSPRG,0,
//        0,0,&TTNeutralMitoSH::getSibProportions,&TTNeutralMitoSH::setSibStats);
//    add("Proportion of paternal half-sib ","prop.phsib",FLAT,OFFSPRG,1,0,0,&TTNeutralMitoSH::getSibProportions,0);
//    add("Proportion of maternal half-sib ","prop.mhsib",FLAT,OFFSPRG,2,0,0, &TTNeutralMitoSH::getSibProportions,0);
//    add("Proportion of non-sib offspring ","prop.nsib",FLAT,OFFSPRG,3,0,0,&TTNeutralMitoSH::getSibProportions,0);
//    add("Coancestry of full-sib offspring","coa.fsib",FLAT,OFFSPRG,0,0,0,&TTNeutralMitoSH::getSibCoaMeans,0);
//    add("Coancestry of paternal half-sib ","coa.phsib",FLAT,OFFSPRG,1,0,0,&TTNeutralMitoSH::getSibCoaMeans,0);
//    add("Coancestry of maternal half-sib ","coa.mhsib",FLAT,OFFSPRG,2,0,0,&TTNeutralMitoSH::getSibCoaMeans,0);
//    add("Coancestry of non-sib offspring ","coa.nsib",FLAT,OFFSPRG,3,0,0,&TTNeutralMitoSH::getSibCoaMeans,0);
//    
//  } else if(token.compare("offsprgfstat") == 0 || token.compare("off.fstat") == 0) {
//	
//    setFstatRecorders(OFFSPRG);
//    
//  } else if(token.compare("adultfstat") == 0 || token.compare("adlt.fstat") == 0) {
//
//    setFstatRecorders(ADULTS);
//    
//  } else if(token.compare("fstat") == 0) {
//    
//    setFstatRecorders(ALL);
//    
//  } else if(token.compare("fstWC") == 0 || token.compare("fstatWC") == 0) {
//    
//    setFstatWCRecorders(ALL);
//     
//  } else if(token.compare("off.fstWC") == 0 || token.compare("off.fstatWC") == 0) {
//    
//    setFstatWCRecorders(OFFSPRG);
//    
//  } else if(token.compare("adlt.fstWC") == 0 || token.compare("adlt.fstatWC") == 0) {
//    
//    setFstatWCRecorders(ADULTS);
//    
//  } else if(token.compare("weighted.fst") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 0);
//    setFstMatrixRecorders(ADULTS, 0);
//    
//  } else if(token.compare("off.weighted.fst") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 0);
//    
//  } else if(token.compare("adlt.weighted.fst") == 0) {
//    
//    setFstMatrixRecorders(ADULTS, 0);
//    
//  } else if(token.compare("weighted.fst.matrix") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 3);
//    setFstMatrixRecorders(ADULTS, 3);
//    
//  } else if(token.compare("off.weighted.fst.matrix") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 3);
//    
//  } else if(token.compare("adlt.weighted.fst.matrix") == 0) {
//    
//    setFstMatrixRecorders(ADULTS, 3);
//    
//  } else if(token.compare("weighted.fst.within") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 1);
//    setFstMatrixRecorders(ADULTS, 1);
//    
//  } else if(token.compare("off.weighted.fst.within") == 0) {
//    
//    setFstMatrixRecorders(OFFSPRG, 1);
//    
//  } else if(token.compare("adlt.weighted.fst.within") == 0) {
//    
//    setFstMatrixRecorders(ADULTS, 1);
//    
//  } else if(token.compare("adlt.NeiDistance") == 0) {
//
//    setNeiGeneticDistanceRecorders(ADULTS, true);
//    
//  } else if(token.compare("off.NeiDistance") == 0) {
//    
//    setNeiGeneticDistanceRecorders(OFFSPRG, true);
//    
//  } else if(token.compare("NeiDistance") == 0) {
//    
//    setNeiGeneticDistanceRecorders(ADULTS, true);
//    setNeiGeneticDistanceRecorders(OFFSPRG, true);
//    
//  } else if(token.compare("mean.NeiDistance") == 0) {
//   
//    setNeiGeneticDistanceRecorders(ADULTS, false);
//    setNeiGeneticDistanceRecorders(OFFSPRG, false);
//    
//  } else if(token.compare("adlt.mean.NeiDistance") == 0) {
//    
//    setNeiGeneticDistanceRecorders(ADULTS, false);
//    
//  } else if(token.compare("off.mean.NeiDistance") == 0) {
//    
//    setNeiGeneticDistanceRecorders(OFFSPRG, false);
//        
//  } else
//	return false;
//  
//  return true;
//}
//// ----------------------------------------------------------------------------------------
//// setCoaMatrixRecorders
//// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::setCoaMatrixRecorders (age_t AGE, unsigned char dim)
//{
//  std::ostringstream name, sub_name;
//  
//  void (TTNeutralMitoSH::* setter) () = (AGE == ADULTS ?
//                                          ( dim & 3 ? &TTNeutralMitoSH::setAdultsCoaMatrix : 
//                                            dim & 2 ? &TTNeutralMitoSH::setAdultsCoaBetween : 
//                                            &TTNeutralMitoSH::setAdultsCoaWithin ) :
//                                          ( dim & 3 ? &TTNeutralMitoSH::setOffsprgCoaMatrix :
//                                            dim & 2 ? &TTNeutralMitoSH::setOffsprgCoaBetween :
//                                            &TTNeutralMitoSH::setOffsprgCoaWithin) );
//  
//  const char *prefix = (AGE == ADULTS ? "adlt." : "off.");
//  
//  unsigned int nbpatch =  _pop->getPatchNbr();
//  unsigned int scale = (unsigned int)pow(10.0, (int)log10((float)nbpatch) + 1);
//  
//  if(dim & 1) {
//    name<<"Wtn Patch Coancestry ("<<prefix<<")";
//    sub_name<< prefix << "theta";
//    add(name.str(), sub_name.str(), FLAT, AGE, 0, &TTNeutralMitoSH::getMeanTheta, 0, 0, setter);
//    name.str("");
//    sub_name.str("");
//  }
//  if(dim & 2){
//    name<<"Btn Patch Coancestry ("<<prefix<<")";
//    sub_name<< prefix << "alpha";
//    add(name.str(), sub_name.str(), FLAT, AGE, 0, &TTNeutralMitoSH::getMeanAlpha, 0, 0, (!(dim&1) ? setter : 0) );
//    name.str("");
//    sub_name.str("");
//  }  
//  if(dim & 1) {
//    for(unsigned int i = 0; i < nbpatch; ++i) {
//      name<<"coancestry "<<i+1<<"."<<i+1;
//      sub_name<< prefix << "coa" << i+1 << "." << i+1;
//      add(name.str(), sub_name.str(), FLAT, AGE, i*scale + i, 0, 0, &TTNeutralMitoSH::getCoa, 0);
//      name.str("");
//      sub_name.str("");
//    }
//  }
//  if(dim & 2){
//    for(unsigned int i = 0; i < nbpatch; ++i) {
//      for(unsigned int j = i+1; j < nbpatch; ++j) {
//        name<<"coancestry "<<i+1<< "." <<j+1;
//        sub_name<< prefix << "coa" << i+1 << "." << j+1;
//        add(name.str(), sub_name.str(), FLAT, AGE, i*scale + j, 0, 0, &TTNeutralMitoSH::getCoa,0);
//        name.str("");
//        sub_name.str("");
//      }
//    }
//  }
//  
//}
//// ----------------------------------------------------------------------------------------
//// setFstatRecorders
//// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::setFstatRecorders (age_t AGE)
//{
//  if(AGE & OFFSPRG) {
//    add("Nbr of Alleles per Locus - local  (offsprg)","off.allnbp",FLAT,OFFSPRG,0,
//        &TTNeutralMitoSH::getNbAllLocal,0,0,&TTNeutralMitoSH::setOffsprgFstat);
//    add("Nbr of Alleles per Locus - global (offsprg)","off.allnb",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getNbAllGlobal,0,0,0);
//    add("Nbr of Fixed Loci per Patch  (offsprg)","off.fixlocp",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFixLocLocal,0,0,0);
//    add("Nbr of Fixed Loci in the Pop (offsprg)","off.fixloc",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFixLocGlobal,0,0,0);
//    add("Ho     (offsprg)","off.ho",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getHo,0,0,0);
//    add("Hs-Nei (offsprg)","off.hsnei",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getHsnei,0,0,0);
//    add("Ht-Nei (offsprg)","off.htnei",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getHtnei,0,0,0);
//    add("Fis    (offsprg)","off.fis",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFis,0,0,0);
//    add("Fst    (offsprg)","off.fst",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFst,0,0,0);
//    add("Fit    (offsprg)","off.fit",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFit,0,0,0);
//  }
//  
//  if(AGE & ADULTS) {
//    add("Nbr of Alleles per Locus - local  (adult)","adlt.allnbp",FLAT,ADULTS,0,
//        &TTNeutralMitoSH::getNbAllLocal,0,0,&TTNeutralMitoSH::setAdultsFstat);
//    add("Nbr of Alleles per Locus - global (adult)","adlt.allnb",FLAT,ADULTS,0,&TTNeutralMitoSH::getNbAllGlobal,0,0,0);
//    add("Nbr of Fixed Loci per Patch  (adult)","adlt.fixlocp",FLAT,ADULTS,0,&TTNeutralMitoSH::getFixLocLocal,0,0,0);
//    add("Nbr of Fixed Loci in the Pop (adult)","adlt.fixloc",FLAT,ADULTS,0,&TTNeutralMitoSH::getFixLocGlobal,0,0,0);
//    add("Ho       (adult)","adlt.ho",FLAT,ADULTS,0,&TTNeutralMitoSH::getHo,0,0,0);
//    add("Hs-Nei   (adult)","adlt.hsnei",FLAT,ADULTS,0,&TTNeutralMitoSH::getHsnei,0,0,0);
//    add("Ht-Nei   (adult)","adlt.htnei",FLAT,ADULTS,0,&TTNeutralMitoSH::getHtnei,0,0,0);
//    add("Fis      (adult)","adlt.fis",FLAT,ADULTS,0,&TTNeutralMitoSH::getFis,0,0,0);
//    add("Fst      (adult)","adlt.fst",FLAT,ADULTS,0,&TTNeutralMitoSH::getFst,0,0,0);
//    add("Fit      (adult)","adlt.fit",FLAT,ADULTS,0,&TTNeutralMitoSH::getFit,0,0,0);
//  }
//}
//// ----------------------------------------------------------------------------------------
//// setFstatWCRecorders
//// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::setFstatWCRecorders (age_t AGE)
//{
//  if(AGE & OFFSPRG) {
//    
//    add("Fst Weir & Cockerham","off.fst.WC",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFstWC,0,0,&TTNeutralMitoSH::setOffspringFstatWeirCockerham);
////    add("Fis Weir & Cockerham","off.fis.WC",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFisWC,0,0,0);
////    add("Fit Weir & Cockerham","off.fit.WC",FLAT,OFFSPRG,0,&TTNeutralMitoSH::getFitWC,0,0,0);
//  }
//  
//  if(AGE & ADULTS) {
//    
//    add("Fst Weir & Cockerham","adlt.fst.WC",FLAT,ADULTS,0,&TTNeutralMitoSH::getFstWC,0,0,&TTNeutralMitoSH::setAdultsFstatWeirCockerham);
////    add("Fis Weir & Cockerham","adlt.fis.WC",FLAT,ADULTS,0,&TTNeutralMitoSH::getFisWC,0,0,0);
////    add("Fit Weir & Cockerham","adlt.fit.WC",FLAT,ADULTS,0,&TTNeutralMitoSH::getFitWC,0,0,0);
//    
//  }
//}
//// ----------------------------------------------------------------------------------------
//// setFstMatrixRecorders
//// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::setFstMatrixRecorders (age_t AGE, unsigned char dim)
//{
//  std::ostringstream name, sub_name;
//  
//  void (TTNeutralMitoSH::* setter) () = (AGE == ADULTS ?
//                                          ( dim & 3 ? &TTNeutralMitoSH::setAdultsFstMatrix : 
//                                            dim & 2 ? &TTNeutralMitoSH::setAdultsFstBetween : 
//                                            &TTNeutralMitoSH::setAdultsFstWithin ) :
//                                          ( dim & 3 ? &TTNeutralMitoSH::setOffsprgFstMatrix :
//                                            dim & 2 ? &TTNeutralMitoSH::setOffsprgFstBetween :
//                                            &TTNeutralMitoSH::setOffsprgFstWithin) );
//  
//  const char *prefix = (AGE == ADULTS ? "adlt." : "off.");
//  
//  unsigned int nbpatch =  _pop->getPatchNbr();
//  unsigned int scale = (unsigned int)pow(10.0, (int)log10((float)nbpatch) + 1);
//  
//  name<<"Weir&Hill weighted Fst ("<<prefix<<")";
//  sub_name<< prefix << "fst.WH";
//  add(name.str(), sub_name.str(), FLAT, AGE, 0, &TTNeutralMitoSH::getWeightedFst, 0, 0, setter);
//  name.str("");
//  sub_name.str("");
//  
//  if(dim & 1) {
//    for(unsigned int i = 0; i < nbpatch; ++i) {
//      name<<"Weighted Fst "<<i+1<<"."<<i+1;
//      sub_name<< prefix << "fst" << i+1 << "." << i+1;
//      add(name.str(), sub_name.str(), FLAT, AGE, i*scale + i, 0, 0, &TTNeutralMitoSH::getFst_ij, 0);
//      name.str("");
//      sub_name.str("");
//    }
//  }
//  if(dim & 2){
//    for(unsigned int i = 0; i < nbpatch; ++i) {
//      for(unsigned int j = i+1; j < nbpatch; ++j) {
//        name<<"Weighted Fst "<<i+1<< "." <<j+1;
//        sub_name<< prefix << "fst" << i+1 << "." << j+1;
//        add(name.str(), sub_name.str(), FLAT, AGE, i*scale + j, 0, 0, &TTNeutralMitoSH::getFst_ij,0);
//        name.str("");
//        sub_name.str("");
//      }
//    }
//  }
//  
//}
//// ----------------------------------------------------------------------------------------
//// setNeiGeneticDistanceRecorders
//// ----------------------------------------------------------------------------------------
//void TTNeutralMitoSH::setNeiGeneticDistanceRecorders (age_t AGE, bool pairwise)
//{
//  std::ostringstream name, sub_name;
//  unsigned int nbpatch =  _pop->getPatchNbr();
//  const char *prefix = (AGE == ADULTS ? "adlt." : "off.");
//  
//  void (TTNeutralMitoSH::* setter) () = (AGE == ADULTS ?
//                                          &TTNeutralMitoSH::setAdltNeiGeneticDistance 
//                                          : &TTNeutralMitoSH::setOffsprgNeiGeneticDistance );
//
//  sub_name<< prefix << "D";
//  add("Average between pop Nei's D", sub_name.str(), FLAT, AGE, 0, 
//      &TTNeutralMitoSH::getMeanNeiGeneticDistance, 0, 0, setter);
//  sub_name.str("");
//  
//  if(pairwise) {
//
//    unsigned int c = 0;
//    
//    for(unsigned int i = 0; i < nbpatch -1; i++){
//      for(unsigned int j = i+1; j < nbpatch; j++) {
//        name<<"Nei's D between pop"<<i+1<<" and pop"<<j+1;
//        sub_name<< prefix << "D" << i+1 << "." << j+1;
//        add(name.str(), sub_name.str(), FLAT, AGE, c++, 0, 0, &TTNeutralMitoSH::getNeiGeneticDistance, 0);
//        name.str("");
//        sub_name.str("");
//      }
//    }
//  }
//}
