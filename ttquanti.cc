/** $Id: ttquanti.cc,v 1.32 2015-05-01 11:29:19 fred Exp $
 *
 *  @file ttquanti.cc
 *  Nemo2
 *
 *   Copyright (C) 2006-2011 Frederic Guillaume
 *   frederic.guillaume@env.ethz.ch
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
 *  created on @date 14.11.2005
 * 
 *  @author fred
 */

#include <sstream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include "ttquanti.h"
#include "filehandler.h"
#include "output.h"
#include "Uniform.h"
#include "tstring.h"
#include "utils.h"
#ifdef HAS_GSL
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#endif

void store_quanti_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                                sex_t SEX, age_idx AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                                unsigned int nTrait, unsigned int TraitIndex);


// ------------------------------------------------------------------------------

//                             TProtoQuanti

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::TProtoQuanti() :
_nb_locus(0),
_nb_traits(0),
_seq_length(0),
_allele_model(0),
_allele_value(0),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0),
_eval(0),
_effects_multivar(0),
_ws(0),
_genomic_mutation_rate(0),
_mutation_correlation(0),
_mutation_sigma(0),
_init_value(0),
_all_chooser(0),
_sizeofLocusType(sizeof(double)),
_eVariance(0),
_stats(0),
_writer(0),
_freqExtractor(0),
_pleio_matx(0)
{
  set_paramset("quantitative_trait", false, this);
  
  add_parameter("quanti_traits",INT,true,false,0,0);
  add_parameter("quanti_loci",INT,true,false,0,0);
  add_parameter("quanti_allele_model",STR,false,false,0,0);
  add_parameter("quanti_allele_value",DBL,false,false,0,0);
  add_parameter("quanti_init_value",MAT,false,false,0,0);
  add_parameter("quanti_init_model",INT,false,true,0,4);
  add_parameter("quanti_environmental_variance",DBL,false,false,0,0);
  add_parameter("quanti_pleio_matrix",INT,false,true,0,1); // allows user to set pleiotropic degree of loci
  
  add_parameter("quanti_mutation_rate",DBL,true,true,0,1, 0);
  add_parameter("quanti_mutation_variance",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_correlation",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_covariance",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_matrix",MAT,false,false,0,0, 0);
  
  //genetic map parameters:
  TTProtoWithMap::addGeneticMapParameters("quanti");
  
  add_parameter("quanti_output",STR,false,false,0,0);
  add_parameter("quanti_logtime",INT,false,false,0,0);
  add_parameter("quanti_dir",STR,false,false,0,0);
  
  add_parameter("quanti_extract_freq",BOOL,false,false,0,0);
  add_parameter("quanti_freq_logtime",INT,false,false,0,0);
  add_parameter("quanti_freq_grain",DBL,false,false,0,0);
}
// ----------------------------------------------------------------------------------------
// copy cstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::TProtoQuanti(const TProtoQuanti& T) : 
_nb_locus(T._nb_locus),
_nb_traits(T._nb_traits),
_seq_length(T._seq_length),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0),
_eval(0),
_effects_multivar(0),
_ws(0),
_genomic_mutation_rate(T._genomic_mutation_rate),
_mutation_correlation(T._mutation_correlation),
_mutation_sigma(0),
_init_value(0),
_allele_value(0),
_all_chooser(0),
_eVariance(0),
_stats(0),
_writer(0),
_freqExtractor(0)
{ 
  _locusByteSize = T._nb_traits * sizeof(double);
  _sizeofLocusType = sizeof(double);
  _paramSet = new ParamSet( *(T._paramSet) ) ;
  _pleio_matx = new TMatrix( *(T._pleio_matx) );
  vector< vector<unsigned int> > _trait_table;
  vector< vector<unsigned int> > _locus_table;
  //vector< vector<double> > _mut_matrix;
}
// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::~TProtoQuanti ()
{
  reset_mutation_pointers();
  
  if(_stats != NULL){delete _stats; _stats = NULL;}
  if(_writer != NULL){delete _writer; _writer = NULL;}
  if(_freqExtractor != NULL){delete _freqExtractor; _freqExtractor = NULL;}
  if(_all_chooser) {delete [] _all_chooser; _all_chooser = NULL;}
  if(_init_value) {delete[] _init_value;_init_value = NULL;}
}
// ----------------------------------------------------------------------------------------
// setParameters
// ----------------------------------------------------------------------------------------
bool TProtoQuanti::setParameters()
{
  //cout << "\nStart of TProtoQuanti::setParameters!\t";

  if(get_parameter("quanti_mutation_covariance")->isSet()) {
    fatal("\"quanti_mutation_covariance\" is deprecated, use \"quanti_mutation_correlation\" instead.\n");
    return false;
  }
  _nb_traits = (unsigned int)get_parameter_value("quanti_traits");
  _nb_locus = (unsigned int)get_parameter_value("quanti_loci");

  //total nbre of values for both traits in a haploid genome
  if(!get_parameter("quanti_pleio_matrix")->isSet()){
	  _seq_length = _nb_traits * _nb_locus; //mutations are pleiotropic if no quanti_pleio_matrix is provided!!!
  }
  _locusByteSize = _nb_traits * sizeof(double); // could move this assignment below pleio_matx reading in order to get correct size
  _sizeofLocusType = sizeof(double);
  
  if(get_parameter("quanti_environmental_variance")->isSet())
    _eVariance = sqrt(get_parameter_value("quanti_environmental_variance"));
  else
    _eVariance = 0;
  
  TMatrix tmp_matx;

  //---------------------------------------------------------------------------------------------
  //initial values:   
  if(_init_value != NULL) {
    delete [] _init_value;
    _init_value = NULL;
  }  
  
  if(get_parameter("quanti_init_value")->isSet()) {
    
    get_parameter("quanti_init_value")->getMatrix(&tmp_matx);
    
    if(tmp_matx.getNbRows() != 1 || tmp_matx.getNbCols() != _nb_traits) {
      fatal("\"quanti_init_value\" must be a vector of length equal to the number of traits!\n");
      return false;
    }
    
    _init_value = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++)
      _init_value[i] = tmp_matx.get(0,i);
    
  }
  else {
    
    _init_value = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++)
      _init_value[i] = 0.0;
  }
  
  if(get_parameter("quanti_init_model")->isSet()) {
    
    _doInitMutation = (unsigned int)get_parameter_value("quanti_init_model");
    
  } else {
    _doInitMutation = 1;
  }
  
  // initializes pleiotropic connectivity matrix -- "quanti_pleio_matx" in init file
  if(!get_parameter("quanti_pleio_matrix")->isSet()) { // if user provides NO pleiotropy matrix
	  _pleio_matx = new TMatrix(_nb_locus, _nb_traits);
	  _pleio_matx->assign(1);

  }	else{
	  if(_pleio_matx == NULL){ _pleio_matx = new TMatrix(); }
	  get_parameter("quanti_pleio_matrix")->getMatrix(_pleio_matx);
	  if(_pleio_matx->getNbRows() != _nb_locus || _pleio_matx->getNbCols() != _nb_traits){
		  return error("\"quanti_pleio-matrix\" must be an m by n matrix of length m equal to the number of loci and n equal to the number of loci!\n");
	  }
	  else{ //return error("\"quanti_pleio_matrix\" IS an m by n matrix. Good work!!!\n");
	  }
  }
  // print _pleio_matrix
  //printf("\n_pleio_matx dimensions: rows(loci) = %i, columns(traits) = %i\n",_pleio_matx->getNbRows(),_pleio_matx->getNbCols());
  //  for(unsigned int i = 0; i < _pleio_matx->getNbRows(); i++) {
  //	  for(unsigned int j = 0; j < _pleio_matx->getNbCols(); j++){
  //		  printf("%.0f ",_pleio_matx->get(i,j));
  //	  	  printf(" | ");
  //	  }
  //	  printf("\n");
  //}

  // create temp_table from _pleiomatx to build _trait_table and _locus_table
  vector< vector<unsigned int> > temp_table(_pleio_matx->getNbCols(),vector<unsigned int>(_pleio_matx->getNbRows()));
  unsigned int pos=1;
  for(unsigned int i=0; i<_pleio_matx->getNbRows(); ++i){
	  for(unsigned int j=0; j<_pleio_matx->getNbCols(); ++j){
		  if(_pleio_matx->get(i,j)==1){
			  temp_table[j][i]=pos;
			  pos++;
		  }
		  else if(_pleio_matx->get(i,j)==0){
			  temp_table[j][i]=0;
		  }
		  else return error("\"quanti_pleio_matrix\" must only contain 0s and 1s \n");
	  }
  }
  // print temp_table (for debugging purposes only)
  /*cout << "temp_table: \n";
  for(unsigned int i=0; i<temp_table.size(); ++i){
	  for(unsigned int j=0; j<temp_table[i].size(); ++j){
		  cout << temp_table[i][j] << " | ";
	  }
  cout << endl;
  }*/

  // build _trait_table from temp_table
  if(_trait_table.size() != 0){
	  _trait_table.clear();
  }
  for(unsigned int i=0; i<temp_table.size(); ++i){
	  vector<unsigned int> temp;
	  for(unsigned int j=0; j<temp_table[i].size(); ++j){
		  if(temp_table[i][j]!=0){
			  temp.push_back(temp_table[i][j]-1);
		  }
	  }
	  _trait_table.push_back(temp);
  }
  // print _trait_table (for debugging purposes only
  /*cout << "trait_table: \n";
  for(unsigned int i=0; i<_trait_table.size(); ++i){
	  for(unsigned int j=0; j<_trait_table[i].size(); ++j){
		  cout << _trait_table[i][j] << " | ";
	  }
	  cout << endl;
  }*/
  //cout << "\n_trait_size: " << _trait_table.size() << endl;

  // build _locus_table from _pleio_matx
  if(_locus_table.size() != 0){
	  _locus_table.clear();
  }
  unsigned int pos4 = 0;
  unsigned int tsize = 0;
  for(unsigned int i=0; i<_pleio_matx->getNbRows(); ++i){
  	vector<unsigned int> temp2;
  	for(unsigned int j=0; j<_pleio_matx->getNbCols(); ++j){
  		if(_pleio_matx->get(i,j)==1){
  			tsize++;
  			pos4++;
  		}
  	}
  	temp2.push_back(pos4-tsize);
  	temp2.push_back(tsize);
  	_locus_table.push_back(temp2);
  	tsize=0;
  }
  // print _locus_table (for debugging purposes only)
  //  cout << "locus_table: \n";
  //for(unsigned int i=0; i<_locus_table.size(); ++i){
  //	  for(unsigned int j=0; j<_locus_table[i].size(); ++j){
  //		  cout << _locus_table[i][j] << " | ";
  //	  }
  //	  cout << endl;
  //}
  //cout << "_locus_size: " << _locus_table.size() << endl;

  // calculate _seq_length from _trait_table
  _seq_length = 0;
  for(unsigned int i=0; i<_trait_table.size(); ++i){
	  for(unsigned int j=0; j<_trait_table[i].size(); ++j){
		  _seq_length++;
	  }
  }
  //cout << "_seq_length: " << _seq_length << endl;
  //return error("\" pleiotropic matrix debug check \" LEGIT");
  //_seq_length = _nb_traits * _nb_locus; // changing back to full pleiotropy for debug check


  //---------------------------------------------------------------------------------------
  if( !setMutationParameters() ) return false; //sets _allele_model
  //bypass the genetic map for free recombination to save significant time:
  //this is done only if the parameter is explicitely set to 0.5 (and not {{0.5}})
  //else, we set the map:
  if ( get_parameter("quanti_recombination_rate")->isSet() && 
      get_parameter_value("quanti_recombination_rate") == 0.5 ) {
  
    _recombRate = 0.5;

    //for use in the inherit_free function:
    if(_all_chooser) delete [] _all_chooser; 
    _all_chooser = new bool[_nb_locus];
    
  } else {
    
    if( !setGeneticMapParameters("quanti") ) return false;

    if(_all_chooser) delete [] _all_chooser; 
    _all_chooser = new bool[_nb_locus]; //may still be used if no map parameters are specified
  
  }

  
  //---------------------------------------------------------------------------------------
  _mutation_func_ptrs.clear();

  for(unsigned int i=0; i < _nb_locus; ++i){
	  if(_locus_table[i][1] == 1) {
    
		  if (_allele_model == 1 || _allele_model == 2) {
			  _mutation_func_ptrs.push_back(&TProtoQuanti::getMutationEffectUnivariateDiallelic);
		  } else {
			  _mutation_func_ptrs.push_back(&TProtoQuanti::getMutationEffectUnivariateGaussian);
		  }
    
	  } else if (_locus_table[i][1] == 2) {

		  if (_allele_model == 1 || _allele_model == 2) {
			  _mutation_func_ptrs.push_back(&TProtoQuanti::getMutationEffectBivariateDiallelic);
		  } else {
			  _mutation_func_ptrs.push_back(&TProtoQuanti::getMutationEffectMultivariateGaussian);
		  }
    
	  } else if (_locus_table[i][1] > 2) {
    
		  if (_allele_model > 2) {
			  _mutation_func_ptrs.push_back(&TProtoQuanti::getMutationEffectMultivariateGaussian);
		  } else {
			  fatal("in \"quanti\" trait, the di-allelic model is only allowed for max. 2 quantitative traits.");
		  }
    
	  }
  }
  //cout << "\nEnd of TProtoQuanti::setParameters!\t";
  return true; 
}
// ----------------------------------------------------------------------------------------
// get_trait_table
// ----------------------------------------------------------------------------------------
vector< vector<unsigned int> >& TProtoQuanti::get_trait_table()
{
	//cout << "\n _trait_table size: " << _trait_table.size() << endl;
	return _trait_table;
}
// ----------------------------------------------------------------------------------------
// get_locus_table
// ----------------------------------------------------------------------------------------
vector< vector<unsigned int> >& TProtoQuanti::get_locus_table()
{
	//cout << "\n _locus_size: " << _locus_table.size() << endl;
	return _locus_table;
}
// ----------------------------------------------------------------------------------------
// setMutationParameters
// ----------------------------------------------------------------------------------------
bool TProtoQuanti::setMutationParameters ()
{
  //cout << "\nStart of TProtoQuanti::setMutationParameters!\n";

  _genomic_mutation_rate = get_parameter_value("quanti_mutation_rate") * 2 * _nb_locus;
  _mutation_correlation = new double [_nb_locus];
  if(get_parameter("quanti_mutation_correlation")->isSet()){
	  for(unsigned int loc = 0; loc < _nb_locus; loc++){
		  _mutation_correlation[loc] = get_parameter_value("quanti_mutation_correlation");
		  //cout << "_mutation_correlation for locus " << loc+1 << ": " << _mutation_correlation[loc] << endl;
	  }
  } else{
	  for(unsigned int loc = 0; loc < _nb_locus; loc++){
		  _mutation_correlation[loc] = 0;
	  }
  }
  reset_mutation_pointers();
  
  //checking allelic model
  if (get_parameter("quanti_allele_model")->isSet()) {
    
    string model = get_parameter("quanti_allele_model")->getArg();
    if (model == "diallelic") {
      
      _allele_model = 1;
      
      return setDiallelicMutationModel ();
      
    } else if (model == "diallelic_HC") {
      
      _allele_model = 2;
      
      return setDiallelicMutationModel ();
      
    } else if (model == "continuous") {
      
      _allele_model = 3;
      
      return setContinuousMutationModel ();
      
    } else if (model == "continuous_HC") {
      
      _allele_model = 4;
      
      return setContinuousMutationModel ();
      
    } else {
      error("\"quanti_allele_model\" has options \"diallelic[_HC]\" or \"continuous[_HC]\" only. \n");
      return false;
    }
    
  } 
  else { //default model
    _allele_model = 3;
    return setContinuousMutationModel ();
  }
  //cout << "\nEnd of TProtoQuanti::setMutationParameters!\t";

  return true;
}
//---------------------------------------------------------------------------------------------
// setDiallelicMutationModel
//---------------------------------------------------------------------------------------------
bool TProtoQuanti::setDiallelicMutationModel ()
{
  //cout << "\nStart of TProtoQuanti::setDiallelicMutationModel!\t";

  if (!get_parameter("quanti_allele_value")->isSet()) {
    error("in \"quanti\" trait, \"quanti_allele_value\" is not set for the diallelic model.\n");
    return false;
    
  } else {
    
    assert(_allele_value == NULL); 
    //should be the case as reset_mutation_pointers has been called in setMutationParameters
    
    
    _allele_value = new double* [_nb_locus];
    
    for (unsigned i = 0; i < _nb_locus; ++i) _allele_value[i] = new double [2];
    
    if (get_parameter("quanti_allele_value")->isMatrix()) { //locus-specific allelic values
      
      TMatrix tmp;
      
      get_parameter("quanti_allele_value")->getMatrix(&tmp);
      
      if (tmp.ncols() != _nb_locus) {
        fatal("\"quanti_allele_value\" must get a matrix with num. columns = num. loci.\n");
        return false;
      }
      
      if (tmp.nrows() == 1) {
        
        for (unsigned i = 0; i < _nb_locus; ++i) {
          _allele_value[i][0] = tmp.get(0, i);
          _allele_value[i][1] = -_allele_value[i][0];
        }
        
      } else if (tmp.nrows() == 2) {
        
        for (unsigned i = 0; i < _nb_locus; ++i) {
          _allele_value[i][0] = tmp.get(0, i);
          _allele_value[i][1] = tmp.get(1, i);
        }
        
      } else {
        fatal("\"quanti_allele_value\" must get a matrix with a max. of 2 rows (and num. columns = num. loci).\n");
        return false;
      }
      
      
    } else { //no locus-specific allelic values
      
      double val = get_parameter_value("quanti_allele_value");
      for (unsigned i = 0; i < _nb_locus; ++i) {
        _allele_value[i][0] = val;
        _allele_value[i][1] = -_allele_value[i][0];
      }
    }
  }
  //cout << "\nEnd of TProtoQuanti::setDiallelicMutationModel!\t";
  return true;
}
//---------------------------------------------------------------------------------------------
// setContinuousMutationModel
//---------------------------------------------------------------------------------------------
bool TProtoQuanti::setContinuousMutationModel ()
{
  //cout << "\nStart of TProtoQuanti::setContinuousMutationModel!\t";

  //_mutation_matrix = new TMatrix;
  _gsl_mutation_matrix = new gsl_matrix* [_nb_locus]();
  _eval = new gsl_vector* [_nb_locus];
  _evect = new gsl_matrix* [_nb_locus];
  _effects_multivar = new gsl_vector* [_nb_locus];
  _ws = new gsl_vector* [_nb_locus];
  _mutation_sigma = new double* [_nb_locus];
  vector< vector<double> > _mut_matrix;

  //assert(_mutation_sigma == NULL);
  if(_mut_matrix.size() != 0){
	  _mut_matrix.clear();
  }
  for(unsigned int loc = 0; loc < _nb_locus; loc++){
	  unsigned int pleio_deg = _locus_table[loc][1];
	  _mutation_sigma[loc] = new double [pleio_deg];

	  //setting the mutation variance-covariance matrix
	  if(get_parameter("quanti_mutation_variance")->isSet()) {
    
		  if(get_parameter("quanti_mutation_matrix")->isSet()) {
			  warning("both \"quanti_mutation_variance\" and \"quanti_mutation_matrix\" are set, using the matrix only!\n");
		  } else {
      
			  //_mutation_sigma[loc] = new double [pleio_deg];
      
			  double sigma = sqrt( get_parameter_value("quanti_mutation_variance") );
      
			  for(unsigned int i = 0; i < pleio_deg; i++)
				  _mutation_sigma[loc][i] = sigma;
      
			  if(pleio_deg > 1) {
				  //setting the mutation matrix
				  _gsl_mutation_matrix[loc] = gsl_matrix_alloc(pleio_deg, pleio_deg);
        
				  double covar, var = get_parameter_value("quanti_mutation_variance");

				  covar = _mutation_correlation[loc] * var;
        
				  for(unsigned int i = 0; i < pleio_deg; i++)
					  gsl_matrix_set(_gsl_mutation_matrix[loc], i, i, var);
        
				  for(unsigned int i = 0; i < pleio_deg - 1; i++)
					  for(unsigned int j = i + 1 ; j < pleio_deg; j++) {
						  gsl_matrix_set(_gsl_mutation_matrix[loc], i, j, covar);
						  gsl_matrix_set(_gsl_mutation_matrix[loc], j, i, covar);
					  }

//				  cout << "_gsl_mutation_matrix for locus " << loc+1 << " \n";
//				  for(unsigned int i=0; i<pleio_deg; ++i){
//					  for(unsigned int j=0; j<pleio_deg; ++j){
//						  cout << gsl_matrix_get(_gsl_mutation_matrix[loc],i,j) << " | ";
//				      }
//				      cout << endl;
//				  }

				  _evect[loc] = gsl_matrix_alloc(pleio_deg, pleio_deg);
				  _eval[loc] = gsl_vector_alloc(pleio_deg);

				  set_mutation_matrix_decomposition(loc,pleio_deg);

			  	  if(pleio_deg > 1) {
		             _effects_multivar[loc] = gsl_vector_alloc(pleio_deg);
					 _ws[loc] = gsl_vector_alloc(pleio_deg);
			      }
			  }
				  //_mutation_matrix->set_from_gsl_matrix(_gsl_mutation_matrix[loc]);
				  //_mutation_matrix->get_dims(&dims[0]);
				  //cout << "_mutation_matrix dims: " << dims[0] << " | " << dims[1] << endl;

				  // TODO create a _mut_matrix from var and covar as below? Not sure if needed.
				  /*unsigned int pos = 0;
				  for(unsigned int i=0; i<_nb_traits-1; i++){
					  _mut_matrix[loc][pos+i]=var;
					  for(unsigned int j=0; j<_nb_traits; j++){
						  _mut_matrix[loc][pos+i+j+1]=covar;
					  }
					  pos += _nb_traits;
				  }
				  _mut_matrix[loc][_nb_traits*_nb_traits-1]=var;*/
		  }
	  } // END if(get_parameter("quanti_mutation_variance")->isSet())
  
	  if(get_parameter("quanti_mutation_matrix")->isSet()) {

		  if(_mut_matrix.size() != 0){
			 _mut_matrix.clear();
	  	  }
		  get_parameter("quanti_mutation_matrix")->getVariableMatrix(&_mut_matrix);

		  // Print out mut_matrix for debugging purposes
//		  cout << "_mut_matrix for locus " << loc+1 << "\n";
//   		  for(unsigned int j=0; j<_mut_matrix[loc].size(); ++j){
//   			  cout << _mut_matrix[loc][j] << " | ";
//   		  }
//   		  cout << endl;

		  if( _mut_matrix.size() != _nb_locus || _mut_matrix[loc].size() != (pleio_deg*pleio_deg)) {
			  error("\"quanti_mutation_matrix\" must be a matrix of size = \"quanti_loci\" BY (\"quanti_traits\" x \"quanti_traits\")\n");
		  	  return false;
		  }

		  //_mutation_matrix[loc].get_gsl_matrix(_gsl_mutation_matrix);
		  _gsl_mutation_matrix[loc] = gsl_matrix_alloc(pleio_deg, pleio_deg);
		  unsigned int mmpos2 = 0;
		  for(unsigned long int i=0; i<pleio_deg; i++){
			  for(unsigned long int j=0; j<pleio_deg; j++){
				  gsl_matrix_set(_gsl_mutation_matrix[loc],i,j,_mut_matrix[loc][mmpos2]);
				  mmpos2++;
			  }
		  }
//		  cout << "\n_gsl_mutation_matrix for locus " << loc+1 << " \n";
//		  for(unsigned int i=0; i<pleio_deg; ++i){
//			  for(unsigned int j=0; j<pleio_deg; ++j){
//				  cout << gsl_matrix_get(_gsl_mutation_matrix[loc],i,j) << " | ";
//		      }
//		      cout << endl;
//		  }

		  _mutation_sigma[loc] = new double [pleio_deg];
    
		  if(pleio_deg == 1){
			  _mutation_sigma[loc][0] = sqrt(_mut_matrix[loc][0]);
			  //cout << "mutation_sigma :" << _mutation_sigma[loc][0] << endl;
		  }
//		  else if(pleio_deg == 2){
//			  _mutation_sigma[loc][0] = sqrt(_mut_matrix[loc][0]);
//			  _mutation_sigma[loc][1] = sqrt(_mut_matrix[loc][3]);
//			  //cout << "_mutation_sigma 1 and 2: " << _mutation_sigma[loc][0] << " | " <<_mutation_sigma[loc][1] << endl;
//			  _mutation_correlation[loc] = _mut_matrix[loc][1] / sqrt(_mut_matrix[loc][0]*_mut_matrix[loc][3]);
//			  //cout << "Mutation correlation for locus " << loc+1 << ": " << _mutation_correlation[loc] << "\n";
//
//			  _evect[loc] = gsl_matrix_alloc(pleio_deg, pleio_deg);
//			  _eval[loc] = gsl_vector_alloc(pleio_deg);
//
//			  set_mutation_matrix_decomposition(loc,pleio_deg);
//		  }
		  else if(pleio_deg > 1){
			  for(unsigned int i = 0; i < pleio_deg; i++){
				  _mutation_sigma[loc][i] = sqrt(gsl_matrix_get(_gsl_mutation_matrix[loc],i,i));
			  } // Would _mutation_sigma ever be needed for > 2 traits?
			  _mutation_correlation[loc] = _mut_matrix[loc][1] / sqrt(_mut_matrix[loc][0]*_mut_matrix[loc][pleio_deg*pleio_deg-1]);
			  //cout << "Mutation correlation for locus " << loc+1 << ": " << _mutation_correlation[loc] << "\n";

			  _evect[loc] = gsl_matrix_alloc(pleio_deg, pleio_deg);
			  _eval[loc] = gsl_vector_alloc(pleio_deg);

			  set_mutation_matrix_decomposition(loc,pleio_deg);

			  _effects_multivar[loc] = gsl_vector_alloc(pleio_deg);
			 _ws[loc] = gsl_vector_alloc(pleio_deg);
		  }
	  } // END if(get_parameter("quanti_mutation_matrix")->isSet())
	  else if(!get_parameter("quanti_mutation_variance")->isSet()) {
		  error("\"quanti_mutation_matrix\" or \"quanti_mutation_variance\" must be specified!\n");
		  return false;
	  }

#ifdef _DEBUG_
    message("-- Mutation matrix:\n");
    _mutation_matrix->show_up();
    message("-- MMatrix decomposition:\n");
    for(unsigned int i = 0; i < _nb_traits; i++)
      cout<<gsl_vector_get(_eval[loc],i)<<" ";
    cout<<endl;
    if(_nb_traits == 2) message("-- mutation correlation: %f\n",_mutation_correlation[loc]);
#endif
//    cout << "Mutation sigmas for locus " << loc+1 << ":  ";
//    for(unsigned int i = 0; i < pleio_deg; i++){
//   		cout <<	_mutation_sigma[loc][i] << " | ";
//   	} cout << "\n";
  }
  //cout << "\nEnd of TProtoQuanti::setContinuousMutationModel!\t";
  return true;
}
// ----------------------------------------------------------------------------------------
// set_mutation_matrix_decomposition
// ----------------------------------------------------------------------------------------
void TProtoQuanti::set_mutation_matrix_decomposition	(unsigned int loc, unsigned int pleio_deg)
{
  //  cout << "\nStart of TProtoQuanti::set_mutation_matrix_decomposition!\t";

	gsl_matrix *E = gsl_matrix_alloc(pleio_deg,pleio_deg);
	gsl_matrix_memcpy(E,_gsl_mutation_matrix[loc]);
  	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (pleio_deg);
  	gsl_eigen_symmv (E, _eval[loc], _evect[loc], w);
  	gsl_eigen_symmv_free (w);
  	gsl_matrix_free(E);
#ifdef _DEBUG_
  message("-- Mutation matrix eigenvalues:\n");
  for(unsigned int i = 0; i < _nb_traits; i++)
    cout<<gsl_vector_get(_eval[loc],i)<<" ";
  cout<<endl;
#endif
  	double eval;
  	//take square root of eigenvalues, will be used in Gaussian as stdev
  	for(unsigned int i = 0; i < pleio_deg; i++) {
  		eval = gsl_vector_get(_eval[loc],i);
  		eval = (eval < 0.000001 ? 0 : eval);
  		gsl_vector_set( _eval[loc], i, sqrt(eval) );
  	}

  	// DEBUG OUTPUT ------------------------------------------------------
/*  	cout << "\n_gsl_mutation_matrix eigen values: ";
  	for(unsigned int i = 0; i < pleio_deg; i++)
  		cout<<gsl_vector_get(_eval[loc],i)<<" ";
  	cout << "\n_gsl_mutation_matrix eigen vectors: \n";
  	for(unsigned int i = 0; i < pleio_deg; i++){
  		for(unsigned int j = 0; j < pleio_deg; j++){
  			cout<<gsl_matrix_get(_evect[loc],i,j)<<" | ";
  		}
  		cout<<endl;
  	}*/
  	// DEBUG OUTPUT ------------------------------------------------------
   //cout << "\nEnd of TProtoQuanti::set_mutation_matrix_decomposition!\t";
}
// ----------------------------------------------------------------------------------------
// getMutationEffectMultivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectMultivariateGaussian (unsigned int loc)
{
  RAND::MultivariateGaussian(_eval[loc], _evect[loc], _ws[loc], _effects_multivar[loc]);
  //cout << "_effects_multivar 1: " << _effects_multivar[loc] << " | " << endl;
  return _effects_multivar[loc]->data;
}
// ----------------------------------------------------------------------------------------
// getMutationEffectBivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectBivariateGaussian   (unsigned int loc)
{
  RAND::BivariateGaussian(_mutation_sigma[loc][0], _mutation_sigma[loc][1], _mutation_correlation[loc],
                          &_effects_bivar[0], &_effects_bivar[1]);
  //cout << "\n_effects_bivar 1 and 2: " << _effects_bivar[0] << " | " << _effects_bivar[1];
  return &_effects_bivar[0];
}

// ----------------------------------------------------------------------------------------
// getMutationEffectUnivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectUnivariateGaussian   (unsigned int loc)
{
  _effects_bivar[0] = RAND::Gaussian(_mutation_sigma[loc][0]);
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// getMutationEffectUnivariateDiallelic
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectUnivariateDiallelic   (unsigned int loc)
{
  _effects_bivar[0] = _allele_value[loc][ RAND::RandBool() ];
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// getMutationEffectBivariateDiallelic
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectBivariateDiallelic   (unsigned int loc)
{
  bool pos = RAND::RandBool();
  _effects_bivar[0] = _allele_value[loc][pos];
  //  _effects_bivar[1] = _allele_value[loc][pos]; // effects on both traits always the same -->> GWAS version!
  _effects_bivar[1] = _allele_value[loc][ (RAND::Uniform() < _mutation_correlation[loc] ?
                                           pos : RAND::RandBool()) ]; // effects on both traits can be different
  //  cout << "_effects_bivar 1 and 2: " << _effects_bivar[0] << " | " << _effects_bivar[1] << endl;
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// reset_mutation_pointers
// ----------------------------------------------------------------------------------------
void TProtoQuanti::reset_mutation_pointers()
{
  if(_allele_value) {
    for(unsigned int i = 0; i < _nb_locus; ++i)
      delete [] _allele_value[i];
    delete [] _allele_value;
    _allele_value = NULL;
  }
  
//	  if(_gsl_mutation_matrix != NULL)
//		  gsl_matrix_free(_gsl_mutation_matrix);

  if(_gsl_mutation_matrix != NULL){
	  delete [] _gsl_mutation_matrix;
	  _gsl_mutation_matrix = NULL;
  }

//  if(_mutation_correlation != NULL){
//    delete [] _mutation_correlation;
//    _mutation_correlation = NULL;
//  }

  if(_mutation_sigma != NULL){
    delete [] _mutation_sigma;
    _mutation_sigma = NULL;
  }
  if(_evect != NULL) {
    delete [] _evect;
    _evect = NULL;
  }
  
  if(_eval != NULL){
  	  delete [] _eval;
  	  _eval = NULL;
  }
  
  if(_effects_multivar != NULL) {
    delete []_effects_multivar;
    _effects_multivar = NULL;
  }
  
  if(_ws) {
    delete [] _ws;
    _ws = 0;
  }
  //cout << "reset_mutation_pointers " << endl;
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoQuanti::inherit_free (sex_t SEX, double* seq, double** parent)
{
  //cout << "\nStart of TProtoQuanti::inherit_free!\t";

  register unsigned int bloc;
  assert(_all_chooser);
  for(unsigned int i = 0; i < _nb_locus; ++i)
    _all_chooser[i] = RAND::RandBool();
  
  for(unsigned int i = 0; i < _nb_locus; ++i) {
    
    //bloc = i*_nb_traits;
    bloc = _locus_table[i][0];

//    memcpy(&seq[bloc], &parent[ _all_chooser[i] ][bloc], _locusByteSize);
    memcpy(&seq[bloc], &parent[ _all_chooser[i] ][bloc], _locus_table[i][1]*_sizeofLocusType);
  }
  //cout << "\nEnd of TProtoQuanti::inherit_free!\t";
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoQuanti::inherit_low (sex_t SEX, double* seq, double** parent)
{
  //cout << "\nStart of TProtoQuanti::inherit_low!\t";
  register unsigned int prevLoc = 0, chrm_bloc = 0, prev_bloc, cpy_bloc;
  register bool flipper;
  
  //the table containing the loci at which x-overs happen
  vector< unsigned int >& recTable = _map.getRecLoci(SEX, _mapIndex);
  
  //the table containing which homologous chromosome copy we start with, for each chromosome
  vector< bool > & firstRecPos = _map.getFirstRecPosition(SEX);
  
  //number of x-overs
  unsigned int nbRec = recTable.size();
  
  //  cout << "TProtoQuanti::inherit; sex="<<SEX<<"; nb Rec = "<<nbRec;//<<endl;
  
  // c is the chromosome number
  // stride is the number of loci considered so-far
  // rec is the number of x-over done so-far
  for(unsigned int c = 0, stride = 0, rec = 0; c < _numChromosome; ++c) {
    
    //the copy of the chromosome with which with start
    flipper = firstRecPos[c];
    
    //number of loci copied so-far
    chrm_bloc = stride + _numLociPerChrmsm[c];
    
    //last locus at which a x-over happened, will be first position on current chromosome
    prevLoc = stride;
    
    //    cout<<"chrm "<<c<<" side="<<firstRecPos[c]<<endl;
    
    // copy blocs of loci between x-over points on current chromosome
    // skip it if locus is not on this chromosome but a latter one
    for(; recTable[rec] < chrm_bloc && rec < nbRec; rec++) {
      
      // start position of the bloc of values to copy (_nb_traits values per locus)
      //prev_bloc = prevLoc * _nb_traits;
      prev_bloc = _locus_table[prevLoc][0]; // for non-universal pleiotropy
      
      // size of the bloc to copy
      //cpy_bloc = (recTable[rec] - prevLoc) * _nb_traits;
      cpy_bloc = 0;
      for (unsigned int i = prevLoc; i < recTable[rec]; ++i){ // for non-universal pleiotropy
    	  cpy_bloc += _locus_table[i][1];
      } // for each locus between current x-over and previous x-over sum their pleiotropic degrees
      //      cout<<"copy seq from "<<prevLoc<<"("<<prev_bloc<<") to "<<recTable[rec]
      //      <<"("<<(recTable[rec] - prevLoc)<<" loc) ("<<cpy_bloc*_locusByteSize<<"B) side "<<flipper<<endl;
      
      //memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);
      memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);

      //update the starting locus to the next recombination point
      prevLoc = recTable[rec];
      
      //switch side for next bloc to copy on this chromosome
      flipper = !flipper;
    }
    
    //prev_bloc = prevLoc * _nb_traits;
    prev_bloc = _locus_table[prevLoc][0]; // for non-universal pleiotropy

    //cpy_bloc = (chrm_bloc - prevLoc) * _nb_traits;
    cpy_bloc = 0;
    for (unsigned int i = prevLoc; i < chrm_bloc; ++i){ // for non-universal pleiotropy
  	  cpy_bloc += _locus_table[i][1];
    } // for each locus between prevLoc and end of chrmsm sum their pleiotropic degrees
    //    cout << "copy end of chrmsm from "<<prevLoc<<" to "<<chrm_bloc
    //         <<"("<<(chrm_bloc - prevLoc)<<" loc) on side "<<flipper<<endl;
    
    //copy what's left between the last x-over point and the end of the chrmsme
    //memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);
    memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);

    stride += _numLociPerChrmsm[c];
  }
  //cout << "\nEnd of TProtoQuanti::inherit_low!\t";
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTQuanti* TProtoQuanti::hatch()
{
  //cout << "\nStart of TProtoQuanti::hatch!\t";

  TTQuanti* kid = new TTQuanti();
  kid->set_proto(this);
  kid->set_nb_locus(_nb_locus);
  kid->set_nb_traits(_nb_traits);
  kid->set_seq_length(_seq_length);
  kid->set_genomic_mutation_rate(_genomic_mutation_rate);
  kid->set_init_value(_init_value, _doInitMutation);
  kid->set_mutation_fptr((_allele_model != 1 && _allele_model != 3));
  if(_recombRate == 0.5) //member of TTProtoWithMap
    kid->set_inherit_fptr(&TProtoQuanti::inherit_free);
  else
    kid->set_inherit_fptr(&TProtoQuanti::inherit_low);
  kid->set_eVariance(_eVariance);
  //cout << "\nEnd of TProtoQuanti::hatch!\t";

  return kid;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TProtoQuanti::loadFileServices  (FileServices* loader)
{ 
  //cout << "\nStart of TProtoQuanti::loadFileServices!\t";

  int logtime = 0;
  //writer
  if(get_parameter("quanti_output")->isSet()) {
    
    if(_writer == NULL) _writer = new TTQuantiFH(this);
    
    _writer->setOutputOption(get_parameter("quanti_output")->getArg());
    
    Param* param = get_parameter("quanti_logtime");
    
    if(param->isMatrix()) {
      
      TMatrix temp;
      param->getMatrix(&temp);
      _writer->set_multi(true, true, 1, &temp, get_parameter("quanti_dir")->getArg());
      
    } else   //  rpl_per, gen_per, rpl_occ, gen_occ, rank (0), path, self-ref      
      _writer->set(true, true, 1, (param->isSet() ? (int)param->getValue() : 0),
                   0, get_parameter("quanti_dir")->getArg(),this);
    
    loader->attach(_writer);
    
  } else if(_writer != NULL) {
    delete _writer;
    _writer = NULL;
  }
  
  //freq extractor
  if(get_parameter("quanti_extract_freq")->isSet()) {
    
    if(_freqExtractor == NULL) _freqExtractor = new TTQFreqExtractor(this);
    
    Param* param = get_parameter("quanti_freq_logtime");
    
    logtime = (param->isSet() ? (int)param->getValue() : 0);
    
    _freqExtractor->set(true,(logtime != 0),1,logtime,0,get_parameter("quanti_dir")->getArg(),this);
    
    param = get_parameter("quanti_freq_grain");
    
    if( ! param->isSet() ) {
      warning(" parameter \"quanti_freq_grain\" is not set, using 0.1\n");
      _freqExtractor->set_granularity(0.1);
    } else 
      _freqExtractor->set_granularity(param->getValue());
    
    loader->attach(_freqExtractor);
    
  } else if(_freqExtractor != NULL) {
    delete _freqExtractor;
    _freqExtractor = NULL;
  }
  //cout << "\nEnd of TProtoQuanti::loadFileServices!\t";
  
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoQuanti::loadStatServices  (StatServices* loader)
{
  //cout << "\nStart of TProtoQuanti::loadStatServices!\t";
  //allocate the stat handler
  if(_stats == NULL)
    _stats = new TTQuantiSH(this);
  
  if(_stats != NULL) {
    loader->attach(_stats);
  }
  //cout << "\nEnd of TProtoQuanti::loadStatServices!\t";

}
// ------------------------------------------------------------------------------

//                             TTQuanti

// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTQuanti& TTQuanti::operator= (const TTrait& T)
{  
  cout << "\nStart of TTQuanti::operator!\t";

  const TTQuanti& TQ = dynamic_cast<const TTQuanti&>(T);
  
  if(this != &TQ) {
    
    _nb_locus = TQ._nb_locus;
    _nb_traits = TQ._nb_traits;
    _seq_length = TQ._seq_length;
    reset();
    init();
    memcpy(_sequence[0],TQ._sequence[0],_seq_length*sizeof(double));
    //cout << "memcpy1\n";
    memcpy(_sequence[1],TQ._sequence[1],_seq_length*sizeof(double));
    //cout << "memcpy2\n";
    set_value();
  }
  cout << "\nEnd of TTQuanti::operator!\t";
  
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTQuanti::operator== (const TTrait& T)
{ 
  if(this->get_type().compare(T.get_type()) != 0) return false;
  
  const TTQuanti& TQ = dynamic_cast<const TTQuanti&>(T);
  
  if(this != &TQ) {
    if(_nb_locus != TQ._nb_locus) return false;
    if(_nb_traits != TQ._nb_traits) return false;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTQuanti::operator!= (const TTrait& T)
{
  if(!((*this) == T) )
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// set_init_value
// ----------------------------------------------------------------------------------------
void TTQuanti::set_init_value             (double* val, unsigned int doInit)
{
  //cout << "\nStart of TTQuanti::set_init_value(double,int)!\t";
  set_init_value(val);
  
  _doInitMutation = doInit;
  //cout << "\nEnd of TTQuanti::set_init_value!\t";

}
// ----------------------------------------------------------------------------------------
// set_init_value
// ----------------------------------------------------------------------------------------
void TTQuanti::set_init_value             (double* val) 
{ 
  //cout << "\nStart of TTQuanti::set_init_value!(double)\t";
  assert(_nb_traits != 0);
  
  if(_init_value) delete [] _init_value;
  
  _init_value = new double [_nb_traits];
  
  for(unsigned int i = 0; i < _nb_traits; ++i) 
    _init_value[i] = val[i];

  //cout << "\nEnd of TTQuanti::set_init_value!(double)\t";
}  
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
inline void TTQuanti::init ()
{
  //cout << "\nStart of TTQuanti::init!\t";
  _sequence = new double*[2];
  _sequence[0] = new double [_seq_length];
  _sequence[1] = new double [_seq_length];
  if(!_phenotypes) _phenotypes = new double [_nb_traits];
  
  if(!_init_value) {
    _init_value = new double [_nb_traits];
    for(unsigned int i = 0; i < _nb_traits; ++i) 
      _init_value[i] = _myProto->get_init_value(i);
  }
  //cout << "\nEnd of TTQuanti::init!\t";
  //cout << _seq_length << sizeof(_sequence[1]) << sizeof(_sequence[0]) << sizeof(_seq_length);
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
inline void TTQuanti::init_sequence ()
{
//  cout << "\nStart of TTQuanti::init_sequence!\t\n";
  //unsigned int pos;
  vector< vector<unsigned int> > ttable = _myProto->get_trait_table();
  vector< vector<unsigned int> > ltable = _myProto->get_locus_table();

  //options:
  //0: no variation, init value = (trait value)/(2*_nb_locus)
  //1: init value = (trait value)/(2*_nb_locus) + 1 mutation/locus
  //2: init value = (trait value)/(2*_nb_locus) + 1 mutation/locus+make parts ==> mean trait value doesn't change
  //3: init value = (trait value + random deviate N(0,sdVm))/(2*_nb_locus) + 1 mutation/locus+make parts
  //4: no variation and initialize with opposite effect alleles at alternating loci.
  
  //decide what initial value to use
  
  //Note: this is kept here as the initial values may have been set individually by LCE_quanti
  //      it wouldn't make sense then to store the init values in the prototype only
//  double my_init[_nb_traits]; // new for model 4?

//  cout << "\n_doInitMutation:\t" << _doInitMutation << "\n";
//  cout << "\n_allele_model:\t" << _myProto->_allele_model << "\n";
  if(_doInitMutation == 3) {
    
    double sdVm;
    
    for(unsigned int j = 0; j < _nb_traits; j++) {
      //      sdVm = sqrt(4*_nb_locus*_mut_rate*_myProto->get_trait_var(j)); //trait variance = 2Vm
      sdVm = 0.25;
//      my_init[j] = (_init_value[j] + RAND::Gaussian(sdVm)) / (2*_nb_locus);
//      _init_value[j] = (_init_value[j] + RAND::Gaussian(sdVm)) / (2*_nb_locus); // used before model 4 was implemented
      _init_value[j] = (_init_value[j] + RAND::Gaussian(sdVm)) / (2*ttable[j].size());
    }
    
  } else {
    
    for(unsigned int j = 0; j < _nb_traits; j++){
//    	cout << "\n\t_init_value before: " << _init_value[j];
//    	my_init[j] = _init_value[j] / (2*_nb_locus);
//    	_init_value[j] /= (2*_nb_locus); // used before model 4 was implemented
    	_init_value[j] /= (2*ttable[j].size()); // used before model 4 was implemented
//    	cout << "\n\t_init_value after: " << _init_value[j];
    }
    
  }
  
  if(_myProto->_allele_model < 3) { //for the di-allelic models

    if (_doInitMutation != 4){
        for(unsigned int i = 0; i < ttable.size(); i++) {
//            pos = i * _nb_traits;
            for(unsigned int j = 0; j < ttable[i].size(); j++) {
                _sequence[0][ttable[i][j]] = _myProto->_allele_value[i][0]; //set with the positive allele value
                _sequence[1][ttable[i][j]] = _myProto->_allele_value[i][0];
//                pos++;
            }
        }

    } else {
        //this is intended for a diallelic model with no initial variance
//    	cout << "Diallelic polarized condition" << endl;
        for(unsigned int i = 0; i < ttable.size(); i++) {
//            pos = i * _nb_traits;
            for(unsigned int j = 0; j < ttable[i].size(); j++) {
                if (j % 2 == 0){
                    _sequence[0][ttable[i][j]] = _myProto->_allele_value[i][0]; //set with the positive allele value
                    _sequence[1][ttable[i][j]] = _myProto->_allele_value[i][0]*-1;
//                    pos++;
                } else {
                    _sequence[0][ttable[i][j]] = _myProto->_allele_value[i][0]*-1; //set with the negative allele value
                    _sequence[1][ttable[i][j]] = _myProto->_allele_value[i][0];
//                    pos++;
                }
            }
        }
        cout << "\t after init_sequence  _doInitMutation=4 (Polarize)\n";
        for(unsigned int i=0; i<2; ++i){
          for(unsigned int j=0; j<_seq_length; ++j){
        	  cout << _sequence[i][j] << " | ";
          }
          cout << endl;
        }
        cout << endl;

    }
  } else {
	  /////////////////////////////////////////////////////////////////
  //set the allele values from the trait value
  //if(get_parameter("quanti_pleio_matrix")->isSet()){
	  for(unsigned int i=0; i<ttable.size(); ++i){
		  double myinit = (_init_value[i] * _nb_locus * 2) / (2* ttable[i].size());
		  for(unsigned int j=0; j<ttable[i].size(); ++j){
			  //cout << "\t_init_value after after: " << myinit << endl;
			  _sequence[0][ttable[i][j]]=myinit;
			  _sequence[1][ttable[i][j]]=myinit;
		  }
	  }
  }
  //print out _sequence (for debugging purposes only
// for(unsigned int i=0; i<2; ++i){
//	  for(unsigned int j=0; j<_seq_length; ++j){
//		  cout << _sequence[i][j] << " | ";
//	  }
//	  cout << endl;
//  }

  // OLD _sequence builder
/*  for(unsigned int i = 0; i < _nb_locus; i++) {
    pos = i * _nb_traits;
    for(unsigned int j = 0; j < _nb_traits; j++) {
      _sequence[0][pos] = _init_value[j];
      _sequence[1][pos] = _init_value[j];
      pos++;
    }
  }*/
  
  //add random effects to allele values
  if(_doInitMutation != 0 && _doInitMutation !=4) {
    double *mut1, *mut2;
    unsigned int L;
    
    for(unsigned int i = 0; i < _nb_locus; i++) {
      
    //  pos = i * _nb_traits;
      
      mut1 = _myProto->getMutationEffects(i);
      mut2 = _myProto->getMutationEffects(i);

      if(_myProto->_allele_model < 3) {

    	  for(unsigned int j=0; j<ltable[i][1]; j++){
        	 //cout << "\tmut1:" << i << ": " << mut1[j] << endl;
             _sequence[0][(ltable[i][0])+j] = mut1[j];
             _sequence[1][(ltable[i][0])+j] = mut2[j];
          }

      } else {


      	  for(unsigned int j=0; j<ltable[i][1]; j++){
      		  //cout << "\tmut1:" << i << ": " << mut1[j] << endl;
      		  _sequence[0][(ltable[i][0])+j] += mut1[j];
      		  _sequence[1][(ltable[i][0])+j] += mut2[j];
      	  }

      	  /* Old mut effects
      	  for(unsigned int j = 0; j < _nb_traits; j++) {
        	_sequence[0][pos] += mut1[j];
        	_sequence[1][pos] += mut2[j];
        	pos++;
      	  } */
      
      	  if(_doInitMutation > 1) { // the make-parts algorithm
        
      		  //select a random locus
      		  do{
      			  L = RAND::Uniform(_nb_locus);
      		  }while(L == i);
        
      		  //subtract the previous random deviates from that locus
      		  for(unsigned int j=0; j<ltable[L][1]; j++){
      			  //cout << "\tunmut1: " << L << ":" << mut1[j] << endl;
      			  _sequence[0][(ltable[L][0])+j] -= mut1[j]; //TODO update to add correct effects to correct part of _sequence
      			  _sequence[1][(ltable[L][0])+j] -= mut2[j];
      		  }

      		  /*        pos = L * _nb_traits;
        	// OLD subtract the previous random deviates from that locus
        	for(unsigned int j = 0; j < _nb_traits; j++) {
          	  _sequence[0][pos] -= mut1[j];
          	  _sequence[1][pos] -= mut2[j];
          	  pos++;
        	}*/
          }
      }
    }
  } 
//  cout << "\t after init_sequence  _doInitMutation \n";
//  for(unsigned int i=0; i<2; ++i){
//	  for(unsigned int j=0; j<_seq_length; ++j){
//		  cout << _sequence[i][j] << " | ";
//  	  }
//  	  cout << endl;
//  }
//  cout << endl;
//  cout << "\nEnd of TTQuanti::init_sequence!\t\n";
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
inline void TTQuanti::reset ()
{
  if(_sequence != NULL) {
    delete [] _sequence[0]; 
    delete [] _sequence[1];
    delete [] _sequence; 
    _sequence = NULL;
  }
  if(_phenotypes != NULL) delete [] _phenotypes;
  _phenotypes = NULL;
  
  if(_init_value) delete [] _init_value;
  _init_value = NULL;
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TTQuanti::inherit (TTrait* mother, TTrait* father)
{

 // cout << "\nStart of TTQuanti::inherit!\t";
/*  cout << " Before sequence: \n";
  for(unsigned int i=0; i<2; ++i){
    for(unsigned int j=0; j<_seq_length; ++j){
	  cout << _sequence[i][j] << " | ";
	}
	cout << endl;
  }*/

  double** mother_seq = (double**)mother->get_sequence();
  double** father_seq = (double**)father->get_sequence();
  
  (_myProto->* _inherit) (FEM, _sequence[FEM], mother_seq);
  
  (_myProto->* _inherit) (MAL, _sequence[MAL], father_seq);

/*  cout << " After sequence: \n";
  for(unsigned int i=0; i<2; ++i){
    for(unsigned int j=0; j<_seq_length; ++j){
	  cout << _sequence[i][j] << " | ";
	}
	cout << endl;
  }*/
 // cout << "\nEnd of TTQuanti::inherit!\t";
 }
// ----------------------------------------------------------------------------------------
// mutate_noHC
// ----------------------------------------------------------------------------------------
inline void TTQuanti::mutate_noHC ()
{
  //cout << "\tStart of TTQuanti::mutate_noHC!\n";
  /*cout << " Before sequence: \n";
  for(unsigned int i=0; i<2; ++i){
	  for(unsigned int j=0; j<_seq_length; ++j){
		  cout << _sequence[i][j] << " | ";
	  }
	  cout << endl;
  }
  cout << endl;*/

  //unsigned int NbMut = 1;
  vector< vector<unsigned int> > ltable = _myProto->get_locus_table();
//  unsigned int NbMut = (unsigned int)RAND::Poisson(_genomic_mutation_rate);
  unsigned int NbMut = (unsigned int)RAND::Binomial(_genomic_mutation_rate/2/_nb_locus, _nb_locus*ltable[0][1]); // Binomial takes per-locus mutation rate and genome size (N.B. Here we have assumed that all loci have the same pleiotropic degree as locus 1)
  unsigned int mut_locus, mut_all;// pos;
  double *effects;

  while(NbMut != 0) {
      mut_locus = RAND::Uniform(_nb_locus);
      effects = _myProto->getMutationEffects(mut_locus);
      mut_all = RAND::RandBool();
      for(unsigned int i=0; i<ltable[mut_locus][1]; i++){
        _sequence[mut_all][(ltable[mut_locus][0])+i] += effects[i];///mutations are added to existing alleles
      }
      NbMut--;
  }

/*  cout << "Locus Mutation Effects: \n";
  for(unsigned int i=0; i<_nb_locus; i++){
	  effects = (_myProto->*_getMutationValues)(i);
	  cout << "\tEffects for Locus " << i << ": ";
	  for(unsigned int j=0; j<_nb_traits; j++){
		  cout << effects[j] << " | ";
	  }
	  cout << endl;
  }
  cout << endl;*/


  // OLD Mutation loop noHC
/*  while(NbMut != 0) {
    mut_locus = RAND::Uniform(_nb_locus);
    effects = (_myProto->*_getMutationValues)(mut_locus);
    mut_all = RAND::RandBool();
    pos = mut_locus*_nb_traits;
    for(unsigned int i = 0; i < _nb_traits; i++)
      _sequence[mut_all][pos + i] += effects[i];///mutations are added to existing alleles
    
    NbMut--; 
  }*/

  //cout << "mutate_noHC Locus table: " << ltable.size() << endl;
  //print out _sequence (for debugging purposes only
  /*cout << "\t mutate_noHC After sequence: \n";
  for(unsigned int i=0; i<2; ++i){
	  for(unsigned int j=0; j<_seq_length; ++j){
		  cout << _sequence[i][j] << " | ";
	  }
	  cout << endl;
  }
  cout << endl;*/
  //cout << "\tEnd of TTQuanti::mutate_noHC!\n";

}
// ----------------------------------------------------------------------------------------
// mutate_HC
// ----------------------------------------------------------------------------------------
inline void TTQuanti::mutate_HC ()
{
  vector< vector<unsigned int> > ltable = _myProto->get_locus_table();

//  cout << "\t mutate_HC Before sequence: \n";
//  for(unsigned int i=0; i<2; ++i){
//	  for(unsigned int j=0; j<_seq_length; ++j){
//		  cout << _sequence[i][j] << " | ";
//	  }
//	  cout << endl;
//  }
//  cout << endl;


//  unsigned int NbMut = (unsigned int)RAND::Poisson(_genomic_mutation_rate); // instead use binomial as done neutral traits
  unsigned int NbMut = (unsigned int)RAND::Binomial(_genomic_mutation_rate/2/_nb_locus, _nb_locus*ltable[0][1]); // Binomial takes per-locus mutation rate and genome size (N.B. Here we have assumed that all loci have the same pleiotropic degree as locus 1)
  unsigned int mut_locus, mut_all;//, pos;
  double *effects;
//  cout << "\nNumber of Mutations: " << NbMut << endl;
  while(NbMut != 0) {
      mut_locus = RAND::Uniform(_nb_locus);
      effects = _myProto->getMutationEffects(mut_locus);
      mut_all = RAND::RandBool();
      for(unsigned int i=0; i<ltable[mut_locus][1]; i++){
        _sequence[mut_all][(ltable[mut_locus][0])+i] = effects[i];///mutations replace existing alleles
      }
      NbMut--;
  }

  // OLD Mutation loop HC
/*  while(NbMut != 0) {
    mut_locus = RAND::Uniform(_nb_locus);
    effects = _myProto->getMutationEffects(mut_locus);
    mut_all = RAND::RandBool();
    pos = mut_locus*_nb_traits;
    for(unsigned int i = 0; i < _nb_traits; i++)
      _sequence[mut_all][pos + i] = effects[i]; ///mutations replace existing alleles
    
    NbMut--; 
  }*/

//  cout << "\t mutate_HC After sequence: \n";
//  for(unsigned int i=0; i<2; ++i){
//	  for(unsigned int j=0; j<_seq_length; ++j){
//		  cout << _sequence[i][j] << " | ";
//	  }
//	  cout << endl;
//  }
//  cout << endl;

} 
// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
inline void TTQuanti::set_value ()
{
  //cout << "\nStart of TTQuanti::set_value!\t";

  //register unsigned int loc;
  vector< vector<unsigned int> > ttable = _myProto->get_trait_table();
  // vector< vector<unsigned int> > ltable = _myProto->get_locus_table();

  // print _sequence and ttable for debugging purposes
//  for(unsigned int i=0; i<2; ++i){
// 	  for(unsigned int j=0; j<_seq_length; ++j){
// 		  cout << _sequence[i][j] << " | ";
// 	  }
// 	  cout << endl;
//  }
//  cout << endl;
/*  cout << "\nttable: ";
  for(unsigned int i=0; i<ttable.size(); ++i){
	  for(unsigned j=0; j<ttable[i].size(); ++j){
		  cout << ttable[i][j] << "|";
	  }
	  cout << endl;
  }*/


  for(unsigned int i = 0; i < _nb_traits; ++i) 
    _phenotypes[i] = 0;
  
  for(unsigned int j = 0; j < ttable.size(); ++j) {
      for(unsigned int i = 0; i < ttable[j].size(); ++i) {
        _phenotypes[j] += (_sequence[0][ttable[j][i]] + _sequence[1][ttable[j][i]]);
        if(abs(_phenotypes[j]) < 0.0000000001){_phenotypes[j]=0;} // min_tolerance checks in all cases where floating point addition of doubles is used
//        cout << "Added sequence[ttable pos]: " << j << "|" << i << ". Values: " << _sequence[0][ttable[j][i]] << "|" << _sequence[1][ttable[j][i]] << ". New Phenotype: " << _phenotypes[j] << endl;
      }
  }

  //print out values
//  for(unsigned int i = 0; i < _nb_traits; ++i)
//	  cout << "set_value Tval" << i << ": "<< _phenotypes[i] << " | ";
//  cout << endl;

/*  // OLD set_value
  for(unsigned int j = 0; j < _nb_locus; ++j) {
    loc = j * _nb_traits;
    for(unsigned int i = 0; i < _nb_traits; ++i) {
      _phenotypes[i] += _sequence[0][loc] + _sequence[1][loc];
      loc++;
    }
  }*/
  
  if(_eVariance != 0)
    for(unsigned int i = 0; i < _nb_traits; ++i) 
      _phenotypes[i] += RAND::Gaussian(_eVariance);

 // cout << "\tEnd of TTQuanti::set_value!\n";
}

// ----------------------------------------------------------------------------------------
// get_genotype
// ----------------------------------------------------------------------------------------
double TTQuanti::get_genotype (unsigned int trait)
{
  //unsigned int loc;
  double genotype = 0;
  vector< vector<unsigned int> > ttable = _myProto->get_trait_table();
  
  for(unsigned int i = 0; i < ttable[trait].size(); ++i) {
	  genotype += _sequence[0][ttable[trait][i]] + _sequence[1][ttable[trait][i]];
      if(abs(genotype) < 0.0000000001){genotype=0;} // min_tolerance checks in all cases where floating point addition of doubles is used
  }
//  cout << "get_genotype Trait" << trait << ": "<< genotype << endl;

  // OLD genotype calculator
/*  for(unsigned int j = 0; j < _nb_locus; ++j) {
    loc = j * _nb_traits + trait;
    genotype += _sequence[0][loc] + _sequence[1][loc];
  } */

  return genotype;
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTQuanti::show_up  ()
{
  message("\
          Trait's type: QUANTI\n\
          traits: %i\n\
          loci: %i\n",_nb_traits,_nb_locus);
  
  for(unsigned int i = 0; i < _nb_traits; i++)
    message("phenotype %i: %f\n",i+1,_phenotypes[i]);
  
  message("_sequence: \n0:");
  unsigned int loc;
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("\nt%i:",i+1);
    for(unsigned int j = 0; (j < _nb_locus); ++j) {
      loc = j * _nb_traits + i;
      message("%.3f,",_sequence[0][loc]);
    }
  }
  message("\n1:");
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("\nt%i:",i+1);
    for(unsigned int j = 0; (j < _nb_locus); ++j) {
      loc = j * _nb_traits + i;
      message("%.3f,",_sequence[1][loc]);
    }
  }
  message("\n");
  
}  
// ------------------------------------------------------------------------------

//                             TTQuantiSH

// ----------------------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------------------
void TTQuantiSH::resetPtrs ( )
{  
  
  if(_G != NULL) gsl_matrix_free(_G);
  if(_eval != NULL) gsl_vector_free(_eval);
  if(_evec != NULL) gsl_matrix_free(_evec);
  if(_ws != NULL)  gsl_eigen_symmv_free (_ws);
  
  if(_meanP != NULL) delete [] _meanP;
  if(_meanG != NULL) delete [] _meanG;
  if(_Va != NULL) delete [] _Va;
  if(_Vb != NULL) delete [] _Vb;
  if(_Vp != NULL) delete [] _Vp;
  if(_covar != NULL) delete [] _covar;
  if(_eigval != NULL) delete [] _eigval;
  
  
  if(_eigvect) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _eigvect[i];
    delete [] _eigvect;
  }
  if(_pVa) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pVa[i];
    delete [] _pVa;
  }
  if(_pVp) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pVp[i];
    delete [] _pVp;
  }
  if(_pmeanP) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pmeanP[i];
    delete [] _pmeanP;
  }
  if(_pmeanG) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pmeanG[i];
    delete [] _pmeanG;
  }
  if(_peigval) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _peigval[i];
    delete [] _peigval;
  }
  
  if(_pcovar != NULL) {  
    for(unsigned int i = 0; i < _patchNbr; i++) delete [] _pcovar[i];
    delete [] _pcovar;
  }
  
  if(_peigvect != NULL) {
    for(unsigned int i = 0; i < _patchNbr; i++) delete [] _peigvect[i];
    delete []  _peigvect;
  }
  
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTQuantiSH::init()
{  
  //cout << "\nStart of TTQuantiSH::init!\t";

  StatHandlerBase::init();
  
  _eVar = (_SHLinkedTrait->get_env_var() != 0);
  
  resetPtrs(); //deallocate anything that was previously allocated
  
  if(_patchNbr != _pop->getPatchNbr())
    _patchNbr = _pop->getPatchNbr();
  
  if(_nb_trait != _SHLinkedTrait->get_nb_traits())
    _nb_trait = _SHLinkedTrait->get_nb_traits();
  
  _G = gsl_matrix_alloc(_nb_trait,_nb_trait);
  _eval = gsl_vector_alloc (_nb_trait);
  _evec = gsl_matrix_alloc (_nb_trait, _nb_trait);
  _ws = gsl_eigen_symmv_alloc (_nb_trait);
  
  _meanP = new double [_nb_trait];
  _meanG = new double [_nb_trait];
  _Va = new double [_nb_trait];
  _Vb = new double [_nb_trait];
  _Vp = new double [_nb_trait];
  _eigval = new double [_nb_trait];
  
  _eigvect = new double* [_nb_trait];
  _pVa = new double* [_nb_trait];
  _pVp = new double* [_nb_trait];
  _pmeanP = new double* [_nb_trait];
  _pmeanG = new double* [_nb_trait];
  _peigval = new double* [_nb_trait];
  
  for(unsigned int i=0; i < _nb_trait; ++i) {
    
    _eigvect[i]  = new double [_nb_trait];
    
    _pVa[i] = new double [_patchNbr];
    _pVp[i] = new double [_patchNbr];
    _pmeanP[i] = new double [_patchNbr];
    _pmeanG[i] = new double [_patchNbr];
    _peigval[i] = new double [_patchNbr];
    
  }
  
  _peigvect = new double* [_patchNbr];
  for(unsigned int i = 0; i < _patchNbr; i++)     
    _peigvect[i] = new double [_nb_trait*_nb_trait]; 
  
  
  if(_nb_trait > 1) {
    
    _covar = new double [_nb_trait*(_nb_trait -1)/2];
    
    _pcovar = new double* [_patchNbr];

    for(unsigned int i = 0; i < _patchNbr; i++) 
      _pcovar[i] = new double [_nb_trait*(_nb_trait - 1)/2];
    
  }
  //cout << "\nEnd of TTQuantiSH::init!\t";
  
}
// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTQuantiSH::setStatRecorders ");
#endif   
  string age_tag = token.substr(0,token.find_first_of("."));
  string sub_token;
  age_t AGE = ALL;
  
  if (age_tag.size() != 0 && age_tag.size() != string::npos) {
    
    if (age_tag == "adlt") AGE = ADULTS;
    
    else if (age_tag == "off") AGE = OFFSPRG;
    
    else age_tag = "";
    
  } else {
    age_tag = "";
  }
  
  if (age_tag.size() != 0) 
    sub_token = token.substr(token.find_first_of(".") + 1, string::npos);
  else
    sub_token = token;
  
  //!!! attention, y a un prob ici quand pas de prefix d'age!!!!
  
  if(sub_token == "quanti") {
    addQuanti(AGE);
  } else if(sub_token == "quanti.eigen") {  //based on Vb; among-patch (D) matrix
    addEigen(AGE);
  } else if(sub_token == "quanti.eigenvalues") {  //based on Vb; among-patch (D) matrix
    addEigenValues(AGE);
  } else if(sub_token == "quanti.eigenvect1") {  //based on Vb; among-patch (D) matrix
    addEigenVect1(AGE);
  } else if(sub_token == "quanti.patch") {
    addQuantiPerPatch(AGE);
  } else if(sub_token == "quanti.mean.patch") {
    addAvgPerPatch(AGE);
  } else if(sub_token == "quanti.var.patch") {
    addVarPerPatch(AGE);
  } else if(sub_token == "quanti.covar.patch") {
    addCovarPerPatch(AGE);
  } else if(sub_token == "quanti.eigen.patch") {
    addEigenPerPatch(AGE);
  } else if(sub_token == "quanti.eigenvalues.patch") {
    addEigenValuesPerPatch(AGE);
  } else if(sub_token == "quanti.eigenvect1.patch") {
    addEigenVect1PerPatch(AGE);
  } else if(sub_token == "quanti.skew.patch") {
    addSkewPerPatch(AGE);
  } else
    return false;
  
  return true;
}// ----------------------------------------------------------------------------------------
// addQuanti
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addQuanti (age_t AGE)
{
  if (AGE == ALL) {
    addQuanti(ADULTS);
    addQuanti(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name = suffix + "q";
  string t1, t2;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("",suffix + "q1",AGE,0,0,0,&TTQuantiSH::getMeanPhenot,0,setter);
  
  for(unsigned int i = 1; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1,AGE,i,0,0,&TTQuantiSH::getMeanPhenot,0,0);
  } 
  
  //Va
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Va",AGE,i,0,0,&TTQuantiSH::getVa,0,0);
  }
  
  //Vb
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Vb",AGE,i,0,0,&TTQuantiSH::getVb,0,0);
  }
  
  //Vp
  if(_eVar) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      t1 = tstring::int2str(i+1);
      add("", name + t1 +".Vp",AGE,i,0,0,&TTQuantiSH::getVp,0,0);
    }
  }
  //Qst
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Qst",AGE,i,0,0,&TTQuantiSH::getQst,0,0);
  }
  
  //cov
  if (_nb_trait > 1) {
    unsigned int c = 0, cov;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1; v < _nb_trait; ++v) {
    	cov = (t+1)*10+(v+1) ;
    	t1 = tstring::int2str(cov);
        add("", name + t1 +".cov",AGE,c++,0,0,&TTQuantiSH::getCovar,0,0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addEigen
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigen (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.eval1",AGE,0,0,0,&TTQuantiSH::getEigenValue, 0, setter); //this one calls the setter
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "q.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiSH::getEigenValue, 0, 0);
  
  for(unsigned int t = 0; t< _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v)
      add("", suffix + "q.evect" + tstring::int2str((t+1)*10+(v+1)),AGE,t,v,0,0,&TTQuantiSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addEigenValues
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenValues (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.eval1",AGE,0,0,0,&TTQuantiSH::getEigenValue, 0, setter);
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "q.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiSH::getEigenValue, 0, 0);
  
  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1 : save only the first eigenvector
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenVect1 (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.evect11",AGE,0,0,0,0,&TTQuantiSH::getEigenVectorElt,setter);
  
  for(unsigned int v = 1; v < _nb_trait; ++v)
    add("", suffix + "q.evect1" + tstring::int2str(v+1),AGE,0,v,0,0,&TTQuantiSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addQuantiPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addQuantiPerPatch (age_t AGE)
{
  
  if (AGE == ALL) {
    addQuantiPerPatch(ADULTS);
    addQuantiPerPatch(OFFSPRG);
    return;
  }
  
  addAvgPerPatch(AGE);
  addVarPerPatch(AGE);
  addCovarPerPatch(AGE);
  addEigenPerPatch(AGE);
  
}
// ----------------------------------------------------------------------------------------
// addAvgPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addAvgPerPatch (age_t AGE)
{
  if (AGE == ALL) {
    addAvgPerPatch(ADULTS);
    addAvgPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("Mean phenotype of trait 1 in patch 1", suffix + "q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getMeanPhenotPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Mean phenotype of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getMeanPhenotPerPatch, 0);
    } 
  }
  
}
// ----------------------------------------------------------------------------------------
// addVarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addVarPerPatch (age_t AGE)
{
  if (AGE == ALL) {
    addVarPerPatch(ADULTS);
    addVarPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats :
                                        &TTQuantiSH::setOffsprgStats);
  
  add("Genetic variance of trait 1 in patch 1", suffix + "Va.q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getVaPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Va." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getVaPerPatch, 0);
    }
  }
  
  if(_eVar) {
    for(unsigned int p = 0; p < patchNbr; p++) {
      for(unsigned int i = 0; i < _nb_trait; i++) {
        name = "Phenotypic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
        t1 = "q" + tstring::int2str(i+1);
        patch = ".p" + tstring::int2str(p+1);
        add(name, suffix + "Vp." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getVpPerPatch, 0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addCovarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addCovarPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording traits covariance with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addCovarPerPatch(ADULTS);
    addCovarPerPatch(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  string cov;
  unsigned int patchNbr = _pop->getPatchNbr();
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("Genetic covariance of trait 1 and trait 2 in patch 1", suffix + "cov.q12.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getCovarPerPatch, setter);
  
  unsigned int c;
  for(unsigned int p = 0; p < patchNbr; p++) {
    patch = ".p" + tstring::int2str(p+1);
    c = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1; v < _nb_trait; ++v){
        if(p==0 && t==0 && v==1) {c++; continue;}
        cov = tstring::int2str((t+1)*10+v+1);
        add("", suffix + "cov.q" + cov + patch,  AGE, p, c++, 0, 0, &TTQuantiSH::getCovarPerPatch, 0);
      }
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  unsigned int pv =0;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvalue in patch 1", suffix + "qeval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "qeval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiSH::getEigenValuePerPatch,0);
    }
  }
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t)
      for(unsigned int v = 0; v < _nb_trait; ++v)
        add("", suffix + "qevect" + tstring::int2str((t+1)*10+v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiSH::getEigenVectorEltPerPatch,0);
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenValuesPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenValuesPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("First G-matrix eigenvalue in patch 1", suffix + "qeval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "qeval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiSH::getEigenValuePerPatch,0);
    }
  }  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1PerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenVect1PerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  unsigned int pv =0;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : 
                                        &TTQuantiSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvector in patch 1", suffix + "qevect11.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenVectorEltPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    //    for(unsigned int t = 0; t < _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v){
      if(p==0 && v==0) {pv++; continue;}
      add("", suffix + "qevect1" + tstring::int2str(v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiSH::getEigenVectorEltPerPatch,0);
    }
  }
}
// ----------------------------------------------------------------------------------------
// addSkewPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addSkewPerPatch(age_t AGE)
{
  if (AGE == ALL) {
    addSkewPerPatch(ADULTS);
    addSkewPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats :
                                        &TTQuantiSH::setOffsprgStats);
  
  add("Genetic skew of trait 1 in patch 1", suffix + "Sk.q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getSkewPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic skew of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Sk." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getSkewPerPatch, 0);
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// getSkewPerPatch
// ----------------------------------------------------------------------------------------
double TTQuantiSH::getSkewPerPatch (unsigned int i, unsigned int p)
{
  double skew = 0;
  
  double *phenot = _phenoTable.getClassWithinGroup(i, p);
  unsigned int patch_size = _phenoTable.size(i, p);  
  
  for(unsigned int k = 0; k < patch_size; ++k)
    skew += pow( phenot[k] - _pmeanP[i][p], 3);  //the mean has been set by setStats()
  
  return skew / patch_size;
}
// ----------------------------------------------------------------------------------------
// setDataTables
// ----------------------------------------------------------------------------------------
void TTQuantiSH::setDataTables(age_t AGE) 
{
  unsigned int **sizes;
  unsigned int nb_patch = _pop->getPatchNbr();
  
  sizes = new unsigned int * [_nb_trait];
  
  for(unsigned int i = 0; i < _nb_trait; ++i) {
    sizes[i] = new unsigned int [nb_patch];
    for(unsigned int j = 0; j < nb_patch; ++j)
      sizes[i][j] = _pop->size(AGE, j);
  }
  
  _phenoTable.update(_nb_trait, nb_patch, sizes);
  _genoTable.update(_nb_trait, nb_patch, sizes);
  
  for(unsigned int i = 0; i < _nb_trait; ++i)
    delete [] sizes[i];
  delete [] sizes;
  
  Patch* patch;
  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);
  
  for(unsigned int i = 0, n; i < nb_patch; i++) {
    
    patch = _pop->getPatch(i);
    
    n=0;
    
    if ((patch->size(MAL, age)+patch->size(FEM, age)) != _phenoTable.size(0,i)) {
      fatal("problem while recording quanti trait values; table size doesn't match patch size.\n");
    }
    store_quanti_trait_values(patch, i, patch->size(MAL, age), &n, MAL, age, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    store_quanti_trait_values(patch, i, patch->size(FEM, age), &n, FEM, age, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    if (n != _phenoTable.size(0,i) || n != _genoTable.size(0,i)) {
      fatal("problem while recording quanti trait values; size counter doesn't match table size.\n");
    }
  }
}
// ----------------------------------------------------------------------------------------
// store_trait_values
// ----------------------------------------------------------------------------------------
void store_quanti_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                         sex_t SEX, age_idx AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                         unsigned int nTrait, unsigned int TraitIndex)
{
  double *phe;
  TTQuanti *trait;
  
  for(unsigned int j = 0; j < size; ++j) {
    
    trait = dynamic_cast<TTQuanti*> (patch->get(SEX, AGE, j)->getTrait( TraitIndex ));
    
    phe = (double*)trait->getValue();
    
    for(unsigned int k = 0; k < nTrait; k++){
      ptable->set( k, patchID, (*cntr), phe[k]);
      gtable->set( k, patchID, (*cntr), trait->get_genotype(k));
    }
    (*cntr)++;
  }
  
}
// ----------------------------------------------------------------------------------------
// setStats
// ----------------------------------------------------------------------------------------
void TTQuantiSH::setStats (age_t AGE)
{  
  if(_table_set_age == AGE 
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate())
    return;
  
  unsigned int pop_size = _pop->size(AGE);
  unsigned int patch_size;
  double *phenot1, *genot1, *genot2;
  
  unsigned int nb_patch = _pop->getPatchNbr();
  
  if(nb_patch < _patchNbr) {
    warning("increase in patch number detected (in Quanti Stat Handler),");
    warning("stats for quanti trait will not be recorded in new patches, patch identity may have changed.\n");
    _patchNbr = nb_patch; // record stat only in remaining patches
  }
    
  if(nb_patch > _patchNbr)  nb_patch = _patchNbr;  //tables have been allocated for _patchNbr patches
    
  setDataTables(AGE);
  
#ifdef HAS_GSL
  unsigned int c = 0; //covariance position counter
  unsigned int pv = 0; //eigenvector position counter
  //within deme stats:
  for(unsigned int j = 0; j < nb_patch; j++) { 
    
    patch_size = _pop->size(AGE, j);
    
    for(unsigned int t=0; t < _nb_trait; t++) {
      
      phenot1 = _phenoTable.getClassWithinGroup(t,j);
      genot1  = _genoTable.getClassWithinGroup(t,j);
      
      _pmeanP[t][j] = my_mean (phenot1, patch_size );
      _pmeanG[t][j] = my_mean (genot1, patch_size );
      _pVp[t][j] = my_variance_with_fixed_mean (phenot1, patch_size, _pmeanP[t][j]);
      _pVa[t][j] = my_variance_with_fixed_mean (genot1,  patch_size, _pmeanG[t][j]);
           
    }
    
    if(_nb_trait > 1) { 
      c = 0;
      //    calculate the covariances and G, need to adjust dimensions in class declaration
      for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
        //      set the diagonal elements of G here
        
        gsl_matrix_set(_G, t1, t1, _pVa[t1][j]);
        
        for(unsigned int t2 = t1 + 1; t2 < _nb_trait; t2++) {
          
          genot1  = _genoTable.getClassWithinGroup(t1,j);
          genot2  = _genoTable.getClassWithinGroup(t2,j);
          
          _pcovar[j][c] = gsl_stats_covariance_m (genot1, 1, genot2, 1, patch_size,
                                                  _pmeanG[t1][j], _pmeanG[t2][j]);
          
          gsl_matrix_set(_G, t1, t2, _pcovar[j][c]);
          gsl_matrix_set(_G, t2, t1, _pcovar[j][c++]);
        }
      }
      
      gsl_eigen_symmv (_G, _eval, _evec, _ws);
      gsl_eigen_symmv_sort (_eval, _evec, GSL_EIGEN_SORT_VAL_DESC);
      
      pv = 0;
      
      for(unsigned int t = 0; t < _nb_trait; t++) {
        _peigval[t][j] = gsl_vector_get (_eval, t);
        for(unsigned int v = 0; v < _nb_trait; v++)
          _peigvect[j][pv++] = gsl_matrix_get (_evec, v, t); //read eigenvectors column-wise
      }
    }
  }
  
  double meanGamong1, meanGamong2;
  c = 0; //reset covariance positioner
  //among demes stats:
  for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
    
    phenot1 = _phenoTable.getGroup(t1);
    genot1  = _genoTable.getGroup(t1);
    
    _meanP[t1] = my_mean (phenot1, pop_size ); //grand mean, pooling all individuals
    _meanG[t1] = my_mean (genot1, pop_size ); //grand mean, pooling all individuals
    
    _Vp[t1] = my_mean_no_nan (_pVp[t1], nb_patch); //mean within patch variance
    _Va[t1] = my_mean_no_nan (_pVa[t1], nb_patch); //mean within patch variance
    
    meanGamong1 = my_mean_no_nan (_pmeanG[t1], nb_patch); //mean of within patch mean genotypic values
    
    _Vb[t1] = my_variance_with_fixed_mean_no_nan (_pmeanG[t1], nb_patch,  meanGamong1); //variance of patch means
    
    gsl_matrix_set(_G, t1, t1, _Vb[t1]); //_G here becomes the D-matrix, the among-deme (Difference) covariance matrix 
    
    for(unsigned int t2 = t1 + 1; t2 < _nb_trait; t2++) {
      
      meanGamong2 = my_mean (_pmeanG[t2], nb_patch);
      
      _covar[c] = gsl_stats_covariance_m (_pmeanG[t1], 1, _pmeanG[t2], 1, nb_patch, meanGamong1, meanGamong2); //covariance of patch means
      
      gsl_matrix_set(_G, t1, t2, _covar[c]);
      gsl_matrix_set(_G, t2, t1, _covar[c++]);
      
    }
  }
  
  if(_nb_trait > 1) { 
    gsl_eigen_symmv (_G, _eval, _evec, _ws);
    gsl_eigen_symmv_sort (_eval, _evec, GSL_EIGEN_SORT_VAL_DESC);
    
    for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
      _eigval[t1] = gsl_vector_get (_eval, t1);
      for(unsigned int t2 = 0; t2 < _nb_trait; t2++) {
        _eigvect[t2][t1] = gsl_matrix_get (_evec, t2, t1);      
      }
    }
  }
  
#else
  fatal("install the GSL library to get the quanti stats!\n");
#endif
  
  _table_set_age = AGE;
  _table_set_gen = _pop->getCurrentGeneration();  
  _table_set_repl = _pop->getCurrentReplicate();
  
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFH::FHwrite()
{
  vector< vector<unsigned int> > ttable = _FHLinkedTrait->get_trait_table();
  Metapop* pop = get_pop_ptr();
  
  if (!pop->isAlive()) return;
  
  std::string filename = get_filename();
  
  std::ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
  
  bool print_gene = (_output_option == "genotypes" || _output_option == "genotype");
  bool print_genotype = (_FHLinkedTrait->get_env_var() != 0);
  
  FILE<<"pop ";
  if(print_gene) {
    // OLD file header creation for full pleiotropy (non-variable)
	/*for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++)
      for(unsigned int l = 0; l < _FHLinkedTrait->get_nb_locus(); l++)
    	          FILE<<"t"<<k+1<<"l"<<l+1<<"1 "<<"t"<<k+1<<"l"<<l+1<<"2 ";*/
    for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++) {
      for(unsigned int l = 0; l < ttable[k].size(); l++) {
    	  	  	  FILE<<"t"<<k+1<<"l"<<l+1<<"1 "<<"t"<<k+1<<"l"<<l+1<<"2 ";
      }
    }
  }
  
  for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++) {
    FILE<<"P"<<k+1<< " "; 
    if(print_genotype) FILE<<"G"<<k+1<< " ";
  }
  
  FILE<<"age sex home ped isMigrant father mother ID\n";
  
  age_t pop_age = pop->getCurrentAge(); //flag telling which age class should contain individuals
  
  //we print anything that is present in the pop:
  if( (pop_age & OFFSPRG) != 0) print(FILE, OFFSx, print_gene, print_genotype);
  if( (pop_age & ADULTS) != 0) print(FILE, ADLTx, print_gene, print_genotype);
  
  FILE.close();  
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFH::print(ofstream& FH, age_idx Ax, bool print_gene, bool print_genotype)
{
  Metapop* pop = get_pop_ptr();
  int patchNbr = pop->getPatchNbr();
  Patch* current_patch;
  Individual* ind;
  TTQuanti* trait;
  double* Tval;
  double **genes;
  unsigned int loc, nb_trait=_FHLinkedTrait->get_nb_traits();//, nb_locus = _FHLinkedTrait->get_nb_locus();
  vector< vector<unsigned int> > ttable =_FHLinkedTrait->get_trait_table();

  for(int i = 0; i < patchNbr; i++) {
    
    current_patch = pop->getPatch(i);
    
    for(unsigned int j = 0, size = current_patch->size(FEM, Ax); j < size; j++) {
      
      ind = current_patch->get(FEM, Ax, j);
      trait = dynamic_cast<TTQuanti*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" ";
      Tval = (double*)trait->getValue();
      
      if(print_gene){
        genes = (double**)trait->get_sequence();
        
        FH.precision(6);
        // OLD gene value output for full pleiotropy (non-variable)
        /*for(unsigned int k = 0; k < nb_trait; k++) {
		  for(unsigned int l = 0; l < nb_locus; l++) {
		    loc = l * nb_trait + k;
		    FH<<genes[0][loc]<<" "<<genes[1][loc]<<" ";
		  }
		}*/
        for(unsigned int k = 0; k < ttable.size(); k++) {
          for(unsigned int l = 0; l < ttable[k].size(); l++) {
        	    FH<<genes[0][ttable[k][l]]<<" "<<genes[1][ttable[k][l]]<<" ";
          }
        }
      }
      
      FH.precision(4);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<Tval[k]<<" ";
        if(print_genotype) FH << trait->get_genotype(k) << " ";
      }
      
      FH<<Ax<<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
      << (ind->getFather() && ind->getMother() ?
          (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)
      <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;
    }
    
    for(unsigned int j = 0, size = current_patch->size(MAL, Ax); j < size; j++) {
      
      ind = current_patch->get(MAL, Ax, j);
      trait = dynamic_cast<TTQuanti*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" ";
      
      Tval = (double*)trait->getValue();
      
      if(print_gene){
        genes = (double**)trait->get_sequence();
        
        FH.precision(6);
        // OLD gene value output for full pleiotropy (non-variable)
        /*for(unsigned int k = 0; k < nb_trait; k++) {
          for(unsigned int l = 0; l < nb_locus; l++) {
            loc = l * nb_trait + k;
            
            FH<<genes[0][loc]<<" "<<genes[1][loc]<<" ";
          }
        }*/
        for(unsigned int k = 0; k < ttable.size(); k++) {
          for(unsigned int l = 0; l < ttable[k].size(); l++) {
        	    FH<<genes[0][ttable[k][l]]<<" "<<genes[1][ttable[k][l]]<<" ";
          }
        }
      }
      
      FH.precision(4);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<Tval[k]<<" ";
        if(print_genotype) FH << trait->get_genotype(k) << " ";
      }
      
      FH<<Ax<<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass()<<" "
      << (ind->getFather() && ind->getMother() ?
          (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)<<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;
    }
  }
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQFreqExtractor::FHwrite()
{
  Metapop* pop = get_pop_ptr();
  int patchNbr = pop->getPatchNbr();
  Patch* current_patch;
  Individual* ind;
  double **seq;
  double **val_t1, **val_t2, *dist_t1[2];//, *dist_t2[2];
  unsigned int t1=0, t2=0, bloc, genome_size=_FHLinkedTrait->get_nb_locus(),
  _nb_trait=_FHLinkedTrait->get_nb_traits(), nb_locus = _FHLinkedTrait->get_nb_locus();
  
  std::string filename = get_filename();
  
  std::ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
  
  FILE.close();
  
  age_t pop_age = pop->getCurrentAge();
  
  if( (pop_age & ADULTS) != 0) {
    
    unsigned int val_size = pop->size(ADULTS) * _nb_trait * 2;
    //    cout<<"allocating the val arrays (size "<<val_size<<")"<<endl;
    val_t1 = new double* [ nb_locus ];
    val_t2 = new double* [ nb_locus ];
    
    for(unsigned int i = 0; i < nb_locus; i++) {
      val_t1[i] = new double [ val_size ];
      val_t2[i] = new double [ val_size ];
    }
    //    cout<<"getting the values"<<endl;
    for(int i = 0; i < patchNbr; i++) {
      
      current_patch = pop->getPatch(i);
      
      for(unsigned int j = 0, size = current_patch->size(FEM, ADLTx); j < size; j++) {
        
        ind = current_patch->get(FEM, ADLTx, j);
        
        seq = (double**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
        
        for(unsigned int k = 0; k < nb_locus; ++k) {
          bloc = k * _nb_trait;
          t2 = t1 = bloc * 2;
          val_t1[k][ t1++ ] = seq[0][bloc];
          val_t1[k][ t1++ ] = seq[1][bloc];
          val_t1[k][ t1++ ] = seq[0][bloc + genome_size];
          val_t1[k][ t1 ] = seq[1][bloc + genome_size];
          bloc++;
          val_t2[k][ t2++ ] = seq[0][bloc];
          val_t2[k][ t2++ ] = seq[1][bloc];
          val_t2[k][ t2++ ] = seq[0][bloc + genome_size];
          val_t2[k][ t2 ] = seq[1][bloc + genome_size];
        }
        
      }
      
      for(unsigned int j = 0, size = current_patch->size(MAL, ADLTx); j < size; j++) {
        ind = current_patch->get(MAL, ADLTx, j);
        seq = (double**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
        
        for(unsigned int k = 0; k < nb_locus; ++k) {
          bloc = k * _nb_trait;
          t2 = t1 = bloc * 2;
          val_t1[k][ t1++ ] = seq[0][bloc];
          val_t1[k][ t1++ ] = seq[1][bloc];
          val_t1[k][ t1++ ] = seq[0][bloc + genome_size];
          val_t1[k][ t1 ] = seq[1][bloc + genome_size];
          bloc++;
          val_t2[k][ t2++ ] = seq[0][bloc];
          val_t2[k][ t2++ ] = seq[1][bloc];
          val_t2[k][ t2++ ] = seq[0][bloc + genome_size];
          val_t2[k][ t2 ] = seq[1][bloc + genome_size];
        }
      }
    }//end for patchNbr
    //dump the values to a text file
    filename += "tot";
    
    FILE.open(filename.c_str(), ios::out);
    
    if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
    
    for(unsigned int j = 0; j < nb_locus; j++)
      FILE<<"t1.l"<<j+1<<" ";
    for(unsigned int j = 0; j < nb_locus; j++)
      FILE<<"t2.l"<<j+1<<" ";
    
    FILE<<endl;
    
    for(unsigned int i = 0; i < val_size; i++) {
      
      for(unsigned int j = 0; j < nb_locus; j++)
        FILE<<val_t1[j][i]<<" ";
      
      for(unsigned int j = 0; j < nb_locus; j++)
        FILE<<val_t2[j][i]<<" ";
      
      FILE<<endl;
    }
    FILE.close();
    
    //find min and max:
    double min1=0, max1=0, range1, min2=0, max2=0, range2;
    //    cout<<"find the bounds"<<endl;
    for(unsigned int i = 0; i < nb_locus; i++) 
      for(unsigned int j = 0; j < val_size; i++) {
        min1 = (min1 < val_t1[i][j] ? min1 : val_t1[i][j]);
        max1 = (max1 > val_t1[i][j] ? max1 : val_t1[i][j]);
        min2 = (min2 < val_t2[i][j] ? min2 : val_t2[i][j]);
        max2 = (max2 > val_t2[i][j] ? max2 : val_t2[i][j]);
      }
    
    range1 = max1 - min1;
    range2 = max2 - min2;
    
    //    cout<<"trait1: max "<<max1<<" min "<<min1<<" range "<<range1<<endl;
    //    cout<<"trait2: max "<<max2<<" min "<<min2<<" range "<<range2<<endl;
    max1 = (max1 > max2 ? max1 : max2);
    min1 = (min1 < min2 ? min1 : min2);
    range1 = max1 - min1;
    //    cout<<"total: max "<<max1<<" min "<<min1<<" range "<<range1<<endl;
    
    unsigned int dist_size1 = (unsigned int)ceil(range1 / _granularity);
    //    cout<<"allocating the dist array size "<<dist_size1<<endl;
    dist_t1[0] = new double [dist_size1];
    dist_t1[1] = new double [dist_size1];
    
    for(unsigned int i = 0; i < dist_size1; ++i) 
      dist_t1[0][i] = dist_t1[1][i] = 0;
    
    //    unsigned int dist_size2 = (unsigned int)ceil(range1 / _granularity);
    //    cout<<"allocating the dist array size "<<dist_size<<endl;
    //    dist_t2[0] = new double [dist_size2];
    //    dist_t2[1] = new double [dist_size2];    
    //    
    //    for(unsigned int i = 0; i < dist_size2; ++i) 
    //      dist_t2[0] = dist_t2[1] = 0;
    
    for(unsigned int i = 0; i < nb_locus; i++)
      for(unsigned int j = 0; j < val_size; j++) {
        val_t1[i][j] -= min1;
        dist_t1[1][ (unsigned int)floor(val_t1[i][j] / _granularity) ]++;
        
        val_t2[i][j] -= min1;
        dist_t1[1][ (unsigned int)floor(val_t2[i][j] / _granularity) ]++;
      }
    
    min1 += _granularity/2.0;
    //    min2 += _granularity/2.0;
    filename = get_filename();
    
    FILE.open(filename.c_str(), ios::out);
    
    if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
    
    FILE<<"val freq"<<endl;
    
    val_size *= 2 * nb_locus;
    for(unsigned int i = 0; i < dist_size1; ++i) {
      dist_t1[1][i] /= val_size;
      dist_t1[0][i] = min1 + i * _granularity;
      FILE<<dist_t1[0][i]<<" "<<dist_t1[1][i]<<endl;
    }
    
    FILE.close();
    
    //    cout<<"trait2 distribution:"<<endl;
    //    for(unsigned int i = 0; i < dist_size2; ++i) {
    //      dist_t2[1][i] /= val_size;
    //      dist_t2[0][i] = min2 + i * _granularity;
    //      cout<<dist_t2[0][i]<<"\t"<<dist_t2[1][i]<<endl;
    //    }
    
    for(unsigned int i = 0; i < nb_locus; i++){
      delete [] val_t1[i];
      delete [] val_t2[i];
    }
    
    delete [] val_t1;
    delete [] val_t2;
    
  } else
    warning("TTQuantiFH::FHwrite:metapop empty or age flag not set, not writing.\n");
}
