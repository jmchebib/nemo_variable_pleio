/**  $Id: ttrait_with_map.h,v 1.12 2016-09-28 15:23:53 fred Exp $
 *
 *  @file ttrait_with_map.h
 *  Nemo2
 *
 *   Copyright (C) 2011 Frederic Guillaume
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
 *  created on @date 15.12.2011
 *  @author fred
 */


#ifndef TTRAIT_WITH_MAP_H
#define TTRAIT_WITH_MAP_H

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "ttrait.h"


// ------------------------------------------------------------------------------
/**
 *  GeneticMap
 */
// ------------------------------------------------------------------------------

class GeneticMap {
  
private:
  
  unsigned long _currentIndividual;
  
  /** Table mapping trait type to its position index in the following tables.*/
  map< trait_t, unsigned int > _traits;
  
  /** Number of traits registered in the map. Length of following arrays.*/
  unsigned int _nTrait;
  
  /** A list of tables that map the map position (cM) to a locus, for each trait. 
   The length of the table is the length of the genetic map, which depends on the 
   map resolution (cM by default).*/
  vector< unsigned int* >  _lociLookupTable; //this is the actual map, for each trait
  /** Vector of number of chromosomes for each trait.*/
  vector< unsigned int >   _numChrsmPerTrait;
  /** Vector of number of loci for each trait.*/
  vector< unsigned int >   _numLociPerTrait;
  /** Vector containing a table of number of loci per chromosome for each trait.*/
  vector< unsigned int* >  _numLociPerChrsmPerTrait;
  /** Vector containing the table of map position for the loci of each trait.*/
  vector< unsigned int* >  _locPositionsPerTrait;
  
  /**Vector of tables containing, for each trait, the locus number at which x-overs happen.
   Updated at each generation. There is matrix per sex/gamete with one trait per row and a variable
   number of elements in each row.*/
  vector< vector < unsigned int > > _recPositions[2];
  
  /**Two vectors holding the starting copy of each chromosome to use when creating the two gametes
   that are used to create a new individual.*/
  vector< bool >  _chrsmFirstRecombPosition[2];
  
  unsigned int  _numChromosome;
  unsigned int  *_perChrsmLength;          //array length is _numChromosome
  unsigned int  *_chrsmFirstLocusPosition; //array length is _numChromosome
  unsigned int  _totalLength;  // total map length
  unsigned int  _recombLength; // this is _totalLength + 1; used to draw the x-over spots
  unsigned int  _totalNumLoci;
  double        _resolution;
  double        _totRecombEventsMean;
  
  
public:
  
  GeneticMap() :_currentIndividual(0), _nTrait(0),
  _numChromosome(0), _perChrsmLength(0), 
  _chrsmFirstLocusPosition(0), _totalLength(0), _totalNumLoci(0), _resolution(1),
  _totRecombEventsMean(0) {}
  
  ~GeneticMap() {reset_tables();}
  
  bool getGeneticMap (trait_t trait, double** table, unsigned int table_length);

  double getResolution ( ) {return _resolution;}
  
  double setResolution (double val) 
  {
    _resolution = (val < _resolution ? val : _resolution); 
    return _resolution;
  }
  
  void   rescaleMap (double val);
  void   reset_tables ();
  void   clear ();
  
  /** Bbuilds the lookup table for each trait. A lookup table maps a chromosomal position to a locus
   so that when a x-over is placed at a given map position, the corresponding locus can be directly found.
   The size of the lookup table depends on the _totalLength and _resolution of the map. The lookup table
   for a 1M map at the 0.01 cM scale will have 10,000 elements.
   @param idx the index of the trait in the table (stored in the trait prototype)*/
  void   setLookupTable (unsigned int idx);
  
  /** Called by TTProtoWithMap::recombine twice to create the two gametes necessary for the creation
   of a new individual. This function sets the location of the x-over on the chromosome map and 
   records the loci positions at which x-over happen for each trait separately. The positions are 
   stored in GeneticMap::_recPositions, which are then passed to the trait-specific inheritance
   functions when setting up the genetics of the traits in new individuals.
   @param SEX the origin of the gamete, whether from the father or the mother
   */
  void   recombine  (sex_t SEX);

  /** Called by TTProtoWithMap::recombine with individual ID passed down from
   Individual::recombine. Returns false when already called by the same individual.
   The function is called for each trait separately but recombination must be computed once
   per individual only. It is so because some functions create new individual's trait separately, 
   for instance in breed_selection when traits under viability selection are created before neutral
   traits.
   @param ID the ID number of the individual calling the recombination function*/
  bool   registerIndForRecombine (unsigned long ID) 
  {
    if (ID == _currentIndividual) return false;
    else  _currentIndividual = ID;
    return true;
  }
  
  /** Returns the table index for the registered trait. */
  unsigned int addTrait (trait_t trait, unsigned int nChrm, unsigned int nLoc, unsigned int* nLocChrm,
                         double resolution, unsigned int* locPositions);
  
  void unregisterTrait (trait_t trait);
  
//  void unregisterTrait_at (unsigned int traitIdx);
  
  /** Returns a vector of the loci where crossing-overs take place.*/
  vector< unsigned int>&  getRecLoci (sex_t SEX, unsigned int trait)
  {
    return _recPositions[SEX][trait];
  }
  
  /** Returns the vector of the first chromosome position for recombination, used for all traits.*/
  vector< bool >& getFirstRecPosition (sex_t SEX) 
  { return _chrsmFirstRecombPosition[SEX]; }
};


// ------------------------------------------------------------------------------
/**
 *  TTProtoWithMap
 */
// ------------------------------------------------------------------------------
class TTProtoWithMap : public TraitPrototype {
  
  string _paramPrefix;
  
  bool _isRegistered;
  
protected:
  
  //recombination:
  unsigned int _mapIndex;
  double       _totRecombEventsMean;
  double       _recombRate;
  double       _mapResolution;
  unsigned int _numChromosome;
  unsigned int _numLoci;
  double       *_recombRatePerChrmsm;
  unsigned int *_numLociPerChrmsm;
  unsigned int *_chrsmLength;
  unsigned int *_lociMapPositions;
  
  friend class TTraitWithMap;
  
public:
  
  
  static GeneticMap   _map;
  
  
  TTProtoWithMap():_mapIndex(0), _isRegistered(0),_totRecombEventsMean(0), 
  _recombRate(0), _mapResolution(1), _numChromosome(0), _numLoci(0),
  _recombRatePerChrmsm(0), _numLociPerChrmsm(0), _chrsmLength(0),
  _lociMapPositions(0) 
  {}
  
  TTProtoWithMap(const TTProtoWithMap& TP);
  
  virtual ~TTProtoWithMap(); 
  
  void setMapIndex (unsigned int idx) {_mapIndex = idx;}
  
  unsigned int getMapIndex () {return _mapIndex;}
  
  bool setGeneticMapParameters (string prefix);

  void addGeneticMapParameters (string prefix);
  
  bool setRecombinationMapRandom ();
  
  bool setRecombinationMapNonRandom (vector< vector<double> >* lociPositions);
  
  bool setRecombinationMapFixed ();
  
  bool setNumLociPerChromosome (string param_name);
  
  void reset_recombination_pointers();
  
  void registerGeneticMap ();
  
  void recordRandomMap ();

  static void recombine (unsigned long indID);
  
  virtual void reset ();
  
};


// ------------------------------------------------------------------------------
/**
 *  TTraitWithMap
 */
// ------------------------------------------------------------------------------
class TTraitWithMap : public TTrait {
  
protected:
  
  TTProtoWithMap* _myProto;
  
public:
  
  TTraitWithMap():_myProto(0) {}
  
  virtual ~TTraitWithMap() {}
  
};

#endif
