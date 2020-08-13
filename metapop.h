/** $Id: metapop.h,v 1.16 2016-09-28 14:57:03 fred Exp $
*
*  @file metapop.h
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
*  Created on @date 08.07.2004 
*
*  @author fred
*/

#ifndef __METAPOP_H
#define __METAPOP_H

#include <list>
#include <deque>
#include <map>
#include <time.h>
#include "types.h"
#include "indfactory.h"
#include "MPStatHandler.h"
#include "binarydataloader.h"
#include "MPImanager.h"
#include "tmatrix.h"
#include "output.h"
#include "filehandler.h"

class Patch;

class LifeCycleEvent;

class MPFileHandler;

//CLASS METAPOP

/**Top class of the metapopulation structure, contains the patches.

 * The basic design for the metapopulation structure is a top-down chain of responsibility where the Metapop class
 * takes care of the patches it contains which are themselves concerned by the management of their individual containers. 
 * The Individual class is only concerned by the management of its traits. Thereby, a metapopulation can be viewed
 * as an interleaving of containers where one container class takes care of its directly contained class only, 
 * without knowledge of the upward container state.

 * The Metapop class thus implements methods used to manage and get information from the patches, to manage the parameters
 * necessary to build a population and to get the population state information. It also implements the methods used to load
 * a population from different source files.
 
 * <b>Population states:</b> given by the number and the position of the individuals (both spatially in demes and temporally
 * in age class containers). The Metapop::_currentAge flag is set according to the age state of the metapopulation.
 * The life cycle events modify that state by moving individuals among individuals containers within the metatpopulation.
 * The Metapop::_currentAge flag can contain the following age class bits as defined in types.h.
 * The OFFSPRNG (=1) age class bit is set whenever the offspring containers are not empty. The ADULTS (=4) age class bit is set
 * whenever the adult containers are not empty. The POSTDISP (=2) age class bit informs about the content of the
 * post-dispersal containers. The ALL (=7) age class is the addition of the previous tags and NONE (=0) is the negation of them.
 * The individual containers are stored in the patches and handled through the Patch class interface. Each age class is represented
 * by two containers, one for the males (index 0) and the other for the females (index 1). These containers are store in a table
 * and are accessed through their age class index as defined by the age_idx enum (see types.h). These indexes are as follows:
 * The OFFSPRNG age class has index OFFSx = 0, the POSTDISP age class has index PDISPx = 1 and the ADULTS age class has index ADLTx = 2.
 *
*/
class Metapop : public StorableComponent, public SimComponent, public IndFactory
{
    friend Metapop *BinaryDataLoader::extractPop(std::string&, unsigned int, SimBuilder*, Metapop*);

private:
  MPImanager *_mpimgr;

  /**The stat handler for the population stats*/
  MPStatHandler _statHandler;

  /**The file handler used to save pedigree info.*/
  MPFileHandler *_writer;

  /**The Patch container*/
  deque< Patch* > _vPatch;
  
  //binary source loader:
  /**A BinaryDataLoader to load a population from a binary data file.*/
  BinaryDataLoader _loader;
  /**A source population as loaded from a binary data file, used to build a population.*/
  Metapop* _source;
  /**Flag to specify the loading mode, true means the source pop is the actual population.*/
  bool _source_preserve;
  /**Flage to specify that the population should be built from the binary loaded source population.*/
  bool _source_load;
  /**The number of source files, each file being a different replicate of the source population. Used to 
     automatically build the source filename.**/
  unsigned int _source_replicates;
  /**Number of digits in the replicate counter filename extension.*/
  unsigned int _source_replicate_digits;
  /**The replicate number to start loading from.*/
  unsigned int _source_start_at_replicate;
  /**The generation to load from the binary file source file*/
  unsigned int _source_generation;
  /**The base filename of source population files*/
  std::string _source_name;
  /**The trait type to load from*/
  std::string _source_filetype;
  /**The age class to fill with the source population.*/
  std::string _source_required_age;
  /**The number of replicates to source from a single source replicate.*/
  unsigned int _source_load_periodicity;
  
  /**The age class flag that is required to fill the population at the beginning of a replicate.*/
  age_t _requiredAge;
  
  //parameters:
  /**Number of patches in the population.*/
  unsigned int _patchNbr;
  /**Patch carrying capacity.*/
  unsigned int _patchK;
  /**Sex specific carrying capacities.*/
  unsigned int _patchKfem, _patchKmal;
  /**Matrix of the deme sizes, row 0 for the males, row 1 for the females.*/
  TMatrix _patchSizes;
//  /**Patch init sizes.*/
//  unsigned int *_patch_init_size;
  
  /**Number of generations to iterate.*/
  unsigned int _generations;
  /**Number of replicates to iterate.*/
  unsigned int _replicates;
  //counters:
  /**The current generation in the generation loop, starts at 1.*/
  unsigned int _currentGeneration;
  /**The current replicate in the replicate loop, starts at 1.*/
  unsigned int _currentReplicate;
//  unsigned int Current_LC_Rank;
  /**The current age class, might be changed by the LCEs.*/
  age_t        _currentAge;
  
public:	
  
  Metapop();
  virtual ~Metapop();
  /**Inits the population parameters from the ParamSet and builds the pop (adds patches), the prototypes and the life cycle.
    Called at the start of each simulation, resets the individual garbage collector.    */
  bool init();
  
  virtual bool setParameters ();
  /**Population's size parameters initializing procedure.*/
  bool setPopulationParameters();
  /**Setter for source population parameters.*/
  bool setSourceParameters();
  /**Called during simulation to change the population's parameters (temporal argument).*/
  bool updatePopulationParameters();
  /**Called to empty the patches, individuals are move to the garbage collector.*/
  void reset();
  /**Called at the end of each simulation, empties the pop and the garbage collector; the Individuals are destroyed.*/
  void clear();

  void setMPImanager(MPImanager *mgr) {_mpimgr = mgr;}
    
  ///@name Population builders
  ///@{
  /**Resets the patch container to the right number of patches as set by _patchNbr. Called at the beginning of each 
     new simulation. Extra patches are destroyed and new ones are added if missing.*/
  void resizePatchArray ();
  /**Builds the new population from parameter values. Supernumerary patches are deleted. All patches are empty.*/
  void buildPatchArray();
  /**Called during simulation to modify the meta-population size. 
     Patch number and capacities are updated, new patches are empty, existing patches containers are untouched.*/
  void updatePatchArray();
  /**Update the patch capacities and patch ID (reset to array position).*/
  void updatePatchState();
  /**Sets the deme capacity matrix from parameter values.*/
  void setPatchCapacities();
  /**Builds the new population from a single matrix of deme sizes. The sex-specific
     patch capacities will be half of the numbers given in the matrix.
   *@param param the name of the parameter to take the matrix argument from.
   **/  
  void setPatchCapacities(string param);
  /**Builds the new population from a matrix of deme sizes but for one sex only. 
   * Sizes for the other sex class(es) must be given seperately.
    *@param SEX the sex class of the given deme sizes.
    *@param param the name of the parameter to take the matrix argument from.
   **/
  void setPatchCapacities(sex_t SEX, string param);
  /**Builds the new population from matrices of deme sizes. 
    *@param paramfem the name of the parameter to take the matrix argument from.
    *@param parammal the name of the parameter to take the matrix argument from.
   **/
  void setPatchCapacities(string paramfem, string parammal);
  /**Loads a population from a soure population.*/
  void loadSourcePopulation ( );
  /**Loads the population from a binary data file when setting the first generation of a replicate. */
  void loadPopFromBinarySource ( string &filename );
  /**Loads a population from a trait's data file (text file).*/
  void loadPopFromTraitFile ( string &filename );
  /**Sets the population for the first generation of each replicates.*/
  void setPopulation (unsigned int currentReplicate, unsigned int replicates);
  void setPopulationFromSourceInPreserveMode ();
  void setPopulationFromSource () ;
  /**Fills the population of the first generation of each replicates with individuals from a population source.
    @param AGE age of the individuals to fetch from the source pop
    @param SEX sex of the individuals to fetch from the source pop
    @param src_pool the container where the source individuals will be placed, they are not removed from the source.
  */
  void fillPopulationFromSource(age_idx AGE, sex_t SEX, deque<Individual*>& src_pool);
  /**Fills a patch from a source patch loaded from a binary file, used when setting the population in preserve mode.
   @param SEX sex of the individuals to fetch from the source pop
   @param src the source patch that will be copied.
   @param patch the local patch to be filled with individuals copied from the source patch.
   @param AGE age class to copy individuals from the source.
   */
  void fillPatchFromSource(sex_t SEX, Patch* src, Patch* patch, age_t AGE);
  ///@}
  
  ///@name Implementations
  ///@{
  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader );
  
  virtual void loadStatServices ( StatServices* loader ) {loader->attach(&_statHandler);}
  
  //StoprableComponent implementation:
  virtual void store_data    ( BinaryStorageBuffer* saver  );
  
  virtual bool retrieve_data ( BinaryStorageBuffer* reader );
  ///@}
  
  /**Iterates through the individuals containers to store the trait data to a binary file.*/
  void store_trait (int trait_idx, BinaryStorageBuffer* saver);
  
  /**Iterates through the individuals containers to retrieve the trait data from a binary file.*/
  void read_trait (int trait_idx, BinaryStorageBuffer* loader);
  
  
  ///@name Getters
  ///@{
  
  /**Patch accessor, return the ith+1 patch in the metapop.*/
  Patch*       getPatch              (unsigned int i) {return (i > _vPatch.size() -1 ? 0 : _vPatch[i]);}
  
  /**A secure version of the getPatch() method.*/
  Patch*       getPatchPtr           (unsigned int patch){
    if(!(patch < _vPatch.size())) 
      fatal("Metapop::getPatchPtr()::_vPatch overflow (id=%i nb=%i)\n", patch, _vPatch.size());
    
    if (_vPatch[patch] == NULL) fatal("Metapop::getPatchPtr()::NULL ptr\n");
    
    return _vPatch[patch];
  }
  
  deque< Patch* >* getPatchArray     ( ) {return &_vPatch;}
  unsigned int getPatchArraySize     ( ) {return _vPatch.size();}
  void getAllIndividuals(age_idx AGE, deque<Individual*>& fem_pool, deque<Individual*>& mal_pool);
  void         setGenerations        (unsigned int gen) {_generations = gen;}
  unsigned int getGenerations        ( ) {return _generations;}
  void         setReplicates         (unsigned int repl) {_replicates = repl;}
  unsigned int getReplicates         ( ) {return _replicates;}
  unsigned int getPatchNbr           ( ) {return _patchNbr;}
  unsigned int getPatchKFem          ( ) {return _patchKfem;}
  unsigned int getPatchKMal          ( ) {return _patchKmal;}
  unsigned int getPatchCapacity      ( ) {return _patchK;}
  unsigned int getPatchCapacity      (sex_t SEX, unsigned int patch) {return (unsigned int)_patchSizes.get(SEX, patch);}
  TMatrix*     getPatchCapacities    ( ) {return &_patchSizes;}
  bool         isSourceLoad          ( ) {return _source_load;}
  string       getSourceName         ( ) {return _source_name;}
  string       getSourceFileType     ( ) {return _source_filetype;}
  unsigned int getSourceReplDigits   ( ) {return _source_replicate_digits;}
  
  ///@}
  
  ///@name Population state interface
  ///@{
  unsigned int getCurrentReplicate   ( ) {return _currentReplicate;}
  unsigned int getCurrentGeneration  ( ) {return _currentGeneration;}
  void setCurrentReplicate   (unsigned int repl) {_currentReplicate = repl;}
  void setCurrentGeneration  (unsigned int gen) {_currentGeneration = gen;}
  age_t        getCurrentAge         ( ) {return _currentAge;}
  
  /**Sets the age flag. 
    @param age the current age. */
  void setCurrentAge                 (age_t age) {_currentAge = age;}
  /**Set the age flag from a LifeCycleEvent object.
    @param LCE the LifeCycleEvent object. */
  void setCurrentAge                 (LifeCycleEvent* LCE) ;

  /**Checks if the population still contains at least one individual in any sex or age class.*/
  bool         isAlive               ( ) {return size() != 0;}

  /**Get the total number of individuals present in the population, all sex and age classes together.*/
  unsigned int size                  ( ) {return size(ALL);}
  
  /**Interface to get the size of a praticular age and sex class(es).
    @param AGE age class flags
    @param SEX sex class
    */
  unsigned int size ( sex_t SEX, age_t AGE );
  
  /**Interface to get the size of a praticular age class and sex class.
    @param SEX sex class
    @param IDX index of age class
    */
  unsigned int size ( sex_t SEX, age_idx IDX );
 
  /**Returns the size of the container for the appropriate age class for both sexes.
   @param IDX the index of the age class
   */
  unsigned int size (age_idx IDX);
  unsigned int size (age_idx IDX, unsigned int deme);
  unsigned int size (sex_t SEX, age_idx IDX, unsigned int deme);
  
  /**Interface to get the size of a praticular age and sex class within a patch.
    @param AGE age class flags
    @param SEX sex class
    @param deme the focal patch
    */
  unsigned int size (sex_t SEX, age_t AGE, unsigned int deme);
    
  /**Simplified interface to get the size of both sexes of the appropriate age class(es) in the whole population.
    @param AGE age class flags
  */
  unsigned int size ( age_t AGE )
  { return size( FEM, AGE ) + size( MAL, AGE );}
  
  /**Simplified interface to get the size of both sexes of the appropriate age class(es) in one patch.
    @param AGE age class flags
    @param deme the focal deme
    */
  unsigned int size ( age_t AGE, unsigned int deme )
  { return size( FEM, AGE, deme ) + size( MAL, AGE, deme );}
  
  /**Returns a pointer to the appropriate individual.
    @param SEX sex class container index
    @param AGE age class container index
    @param at the index of the individual in its container
    @param deme the patch where to grab the individual*/
  Individual* get (sex_t SEX, age_idx AGE, unsigned int at, unsigned int deme);
  
  /**Moves an individual from a deme to an other one, both demes sizes are modified.
    @param SEX sex class container index
    @param from_age age class container index in the deme of origin
    @param from_deme index of the deme of origin
    @param to_age age class container index in the destination deme
    @param to_deme index of the destination deme
    @param at index of the focal individual in the 'from' deme
    */
  void move (sex_t SEX, age_idx from_age, unsigned int from_deme, 
             age_idx to_age, unsigned int to_deme, unsigned int at);
  
  /**Removes all individual pointers and flush them into the recycling pool.
    Container sizes are reset to null values.
    @see Patch::flush()
    */  
  void flush ();
  /**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
    Container sizes are reset to null values.
    @see Patch::flush()
    @param SEX the sex class of the individual
    @param AGE the index of the age class to be flushed
    */  
  void flush (sex_t SEX, age_idx AGE);
  /**Removes all individual pointers of both sexes and specified age class and flush them into the recycling pool.
    Container sizes are reset to null values.
    @see Patch::flush()
    @param AGE the index of the age class to be flushed
    */  
  void flush (age_idx AGE);
  /**Removes all individual pointers of both sexes and specified age class(es) and flush them into the recycling pool.
    Container sizes are reset to null values.
    @see Patch::flush()
    @param AGE the age class(es) to be flushed
    */    
  void flush (age_t AGE);
  /**Removes a patch from the patch array and returns it pointer. The patch and its content are NOT deleted. 
     The IDs of the remaining patches are *NOT* updated. 
     @param i the index of the patch to remove.*/
  Patch* removePatch (unsigned int i);
  /**Removes a patch from the patch array and deletes it and its content.  
     The IDs of the remaining patches are updated. 
     @param i the index of the patch to remove.*/
  void deletePatch (unsigned int i);
  /**Adds a patch to the population. The patch is added at the end of the array 
     and it is empty. */
  void addPatch (Patch* patch);
  /**Adds num patches to the population.*/
  void addPatch (unsigned int num);
  ///@}
  
  void show_up();
};
//------------------------------------------------------------------------------------------
//
//                                        CLASS PATCH
//
//------------------------------------------------------------------------------------------

/**Second class in the metapopulation design structure, between the Metapop and Individual classes.
 * The Patch class is an abstraction of a sub-population or patch concept (also called a deme) in a metapopulation context.
 * It contains the individual containers for the different age classes. Three main age classes are currently implemented,
 * the offspring, post-dispersal and adult classes (see the age_t enum) which are all subdivided into male and female individuals.
 * These containers are accessed using the interface defined here or through the Metapop class interface.
 * The different LCEs will use these interfaces to handle the individuals. They are also responsible to change the age flag of the population (see Metapop).
 *
 * The individuals are accessed using their age index value and not the age flag value (e.g., using the Patch::get() method). 
 * These indexes are defined in the age_idx enum (see type.h) and differ from the age class flag values. For instance the adults'
 * containers have the index 2 (ADLTx = 2) whereas their age flag is 4 (ADULTS = 4). This might be confusing but it saves a lot
 * of checks at runtime! It also allows to give a flag containing several class bits set instead of a unique index value when needed
 * (see the Patch::size() and Metapop::size() suite of functions).  
*/
class Patch
{
  /**Patch ID is equal to its position in the metapop patch array.*/
  unsigned int _ID;
  /**Carrying capacity for males and females*/
  unsigned int _K;
  /**Sex specific carrying capacity*/
  unsigned int _KFem, _KMal;
  /**Extinction flag*/
  bool _isExtinct;
  /**age since last extinction.*/
  unsigned int _age;
  /**Number of age classes present.*/
  unsigned int _nb_age_class;
  /**Containers size counters, sex X age.*/
  unsigned int _sizes[2][3];
  /**Total size of the containers, amount of allocated memory.*/
  unsigned int _capacities[2][3];
  /**Individuals containers, sex X age.*/
  deque <Individual*> _containers[2][3];
  
 public:
//counters:
  unsigned short nbEmigrant, nbImigrant, nbPhilopat;
  short nbKolonisers;

//at construction, the capacities should be at least 1 to allow a patch to be filled (see add)
  Patch() : _ID(0), _K(1), _KFem(1), _KMal(1), _isExtinct(0), _age(0), _nb_age_class(3) 
  { for(unsigned int i = 0; i < _nb_age_class; i++) {
    _sizes[MAL][i] = 0;
    _sizes[FEM][i] = 0;
    _capacities[MAL][i] = 0;
    _capacities[FEM][i] = 0;}
  }
  ~Patch();
  Patch*        init                      (unsigned int nbfem, unsigned int nbmal, unsigned int id);
  ///@name Setters
  ///@{
  void          setID                     (unsigned int i) {_ID = i;}
  void          set_K                     (unsigned int k) {_K = k;}
  void          set_KFem                  (unsigned int k) {_KFem = k;}
  void          set_KMal                  (unsigned int k) {_KMal = k;}
  void          set_isExtinct             (bool status)    {_isExtinct = status;}
  void          set_age                   (unsigned int a) {_age = a;}
//  void          set_growth_rate           (double r)       {_r = r;}
  ///@}
  ///@name Getters
  ///@{
  unsigned int  getID                     ()               {return _ID;}
  unsigned int  get_K                     ()               {return _K;}
  unsigned int  get_K                     (sex_t SEX)      {return (SEX ? _KFem : _KMal);}
  unsigned int  get_KFem                  ()               {return _KFem;}
  unsigned int  get_KMal                  ()               {return _KMal;}
//  double        get_growth_rate           ()               {return _r;}
  bool          get_isExtinct             ()               {return _isExtinct;}
  unsigned int  get_age                   ()               {return _age;}
  bool          isEmpty                   ()               {return (size(ALL) == 0);}
  unsigned int  getAdultsNumber           ()               {return size(ADLTx);}
  double        getDensity                (age_idx age)    {return (double)size(age)/_K;}
  ///@}
  
  ///@name State getter and modifier functions
  ///@{
  /**Returns the size of the container of the appropriate age class(es) for both sexes.
    @param AGE the flag value of the age class
    */
  unsigned int size       (age_t AGE)
  { return size(MAL,AGE) + size(FEM,AGE); }
  
  /**Returns the size of the container for the appropriate sex and age classes present in the age flag.
    @param SEX the sex class
    @param AGE the flag value of the age class
    */
  unsigned int size       (sex_t SEX, age_t AGE)
  { 
    unsigned int mask = 1, s = 0;
    for(unsigned int i = 0; i < _nb_age_class; i++) {
      if( (mask & AGE) != 0) s +=  _sizes[SEX][i];
      mask <<= 1;
    }
    return s;
  }
  
  /**Returns the size of the container for the appropriate sex and age class.
    @param SEX the sex class
    @param AGE the index of the age class
  */
  unsigned int size       (sex_t SEX, age_idx AGE)
  { return _sizes[SEX][AGE]; }
  
  /**Returns the size of the container for the appropriate age class for both sexes.
    @param AGE the index of the age class
    */
  unsigned int size       (age_idx AGE)
  { return _sizes[0][AGE] + _sizes[1][AGE]; }
  
  /**Returns a pointer to the individual sitting at the index passed.
    \b Note: the get operations are unchecked! It's up to the user to check for overflows.
    @param SEX the sex class of the individual
    @param AGE the index of the age class
    @param at the index of the individual in the container
    */
  Individual*  get        (sex_t SEX, age_idx AGE, unsigned int at)
  { return _containers[SEX][AGE][at]; }
  
  /**Modifies the appropriate container with value of the pointer given.
    @param SEX the sex class of the individual
    @param AGE the index of the age class
    @param at the index of the individual in the container
    @param ind the pointer to the individual
    */
  void         set        (sex_t SEX, age_idx AGE, unsigned int at, Individual* ind)
  { _containers[SEX][AGE][at] = ind; }
  
  /**Adds an individual to the appropriate container, increments its size, eventually resizing it.
    @param SEX the sex class of the individual
    @param AGE the index of the age class
    @param ind the pointer to the individual
    */
  void         add        (sex_t SEX, age_idx AGE, Individual* ind)
  {
    
    if( _sizes[SEX][AGE] + 1 > _capacities[SEX][AGE] ) {
      _containers[SEX][AGE].resize( _capacities[SEX][AGE] + (_K + 1) );
      _capacities[SEX][AGE] += (_K + 1); //the +1 is here to avoid seg faults when K=0
    }
    
    _containers[SEX][AGE][ _sizes[SEX][AGE]++ ] = ind;
  }
  
  /**Assigns a new container of given size for the sex and age class passed, sets all values to NULL.*/
  void         assign     (sex_t SEX, age_idx AGE, unsigned int n)
  { _containers[SEX][AGE].assign(n,0);
   _sizes[SEX][AGE] = 0;
   _capacities[SEX][AGE] = n;
  }
  
  /**Removes the individual sitting at the given index in the appropriate container.
    @param SEX the sex class of the individual
    @param AGE the index of the age class
    @param at the index of the individual in the container
    @return pointer to the individual that has been removed
    */
  Individual*   remove     (sex_t SEX, age_idx AGE, unsigned int at)
  {
    if(_sizes[SEX][AGE] == 0) {
      error("Patch::remove:: container already empty!!");
      return NULL;
    }
    unsigned int last = _sizes[SEX][AGE] - 1;
    Individual* ind = _containers[SEX][AGE][at];
    _containers[SEX][AGE][at] = _containers[SEX][AGE][ last ];
    _containers[SEX][AGE][ last ] = 0;
    _sizes[SEX][AGE]--;
    return ind;
  }
  
  /**Moves an individual from an age class to an other one.
    \b Note: both containers are transformed by this operation. The 'from'
             container size is reduced by one while the 'to' container size
             is increased by one.
    @param SEX the sex class of the individual
    @param from the original age class of the individual
    @param to the destination age class of the individual
    @param at the index of the individual in the container
    */
  void         move       (sex_t SEX, age_idx from, age_idx to, unsigned int at)
  {
    add( SEX, to, _containers[SEX][from][at] );
    remove( SEX, from, at );
  }
  
  /**Copies all elements in the 'from' age-class container to the 'to' age-class container of the same sex.
    The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
    before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
    @param SEX the sex class of the individual
    @param from the original age class of the individual
    @param to the destination age class of the individual
  */
  void         swap       (sex_t SEX, age_idx from, age_idx to)
  {
    if( _sizes[SEX][from] > _capacities[SEX][to] ) {
      _containers[SEX][to].resize( _sizes[SEX][from] );
      _capacities[SEX][to] = _sizes[SEX][from];
    }
    
    for(unsigned int i = 0; i < _sizes[SEX][from]; ++i)
      _containers[SEX][to][i] = _containers[SEX][from][i];
    
    _sizes[SEX][to] = _sizes[SEX][from];
    clear(SEX, from);
  }
    
  /**Sets the size of the appropriate container to zero.
    \b Note: no memory operation is performed, the capacity of the container is thus not affected.
            The individual pointers are not flushed to the recycling pool, they will be overwritten by
            subsequent operations. It is thus a good idea to consider using Patch::flush to be sure no 
            pointers remained in the container.
    @see flush()
    @param SEX the sex class
    @param AGE the index of the age class
    */
  void         clear      (sex_t SEX, age_idx AGE) { _sizes[SEX][AGE] = 0;}
  void         clear      () { for(int i = 0; i < 3; i++) {_sizes[0][i] = 0;_sizes[1][i] = 0;}}

  /**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
    Container sizes are reset to null values.
    \b Note: not memory operation is performed, the total amount of memory allocated is
    left untouched.
    @param SEX the sex class of the individual
    @param AGE the index of the age class
    @param pop the pointer to the metapop for access to the recycling pool
  */    
  void         flush      (sex_t SEX, age_idx AGE, Metapop* pop)
  {
    for (unsigned int i = 0; i < _sizes[SEX][AGE]; ++i) {
      pop->recycle(_containers[SEX][AGE][i]);
      _containers[SEX][AGE][i] = 0;
    }
    _sizes[SEX][AGE] = 0;
  }
  
  void         flush      (age_idx AGE, Metapop* pop)
  { flush(FEM, AGE, pop); flush(MAL, AGE, pop); }
  
  /**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
    @param AGE an unsigned int containing the flags of the age classes to flush
    @param pop the pointer to the metapop for access to the recycling pool
    @see flush()
    */   
  void         flush      (age_t AGE, Metapop* pop)
  {
    unsigned int mask = 1;
    
    for(unsigned int i = 0; i < _nb_age_class; i++) {
      if( (mask & AGE) != 0) {
        flush(MAL, static_cast<age_idx>(i), pop);
        flush(FEM, static_cast<age_idx>(i), pop);
      }
      mask <<= 1;
    }
  }
  
  /**Removes all individual pointers of all sex and age classes and flush them into the recycling pool.*/
  void         flush      (Metapop* pop)
  {
    for(unsigned int i = 0; i < _nb_age_class; i++) {
      flush(MAL, static_cast<age_idx>(i), pop);
      flush(FEM, static_cast<age_idx>(i), pop);
    }
  }
  
  void         getCopy   (sex_t SEX, age_idx AGE, deque< Individual* >& to)
  {
    for (unsigned int i = 0; i < _sizes[SEX][AGE]; ++i) {
      to.push_back(_containers[SEX][AGE][i]);
    }
  }
  
  void       copy2patch (sex_t from_sex, sex_t to_sex, age_idx from_age, age_idx to_age, Patch* to_patch)
  {
    for (unsigned int i = 0; i < _sizes[from_sex][from_age]; ++i)
      to_patch->add(to_sex, to_age, _containers[from_sex][from_age][i] );
  }
  
  void       copy2patch (sex_t SEX, age_idx AGE, Patch* patch)
  {
    for (unsigned int i = 0; i < _sizes[SEX][AGE]; ++i)
      patch->add(SEX, AGE, get(SEX, AGE, i) );
  }
    
  void       copy2patch (age_idx AGE, Patch* patch)
  {
    copy2patch(FEM, AGE, patch);
    copy2patch(MAL, AGE, patch);
  }
  
  void       copy2patch (Patch* patch)
  {
    for (unsigned int i = 0; i < _nb_age_class; ++i){
      copy2patch(FEM, static_cast<age_idx> (i), patch);
      copy2patch(MAL, static_cast<age_idx> (i), patch);
    }
  }
  ///@}
  
  void reset_counters();
  void reset_containers();
  /**Fills the patch containers corresponding to the age flags passed, for both sexes.*/
  void setNewGeneration(age_t AGE, Metapop* pop);
  /**Fills the patch container corresponding to the age class index passed, for both sexes.*/
  void setNewGeneration(age_idx AGE, Metapop* pop);
  
  void show_up();
  
};

//------------------------------------------------------------------------------------------
//
//                                       MPFileHandler
//
//------------------------------------------------------------------------------------------
class MPFileHandler : public FileHandler {

	int _patch_sample_size;

public:
	MPFileHandler () : FileHandler(".ped"), _patch_sample_size(0) {}

	virtual ~MPFileHandler() {}

	void setOption(int size){_patch_sample_size = size;}

	virtual void FHwrite();
	virtual void FHread (string& filename){}

	void createAndPrintSample (age_idx AGE, Patch* patch, ofstream& FH);
	void printNoSample (sex_t SEX, age_idx AGE, Patch* patch, ofstream& FH);
};

//------------------------------------------------------------------------------------------
//
//                                  Metapop inline functions
//
//------------------------------------------------------------------------------------------
inline unsigned int Metapop::size (age_idx IDX)
{ 
  return size(FEM, IDX) + size(MAL, IDX); 
}

inline unsigned int Metapop::size (age_idx IDX, unsigned int deme)
{ 
  return size(FEM, IDX, deme) + size(MAL, IDX, deme); 
}

inline unsigned int Metapop::size ( sex_t SEX, age_idx IDX )
{
  unsigned int s = 0;
  for(unsigned int i = 0; i < _patchNbr; i++)
    s += _vPatch[i]->size(SEX, IDX);
  return s;
}

inline unsigned int Metapop::size ( sex_t SEX, age_idx IDX, unsigned int deme)
{
  return _vPatch[deme]->size(SEX, IDX);
}

inline unsigned int Metapop::size ( sex_t SEX, age_t AGE )
{
  unsigned int s = 0;
  for(unsigned int i = 0; i < _vPatch.size(); i++)
    s += _vPatch[i]->size(SEX, AGE);
  return s;
}

inline unsigned int Metapop::size (sex_t SEX, age_t AGE, unsigned int deme)
{
  Patch* patch = getPatch(deme);
  return (patch!=0? patch->size(SEX, AGE) : 0);
}  

inline Individual* Metapop::get (sex_t SEX, age_idx AGE, unsigned int at, unsigned int deme)
{ return getPatchPtr(deme)->get(SEX, AGE, at); }

inline void Metapop::move (sex_t SEX, age_idx from_age, unsigned int from_deme, 
                           age_idx to_age, unsigned int to_deme, unsigned int at)
{//cout << "  add "<<get(SEX, from_age, at, from_deme)->getID()<<" to "<<to_deme<<endl;
  _vPatch[to_deme]->add( SEX, to_age, get(SEX, from_age, at, from_deme));
  //cout << "  remove "<<get(SEX, from_age, at, from_deme)->getID()<<" from "<<from_deme<<endl;
  _vPatch[from_deme]->remove(SEX, from_age, at);
  //cout << "  sizes: "<<from_deme<<": "<<_vPatch[from_deme]->size(SEX, from_age)<<", "
  //     << to_deme << ": " <<_vPatch[to_deme]->size( SEX, to_age)<<endl;
}

inline void Metapop::flush() 
{
  for(unsigned int i = 0; i < _patchNbr; i++) _vPatch[i]->flush(this);
}

inline void Metapop::flush(sex_t SEX, age_idx AGE) 
{
  for(unsigned int i = 0; i < _patchNbr; i++) _vPatch[i]->flush(SEX, AGE, this);
}

inline void Metapop::flush(age_idx AGE) 
{
  for(unsigned int i = 0; i < _patchNbr; i++) _vPatch[i]->flush(AGE, this);
}

inline void Metapop::flush(age_t AGE) 
{
  for(unsigned int i = 0; i < _patchNbr; i++) _vPatch[i]->flush(AGE, this);
}

inline Patch* Metapop::removePatch (unsigned int i)
{
  Patch* patch = getPatchPtr(i);
  _vPatch.erase(_vPatch.begin() + i);
  return patch;
}

inline void Metapop::deletePatch (unsigned int i)
{
  delete _vPatch[i];
  for (unsigned int k = i; k < _vPatch.size() -1; k++) {
    _vPatch[k] = _vPatch[k + 1];
    _vPatch[k]->setID(k);
  }
  _vPatch.pop_back();
}

inline void Metapop::addPatch (Patch* patch)
{
  _vPatch.push_back(patch);
}

inline void Metapop::addPatch (unsigned int num)
{
  for (unsigned int i = 0; i < num; i++)
    _vPatch.push_back(new Patch());
}




#endif
