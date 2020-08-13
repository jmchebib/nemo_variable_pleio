/** $Id: binarydataloader.h,v 1.8 2015-07-13 08:52:58 fred Exp $
*
*  @file binarydataloader.h
*  Nemo2
*
*  Copyright (C) 2006-2015 Frederic Guillaume
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
*  Created on @date 19.10.2005
*  @author fred
*/
#ifndef BINARYDATALOADER_H
#define BINARYDATALOADER_H

#include <string>
#include <map>
#include "binarystoragebuffer.h"
#include "fileparser.h"
#include "MPImanager.h"

class SimBuilder;
class Metapop;

/**A class to load a whole population from a binary data file.*/
class BinaryDataLoader {
  
  SimBuilder* _current_sim;
  
  BinaryStorageBuffer _buff;
  
  Metapop *_in_pop;
  
  SimBuilder* _new_sim;
  
  unsigned int _gen;
  
  std::string _filename;
  
  std::map<unsigned int, unsigned int> _offset_table;

  BinaryFileParser _pExtractor;
  
public:
    
  BinaryDataLoader() : _current_sim(0),_buff(),_in_pop(0),
                       _new_sim(0),_gen(0),_filename(),_pExtractor(0)
    { }
  
  ~BinaryDataLoader();
 
  int extractOffsetTable (int FD);
  
  Metapop* extractPop (string& filename, unsigned int generation, SimBuilder* sim, Metapop* popPtr);
    
  Metapop* getPop ( ) const {return _in_pop;}
  
  unsigned int getGeneration ( ) const {return _gen;}
  
  void setGeneration (unsigned int gen) {_gen = gen;}
  
  const BinaryStorageBuffer* getBuffer ( ) const {return &_buff;}
  
  void clear ( );
};
#endif
