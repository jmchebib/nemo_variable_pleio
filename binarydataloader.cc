/** $Id: binarydataloader.cc,v 1.11 2015-07-13 08:52:59 fred Exp $
*
*  @file binarydataloader.cc
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
*  Created on @date 19.10.2005
*  @author fred
*/

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <functional>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/mman.h>
#include <assert.h>
#include "tstring.h"
#include "Uniform.h"
#include "binarydataloader.h"
#include "basicsimulation.h"
#include "metapop.h"

BinaryDataLoader::~BinaryDataLoader()
{
  if(_in_pop != NULL) delete _in_pop;
  if(_new_sim != NULL) delete _new_sim;
}
void BinaryDataLoader::clear ( )
{
  if(_in_pop != NULL) delete _in_pop;
  if(_new_sim != NULL) delete _new_sim;
  _in_pop = NULL;
  _new_sim = NULL;
  _buff.clear();
  
}
int BinaryDataLoader::extractOffsetTable (int FD) 
{
  int rout, nb_recgen = 0, off_table = 0;
  lseek(FD,0,SEEK_END);
  off_t current_pos = 0L;
  
  current_pos = lseek(FD,0,SEEK_CUR);
//#ifdef _DEBUG_ 
//  message("eof at: %i\n",current_pos);
//  message("position at the end of the file, at the position to read table info: ");
//#endif
  if((current_pos = lseek(FD,current_pos -2*sizeof(int),SEEK_SET)) == -1) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::lseek failed on %s\n",strerror(errno));
  
//#ifdef _DEBUG_ 
//  message("now at: %i\n",current_pos);
//  message("get the number of generations stored in the file: "); 
//#endif
  unsigned int tot=0;
  do{
  if((rout = read(FD,&nb_recgen,sizeof(int))) == -1)// || rout != sizeof(int)) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::read failed on %s\n",strerror(errno));
    tot+=rout;
    }while(tot < sizeof(int) && rout != 0);
//#ifdef _DEBUG_ 
//  current_pos = lseek(FD,0,SEEK_CUR);
//  message("%i, bytes read:%i, current offset %i\n",nb_recgen, rout, current_pos);
//  message("get offset of the beginning of the offset table: "); 
//#endif
  if((rout = read(FD,&off_table,sizeof(int))) == -1)// || rout != sizeof(int)) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::read failed on %s\n",strerror(errno));
//#ifdef _DEBUG_ 
//  current_pos = lseek(FD,0,SEEK_CUR);
//  message("%i, bytes read:%i, current offset %i\n",off_table, rout, current_pos);
//  message("position at the beginning of the offset table:\n"); 
//#endif
  if((current_pos = lseek(FD,off_table,SEEK_SET)) == -1) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::lseek failed on %s",strerror(errno));
 
//#ifdef _DEBUG_ 
//  message("now at: %i\n",current_pos);
//  message("check if we are correctly positionned, get the separator: "); 
//#endif
  unsigned char separator[3];
  if((rout = read(FD, &separator, 3)) == -1)// || rout != 3) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::read failed on %s\n",strerror(errno));
  
  if( separator[0] != '@' || separator[1] != 'O' || separator[2] != 'T')
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable:: wrong separator\n");
//#ifdef _DEBUG_ 
//  message(" ok\n read the table elements:\n"); 
//#endif
  int table_elt[2];
  _offset_table.clear();
  for(int i = 0; i < nb_recgen; ++i) {
    
    if((rout = read(FD,&table_elt,2*sizeof(int))) == -1 || rout != 2*sizeof(int)) 
      fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractOffsetTable::read failed on %s\n",strerror(errno));
    
    _offset_table[table_elt[0]] = table_elt[1];
//#ifdef _DEBUG_ 
//    message("gen %i at %i\n",table_elt[0], table_elt[1]); 
//#endif
  }
//#ifdef _DEBUG_ 
//    message("returning %i\n",off_table); 
//#endif
  return off_table;
}

Metapop* BinaryDataLoader::extractPop (string& file, unsigned int generation, SimBuilder* sim, Metapop* popPtr)
{
  int FD, gen_length, gen_offset=0, otab_offset;
  string cmd, magic_name, old_name;
  bool do_compress = 0;
  _filename = file;
  _gen = generation;
  _pExtractor.setName(file.c_str());
  
  if(_in_pop != NULL) delete _in_pop;
  _in_pop = new Metapop();
  _in_pop->setMPImanager( popPtr->_mpimgr );
  
  //create new sim builder from copy of existing one, will copy templates
  if(_new_sim != 0) delete _new_sim;
  _new_sim = new SimBuilder(*sim);
  _current_sim = sim;

  //add the current metapop paramSet:
  _new_sim->add_paramset(_in_pop->get_paramset());

  if((FD = open(_filename.c_str(),O_RDONLY)) == -1){
    close(FD);
    //check if we have to uncompress it:
    string alt_name = _filename + ".bz2";
    magic_name = _filename + tstring::int2str(RAND::Uniform(5477871));
    int status;
    if((FD = open(alt_name.c_str(),O_RDONLY)) == -1){
      close(FD);
      //check if gziped:
      alt_name = _filename + ".gz";
      if((FD = open(alt_name.c_str(),O_RDONLY)) == -1){
        close(FD);
        fatal("BinaryDataLoader::extractPop::open failed on %s: %s\n",_filename.c_str(),strerror(errno));
      } else {
        //try to ungzip it:
        cmd = "gzip -cd " + alt_name + " > " + magic_name;

      }
    }
    close(FD);
    //try to bunzip2:
    cmd = "bunzip2 -ck " + alt_name + " > " + magic_name;

    status = system(cmd.c_str());
    if(status != 0)
      fatal("BinaryDataLoader::extractPop::bunzip2 failed on %s\n", alt_name.c_str());
    
    if((FD = open(magic_name.c_str(),O_RDONLY)) == -1) {//this should work now...
      close(FD);
      fatal("BinaryDataLoader::extractPop::open failed on %s: %s\n", _filename.c_str(), strerror(errno));
    }
    old_name = _filename;
    _filename = magic_name;
    _pExtractor.setName(magic_name.c_str());
    do_compress = 1; //to remove the duplicated source file
  }
  
  //1. read the offset table:
  otab_offset = extractOffsetTable(FD);

  //2. read the params and build the prototypes
  //extract params from binary file and set the sim params
  // !!! the param values in init part might not reflect pop state, esp if temporal arguments used !!!
  if( !(_new_sim->build_currentParams(_pExtractor.getParameters(NULL))) ){
    error("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop::could not set parameters from binary file\n");
    return NULL;
  }
  //set Metapop params
  //build the list of the selected trait templates, will init() the prototypes
  //and load them into the pop, build the individual prototype
//  if( !(_in_pop->init(_new_sim->build_currentTraits(),_new_sim->build_currentLifeCycle())) )
//    return NULL;
  if(!_in_pop->setParameters()) return NULL;  
  _in_pop->buildPatchArray();
  _in_pop->makePrototype(_new_sim->build_currentTraits());
  
  //3. check prototypes coherence with existing pop?
  /* check prototypes coherence with existing pop, might have to add or remove traits? */
  
  //4. load the pop
  //compute the number of bytes to read:
  std::map<unsigned int, unsigned int>::iterator OT = _offset_table.end(); 
  //starting point:
  if(generation == 0) {
    OT--;
    gen_offset = OT->second; //load last generation recorded
    _gen = OT->first;
  } else {
    if( (OT = _offset_table.find(generation)) == _offset_table.end() )
      fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop::generation %i not found in file \"%s\"\n",generation,_filename.c_str());
    else
       gen_offset = OT->second;
  }
//#ifdef _DEBUG_ 
//  message("\ngeneration offset is: %i\notab offset is: %i\n",gen_offset,otab_offset); 
//#endif
  //length
  if(_offset_table.size() == 1 || (++OT) == _offset_table.end())
    gen_length = otab_offset - gen_offset;
  else
    gen_length = OT->second - gen_offset;
  
  assert(gen_length > 0);

  //+ 2 bytes to get the next separator:
  gen_length += 2;

//#ifdef _DEBUG_ 
//  message("nb bytes to mmap are: %i\n",gen_length); 
//#endif
  off_t current_pos = lseek(FD,gen_offset,SEEK_SET);
  
  if(current_pos == -1) 
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop::lseek failed on %s\n",strerror(errno));  
  
  int rout = -1, rcount = 0, rest = gen_length;
  char *data = new char [gen_length];
  
  while(rest != 0 && rout != 0) {
    if( (rout = read(FD,&data[rcount],rest)) == -1) fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop::read data %s\n",strerror(errno));
    rcount += rout;
    rest -= rout;
#ifdef _DEBUG_ 
    message("BinaryDataLoader::extractPop:read %i bytes from file\n",rout);
#endif
  }
  
  close(FD);
  
  _buff.set_buff(data, gen_length);
  
  delete [] data;
  
  unsigned char separator[2];
  
  _buff.BSBread(&separator, 2 * sizeof(unsigned char));
  
  if( separator[0] != '@' || separator[1] != 'G')
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop:: wrong generation separator\n");
  
  _buff.BSBread(&generation, sizeof(unsigned int));
  if(generation != _gen)
    fatal("Binary file appears corrupted:\n >>>> BinaryDataLoader::extractPop:: wrong generation in file\n");
  
  bool status = true;
  //retrieve data from the buffer:
  status = _in_pop->retrieve_data(&_buff);
  
  
  if(do_compress){
    cmd = "rm -f " + magic_name;
    if( system(cmd.c_str()) == 1)
      warning("BinaryDataLoader::extractPop:: deleting duplicated source file failed on %s\n", _filename.c_str());
  }
  
  delete _new_sim;
  _new_sim = NULL;
  
  if( !status ) return NULL;

  else return _in_pop;

}
