/** $Id: binarydatasaver.cc,v 1.11 2015-07-13 08:52:56 fred Exp $
 *
 *  @file binarydatasaver.cc
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
 *  Created on @date 18.10.2005
 *  @author fred
 */

#include <dirent.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <string.h>
#include <fcntl.h>
#include <sstream>
#include "metapop.h"
#include "individual.h"
#include "ttrait.h"
#include "binarydatasaver.h"
#include "output.h"
#include "version.h"
#include "simenv.h"

extern MPIenv *_myenv;
pid_t BinaryDataSaver::PID = 0;

// ----------------------------------------------------------------------------------------
// BinaryDataSaver::cstor
// ----------------------------------------------------------------------------------------
BinaryDataSaver::BinaryDataSaver() 
: LifeCycleEvent("store",""), FileHandler(".bin"), _isPeriodic(0)
{
  ParamUpdater< BinaryDataSaver > * updater =
  new ParamUpdater< BinaryDataSaver > (&BinaryDataSaver::setParameters);
  
  add_parameter("store_generation",INT,true,false,0,0,updater);
  add_parameter("store_recursive",BOOL,false,false,0,0,updater);
  add_parameter("store_dir",STR,false,false,0,0,updater);
  add_parameter("store_nocompress",BOOL,false,false,0,0,updater);
  add_parameter("store_noarchive",BOOL,false,false,0,0,updater);
  add_parameter("store_compress_cmde",STR,false,false,0,0,updater);
  add_parameter("store_compress_extension",STR,false,false,0,0,updater);
  add_parameter("store_uncompress_cmde",STR,false,false,0,0,updater);
  add_parameter("store_archive_cmde",STR,false,false,0,0,updater);
  add_parameter("store_archive_extension",STR,false,false,0,0,updater);
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::setParameters
// ----------------------------------------------------------------------------------------
bool BinaryDataSaver::setParameters()
{
  Param* param;
  
  _generation = (unsigned int)get_parameter_value("store_generation");
  
  //make sure we at least save the last generation:
  if(_generation > _popPtr->getGenerations()) _generation = _popPtr->getGenerations();
  
  param = get_parameter("store_dir");
  if(param->isSet())
    _dir = param->getArg();
  else
    _dir = "";
  
  param = get_parameter("store_recursive");
  if(param->isSet())
    _isPeriodic = true;
  else
    _isPeriodic = false;
  
  param = get_parameter("store_compress_cmde");
  if(param->isSet())
    _comp_cmd = param->getArg();
  else
    _comp_cmd = "bzip2"; //default compressor used
  
  param = get_parameter("store_compress_extension");
  if(param->isSet())
    _comp_ext = param->getArg();
  else
    _comp_ext = ".bz2"; //default
  
  param = get_parameter("store_uncompress_cmde");
  if(param->isSet())
    _uncomp_cmd = param->getArg();
  else
    _uncomp_cmd = "unbzip2"; //default
  
  param = get_parameter("store_archive_cmde");
  if(param->isSet())
    _tar_cmd = param->getArg();
  else
    _tar_cmd = "tar"; //default
  
  param = get_parameter("store_archive_extension");
  if(param->isSet())
    _tar_ext = param->getArg();
  else
    _tar_ext = ".tar"; //default
  
  param = get_parameter("store_nocompress");
  if(param->isSet()) {_comp_cmd = ""; _comp_ext ="";} //compression process desabled
  
  param = get_parameter("store_noarchive");
  if(param->isSet()) {_tar_cmd = ""; _tar_ext = "";} //archiving processed desabled
  
  //set the file handler, will be called by FileServices, will output data at last gen
  FileHandler::set(true, false, 1, _popPtr->getGenerations(), get_rank(), _dir);
  
  _offset_table.clear();
    
  return true;
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::execute
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::execute()
{ 
  //set the buffer ready to be filled
  if(_popPtr->getCurrentGeneration() == 1)
    _buff.set_buff();
  
  //if the generation to store is > than the total nbr of generations, correct this and store the last one
  if(!_isPeriodic && _generation == _popPtr->getCurrentGeneration()) {
    storeData();
  } else
  //if periodic, and generation match, store the data in the buffer
  if( _isPeriodic && !(_popPtr->getCurrentGeneration() % _generation)){
    storeData();
  }  
  
  //if reached last generation, write data to file, 
  //this empties the storage buffer with a call to printData() within FHwrite()
  //thus prevents further calls to FHwrite to rewrite data to disc
  if (_popPtr->getCurrentGeneration() == _popPtr->getGenerations()) {
    //make sure the last generation is always stored before writing data to disc:
    //catch if tot num gene is not a multiple of _generation
    if( _popPtr->getCurrentGeneration() % _generation || 
       //catch if tot num gen is a multiple but _generation < num gen and non periodic
       (!_isPeriodic && _generation != _popPtr->getGenerations() 
        && !(_popPtr->getCurrentGeneration() % _generation))) 
    {
      storeData();
    }
    printHeader();
    printData(); //this calls _buff.clear()
    finish();
  }
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::FHwrite
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::FHwrite()
{
  //need to make sure data is saved even though last generation has not been reached
  //in case of population extinction, the file manager issues a notify() call after
  //the simulation exited from the generations loop, we try to catch this call here

  //keep in mind that this function is always called at the last generation

  //if data is already written to disc, the data buffer is empty
  //if the 'store' LCE is after the 'save_files' LCE in the life cycle, we don't write
  //--> issues two calls to printData, this one would be before storing last generation
  //if 'store' is before 'save_files', the data buffer should be empty by now
  
  if(_buff.getByteLength() > 0 && SIMenv::getCurrentRankInLifeCycle() >= LifeCycleEvent::get_rank()) {
    printHeader();
    printData(); //this calls _buff.clear()
    finish();
  }

}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printHeader
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::printHeader()
{
  //first, record the parameters in text mode:
  ofstream FILE (get_filename().c_str(), ios::out);
  
  if(!FILE) fatal("could not open Binary output file!!\n");
  
  list< ParamSet* > current_params = get_service()->get_params();
  list< ParamSet* >::iterator Pit = current_params.begin();
  
  FILE<<"#NEMO "<<MAIN_VERSION<<" "<<MINOR_VERSION<<" "<<REVISION<<" "<<RELEASE<<" "<<VERSION_DATE<<endl;
  
#ifdef _DEBUG_
  message("BinaryDataSaver::printHeader:storing parameters\n");
#endif
  while(Pit != current_params.end()) {
    (*Pit)->print(FILE);
    Pit++;
  }
  FILE.close();
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::storeData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::storeData()
{
#ifdef _DEBUG_
  message("BinaryDataSaver::storeData\n");
  unsigned int byte_count = _buff.getByteLength();
#endif
  unsigned int  generation = _popPtr->getCurrentGeneration();
  unsigned char separator[2] = {'@','G'}; //generation separator = '@G'
  //store the position in the buffer to the offset table:
  _offset_table[generation] = _buff.getByteLength();
  //store the data, begin with generation separator and number:
  _buff.store(&separator, 2 * sizeof(unsigned char));
  
  _buff.store(&generation, sizeof(unsigned int));
  
  //Metapop is a StorableComponent, call store_data:
  //stores all individual info, including trait sequences
  _popPtr->store_data(&_buff);
  
#ifdef _DEBUG_
  message("BinaryDataSaver::storeData::stored %ikB\n",(_buff.getByteLength()-byte_count)/1024);
  byte_count = _buff.getByteLength();
#endif
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::printData()
{
  
  mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  int flag = O_WRONLY | O_APPEND;
  //open in write and append mode with permission flag set to rw-r--r--
  int fdesc = open(get_filename().c_str(), flag, mode);
  
  if(fdesc == -1) fatal("BinaryDataSaver::printData::open %s:%s\n",get_filename().c_str(),strerror(errno));
  //get the offset position of the EOF, will be the nbr of bytes to add to the generation offsets
  off_t pos = lseek(fdesc,0,SEEK_END);
  
  if(pos == -1) fatal("BinaryDataSaver::printData::lseek %s\n",strerror(errno));
  
  //offset of the offset table:
  int off_table = _buff.getByteLength() + (int)pos;
  //add the offset info at the end of the _buff:
  unsigned char separator[3] = {'@','O','T'}; 
  _buff.store(&separator, 3);
  
  std::map<unsigned int, unsigned int>::iterator IT = _offset_table.begin();
  while(IT != _offset_table.end()) {
    IT->second += (int)pos;
    _buff.store((int*)&IT->first, sizeof(int));
    _buff.store((int*)&IT->second, sizeof(int));
    IT++;
  }
  //finally, record the number of generations recorded and the offset of start of the offset table
  int nb_recgen = _offset_table.size();
  _buff.store(&nb_recgen, sizeof(int));
  _buff.store(&off_table, sizeof(int));
  
  //now write the data to the file:
#ifdef _DEBUG_
  message("BinaryDataSaver::printData:writing %ikB of data ",_buff.getByteLength()/1024);
#endif
  
  if((write(fdesc,_buff.getBuffer(),_buff.getByteLength())) == -1)
    fatal("BinaryDataSaver::printData::write %s\n",strerror(errno));
  
  if((close(fdesc)) == -1)
    error("BinaryDataSaver::printData::close %s\n",strerror(errno));

  //empty the buffer, get ready for next record
  _buff.clear();
  
#ifdef _DEBUG_
  message(" [ok]\n");
#endif  
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::finish
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::finish ()
{
  bool doComp = (_comp_cmd.size() > 0), doTar = (_tar_cmd.size() > 0),
  first = (_popPtr->getCurrentReplicate() 
           == max((unsigned)1,_myenv->slaveRank()));
  
  if (!(doComp || doTar)) return;
  
  stringstream sysCmd;
  
  if (doComp) 
    sysCmd << _comp_cmd << " " << get_filename() << ";";
  
  if (doTar) {
    sysCmd << _tar_cmd;
    if (first) sysCmd << " c"; //first replicate, create the tar archive:
    else 	   sysCmd << " r"; //next replicates, append files to it:
    sysCmd << "f " << get_path() << get_service()->getBaseFileName();
    if (!_myenv->isMaster()) sysCmd << _myenv->slaveRank();
    sysCmd << _tar_ext << " --remove-files " << get_filename() << _comp_ext;
  }
  

  //wait for last child:
  if (!first) waitpid(PID,NULL,0);
  
  pid_t tPID = fork();
  
  if (tPID == 0) {
    
    if (system(sysCmd.str().c_str()) < 0)
      error("BinaryDataSaver::finish system cmd \"%s\" failed: %s\n",
            sysCmd.str().c_str(),strerror(errno));
  
    _exit(EXIT_FAILURE);
  
  } else if (tPID < 0) //fork failed :
    error("BinaryDataSaver::finish::could not fork new process: %s\n",
          strerror(errno));
  
  PID = tPID;
}
