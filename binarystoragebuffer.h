/** $Id: binarystoragebuffer.h,v 1.6 2015-07-13 08:52:57 fred Exp $
*
*  @file binarystoragebuffer.h
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
*  Created on @date 20.10.2005
*
*  @author fred
*/
#ifndef BINARYSTORAGEBUFFER_H
#define BINARYSTORAGEBUFFER_H

#include <iostream>
#include <string.h>
#include "output.h"

#define MAX_BUFF 1024000

/**A class to store any kind of data in a char buffer before unloading it in a binary data file.*/

class BinaryStorageBuffer {
  
private:
  char* _buff;
  unsigned int _len, _bytes_in, _bytes_out;
  
public:
	
  BinaryStorageBuffer() : _buff(NULL),_len(0),_bytes_in(0),_bytes_out(0) {}
  ~BinaryStorageBuffer() {if(_buff != NULL) delete[]_buff;}
  
  char*        getBuffer      ( )  const  {return _buff;}
  unsigned int getByteLength  ( )  const  {return _bytes_in;} 
  unsigned int getBytesOut    ( )  const  {return _bytes_out;} 
  void         clear          ( )         
  {
    if(_buff != NULL) delete [] _buff;
    _buff = NULL;
    _len = _bytes_in = _bytes_out = 0;
  }
  // ----------------------------------------------------------------------------------------
  // set_buff
  // ----------------------------------------------------------------------------------------
  inline void set_buff()
  {
#ifdef __DEBUG__
	std::cout<<"BinaryStorageBuffer::set_buff";
#endif
	
	if(_buff != NULL) delete [] _buff;
	
	_buff = new char[MAX_BUFF];
	
	if(_buff == NULL) fatal("BinaryStorageBuffer::set_buff::memory exhausted !!\n");
	
	_len = MAX_BUFF;
	_bytes_in = 0;
	memset(_buff,'\0',MAX_BUFF);
	
#ifdef _DEBUG_
    std::cout<<"[ok]"<<std::endl;
#endif
  }
  // ----------------------------------------------------------------------------------------
  // set_buff
  // ----------------------------------------------------------------------------------------
  inline void set_buff(void* zone, unsigned int length)
  {
    if(_buff != NULL) delete [] _buff;
    
    _buff = new char [length];
    
    memcpy(_buff, zone, length);
    
    _bytes_in = _len = length ;
    
    _bytes_out = 0;
  }
  // ----------------------------------------------------------------------------------------
  // extend_buff
  // ----------------------------------------------------------------------------------------
  inline void extend_buff()
  {
#ifdef _DEBUG_
	std::cout<<"BinaryStorageBuffer::extend_buff"<<std::flush;
#endif
	
	char *old_buff, *new_buff;
	
	old_buff = _buff;
	
	new_buff = new char[_len + MAX_BUFF];
	
	if(new_buff == NULL) fatal("BinaryStorageBuffer::extend_buff::memory exhausted !!\n");
	
	memcpy(new_buff,_buff,_len);
	
	_buff = new_buff;
	
	_len += MAX_BUFF;
	
	delete [] old_buff;
	
#ifdef _DEBUG_
	std::cout<<"["<<_len/1024<<"kB]"<<std::endl;
#endif
  }
  // ----------------------------------------------------------------------------------------
  // store
  // ----------------------------------------------------------------------------------------
  inline void store (void* stream, unsigned int nb_bytes)
  {
    
	while( !((_bytes_in + nb_bytes) < _len) ) extend_buff();
	
	if(nb_bytes == 1) {
	  _buff[_bytes_in] = *(char*)stream;
	} else {
	  char *tab = (char*)stream;
	  for(unsigned int i = 0; i < nb_bytes; i++)
		_buff[_bytes_in + i] = tab[i];
	}
	_bytes_in += nb_bytes;
  }
  // ----------------------------------------------------------------------------------------
  // read
  // ----------------------------------------------------------------------------------------
  inline void read (void *out, unsigned int nb_bytes)
  {
    
	if(((_bytes_out + nb_bytes) > _bytes_in) ) 
	  fatal("BinaryStorageBuffer::read::attempt to read beyond buffer length (asked %i bytes)\n",nb_bytes);
	else {
	  char *tab = (char*)out;
	  if(nb_bytes == 1) {
        (*tab) = _buff[_bytes_out];
	  } else {
        for(unsigned int i = 0; i < nb_bytes; ++i)
          tab[i] = _buff[_bytes_out + i];
	  }
	}
	_bytes_out += nb_bytes;
  }
  
  void BSBread(void *out, unsigned int nb_bytes)
  { read(out,nb_bytes); }
};

#endif
