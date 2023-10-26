//
// Created by mdupont on 12/03/18.
//

#include "TreeFile.hh"

#include <iostream>

TreeFile::TreeFile()
{


  this->add_size<double>();
  this->add_size<float>();
//  this->add_size<int>();
  this->add_size<uint8_t>();
  this->add_size<uint16_t>();
  this->add_size<uint32_t>();
  this->add_size<uint64_t>();

  this->add_size<int8_t>();
  this->add_size<int16_t>();
  this->add_size<int32_t>();
  this->add_size<int64_t>();

  this->add_size<bool>();
  this->add_size<char>();


}
