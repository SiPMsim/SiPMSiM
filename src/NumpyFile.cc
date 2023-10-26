//
// Created by mdupont on 30/08/17.
//

#include "NumpyFile.hh"
#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cstring>

#include <stdexcept>

using namespace std;


int is_big_endian()
{
  union {
    long int l;
    char c[sizeof (long int)];
  } u;

  u.l = 1;

  if (u.c[sizeof(long int)-1] == 1)
  {
    return 1;
  }
  else
    return 0;
}

string get_endianness_carater()
{
  if (is_big_endian())
    return ">";
  return "<";
}


NumpyFile::NumpyFile() : m_nb_elements(0)
{
  //cout << "NumpyFile::NumpyFile()\n";
  //cout << "NumpyFile::NumpyFile() is_big_endian = " << is_big_endian() << "\n";

  m_tmapOfDefinition[typeid(double)] = get_endianness_carater() + "f8";
  m_tmapOfDefinition[typeid(float)] = get_endianness_carater() + "f4";

  m_tmapOfDefinition[typeid(uint8_t)] = get_endianness_carater() + "u1";
  m_tmapOfDefinition[typeid(uint16_t)] = get_endianness_carater() + "u2";
  m_tmapOfDefinition[typeid(uint32_t)] = get_endianness_carater() + "u4";
  m_tmapOfDefinition[typeid(uint64_t)] = get_endianness_carater() + "u8";

  m_tmapOfDefinition[typeid(int8_t)] = get_endianness_carater() + "i1";
  m_tmapOfDefinition[typeid(int16_t)] = get_endianness_carater() + "i2";
  m_tmapOfDefinition[typeid(int32_t)] = get_endianness_carater() + "i4";
  m_tmapOfDefinition[typeid(int32_t)] = get_endianness_carater() + "i4";
  m_tmapOfDefinition[typeid(int64_t)] = get_endianness_carater() + "i8";
  m_tmapOfDefinition[typeid(bool)] = get_endianness_carater() + "?";
  m_tmapOfDefinition[typeid(char)] = get_endianness_carater() + "b";

//  m_file = ofstream();

//  m_file.open("/tmp/z.npy",std::ofstream::binary);


}


void NumpyFile::writeHeader()
{
  //cout << "NumpyFile::writeHeader()\n";

  m_file << magic_prefix.c_str();


  uint32_t  magic_len = magic_prefix.size() + 2;

  //cout << "magic_len=" <<  magic_len <<'\n';


//  string dico = "{'descr': '<f8', 'fortran_order': False, 'shape': (1,), }";


  std::stringstream ss_dico_before_shape, ss_dico_after_shape;
  string dico_before_shape, dico_after_shape;

  ss_dico_before_shape << "{'descr': [";
//  for (auto&& d : m_vector_names) // access by const reference
//  {
//    ss_dico_before_shape << d;
////    if 'd != m_vector_names.end();
////    ss_dico << ",";
//  }

  for ( auto it = m_vector_of_pointer_to_data.begin(); it != m_vector_of_pointer_to_data.end(); ++it )
  {
    if(it != m_vector_of_pointer_to_data.begin())
      ss_dico_before_shape << ", ";
    ss_dico_before_shape << (*it).m_numpy_description;
  }

  ss_dico_before_shape << "], 'fortran_order': False, 'shape': (";

  string shape = "######";

  stringstream ss_shape;
  ss_shape << std::setw(20) << std::setfill(' ') << 0;
  shape = ss_shape.str();

  ss_dico_after_shape << ",),}";


  dico_before_shape = ss_dico_before_shape.str();
  dico_after_shape = ss_dico_after_shape.str();


  //cout << "dico_before_shape = " << dico_before_shape << "size = " << dico_before_shape.size() << '\n';
  //cout << "ss_dico_after_shape = " << dico_after_shape << "size = " << dico_after_shape.size() << '\n';


  uint32_t  current_header_len = magic_len + 2  + dico_before_shape.size() + shape.size() + dico_after_shape.size() + 1;  // 1 for the newline



  uint32_t  topad = 16 - (current_header_len % 16);

  //cout << "current_header_len = " << current_header_len << " topad = " << topad;

  for(int i = 0; i < topad; ++i)
    dico_after_shape += ' ';

  dico_after_shape += '\n';

  uint16_t hlen = dico_before_shape.size() + shape.size() + dico_after_shape.size();

  //cout << "current_header_len = " << current_header_len << " topad = " << topad << "hlen = " << hlen << "\n";

  uint8_t major = 1;
  uint8_t minor = 0;

  m_file.write((char*)&major, sizeof(major));
  m_file.write((char*)&minor, sizeof(minor));
  m_file.write((char*)&hlen, sizeof(hlen));


  m_file << dico_before_shape.c_str();
  m_position_before_shape = m_file.tellp();
//  m_file << shape.c_str();

  m_file.write(shape.c_str(), shape.size());

  m_position_after_shape = m_file.tellp();
  m_file << dico_after_shape.c_str();

//  cout << "position_before_shape = " << m_position_before_shape << " position_after_shape =  " << m_position_after_shape;

  //cout << endl;




}

ostream &operator<<(ostream &os, const NumpyFile &header)
{


  return os;
}


//void NumpyFile::register_variable(std::string name, void *p, size_t size, string numpy_format)

void NumpyFile::register_variable(std::string name, const void *p, const size_t size, std::string numpy_format, std::type_index t_index)
{
  //cout << "NumpyFile::register_variable() " << " name = " << name << " type = " <<  t_index.name() << " size = " << size << " numpy_format = " << numpy_format << "\n";

  NumpyData d(p, size, name, numpy_format, t_index);

  for (auto&& d : m_vector_of_pointer_to_data) // access by const reference
  {
      if(d.name() == name)
      {
        string s("Error: Key '");
        s+= name;
        s += "' already used !";
        cerr << s << endl;
        throw new runtime_error( s );
      }
  }

  m_vector_of_pointer_to_data.push_back(d);
}

void NumpyFile::register_variable(const std::string name, const char *p, size_t nb_char)
{
  stringstream ss;
  ss << "|S" << nb_char;
  register_variable(name, p, sizeof(char)*nb_char, ss.str(), typeid(char*));
  auto&& data = m_vector_of_pointer_to_data.back();
  data.m_nb_characters = nb_char;
}


void NumpyFile::register_variable(const std::string name, const std::string *p, size_t nb_char)
{
  stringstream ss;
  ss << "|S" << nb_char;
  register_variable(name, p, sizeof(char)*nb_char, ss.str(), typeid(string));
  auto&& data = m_vector_of_pointer_to_data.back();
  data.m_nb_characters = nb_char;
}

//void NumpyFile::register_variable(const std::string name, const char array[])
//{
//    cout << "NumpyFile::register_variable(const std::string name, const char *p)";
//    cout << "name = " << name << "\n";
//    std::cout << "Length of array = " << (sizeof(char)/sizeof(*array)) << std::endl;
//
//}


void NumpyFile::fill()
{
//  cout << "NumpyFile::fill()";
  for (auto&& d : m_vector_of_pointer_to_data) // access by const reference
  {
    if(d.m_nb_characters == 0)
      m_file.write((const char*)d.m_pointer_to_data, d.m_size_of_data);
    else
    {
      if(d.m_type_index == typeid(char*))
      {
        char *p_data = (char*)d.m_pointer_to_data;
        int current_nb_characters = strlen(p_data);
//        cout << " p_data = " << p_data;
//        cout << " current_nb_characters = " << current_nb_characters;
        for(size_t i = current_nb_characters; i < d.m_nb_characters; ++i)
          p_data[i] = '\0';



//        cout << " p_data = " << p_data;

        m_file.write(p_data, d.m_size_of_data);
      } else if (d.m_type_index == typeid(string))
      {
        const string *p_s = (const string*) d.m_pointer_to_data;
        if( p_s->size() > d.m_nb_characters)
        {
          string m;
          m += "NumpyFile::fill ";
          m += "lenght(" + *p_s + ") = (" + std::to_string(p_s->size()) +   ") > " + std::to_string(d.m_nb_characters);
          throw std::length_error(m);
        }

        string s(*p_s); // copy of the data, :-(
        s.resize(d.m_nb_characters, '\0');


        m_file.write(s.c_str(), d.m_size_of_data);
      }
    }
  }

  m_nb_elements++;
//  cout << "\n";
}

void NumpyFile::close()
{


//  cout << "current position = " << m_file.tellp() << "\n";

  m_file.seekp(m_position_before_shape);
  //cout << "current position = " << m_file.tellp() << "\n";
  stringstream ss_shape;
  ss_shape << std::setw(20) << std::setfill(' ') << m_nb_elements;
  string shape = ss_shape.str();
  m_file.write(shape.c_str(), shape.size());
//  cout << "current position = " << m_file.tellp() << "\n";
//
//  cout << "close.." << endl;


  m_file.close();
}

void NumpyFile::open(const char *s)
{
  m_file.open(s, std::ofstream::binary);

}

void NumpyFile::register_variable(std::string name, const void *p, std::type_index t_index)
{
  //std::cout << "NumpyFile::register_variable name = " << name << " t_index = " << t_index.name() << "\n";
  this->register_variable(name, p, m_tmapOfSize.at(t_index), m_tmapOfDefinition.at(t_index), t_index);


}
