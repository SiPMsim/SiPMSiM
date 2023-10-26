//
// Created by mdupont on 30/08/17.
//

#ifndef NUMPY_FROM_CPP_NUMPYHEADER_HH
#define NUMPY_FROM_CPP_NUMPYHEADER_HH

#include <string>
#include <ostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>

#include "TreeFile.hh"

class NumpyData : public Data
{
public:
  NumpyData(const void * pointer_to_data,
       const size_t size_of_data,
       const std::string name,
       const std::string numpy_format,
       const std::type_index type_index

  ) : Data(pointer_to_data, name, type_index),
      m_size_of_data(size_of_data),
      m_numpy_format(numpy_format),
      m_nb_characters(0)
  {
    std::stringstream descr;
    descr << "('" << name << "', '" << numpy_format << "')";
    m_numpy_description = descr.str();
  }

  const std::string &numpy_description() const
  {
    return m_numpy_description;
  }

public:
  const size_t m_size_of_data;
  const std::string m_numpy_format;
  size_t m_nb_characters;

  std::string m_numpy_description;
};

class NumpyFile : public TreeFile
{
public:
  NumpyFile();
  void  open(const char* s) override ;

  friend std::ostream &operator<<(std::ostream &os, const NumpyFile &header);
  void fill() override ;


  void register_variable(std::string name, const void *p, std::type_index t_index) override;
  void register_variable(const std::string name, const char *p, size_t nb_char) override;
  void register_variable(const std::string name, const std::string *p, size_t nb_char) override;
  void close() override ;
  void  writeHeader() override ;


  template<typename T>
  void register_variable(const std::string name, const T *p)
  {
    register_variable(name, p, typeid(T));
  }








private:
  void register_variable(std::string name, const void *p, const size_t size, std::string numpy_format, std::type_index t_index);
  const std::string magic_prefix = "\x93NUMPY";
  uint64_t m_nb_elements;

  std::vector<NumpyData> m_vector_of_pointer_to_data;
  uint64_t m_position_before_shape;
  uint64_t m_position_after_shape;

  std::unordered_map<std::type_index, std::string> m_tmapOfDefinition;


public:
  std::ofstream m_file;


};


#endif //NUMPY_FROM_CPP_NUMPYHEADER_HH
