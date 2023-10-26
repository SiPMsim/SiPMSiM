//
// Created by mdupont on 12/03/18.
//

#ifndef NUMPY_FROM_CPP_TREEFILE_HH
#define NUMPY_FROM_CPP_TREEFILE_HH

#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

class Data
{
public:
  Data(const void * pointer_to_data,
       const std::string name,
       const std::type_index type_index

  ) : m_pointer_to_data(pointer_to_data),
      m_name(name),
      m_type_index(type_index)
  {

  }

  const std::string &name() const
  {
    return m_name;
  }

public:
  const void *m_pointer_to_data;
//  const size_t m_size_of_data;
  const std::string m_name;
  const std::type_index m_type_index;
};



class TreeFile
{
public:
  TreeFile();
public:
  virtual void register_variable(std::string name, const void *p, std::type_index t_index) = 0;
  virtual void register_variable(const std::string name, const std::string *p, size_t nb_char) = 0;
  virtual void register_variable(const std::string name, const char *p, size_t nb_char) = 0;

  virtual void open(const char* s) = 0;
  virtual void close() = 0;
  virtual void fill() = 0;
  virtual void writeHeader() = 0;

private:
  template<typename T>
  void add_size()
  {
    m_tmapOfSize[typeid(T)] = sizeof(T);
  }

protected:
  std::unordered_map<std::type_index, std::size_t> m_tmapOfSize;
};


#endif //NUMPY_FROM_CPP_TREEFILE_HH
