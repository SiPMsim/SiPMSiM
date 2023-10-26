
#include "GateXMLDocument.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define LIBXML_SAX1_ENABLED


//!@name Constructors and destructors
//@{

//! Default constructor.
/**
  * Opens the xml file whose name is given by filename.
  *
  * Use the Ok() method to check if the opening went ok.
  * */
GateXMLDocument::GateXMLDocument(const G4String& filename) :
  m_ok(false), m_reset(true)

//
// SJ COMMENTS## : read the file by using a messenger mechanism
//
{
  m_doc = xmlParseFile(filename.c_str());

  if (m_doc)
  {
    m_cur = xmlDocGetRootElement(m_doc);
    if (m_cur==0) xmlFreeDoc(m_doc);
    else m_ok = true;
  }
  else
  {
    std::cout << "I/O warning: Discard the previous warning if your simulation doesn't transport OPTICAL PHOTONS. \n";
    std::cout << "Otherwise, please copy the "<< filename << " file from the gate-source directory in the directory where you run your main macro.\n";
  }
}


//! Destructor.
GateXMLDocument::~GateXMLDocument()
{
  if (m_doc) xmlFreeDoc(m_doc);
}
//@}

//! Check if file is opened correctly and is ready for read operations
/**
  * Use this check after calling the constructor to check if there were no errors.
  * */
bool GateXMLDocument::Ok() const
{ return m_ok;}

//!@name Getting values
//@{

//! Gets the name of the current node.
/**
  * For example, returns "foo" for the node <foo prop="bar">flierp</foo>
  * */
G4String GateXMLDocument::GetName()
{
  return (char *)m_cur->name;
}

//! Gets the value of the property.
/**
  * For example, returns "bar" for the node <foo prop="bar">flierp</foo> when called as GetProperty("prop").
  * Returns an empty string when the property does not exist.
  * */
G4String GateXMLDocument::GetProperty(const G4String& property)
{
  xmlChar*    key = xmlGetProp(m_cur,(const xmlChar*)property.c_str());
  G4String str = key ? (char *)key : "";
#ifdef WIN32
  xmlFreeFunc((void*)key);
#else
  xmlFree(key);
#endif
  return str;
}

//! Checks if the current node has a certain property.
/**
  * For example, returns true for the node <foo prop="bar">flierp</foo> when called as HasProperty("prop").
  * */
G4bool GateXMLDocument::HasProperty(const G4String& property)
{
  xmlChar* key = xmlGetProp(m_cur,(const xmlChar*)property.c_str());
  G4bool     prop = key ? true : false;
#ifdef WIN32
  xmlFreeFunc((void*)key);
#else
  xmlFree(key);
#endif
  return prop;
}

//! Gets the content of the node.
/**
  * For example, returns "flierp" for the node <foo prop="bar">flierp</foo>.
  * */
G4String GateXMLDocument::GetValue()
{
  xmlChar*    key = xmlNodeListGetString(m_doc, m_cur->xmlChildrenNode, 1);
  G4String str = key ? (char *)key : "";
#ifdef WIN32
  xmlFreeFunc((void*)key);
#else
  xmlFree(key);
#endif
  return str;
}
//@}

//!@name Navigation routines
//@{
//! Returns to the root node.
void GateXMLDocument::Reset()
{
  m_cur = xmlDocGetRootElement(m_doc);
  m_reset = true;
}

//! Goes to the first daughter of the current node.
/**
  * Returns false when the node does not contain a daughter.
  * */
G4bool GateXMLDocument::Enter()
{
  if (m_cur->xmlChildrenNode!=0)
  {
    m_cur   = m_cur->xmlChildrenNode;
    m_reset = true;
    return true;
  }
  else return false;
}

//! Goes to the mother of the current node.
void GateXMLDocument::Leave()
{
  m_cur   = m_cur->parent;
  m_reset = false;
}

//! Goes to the next node.
/**
  * Returns false when there is no more node, true otherwise. The method can therefore be used
  * in a while loop:
  * \code
  * while (xmldoc.Next())
  * {
  *   // do something with the node
  * }
  * \endcode
*/
G4bool GateXMLDocument::Next()
{
  if (m_cur->next!=0)
  {
    m_cur   = m_cur->next;
    m_reset = false;
    return true;
  }
  else
  {
    m_reset = false;
    return false;
  }
}

//! Goes to the previous node.
/**
  * Returns false when there is no more node, true otherwise. The method can therefore be used
  * in a while loop:
  * \code
  * while (xmldoc.Previous())
  * {
  *   // do something with the node
  * }
  * \endcode
  * */
G4bool GateXMLDocument::Previous()
{
  if (m_cur->prev!=0)
  {
    m_cur = m_cur->prev;
    return true;
  }
  else
  {
    m_reset = true;
    return false;
  }
}

//! Goes to the first node.
void GateXMLDocument::First()
{
  while (Previous()) ;
  m_reset = true;
}
//@}

//!@name Finding nodes
//@{

//! Finds the next node the name given by 'tag'.
/**
  * Returns true when found false otherwise.
  * Find only looks at the current depth.
  * */
G4bool GateXMLDocument::Find(const G4String& tag)
{
  if (!m_reset)
  { if (!Next()) return false;}

  do
  {
    if (GetName()==tag)
    {
      m_reset = false;
      return true;
    }
  }
  while (Next());

  return false;
}

//! Finds the next node the name given by 'tag' and the property name equal to 'value'.
/**
  * Returns true when found false otherwise.
  * Find only looks at the current depth.
  * */
G4bool GateXMLDocument::Find(const G4String& tag, const G4String& name)
{
  return Find(tag,"name",name);
}

//! Finds the next node the name given by 'tag' and the property 'property' equal to 'value'.
/**
  * Returns true when found false otherwise.
  * Find only looks at the current depth.
  * */
G4bool GateXMLDocument::Find(const G4String& tag, const G4String& property, const G4String& value)
{
  if (!m_reset)
  { if (!Next()) return false;}

  do
  {
    if (GetName()==tag)
    {
      if (GetProperty(property)==value)
      {
    m_reset = false;
    return true;
      }
    }
  }
  while (Next());

  return false;
}
//@}

//! Gets the current position in the document.
/**
  * Can be used in combination with SetState() to return to the previous position in the document.
  * */
GateXMLDocumentState GateXMLDocument::GetState()
{
  GateXMLDocumentState state;
  state.cur   = m_cur;
  state.reset = m_reset;
  return state;
}

//! Sets the position in the document back to the one given by state.
void GateXMLDocument::SetState(GateXMLDocumentState state)
{
  m_cur   = state.cur;
  m_reset = state.reset;
}


