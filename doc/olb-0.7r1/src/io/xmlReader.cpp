/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010 Jonas Latt, Jonas Fietz, Mathias Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * Input/Output in XML format -- non-generic code.
 */

#include "io/xmlReader.h"
#include <algorithm>
#include <cctype>
#include <iostream>

namespace olb {

XMLreader XMLreader::notFound;

XMLreader::XMLreader()
  : clout(std::cout,"XMLreader") {
  name = "XML node not found";
}

XMLreader::XMLreader( TiXmlNode* pParent )
  : clout(std::cout,"XMLreader") {

  if (singleton::mpi().isMainProcessor()) {
    mainProcessorIni(pParent);
  }
  else {
    slaveProcessorIni();
  }
}

XMLreader::XMLreader(const std::string& fName )
  : clout(std::cout,"XMLreader") {
  TiXmlDocument* doc = 0;
  int loadOK = false;
  if (singleton::mpi().isMainProcessor()) {
    doc = new TiXmlDocument(fName.c_str());
    loadOK = doc->LoadFile();
    if(!loadOK) {
      clout << std::string("Problem processing input XML file ") << fName << std::endl;
    }

  }


  if (singleton::mpi().isMainProcessor()) {
    mainProcessorIni(doc);
    delete doc;
  }
  else {
    slaveProcessorIni();
  }
}

XMLreader::~XMLreader() {
  for (unsigned int iNode=0; iNode<children.size(); ++iNode) {
    delete children[iNode];
  }
}

void XMLreader::mainProcessorIni( TiXmlNode* pParent ) {
  assert (pParent->Type()==TiXmlNode::TINYXML_DOCUMENT || pParent->Type()==TiXmlNode::TINYXML_ELEMENT );

  if (pParent->Type() == TiXmlNode::TINYXML_DOCUMENT) {
    // ignore the surrounding PARAM-block
    pParent = pParent->FirstChild();
  }

  name = pParent->ValueStr();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  singleton::mpi().bCast(&name,1);
#endif

  TiXmlNode * pChild;

  int type = 0;
  for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling())
  {
    type = pChild->Type();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
    singleton::mpi().bCast(&type, 1);
#endif
    if ( type==TiXmlNode::TINYXML_ELEMENT ) {
      children.push_back( new XMLreader( pChild ) );
    }
    else if ( type==TiXmlNode::TINYXML_TEXT ) {
      text = pChild->ToText()->ValueStr();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
      singleton::mpi().bCast(&text,1);
#endif
    }
  }
  type = TiXmlNode::TINYXML_UNKNOWN;
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  singleton::mpi().bCast(&type, 1);
#endif
}

void XMLreader::slaveProcessorIni()
{
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  singleton::mpi().bCast(&name,1);
#endif

  int type=0;
  do {
#ifdef PARALLEL_MODE_MPI  // parallel program execution
    singleton::mpi().bCast(&type, 1);
#endif
    if ( type==TiXmlNode::TINYXML_ELEMENT ) {
      children.push_back( new XMLreader( 0 ) );
    }
    else if ( type==TiXmlNode::TINYXML_TEXT ) {
#ifdef PARALLEL_MODE_MPI  // parallel program execution
      singleton::mpi().bCast(&text,1);
#endif
    }
  }
  while (type != TiXmlNode::TINYXML_UNKNOWN);
}

void XMLreader::print(int indent) const {
  std::string indentStr(indent, ' ');
  clout << indentStr << "[" << name << "]" << std::endl;
  if (!text.empty()) {
    clout << indentStr << "  " << text << std::endl;
  }
  for (unsigned int iNode=0; iNode<children.size(); ++iNode) {
    children[iNode]->print(indent+2);
  }
}

XMLreader const& XMLreader::operator[] (std::string name) const
{
  for (unsigned int iNode=0; iNode<children.size(); ++iNode) {
    if (children[iNode]->name == name) {
      return *children[iNode];
    }
  }
  clout << "Warning: cannot read value from node \"" << name << "\"" << std::endl;
  return notFound;
}

std::vector<XMLreader*>::const_iterator XMLreader::begin() const {
  return children.begin();
}

std::vector<XMLreader*>::const_iterator XMLreader::end() const {
  return children.end();
}

std::string XMLreader::getName() const {
  return name;
}




template <>
bool XMLreader::read(bool& value) const {
  std::stringstream valueStr(text);
  std::string word;
  valueStr >> word;
  // Transform to lower-case, so that "true" and "false" are case-insensitive.
  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  if ((word=="true") || (word=="1")) {
    value = true;
    return true;
  }
  else if ((word=="false") || (word=="0")) {
    value=false;
    return true;
  }
  else {
    clout << std::string("Cannot read boolean value from XML element ") << name << std::endl;
  }
  return false;
}

template <>
bool XMLreader::read(std::string& entry) const {
  if(name == "XML node not found") {
    return false;
  }

  entry = text;
  return true;
}

}  // namespace olb
