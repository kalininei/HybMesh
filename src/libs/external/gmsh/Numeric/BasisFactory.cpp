// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-B-> Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#include "BasisFactory.h"
#include "GmshDefines.h"
#include "polynomialBasis.h"
#include "pyramidalBasis.h"
#include "miniBasis.h"
#include "MetricBasis.h"
#include "CondNumBasis.h"
#include <map>
#include <cstddef>

std::map<int, nodalBasis*> BasisFactory::fs;
std::map<int, MetricBasis*> BasisFactory::ms;
std::map<int, CondNumBasis*> BasisFactory::cs;
std::map<FuncSpaceData, JacobianBasis*> BasisFactory::js;
std::map<FuncSpaceData, bezierBasis*> BasisFactory::bs;
std::map<FuncSpaceData, GradientBasis*> BasisFactory::gs;

const nodalBasis* BasisFactory::getNodalBasis(int tag)
{
  // If the Basis has already been built, return it.
  std::map<int, nodalBasis*>::const_iterator it = fs.find(tag);
  if (it != fs.end()) {
    return it->second;
  }
  // Get the parent type to see which kind of basis
  // we want to create
  nodalBasis* F = NULL;
  if (tag == MSH_TRI_MINI)
    F = new miniBasisTri();
  else if (tag == MSH_TET_MINI)
    F = new miniBasisTet();
  else {
    int parentType = ElementType::ParentTypeFromTag(tag);
    switch(parentType) {
      case(TYPE_PNT):
      case(TYPE_LIN):
      case(TYPE_TRI):
      case(TYPE_QUA):
      case(TYPE_PRI):
      case(TYPE_TET):
      case(TYPE_HEX):
        F = new polynomialBasis(tag);
        break;
      case(TYPE_PYR):
        F = new pyramidalBasis(tag);
        break;
      default:
        Msg::Error("Unknown type of element %d (in BasisFactory)", tag);
        return NULL;
    }
  }

  std::pair<std::map<int, nodalBasis*>::const_iterator, bool> inserted;

#if defined(_OPENMP)
  #pragma omp critical
#endif
    {
      inserted = fs.insert(std::make_pair(tag, F));

      if (!inserted.second)
        delete F;
    }

  return inserted.first->second;
}

const JacobianBasis* BasisFactory::getJacobianBasis(FuncSpaceData fsd)
{
  FuncSpaceData data = fsd.getForNonSerendipitySpace();

  std::map<FuncSpaceData, JacobianBasis*>::const_iterator it = js.find(data);
  if (it != js.end()) return it->second;

  JacobianBasis* J = new JacobianBasis(data);
  js.insert(std::make_pair(data, J));
  return J;
}

const MetricBasis* BasisFactory::getMetricBasis(int tag)
{
  std::map<int, MetricBasis*>::const_iterator it = ms.find(tag);
  if (it != ms.end()) return it->second;

  MetricBasis* M = new MetricBasis(tag);
  ms.insert(std::make_pair(tag, M));
  return M;
}

const CondNumBasis* BasisFactory::getCondNumBasis(int tag, int cnOrder)
{
  std::map<int, CondNumBasis*>::const_iterator it = cs.find(tag);
  if (it != cs.end()) return it->second;

  CondNumBasis* M = new CondNumBasis(tag, cnOrder);
  cs.insert(std::make_pair(tag, M));
  return M;
}

const GradientBasis* BasisFactory::getGradientBasis(FuncSpaceData data)
{
  std::map<FuncSpaceData, GradientBasis*>::const_iterator it = gs.find(data);
  if (it != gs.end()) return it->second;

  GradientBasis* G = new GradientBasis(data);
  gs.insert(std::make_pair(data, G));
  return G;
}

const bezierBasis* BasisFactory::getBezierBasis(FuncSpaceData fsd)
{
  FuncSpaceData data = fsd.getForPrimaryElement();

  std::map<FuncSpaceData, bezierBasis*>::const_iterator it = bs.find(data);
  if (it != bs.end()) return it->second;

  bezierBasis* B = new bezierBasis(data);
  bs.insert(std::make_pair(data, B));
  return B;
}

void BasisFactory::clearAll()
{
  std::map<int, nodalBasis*>::iterator itF = fs.begin();
  while (itF != fs.end()) {
    delete itF->second;
    itF++;
  }
  fs.clear();

  std::map<int, MetricBasis*>::iterator itM = ms.begin();
  while (itM != ms.end()) {
    delete itM->second;
    itM++;
  }
  ms.clear();

  std::map<FuncSpaceData, JacobianBasis*>::iterator itJ = js.begin();
  while (itJ != js.end()) {
    delete itJ->second;
    itJ++;
  }
  js.clear();

  std::map<FuncSpaceData, GradientBasis*>::iterator itG = gs.begin();
  while (itG != gs.end()) {
    delete itG->second;
    itG++;
  }
  gs.clear();

  std::map<FuncSpaceData, bezierBasis*>::iterator itB = bs.begin();
  while (itB != bs.end()) {
    delete itB->second;
    itB++;
  }
  bs.clear();
}
