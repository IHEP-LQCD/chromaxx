// -*- C++ -*-
/*! \file
 *  \brief Write object function std::map
 */

#ifndef __uky_write_obj_funcmap_h__
#define __uky_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function std::map
  /*! \ingroup inlineio */
  namespace UKYWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function std::map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_5(const std::string&,
			     const std::string&, 
			     int, int, int),
		  void (*)(const std::string& buffer_id,
			   const std::string& filename),
		  StringFunctionMapError> >
    TheUKYWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
