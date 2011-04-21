// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
#ifndef __DISCRETE_SYSTEM_MODEL__
#define __DISCRETE_SYSTEM_MODEL__

#include "systemmodel.h"
#include "../pdf/discreteconditionalpdf.h"

namespace BFL
{
  /// Class for discrete System Models
  /** Class representing discrete System Models, ie.  System Models for
      which _BOTH_ states and inputs are discrete variables!
  */
  class DiscreteSystemModel: public SystemModel<int>
    {
    public:
      /// Constructor
      /** @param systempdf ConditionalPdf<int> representing
	  P(X_k | X_{k-1}, U_{k})
	  @see SystemPdf
      */
      DiscreteSystemModel(DiscreteConditionalPdf * systempdf = NULL);
      /// Destructor
      virtual ~DiscreteSystemModel();
      /// Copy constructor
      DiscreteSystemModel(const DiscreteSystemModel &);
      /// Get the number of discrete states
      unsigned int NumStatesGet()const;
    };

} // End namespace BFL
#endif // __DISCRETE_SYSTEM_MODEL__
