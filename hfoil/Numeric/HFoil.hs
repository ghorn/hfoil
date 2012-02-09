{- |
   Module      : Numeric.Hfoil
   Description : Top level module

   This is the top level module which exports the API
 -}

{-# OPTIONS_GHC -Wall #-}

module Numeric.HFoil( -- * Airfoils
                      module Numeric.HFoil.Foil
                      -- * Flow solution
                    , module Numeric.HFoil.Flow
                    -- * Naca 4 utilities
                    , module Numeric.HFoil.Naca4
                    -- * Drawing utilities (gloss backend)
                    , module Numeric.HFoil.Drawing
                    -- * Interactive command line application
                    , module Numeric.HFoil.Repl
                    ) where

import Numeric.HFoil.Naca4
import Numeric.HFoil.Foil
import Numeric.HFoil.Drawing
import Numeric.HFoil.Flow
import Numeric.HFoil.Repl
