{- |
   Module      : Hfoil
   Description : Top level module

   This is the top level module which exports the API
 -}

{-# OPTIONS_GHC -Wall #-}

module HFoil
       ( -- * Airfoils
         module HFoil.Foil
         -- * Flow solution
       , module HFoil.Flow
       -- * Naca 4 utilities
       , module HFoil.Naca4
       -- * Drawing utilities (gloss backend)
       , module HFoil.Drawing
       -- * Interactive command line application
       , module HFoil.Repl
       ) where

import HFoil.Naca4
import HFoil.Foil
import HFoil.Drawing
import HFoil.Flow
import HFoil.Repl
