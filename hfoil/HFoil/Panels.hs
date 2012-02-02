{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Panels( toNodes 
                   , midpoints
                   ) where

import qualified HFoil.Naca4 as Naca4
import Numeric.LinearAlgebra
import Data.Tuple.Utils(fst3)

toNodes :: (Enum a, Floating (Vector a), Floating a, Ord a, Field a) => Naca4.Naca4 a -> Int -> [(a,a)]
toNodes foil nPanels = [(1,0)]++reverse (zip xcs (map negate ycs))++[(0,0)]++zip xcs ycs++[(1,0)]
  where
    xcs = toList $ fst3 $ bunchPanels (Naca4.yt foil) (Naca4.dyt foil) xcs0 0 0
    ycs = map (Naca4.yt foil) xcs
    xcs0 = fromList $ init $ tail $ toList $ linspace nXcs (0,1)
    nXcs = (nPanels + (nPanels `mod` 2)) `div` 2 + 1


midpoints :: Fractional a => [(a,a)] -> ([a],[a])
midpoints panels = unzip $ zipWith (\(x1, y1) (x0, y0) -> (0.5*(x1+x0),0.5*(y1+y0))) (tail panels) (init panels)



bunchPanels :: (Enum a, Floating (Vector a), Floating a, Ord a, Field a) =>
               (a -> a) -> (a -> a) -> Vector a -> Int -> Int -> (Vector a, Int, Int)
bunchPanels yt dyt xcs nIter nBadSteps 
  | nIter     > 300  = error "panel buncher exceeded 300 iterations"
  | nBadSteps > 1000 = error "panel buncher exceeded 1000 bad steps"
  | sum (toList (abs deltaXcs)) < 1e-12 = (xcs, nIter, nBadSteps)
  | otherwise                           = bunchPanels yt dyt goodXcs (nIter+1) (nBadSteps + length badOnes)
  where
    (badOnes, goodXcs:_) = break (\xs -> all (>0) (toList xs)) nextXcs
    nextXcs = map (\alpha -> xcs + (scale alpha deltaXcs)) (map (2.0**) [0,-1..])
    deltaXcs = xcsStep yt dyt xcs

xcsStep :: (Enum a, Floating (Vector a), Field a) => (a -> a) -> (a -> a) -> Vector a -> Vector a
xcsStep yt dyt xcs = flatten $ -(linearSolveLS mat2 (asColumn rs))
  where
    n = dim xcs
    
    xs = join [fromList [0],              xcs, fromList[1]]
    ys = join [fromList [0], mapVector yt xcs, fromList[0]]
    dxs = (subVector 1 (n+1) xs) - (subVector 0 (n+1) xs)
    dys = (subVector 1 (n+1) ys) - (subVector 0 (n+1) ys)
    deltas = sqrt (dxs*dxs + dys*dys)
    dydxs = mapVector dyt xcs
    
    mat = r1 - r0
      where
        r0 = fromBlocks [[(1><n)[0,0..]],[diagRect 0 diag0 n n]]
        r1 = fromBlocks [                [diagRect 0 diag1 n n], [(1><n)[0,0..]]]
        
        (diag0, diag1) = (subVector 1 n d0, subVector 0 n d1)
          where
            d0 = (dxs + dys*dy0dxs)/deltas
            d1 = (dxs + dys*dy1dxs)/deltas
            dy0dxs = join [fromList [0], dydxs]
            dy1dxs = join [dydxs, fromList [0]]
    
    zeros = (n><1)[0,0..]
    eye = ident n
    frontBunchingParam = 2.0
    diff = (fromBlocks [[zeros, eye]]) - (fromBlocks [[eye*(1+frontBunchingParam/(fromIntegral n)), zeros]])

    mat2 = diff <> mat
    rs = diff <> deltas
