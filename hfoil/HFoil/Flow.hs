{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Flow( solveFlow
                 ) where

import Data.Packed.ST(newMatrix, writeMatrix, freezeMatrix)
import Control.Monad(forM_)
import Control.Monad.ST(runST)
import Numeric.LinearAlgebra hiding(i)
import Foreign.Storable

import HFoil.Foil

solveFlow :: (Num (Vector a), RealFloat a, Field a) => Foil a -> a -> (Vector a, Vector a)
solveFlow panels alpha = ((mapVector (\q -> cos(q - alpha)) angles) + (mV <> qg), qg)
  where
    angles = pAngles panels
    (mA, mV) = getA panels
    b = getB panels alpha
    qg = flatten $ linearSolve mA b

getB :: (Floating a, Storable a) => Foil a -> a -> Matrix a
getB panels alpha = asColumn $ fromList $ (map (\i -> sin $ (angles @> i) - alpha) [0..n-1]) ++
                    [-cos((angles @> 0) - alpha) - cos((angles @> (n-1)) - alpha)]
  where
    angles = pAngles panels
    n = (dim (fst (pNodes panels)))-1    

lnrr :: (Floating a, Storable a) =>
        Foil a -> (Int, Int) -> a
lnrr panels (i,j)
  | i == j = 0
  | otherwise = log(r1/r0)
  where
    (xms, yms) = pMidpoints panels
    (xns, yns) = pNodes panels
    r1 = distance (xms @> i, yms @> i) ( xns @> (j+1), yns @> (j+1))
    r0 = distance (xms @> i, yms @> i) ( xns @>  j   , yns @>  j   )
    distance (x1,y1) (x0,y0) = sqrt $ dx*dx + dy*dy
      where
        dx = x1 - x0
        dy = y1 - y0

beta :: (RealFloat a, Storable a) =>
        Foil a -> (Int, Int) -> a
beta panels (i,j)
  | i == j    = pi
  | otherwise = atan2 (dyjp*dxj - dxjp*dyj)
                      (dxjp*dxj + dyjp*dyj)
  where
    dyj  = (yms @> i) - (yns @>  j   )
    dxj  = (xms @> i) - (xns @>  j   )
    dyjp = (yms @> i) - (yns @> (j+1))
    dxjp = (xms @> i) - (xns @> (j+1))
    (xms, yms) = pMidpoints panels
    (xns, yns) = pNodes panels

getA :: (Num (Vector a), RealFloat a, Container Vector a) =>
        Foil a -> (Matrix a, Matrix a)
getA panels = ( scale (1/(2*pi)) $ fromBlocks [ [ mS, asColumn $ fromList (map (sum . toList) (toRows mV))]
                                              , [ mK, fromLists [[mk]]]]
              , scale (1/(2*pi)) $ fromBlocks [ [-mV, asColumn $ fromList (map (sum . toList) (toRows mS))]]
              )
  where
    n = (dim (fst (pNodes panels)))-1
--    xns = fromList $ fst (unzip panels)
--    yns = fromList $ snd (unzip panels)
--    (xms,yms) = (\(x,y) -> (fromList x, fromList y)) $ midpoints panels
    
    -- all the indices ready to be mapped over
    ijs :: [[(Int,Int)]]
    ijs = map (\i -> map (\j -> (i,j)) [0..n-1]) [0..n-1]

    (lnrrs, betas) = runST $ do
      lnrrsST <- newMatrix 0 n n
      betasST <- newMatrix 0 n n
      
      forM_ [0..n-1] $ \i -> forM_ [0..n-1] $ \j -> do
        _ <- writeMatrix lnrrsST i j (lnrr panels (i,j))
        writeMatrix betasST i j (beta panels (i,j))
      l <- freezeMatrix lnrrsST
      b <- freezeMatrix betasST
      return (l,b)
      
    angles = pAngles panels
    mS = fromLists $ map (map panelVortex) ijs
      where
        panelVortex ij@(i,j) = sin( qi - qj )*(lnrrs @@> ij) + cos( qi - qj )*(betas @@> ij)
          where
            qi = angles @> i
            qj = angles @> j

    mV = fromLists $ map (map panelSource) ijs
      where
        panelSource ij@(i,j) = cos( qi - qj )*(lnrrs @@> ij) - sin( qi - qj )*(betas @@> ij)
          where
            qi = angles @> i
            qj = angles @> j

    mK = asRow $ fromList $ map sum $ map (map kuttaVortex) $ map (\j -> map (\k -> (k,j)) [0,n-1]) [0..n-1]
      where
        kuttaVortex kj@(k,j) = sin( qk - qj )*(betas @@> kj) - cos( qk - qj )*(lnrrs @@> kj)
          where
            qk = angles @> k
            qj = angles @> j

    mk = sum $ concatMap (map kuttaVortex) $ map (\j -> map (\k -> (k,j)) [0,n-1]) [0..n-1]
      where
        kuttaVortex kj@(k,j) = sin( qk - qj )*(lnrrs @@> kj) + cos( qk - qj )*(betas @@> kj)
          where
            qk = angles @> k
            qj = angles @> j
