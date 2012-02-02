{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Flow( solveFlow
                 ) where

import Numeric.LinearAlgebra hiding(i)
import Foreign.Storable

import HFoil.Panels

solveFlow :: (Num (Vector a), RealFloat a, Field a) => [(a, a)] -> a -> (Vector a, Vector a)
solveFlow panels alpha = ((mapVector (\q -> cos(q - alpha)) angles) + (mV <> qg), qg)
  where
    (mA, mV, angles) = getA panels
    b = getB panels alpha
    qg = flatten $ linearSolve mA b

getB :: (RealFloat a, Storable a) => [(a, a)] -> a -> Matrix a
getB panels alpha = asColumn $ fromList $ (map (\i -> sin $ (angles @> i) - alpha) [0..n-1]) ++
                    [-cos((angles @> 0) - alpha) - cos((angles @> (n-1)) - alpha)]
  where
    n = (length panels)-1
    angles = fromList $ zipWith (\(x1,y1) (x0,y0) -> atan2 (y1 - y0) (x1 - x0)) (tail panels) (init panels)

getA :: (Num (Vector a), RealFloat a, Container Vector a) =>
        [(a, a)] -> (Matrix a, Matrix a, Vector a)
getA panels = ( scale (1/(2*pi)) $ fromBlocks [ [ mS, asColumn $ fromList (map (sum . toList) (toRows mV))]
                                              , [ mK, fromLists [[mk]]]]
              , scale (1/(2*pi)) $ fromBlocks [ [-mV, asColumn $ fromList (map (sum . toList) (toRows mS))]]
              , angles
              )
  where
    n = (length panels)-1
    xns = fromList $ fst (unzip panels)
    yns = fromList $ snd (unzip panels)
    (xms,yms) = (\(x,y) -> (fromList x, fromList y)) $ midpoints panels
    angles = fromList $ zipWith (\(x1,y1) (x0,y0) -> atan2 (y1 - y0) (x1 - x0)) (tail panels) (init panels)
    
    -- all the indices ready to be mapped over
    ijs :: [[(Int,Int)]]
    ijs = map (\i -> map (\j -> (i,j)) [0..n-1]) [0..n-1]

    lnrrs = fromLists $ map (map lnrr) ijs
      where
        lnrr (i,j)
          | i == j    = 0
          | otherwise = log(r1/r0)
          where
            r1 = distance (xms @> i, yms @> i) ( xns @> (j+1), yns @> (j+1))
            r0 = distance (xms @> i, yms @> i) ( xns @>  j   , yns @>  j   )

            distance (x1,y1) (x0,y0) = sqrt $ dx*dx + dy*dy
              where
                dx = x1 - x0
                dy = y1 - y0

    betas = fromLists $ map (map beta) ijs
      where
        beta (i,j)
          | i == j    = pi
          | otherwise = v1 - v0
          where
            v1 = atan2 ((yms @> i) - (yns @> (j+1))) ((xms @> i) - (xns @> (j+1)))
            v0 = atan2 ((yms @> i) - (yns @>  j   )) ((xms @> i) - (xns @>  j   ))
    
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
