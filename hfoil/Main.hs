{-# OPTIONS_GHC -Wall #-}

module Main where

import Graphics.Gloss hiding(scale)
import Numeric.LinearAlgebra hiding(i)

import HFoil.Panels
import HFoil.Naca4
import HFoil.Drawing
import HFoil.Flow

normalLength :: Double
normalLength = 0.01

alphaDeg :: Double
alphaDeg = 4*pi/180

nPanels :: Int
nPanels = 200

cpScale :: Double
cpScale = -0.3

nodes :: [(Double, Double)]
nodes = toNodes (naca4 "2412") nPanels

main :: IO ()
main = do
  draw $ (toPic white nodes):(map (toPic green) normals)++[toPic red (zip xs mcps)]
    where
      (xs, _) = midpoints nodes
      mcps = toList (scale cpScale (1 - mV*mV))
        where
          (mV, _) = solveFlow nodes alphaDeg

normals :: [[(Double, Double)]]
normals = zipWith f (tail nodes) (init nodes)
  where
    f (x1,y1) (x0,y0) = [(mx, my), (mx - dy/len, my + dx/len)]
      where
        mx = 0.5*(x1+x0)
        my = 0.5*(y1+y0)
        dx = x1 - x0
        dy = y1 - y0
        len = sqrt(dx*dx + dy*dy)/normalLength
