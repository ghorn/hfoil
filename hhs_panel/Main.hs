{-# OPTIONS_GHC -Wall #-}

module Main where
import HFoil.Panels
import HFoil.Naca4
import HFoil.Drawing

normalLength :: Double
normalLength = 0.01

nPanels :: Int
nPanels = 200

main :: IO ()
main = draw $ map toPicture (nodes:normals)

nodes :: [(Double, Double)]
nodes = toNodes (naca4 "0012") nPanels

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
