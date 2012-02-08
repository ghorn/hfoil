{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Drawing( drawLine
                    , drawLineV
                    , drawFlow
                    , drawFoil
                    , drawOnce
                    , drawNormals
                    ) where

import Graphics.Gloss hiding(Vector)
import Numeric.LinearAlgebra hiding(scale,i)
import Foreign.Storable(Storable)
import qualified Numeric.LinearAlgebra as LA

import HFoil.Foil

xSize, ySize :: Int
xSize = 800
ySize = 500

cpScale :: Double
cpScale = -0.3

normalLengths :: Double
normalLengths = 0.01

drawLine :: Real a => Color -> [(a,a)] -> Picture
drawLine col coords = scale (0.8*(fromIntegral xSize)) (0.8*(fromIntegral xSize))
                      $ translate (-0.5) 0
                      $ color col
                      $ line $ map (\(x,y) -> (realToFrac x, realToFrac y)) coords

drawLineV :: (Real a, Storable a) => Color -> (Vector a, Vector a) -> Picture
drawLineV col (vx, vy) = drawLine col $ zip (toList vx) (toList vy)

drawFoil :: (Real a, Storable a) => Foil a -> Picture
drawFoil foil = drawLineV white (pNodes foil)

drawNormals :: Foil Double -> Picture
drawNormals foil = pictures $ map (\(xy0, xy1) -> drawLine green [xy0, xy1]) (zip xy0s xy1s)
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + (LA.scale normalLengths xUnitNormal))) (toList (ym + (LA.scale normalLengths yUnitNormal)))
    
    (xUnitNormal, yUnitNormal) = pUnitNormals foil
    (xm, ym) = pMidpoints foil
    
drawFlow :: Foil Double -> (Vector Double, b) -> [Picture]
drawFlow foil flowSolution = [drawFoil foil, drawNormals foil, drawLineV red (xs, mcps)]
  where
    mV = fst $ flowSolution
    (xs, _) = pMidpoints foil
    mcps = LA.scale cpScale (1 - (mV*mV))

drawOnce :: [Picture] -> IO ()
drawOnce pics = do
  display 
    (InWindow
     "hfoil"             -- window title
     (xSize, ySize)      -- window size
     (10, 650))          -- window position
    black                -- background color
    (pictures pics)      -- picture to display

--  let line = plot_lines_values ^= [[ (xc, yt (naca4 "0012") xc)
--                                   | xc <- [0,0.01..0.99::Double]]]
--             $ plot_lines_title ^= "naca 0012"
--             $ defaultPlotLines
--  
--      chart = layout1_title ^= "naca yo"
--              $ layout1_plots ^= [Left (toPlot line)]
--              $ defaultLayout1
--  
--  renderableToWindow (toRenderable chart) 640 480
--  _ <- renderableToPNGFile (toRenderable chart) 640 480 "mDiv_vs_tc.png"
--  return ()
