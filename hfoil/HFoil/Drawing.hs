{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Drawing( drawLine
                    , drawLineV
                    , drawFlow
                    , drawFoil
                    , drawOnce
                    , drawNormals
                    , drawForces
                    , drawCps
                    ) where

import Graphics.Gloss hiding(Vector)
import Numeric.LinearAlgebra hiding(scale,i)
import Foreign.Storable(Storable)
import qualified Numeric.LinearAlgebra as LA
import Text.Printf

import HFoil.Flow
import HFoil.Foil

xSize, ySize :: Int
xSize = 800
ySize = 500

cpScale :: Double
cpScale = -0.25

border :: Float
border = 0.7

normalLengths :: Double
normalLengths = 0.01

drawText :: Color -> (Float, Float) -> Float -> String -> Picture
drawText col (x,y) size str = translate (0.5*x*(fromIntegral xSize)) (0.5*y*(fromIntegral ySize))
                              $ scale size size
                              $ color col
                              $ Text str

drawLine :: Real a => Color -> [(a,a)] -> Picture
drawLine col coords = scale (border*(fromIntegral xSize)) (border*(fromIntegral xSize))
                      $ translate (-0.5) 0
                      $ color col
                      $ line $ map (\(x,y) -> (realToFrac x, realToFrac y)) coords

drawLineV :: (Real a, Storable a) => Color -> (Vector a, Vector a) -> Picture
drawLineV col (vx, vy) = drawLine col $ zip (toList vx) (toList vy)

drawFoil :: (Real a, Storable a, Show a) => Foil a -> Picture
drawFoil foil = drawLineV white (pNodes foil)

drawNormals :: Foil Double -> Picture
drawNormals foil = pictures $ map (\(xy0, xy1) -> drawLine green [xy0, xy1]) (zip xy0s xy1s)
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + (LA.scale normalLengths xUnitNormal))) (toList (ym + (LA.scale normalLengths yUnitNormal)))
    
    (xUnitNormal, yUnitNormal) = pUnitNormals foil
    (xm, ym) = pMidpoints foil
    
drawForces :: FlowSol Double -> Picture
drawForces flow = pictures $ map (\(xy0, xy1) -> drawLine blue [xy0, xy1]) (zip xy0s xy1s)
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + cps*xUnitNormal)) (toList (ym + (cps*yUnitNormal)))
    cps = LA.scale (0.5*cpScale) (fsCps flow)
    (xUnitNormal, yUnitNormal) = pUnitNormals $ fsFoil flow
    (xm, ym) = pMidpoints $ fsFoil flow

drawFlow :: FlowSol Double -> Picture
drawFlow flow = pictures [drawFoil foil, drawCps flow, drawForces flow]
  where
    foil = fsFoil flow

drawCps :: FlowSol Double -> Picture
drawCps flow = pictures $ [ drawLineV red (xs, mcps)
                          , drawText white (0.45, 0.8) 0.15 m0
                          , drawText white (0.45, 0.65) 0.15 m1
                          , drawText white (0.45, 0.5) 0.15 m2
                          , drawText white (0.45, 0.35) 0.15 m3
                          ]
  where
    foil = fsFoil flow
    cps = fsCps flow
    
    (xs, _) = pMidpoints foil
    mcps = LA.scale cpScale cps
    
    m0 = pName foil
    m1 = printf ("alpha: %.6f") (q*180/pi)
    m2 = printf ("Cl: %.6f") cl
    m3 = printf ("Cd: %.6f") cd
    
    xForces = -cps*(fst $ pNormals foil)
    yForces = -cps*(snd $ pNormals foil)
    xf = sumElements xForces
    yf = sumElements yForces
    
    q = fsAlpha flow
    cd =  xf*(cos q) + yf*(sin q)
    cl = -xf*(sin q) + yf*(cos q)
        

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
