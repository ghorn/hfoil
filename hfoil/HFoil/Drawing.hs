{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Drawing( drawLine
                    , drawLineV
                    , drawSolution
                    , drawFoil
                    , drawOnce
                    , drawNormals
                    , drawForces
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

drawCircle :: Real a => Color -> (a, a) -> Float -> Picture
drawCircle col (x,y) size = scale (border*(fromIntegral xSize)) (border*(fromIntegral xSize))
                            $ translate (-0.5 + realToFrac x) (realToFrac y)
                            $ color col
                            $ Circle size

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

colorFun :: (Fractional a, Real a) => a -> a -> a -> Color
colorFun min' max' x' = makeColor (1-x) (1-x) x 1
  where
    x = realToFrac $ (x' - min')/(max'-min')

drawForces :: FlowSol Double -> Picture
drawForces flow = pictures $ map (\(xy0, xy1, cp) -> drawLine (colorFun minCp maxCp cp) [xy0, xy1])
                  $ zip3 xy0s xy1s (toList (fsCps flow))
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + xPressures)) (toList (ym + yPressures))
    (xPressures, yPressures) = (\(x,y) -> (LA.scale c x/lengths, LA.scale c y/lengths)) (fsForces flow)
    lengths = pLengths $ fsFoil flow
    (xm, ym) = pMidpoints $ fsFoil flow
    
    c = 0.1
    
    maxCp = maxElement (fsCps flow)
    minCp = minElement (fsCps flow)

drawColoredFoil :: [Color] -> Foil Double -> Picture
drawColoredFoil colors foil = pictures $ map (\(xy0, xy1, col) -> drawLine col [xy0, xy1]) (zip3 xy0s xy1s colors)
  where
    xys = (\(x,y) -> zip (toList x) (toList y)) $ pNodes foil
    xy0s = tail xys
    xy1s = init xys

drawSolution :: FlowSol Double -> Picture
drawSolution flow = pictures [ drawText white (0.45, 0.8) 0.15 m0
                             , drawText white (0.45, 0.65) 0.15 m1
                             , drawText white (0.45, 0.5) 0.15 m2
                             , drawText white (0.45, 0.35) 0.15 m3
                             , drawForces flow
                             , drawColoredFoil colors foil
                             , drawLineV red (xs, LA.scale cpScale cps) -- cp graph
                             , drawCircle white (fst $ fsCenterPressure flow, snd $ fsCenterPressure flow) 0.006
                             , drawCircle white (fst $ fsCenterPressure flow, 0) 0.006
                             ]
  where
    foil = fsFoil flow
    cps = fsCps flow
    
    (xs, _) = pMidpoints foil
    
    [m0,m1,m2,m3] = [ pName foil
                    , printf ("alpha: %.6f") ((fsAlpha flow)*180/pi)
                    , printf ("Cl: %.6f") (fsCl flow)
                    , printf ("Cd: %.6f") (fsCd flow)
                    ]
        
    colors = map (colorFun (minElement cps) (maxElement cps)) (toList cps)



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
