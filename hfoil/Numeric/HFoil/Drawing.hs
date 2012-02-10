{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module Numeric.HFoil.Drawing( drawLine
                            , drawLineV
                            , drawSolution
                            , drawFoil
                            , drawOnce
                            , drawNormals
                            , drawForces
                            ) where

import Graphics.Gloss hiding(Vector,dim)
import Numeric.LinearAlgebra hiding(Element, scale,i)
import Foreign.Storable(Storable)
import qualified Numeric.LinearAlgebra as LA
import Text.Printf

import Numeric.HFoil.Flow
import Numeric.HFoil.Foil

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
drawFoil (Foil elements _) = pictures $ map drawElement elements

drawElement :: (Real a, Storable a) => Element a -> Picture
drawElement element = drawLineV white (fNodes element)

drawNormals :: Foil Double -> Picture
drawNormals (Foil elements _) = pictures $ map (\(xy0, xy1) -> drawLine green [xy0, xy1]) (zip xy0s xy1s)
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + (LA.scale normalLengths xUnitNormal))) (toList (ym + (LA.scale normalLengths yUnitNormal)))
    
    (xUnitNormal, yUnitNormal) = (\(x,y) -> (join x, join y)) $ unzip $ map fUnitNormals elements
    (xm, ym) = (\(x,y) -> (join x, join y)) $ unzip $ map fMidpoints elements

colorFun :: (Fractional a, Real a) => a -> a -> a -> Color
colorFun min' max' x' = makeColor (1-x) (1-x) x 1
  where
    x = realToFrac $ (x' - min')/(max'-min')

drawForces :: FlowSol Double -> Picture
drawForces flow = pictures $ map (\(xy0, xy1, cp) -> drawLine (colorFun minCp maxCp cp) [xy0, xy1])
                  $ zip3 xy0s xy1s (toList (solCps flow))
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + xPressures)) (toList (ym + yPressures))
    (xPressures, yPressures) = (\(x,y) -> (LA.scale c x/lengths, LA.scale c y/lengths)) (solForces flow)
    lengths = join $ map fLengths $ (\(Foil x _) -> x) $ solFoil flow
    (xm, ym) = (\(x,y) -> (join x, join y)) $ unzip $ map fMidpoints $ (\(Foil x _) -> x) $ solFoil flow
    
    c = 0.1
    
    maxCp = maxElement (solCps flow)
    minCp = minElement (solCps flow)

drawColoredFoil :: [Color] -> Foil Double -> Picture
drawColoredFoil colors foil@(Foil elements _) = pictures $ zipWith drawColoredElement colors' elements
  where
    colors' = groupSomethingByFoil foil colors

drawColoredElement :: [Color] -> Element Double -> Picture
drawColoredElement colors element = pictures $ map (\(xy0, xy1, col) -> drawLine col [xy0, xy1]) (zip3 xy0s xy1s colors)
  where
    xys = (\(x,y) -> zip (toList x) (toList y)) $ fNodes element
    xy0s = tail xys
    xy1s = init xys

groupSomethingByFoil :: Storable a => Foil a -> [b] -> [[b]]
groupSomethingByFoil (Foil elements _) somethings = f somethings (map (dim . fAngles) elements)
  where
    f xs (n:ns) = (take n xs):(f (drop n xs) ns)
    f [] []= []
    f _ _ = error "uh oh (groupSomethingByFoil)"

drawSolution :: FlowSol Double -> Picture
drawSolution flow = pictures $ onscreenText ++
                               [ drawForces flow
                               , drawColoredFoil colors foil
                               , drawCircle white (fst $ solCenterPressure flow, snd $ solCenterPressure flow) 0.006
                               , drawCircle white (fst $ solCenterPressure flow, 0) 0.006
                               ] ++ zipWith (\x y -> drawLineV red (x, y)) xs
                                    (takesV (map dim xs) (LA.scale cpScale cps)) -- cp graph
  where
    foil@(Foil elements name) = solFoil flow
    cps = solCps flow
    
    xs = map (fst . fMidpoints) elements
    
    onscreenText = zipWith (\m y -> drawText white (0.45,y) 0.15 m) msgs
                   $ take (length msgs) [0.8,0.65..]
    msgs = [ name
           , printf ("alpha: %.6f") ((solAlpha flow)*180/pi)
           , printf ("Cl: %.6f") (solCl flow)
           , printf ("Cd: %.6f") (solCd flow)
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
