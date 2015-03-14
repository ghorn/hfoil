{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module Numeric.HFoil.Drawing
       ( drawSolution
       , drawFoil
       , drawNormals
       , drawForces
       , drawKuttas
       , drawOnce
       ) where

import Numeric.LinearAlgebra hiding( Element, scale, i )
import Foreign.Storable ( Storable )
import qualified Numeric.LinearAlgebra as LA
import Text.Printf ( printf )
import Linear ( V3(..) )

import Vis

import Numeric.HFoil.Flow
import Numeric.HFoil.Foil

cpScale :: Fractional a => a
cpScale = -0.25

normalLengths :: Fractional a => a
normalLengths = 0.01

drawLine :: Num a => Color -> [(a,a)] -> VisObject a
drawLine col coords = Line (map (\(x,y) -> V3 x y 0) coords) col

drawCircle :: Num a => Color -> (a, a) -> a -> VisObject a
drawCircle col (x,y) size = Trans (V3 x y 0) $ Sphere size Solid col
--drawCircle col (x,y) size = scale (border*(fromIntegral xSize)) (border*(fromIntegral xSize))
--                            $ translate (-0.5 + realToFrac x) (realToFrac y)
--                            $ color col
--                            $ Circle size

drawLineV :: (Num a, Storable a) => Color -> (Vector a, Vector a) -> VisObject a
drawLineV col (vx, vy) = Line (zipWith (\x y -> V3 x y 0) (toList vx) (toList vy)) col

drawFoil :: (Num a, Storable a) => Foil a -> VisObject a
drawFoil (Foil elements _) = VisObjects $ map drawElement elements

drawElement :: (Num a, Storable a) => Element a -> VisObject a
drawElement element = drawLineV white (fNodes element)

drawNormals :: Foil Double -> VisObject Double
drawNormals (Foil elements _) = VisObjects $ map (\(xy0, xy1) -> drawLine green [xy0, xy1]) (zip xy0s xy1s)
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + (LA.scale normalLengths xUnitNormal))) (toList (ym + (LA.scale normalLengths yUnitNormal)))

    (xUnitNormal, yUnitNormal) = (\(x,y) -> (vjoin x, vjoin y)) $ unzip $ map fUnitNormals elements
    (xm, ym) = (\(x,y) -> (vjoin x, vjoin y)) $ unzip $ map fMidpoints elements

colorFun :: (Fractional a, Real a) => a -> a -> a -> Color
colorFun min' max' x' = makeColor (1-x) (1-x) x 1
  where
    x = realToFrac $ (x' - min') / (max'-min')

drawKuttas :: (Real a, Fractional a, Storable a) => FlowSol a -> VisObject a
drawKuttas flow = VisObjects $ concatMap (\(k0,k1) -> [circ k0, circ k1]) kis
  where
    kis = solKuttaIndices flow
    (xs',ys') = unzip $ map fMidpoints $ (\(Foil els _) -> els) (solFoil flow)
    xs = vjoin xs'
    ys = vjoin ys'
    circ k = drawCircle yellow (xs @> k, ys @> k) 0.006

drawForces :: FlowSol Double -> VisObject Double
drawForces flow = VisObjects $ map (\(xy0, xy1, cp) -> drawLine (colorFun minCp maxCp cp) [xy0, xy1])
                  $ zip3 xy0s xy1s (toList (solCps flow))
  where
    xy0s = zip (toList xm) (toList ym)
    xy1s = zip (toList (xm + xPressures)) (toList (ym + yPressures))
    (xPressures, yPressures) = (\(x,y) -> (LA.scale c x/lengths, LA.scale c y/lengths)) (solForces flow)
    lengths = vjoin $ map fLengths $ (\(Foil x _) -> x) $ solFoil flow
    (xm, ym) = (\(x,y) -> (vjoin x, vjoin y)) $ unzip $ map fMidpoints $ (\(Foil x _) -> x) $ solFoil flow

    c = 0.1

    maxCp = maxElement (solCps flow)
    minCp = minElement (solCps flow)


drawColoredFoil :: (Num a, Storable a) => [Color] -> Foil a -> VisObject a
drawColoredFoil colors foil@(Foil elements _) = VisObjects $ zipWith drawColoredElement colors' elements
  where
    colors' = groupSomethingByFoil foil colors

drawColoredElement :: (Num a, Storable a) => [Color] -> Element a -> VisObject a
drawColoredElement colors element = VisObjects $ map (\(xy0, xy1, col) -> drawLine col [xy0, xy1]) (zip3 xy0s xy1s colors)
  where
    xys = (\(x,y) -> zip (toList x) (toList y)) $ fNodes element
    xy0s = tail xys
    xy1s = init xys

groupSomethingByFoil :: Storable a => Foil a -> [b] -> [[b]]
groupSomethingByFoil (Foil elements _) somethings = f somethings (map (LA.dim . fAngles) elements)
  where
    f xs (n:ns) = (take n xs):(f (drop n xs) ns)
    f [] []= []
    f _ _ = error "uh oh (groupSomethingByFoil)"

drawSolution :: FlowSol Double -> VisObject Double
drawSolution flow = VisObjects $ onscreenText ++
                               [ drawColoredFoil colors foil
                               , drawCircle white (fst $ solCenterPressure flow, snd $ solCenterPressure flow) 0.006
                               , drawCircle white (fst $ solCenterPressure flow, 0) 0.006
                               , drawCircle green (0.25,0) 0.006
                               ] ++ zipWith (\x y -> drawLineV red (x, y)) xs
                                    (takesV (map LA.dim xs) (LA.scale cpScale cps)) -- cp graph
  where
    foil@(Foil elements name) = solFoil flow
    cps = solCps flow

    xs = map (fst . fMidpoints) elements

    onscreenText =
      zipWith (\s k -> Text2d s (30,fromIntegral $ 30*k) Fixed9By15 (makeColor 1 1 1 1))
      msgs (reverse [1..length msgs])

    msgs = [ name
           , printf ("alpha: %.6f deg") ((solAlpha flow)*180/pi)
           , printf ("Cl: %.6f") (solCl flow)
           , printf ("Cd: %.6f") (solCd flow)
           , printf ("Cm: %.6f (c/4, 0)") (solCm flow)
           ]

    colors = map (colorFun (minElement cps) (maxElement cps)) (toList cps)


drawOnce :: Real a => [VisObject a] -> IO ()
drawOnce pics = display (defaultOpts {optWindowName = "hfoil"}) (VisObjects pics)

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
