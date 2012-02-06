{-# OPTIONS_GHC -Wall #-}

module Main where

import System.Console.Haskeline hiding(display)
import Graphics.Gloss.Interface.IO.Animate hiding(scale, Vector)
import Numeric.LinearAlgebra hiding(i)
import Control.Monad.IO.Class
import Control.Concurrent(forkIO)
import Control.Concurrent.MVar(newMVar, readMVar, swapMVar)

import HFoil.Panels
import HFoil.Naca4
import HFoil.Drawing(toPic)
import HFoil.Flow

-- configuration
normalLengths :: Double
normalLengths = 0.01

nPanels :: Int
nPanels = 200

cpScale :: Double
cpScale = -0.3

xSize, ySize :: Int
xSize = 800
ySize = 500

main :: IO ()
main = do
  let naca0 = "2412"
  mpics <- newMVar $ drawFlow (Foil (toNodes (naca4 naca0) nPanels) naca0) (pi/180*4)
  
  putStrLn "Welcome to hfoil\n"
  
  _ <- forkIO $ runInputT defaultSettings (topLoop (\pics -> swapMVar mpics pics >>= (\_ -> return ())))

  animateIO
    (InWindow
     "hfoil"             -- window title
     (xSize, ySize)      -- window size
     (10, 650))          -- window position
    black                -- background color
    (\_ -> readMVar mpics >>= return . pictures) -- draw function

data Foil a = Foil [(a,a)] String

drawFoil :: Real a => Foil a -> Picture
drawFoil (Foil nodes _) = toPic white nodes

drawFlow :: Foil Double -> Double -> [Picture]
drawFlow foil@(Foil nodes _) alphaDeg = (drawFoil foil):(toPic white nodes):
                                        (map (toPic green) (unitNormals normalLengths nodes))
                                        ++ [toPic red (zip xs mcps)]
  where
    mV = fst $ solveFlow nodes (pi/180*alphaDeg)
    (xs, _) = midpoints nodes
    mcps = toList (scale cpScale (1 - mV*mV))


foilLoop :: ([Picture] -> IO ()) -> Foil Double -> InputT IO ()
foilLoop draw foil@(Foil _ name) = do
  minput <- getInputLine $ "\ESC[1;32m\STXhhfoil."++name++">> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> return ()
    Just ('a':'l':'f':'a':' ':alphaDeg) -> do liftIO $ draw $ drawFlow foil (read alphaDeg)
                                              foilLoop draw foil
    Just "" -> return ()
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     foilLoop draw foil

topLoop :: ([Picture] -> IO ()) -> InputT IO ()
topLoop draw = do
  minput <- getInputLine "\ESC[1;32m\STXhfoil>> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> return ()
    Just ('n':'a':'c':'a':' ':spec) -> do parseNaca draw spec
                                          topLoop draw
    Just "" -> topLoop draw
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     topLoop draw

parseNaca :: ([Picture] -> IO ()) -> String -> InputT IO ()
parseNaca draw str 
  | length str == 4 = do let foil = Foil (toNodes (naca4 str :: Naca4 Double) nPanels) str
                         liftIO $ draw [drawFoil foil]
                         foilLoop draw foil
  | otherwise = do outputStrLn $ "Not 4 digits"
                   return ()
