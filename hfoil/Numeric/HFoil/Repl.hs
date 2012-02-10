{-# OPTIONS_GHC -Wall #-}

module Numeric.HFoil.Repl( run
                         ) where

import System.Console.Haskeline hiding(display)
import Graphics.Gloss.Interface.IO.Animate hiding(scale, Vector)
import Control.Monad.IO.Class
import Control.Concurrent(forkIO)
import Control.Concurrent.MVar(newMVar, readMVar, swapMVar)

import Numeric.HFoil.Foil
import Numeric.HFoil.Naca4
import Numeric.HFoil.Drawing
import Numeric.HFoil.Flow

---- configuration
nPanels :: Int
nPanels = 200

xSize, ySize :: Int
xSize = 800
ySize = 500

run :: IO ()
run = do
  let naca0 = "2412"
      alfaDeg0 = 4
      flow0 = solveFlow (panelizeNaca4 (naca4 naca0) nPanels) (pi/180*alfaDeg0)
  mpics <- newMVar $ [drawSolution flow0]
  
  putStrLn "Welcome to hfoil\n"
  
  _ <- forkIO $ runInputT defaultSettings (topLoop (\pics -> swapMVar mpics pics >>= (\_ -> return ())))

  animateIO
    (InWindow
     "hfoil"             -- window title
     (xSize, ySize)      -- window size
     (10, 650))          -- window position
    black                -- background color
    (\_ -> readMVar mpics >>= return . pictures) -- draw function

foilLoop :: ([Picture] -> IO ()) -> Foil Double -> InputT IO ()
foilLoop draw foil@(Foil _ name) = do
  minput <- getInputLine $ "\ESC[1;32m\STXhfoil."++name++">> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> do outputStrLn "gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
                      foilLoop draw foil
    Just ('a':'l':'f':'a':' ':[]) -> do outputStrLn $ "unrecognized command"
                                        foilLoop draw foil
    Just ('a':'l':'f':'a':' ':alphaDeg) -> do let flow = solveFlow foil (pi/180*(read alphaDeg))
                                              liftIO $ draw $ [drawSolution flow]
                                              foilLoop draw foil
    Just "" -> return ()
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     foilLoop draw foil

topLoop :: ([Picture] -> IO ()) -> InputT IO ()
topLoop draw = do
  minput <- getInputLine "\ESC[1;32m\STXhfoil>> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> do outputStrLn "gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
                      topLoop draw
    Just ('n':'a':'c':'a':' ':spec) -> do parseNaca draw spec
                                          topLoop draw
    Just ('l':'o':'a':'d':' ':name) -> do
      foil <- liftIO (loadFoil name)
      case foil of Left errMsg -> outputStrLn errMsg
                   Right foil' -> do liftIO $ draw [drawFoil foil', drawNormals foil']
                                     foilLoop draw foil'
      topLoop draw
    Just ('u':'i':'u':'c':' ':name) -> do
      efoil <- liftIO (getUIUCFoil name)
      case efoil of Left errMsg -> outputStrLn errMsg
                    Right foil -> do liftIO $ draw [drawFoil foil, drawNormals foil]
                                     foilLoop draw foil
      topLoop draw

    Just "" -> topLoop draw
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     topLoop draw

parseNaca :: ([Picture] -> IO ()) -> String -> InputT IO ()
parseNaca draw str 
  | length str == 4 = do let foil = panelizeNaca4 (naca4 str :: Naca4 Double) nPanels
                         liftIO $ draw [drawFoil foil, drawNormals foil]
                         foilLoop draw foil
  | otherwise = do outputStrLn $ "Not 4 digits"
                   return ()