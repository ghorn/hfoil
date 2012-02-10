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

data Config = Config { confForces :: Bool
                     , confKuttas :: Bool
                     , confNormals :: Bool
                     }

defaultConfig :: Config
defaultConfig = Config { confForces = False
                       , confKuttas = False
                       , confNormals = False
                       }

run :: IO ()
run = do
  let naca0 = "2412"
      alfaDeg0 = 4
      flow0 = solveFlow (panelizeNaca4 (naca4 naca0) nPanels) (pi/180*alfaDeg0)
  mpics <- newMVar $ [drawSolution flow0]
  
  putStrLn "Welcome to hfoil\n"
  
  _ <- forkIO $ runInputT defaultSettings
       $ topLoop (\pics -> swapMVar mpics pics >>= (\_ -> return ())) defaultConfig

  animateIO
    (InWindow
     "hfoil"             -- window title
     (xSize, ySize)      -- window size
     (10, 650))          -- window position
    black                -- background color
    (\_ -> readMVar mpics >>= return . pictures) -- draw function

foilLoop :: ([Picture] -> IO ()) -> Config -> Foil Double -> InputT IO ()
foilLoop draw conf foil@(Foil _ name) = do
  minput <- getInputLine $ "\ESC[1;32m\STXhfoil."++name++">> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> do outputStrLn "gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
                      foilLoop draw conf foil
    Just ('a':'l':'f':'a':' ':[]) -> do outputStrLn $ "unrecognized command"
                                        foilLoop draw conf foil
    Just ('a':'l':'f':'a':' ':alphaDeg) -> do let flow = solveFlow foil (pi/180*(read alphaDeg))
                                                  forces = case (confForces conf) of
                                                    True -> [drawForces flow]
                                                    False -> []
                                                  kuttas = case (confKuttas conf) of
                                                    True -> [drawKuttas flow]
                                                    False -> []
                                                  normals = case (confNormals conf) of
                                                    True -> [drawNormals (solFoil flow)]
                                                    False -> []
                                              liftIO $ draw $ forces++kuttas++normals++[drawSolution flow]
                                              foilLoop draw conf foil
    Just ('f':'o':'r':'c':'e':'s':[]) -> do
      let newConf = conf {confForces = not (confForces conf)}
      outputStrLn $ "force drawing set to "++ show (not (confForces conf))
      foilLoop draw newConf foil
    Just ('k':'u':'t':'t':'a':'s':[]) -> do 
      let newConf = conf {confKuttas = not (confKuttas conf)}
      outputStrLn $ "kutta drawing set to "++ show (not (confKuttas conf))
      foilLoop draw newConf foil
    Just ('n':'o':'r':'m':'a':'l':'s':[]) -> do
      let newConf = conf {confNormals = not (confNormals conf)}
      outputStrLn $ "normals drawing set to "++ show (not (confNormals conf))
      foilLoop draw newConf foil
    Just "" -> return ()
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     foilLoop draw conf foil

topLoop :: ([Picture] -> IO ()) -> Config -> InputT IO ()
topLoop draw conf = do
  minput <- getInputLine "\ESC[1;32m\STXhfoil>> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just "quit" -> do outputStrLn "gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
                      topLoop draw conf
    Just ('n':'a':'c':'a':' ':spec) -> do parseNaca draw conf spec
                                          topLoop draw conf
    Just ('l':'o':'a':'d':' ':name) -> do
      foil <- liftIO (loadFoil name)
      case foil of Left errMsg -> outputStrLn errMsg
                   Right foil' -> do liftIO $ draw [drawFoil foil', drawNormals foil']
                                     foilLoop draw conf foil'
      topLoop draw conf
    Just ('u':'i':'u':'c':' ':name) -> do
      efoil <- liftIO (getUIUCFoil name)
      case efoil of Left errMsg -> outputStrLn errMsg
                    Right foil -> do liftIO $ draw [drawFoil foil, drawNormals foil]
                                     foilLoop draw conf foil
      topLoop draw conf

    Just "" -> topLoop draw conf
    Just input -> do outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     topLoop draw conf

parseNaca :: ([Picture] -> IO ()) -> Config -> String -> InputT IO ()
parseNaca draw conf str 
  | length str == 4 = do let foil = panelizeNaca4 (naca4 str :: Naca4 Double) nPanels
                         liftIO $ draw [drawFoil foil, drawNormals foil]
                         foilLoop draw conf foil
  | otherwise = do outputStrLn $ "Not 4 digits"
                   return ()
