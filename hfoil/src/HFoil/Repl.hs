{-# OPTIONS_GHC -Wall #-}
-- {-# Language FlexibleContexts #-}

module HFoil.Repl
       ( run
       ) where

import Control.Concurrent ( forkIO )
import Control.Concurrent.MVar ( newMVar, readMVar, swapMVar )
import Control.Monad.IO.Class ( MonadIO, liftIO )
import Control.Monad.Trans.Class ( lift )
import Control.Monad.Trans.State.Strict ( StateT, evalStateT, get, modify )
import Data.List ( isPrefixOf )
import Linear ( Quaternion(..), V3(..) )
import System.Console.Haskeline ( InputT, runInputT, defaultSettings, getInputLine, outputStrLn, setComplete )
import System.Console.Haskeline.Completion ( CompletionFunc, Completion, completeWord, simpleCompletion )
import Text.Read ( readMaybe )

import Vis

import HFoil.Foil
import HFoil.Naca4
import HFoil.Drawing
import HFoil.Flow

nPanels :: Int
nPanels = 200

-- configuration
data Config = Config { confForces :: Bool
                     , confKuttas :: Bool
                     , confNormals :: Bool
                     }

defaultConfig :: Config
defaultConfig = Config { confForces = False
                       , confKuttas = False
                       , confNormals = False
                       }

data Mode = TopMode | FoilMode

foilCommands :: [(String, String)]
foilCommands =
  [ ("alfa", "alfa [#]")
  , ("forces", "forces")
  , ("kuttas", "kuttas")
  , ("normals", "normals")
  , ("help", "help")
  ]

topCommands :: [(String, String)]
topCommands =
  [ ("naca",         "naca xxxx")
  , ("load",   "load [filename]")
  , ("uiuc",  "uiuc [foil name]")
  ]

topHelp :: InputT (StateT Mode IO) ()
topHelp = mapM_ (outputStrLn . snd) topCommands

foilHelp :: InputT (StateT Mode IO) ()
foilHelp = mapM_ (outputStrLn . snd) foilCommands

comp :: CompletionFunc (StateT Mode IO)
comp = completeWord Nothing " \t" searchFunc
  where
    searchFunc :: String -> (StateT Mode IO) [Completion]
    searchFunc str = do
      mode <- get
      let wordList = case mode of
            TopMode  -> map fst topCommands
            FoilMode -> map fst foilCommands
      return $ map simpleCompletion $ filter (str `isPrefixOf`) wordList

run :: IO ()
run = do
  mpics <- newMVar $ []

  putStrLn "Welcome to hfoil\n"

  let go :: InputT (StateT Mode IO) ()
      go = topLoop (\pics -> swapMVar mpics pics >>= (\_ -> return ())) defaultConfig

      settings = setComplete comp defaultSettings
  _ <- forkIO $ flip evalStateT TopMode $ runInputT settings go

  let toScreen xs =
        RotQuat (Quaternion 0 (V3 1 0 0))
        $ Trans (V3 (-0.5) 0 0)
        $ VisObjects xs
      cam0 =
        Camera0
        { phi0 = 90
        , theta0 = 90
        , rho0 = 2
        }

  animateIO
    (defaultOpts {optWindowName = "hfoil", optInitialCamera = Just cam0})
    (\_ -> fmap toScreen (readMVar mpics))

data FoilState =
  FoilState
  { fsFlowSol :: Maybe (FlowSol Double)
  , fsConf :: Config
  , fsFoil :: Foil Double
  }

drawPicture :: MonadIO m => ([VisObject Double] -> IO ()) -> StateT FoilState m ()
drawPicture draw = do
  fs <- get
  let conf = fsConf fs
      foil = fsFoil fs
      normals = case confNormals conf of
        True -> [drawNormals foil]
        False -> []

  case fsFlowSol fs of
    Nothing -> liftIO $ draw (drawFoil foil: normals)
    Just flow -> do
      let forces = case (confForces conf) of
            True -> [drawForces flow]
            False -> []
          kuttas = case (confKuttas conf) of
            True -> [drawKuttas flow]
            False -> []
      liftIO $ draw $ forces++kuttas++normals++[drawSolution flow]

strip :: String -> String
strip = rstrip . lstrip
  where
    lstrip (' ':xs) = lstrip xs
    lstrip x = x

    rstrip = reverse . lstrip . reverse

foilLoop :: ([VisObject Double] -> IO ()) -> StateT FoilState (InputT (StateT Mode IO)) ()
foilLoop draw = do
  lift (lift (modify (const FoilMode)))
  fs <- get
  let foil@(Foil _ name) = fsFoil fs
      conf = fsConf fs
  drawPicture draw
  minput <- lift $ getInputLine $ "\ESC[1;32m\STXhfoil."++name++">> \ESC[0m\STX"

  case fmap strip minput of
    Nothing -> return ()
    Just "quit" -> do lift $ outputStrLn "not-gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
                      foilLoop draw
    Just ('a':'l':'f':'a':' ':alphaDeg') -> do
      case readMaybe alphaDeg' of
        Nothing -> do lift $ outputStrLn $ "parse fail on " ++ show alphaDeg'
                      foilLoop draw
        Just alphaDeg -> do let flow :: FlowSol Double
                                flow = solveFlow foil (pi/180*alphaDeg)
                            modify (\fs' -> fs' {fsFlowSol = Just flow})
                            foilLoop draw
    Just "forces" -> do
      let newConf = conf {confForces = not (confForces conf)}
      lift $ outputStrLn $ "force drawing set to "++ show (not (confForces conf))
      modify (\fs' -> fs' {fsConf = newConf})
      foilLoop draw
    Just "kuttas" -> do
      let newConf = conf {confKuttas = not (confKuttas conf)}
      lift $ outputStrLn $ "kutta drawing set to "++ show (not (confKuttas conf))
      modify (\fs' -> fs' {fsConf = newConf})
      foilLoop draw
    Just "normals" -> do
      let newConf = conf {confNormals = not (confNormals conf)}
      lift $ outputStrLn $ "normals drawing set to "++ show (not (confNormals conf))
      modify (\fs' -> fs' {fsConf = newConf})
      foilLoop draw
    Just "help" -> lift foilHelp >> foilLoop draw
    Just "h"    -> lift foilHelp >> foilLoop draw
    Just "?"    -> lift foilHelp >> foilLoop draw
    Just "" -> return ()
    Just input -> do lift $ outputStrLn $ "unrecognized command \"" ++ input ++ "\""
                     foilLoop draw


topLoop :: ([VisObject Double] -> IO ()) -> Config -> InputT (StateT Mode IO) ()
topLoop draw conf = do
  lift (modify (const TopMode))
  minput <- getInputLine "\ESC[1;32m\STXhfoil>> \ESC[0m\STX"
  case minput of
    Nothing -> return ()
    Just msg -> do runTop draw conf msg
                   topLoop draw conf

runTop :: ([VisObject Double] -> IO ()) -> Config -> String -> InputT (StateT Mode IO) ()
runTop draw conf msg = case strip msg of
  "quit" -> outputStrLn "not-gloss won't let you quit :(\ntry ctrl-c or hit ESC in drawing window"
  ('n':'a':'c':'a':' ':spec) -> do
    case naca4 spec :: Maybe (Naca4 Double) of
      Nothing -> outputStrLn "not a valid naca4"
      Just n4 -> runFoil draw conf (panelizeNaca4 n4 nPanels)
  ('l':'o':'a':'d':' ':name) -> do
    mfoil <- liftIO (loadFoil name)
    case mfoil of Left errMsg -> outputStrLn errMsg
                  Right foil -> do runFoil draw conf foil
  ('u':'i':'u':'c':' ':name) -> do
    efoil <- liftIO (getUIUCFoil name)
    case efoil of Left errMsg -> outputStrLn errMsg
                  Right foil -> do let Foil els _ = foil
                                   outputStrLn $ "got " ++ show (length els) ++ " elements"
                                   runFoil draw conf foil
  "help" -> topHelp
  "h"    -> topHelp
  "?"    -> topHelp
  "" -> return ()
  other -> outputStrLn $ "unrecognized command \"" ++ other ++ "\""


runFoil :: ([VisObject Double] -> IO ()) -> Config -> Foil Double -> InputT (StateT Mode IO) ()
runFoil draw conf foil = do
  let state0 =
        FoilState
        { fsFlowSol = Nothing
        , fsConf = conf
        , fsFoil = foil
        }
  let go :: StateT FoilState (InputT (StateT Mode IO)) ()
      go = foilLoop draw
  flip evalStateT state0 go
