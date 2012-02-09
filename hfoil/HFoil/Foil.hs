{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Foil( Foil(..)
                 , Element(..)
                 , panelizeNaca4
                 , loadFoil
                 , getUIUCFoil
                 ) where

import System.Directory(doesFileExist)
import Numeric.LinearAlgebra hiding (Element)
import Data.Tuple.Utils(fst3)
import Network.HTTP(simpleHTTP, getRequest, getResponseBody)
import Foreign.Storable(Storable)

import qualified HFoil.Naca4 as Naca4

data Foil a = Foil [Element a] String

instance (Storable a) => Show (Foil a) where
  show (Foil [el] name) = "{"++name++": " ++ show (1 + dim (fLengths el)) ++ " nodes}"
  show (Foil els  name) = "{"++name++": " ++ show (length els) ++ " elements, " ++
              show nodesPerEl ++ " nodes == "++show (sum nodesPerEl)++" total nodes}"
    where
      nodesPerEl = map (\x -> 1 + dim (fLengths x)) els

data Element a = Element { fNodes :: (Vector a, Vector a)
                         , fLengths :: Vector a
                         , fAngles :: Vector a
                         , fMidpoints :: (Vector a, Vector a)
                         , fTangents :: (Vector a, Vector a)
                         , fNormals :: (Vector a, Vector a)
                         , fUnitNormals :: (Vector a, Vector a)
                         , fInits :: (Vector a, Vector a)
                         , fTails :: (Vector a, Vector a)
                         }

toElement :: (Num (Vector a), RealFloat a, Container Vector a) =>
             [(a, a)] -> Element a
toElement xynodes = Element { fNodes = (xNodes, yNodes)
                            , fLengths = lengths
                            , fAngles = zipVectorWith atan2 yTangents xTangents
                            , fMidpoints = (xMids, yMids)
                            , fTangents = (xTangents, yTangents)
                            , fNormals = (xNormals, yNormals)
                            , fUnitNormals = (xUnitNormals, yUnitNormals)
                            , fInits = (xInits, yInits)
                            , fTails = (xTails, yTails)
                            }
  where
    n = (dim xNodes) - 1
    (xNodes, yNodes) = (\(xs,ys) -> (fromList xs, fromList ys)) $ unzip xynodes
    (xInits, yInits) = (subVector 0 n xNodes, subVector 0 n yNodes)
    (xTails, yTails) = (subVector 1 n xNodes, subVector 1 n yNodes)
    (xTangents, yTangents) = (xTails - xInits, yTails - yInits)
    (xMids, yMids) = (0.5*(xInits + xTails), 0.5*(yInits + yTails))
    (xNormals, yNormals) = (-yTangents, xTangents)
    lengths = mapVector sqrt $ xTangents*xTangents + yTangents*yTangents
    (xUnitNormals, yUnitNormals) = (xNormals/lengths, yNormals/lengths)


getUIUCFoil :: (Num (Vector a), Read a, RealFloat a, Container Vector a) =>
               String -> IO (Foil a)
getUIUCFoil name = do
  let file = "http://www.ae.illinois.edu/m-selig/ads/coord/" ++ name ++ ".dat"
  dl <- simpleHTTP (getRequest file) >>= getResponseBody
  return (parseRawFoil dl name)

parseRawFoil :: (Num (Vector a), Read a, RealFloat a, Container Vector a) =>
                String -> String -> Foil a
parseRawFoil raw
  -- bad data
  | any (\x -> 2 /= length x) (concat groupsOfLines) = error $ "parseRawFoil fail, bad foil data?" ++ show groupsOfLines
  -- single element
  | length elements > 0 = Foil (map toElement elements)
  | otherwise =  error "parseRawFoil fail, bad foil date?"
  where
    elements = map (map (\(x:y:[]) -> (read x, read y))) groupsOfLines
    
    -- group lines by splitting at empty lines
    groupsOfLines :: [[[String]]]
    groupsOfLines = f (lines raw)
      where
        f [] = []
        f ([]:xs) = f xs
        f xs = (map words fst'):(f snd')
          where
            (fst', snd') = break (\x -> 0 == length x) xs
  
loadFoil :: FilePath -> IO (Maybe (Foil Double))
loadFoil filename = do
  exists <- doesFileExist filename
  if exists
    then do rawData <- readFile filename
            return (Just (parseRawFoil rawData filename)) -- use filename as name
    else do putStrLn $ "ERROR: file \"" ++ filename ++ "\" couldn't be found"
            return Nothing


panelizeNaca4 :: (Enum a, Floating (Vector a), RealFloat a, Field a) =>
                Naca4.Naca4 a -> Int -> Foil a
panelizeNaca4 foil nPanels = Foil [toElement $ [(1,0)]++reverse lower++[(0,0)]++upper++[(1,0)]]
                             (Naca4.naca4_name foil)
  where
    (upper, lower) = unzip $ map (Naca4.coords foil) xcs
    xcs = toList $ fst3 $ bunchPanels (Naca4.yt foil) (Naca4.dyt foil) xcs0 0 0
    xcs0 = fromList $ init $ tail $ toList $ linspace nXcs (0,1)
    nXcs = (nPanels + (nPanels `mod` 2)) `div` 2 + 1
    
bunchPanels :: (Enum a, Floating (Vector a), Floating a, Ord a, Field a) =>
               (a -> a) -> (a -> a) -> Vector a -> Int -> Int -> (Vector a, Int, Int)
bunchPanels yt dyt xcs nIter nBadSteps 
  | nIter     > 300  = error "panel buncher exceeded 300 iterations"
  | nBadSteps > 1000 = error "panel buncher exceeded 1000 bad steps"
  | sum (toList (abs deltaXcs)) < 1e-12 = (xcs, nIter, nBadSteps)
  | otherwise                           = bunchPanels yt dyt goodXcs (nIter+1) (nBadSteps + length badOnes)
  where
    (badOnes, goodXcs:_) = break (\xs -> all (>0) (toList xs)) nextXcs
    nextXcs = map (\alpha -> xcs + (scale alpha deltaXcs)) (map (2.0**) [0,-1..])
    deltaXcs = xcsStep yt dyt xcs

xcsStep :: (Enum a, Floating (Vector a), Field a) => (a -> a) -> (a -> a) -> Vector a -> Vector a
xcsStep yt dyt xcs = flatten $ -(linearSolveLS mat2 (asColumn rs))
  where
    n = dim xcs
    
    xs = join [fromList [0],              xcs, fromList[1]]
    ys = join [fromList [0], mapVector yt xcs, fromList[0]]
    dxs = (subVector 1 (n+1) xs) - (subVector 0 (n+1) xs)
    dys = (subVector 1 (n+1) ys) - (subVector 0 (n+1) ys)
    deltas = sqrt (dxs*dxs + dys*dys)
    dydxs = mapVector dyt xcs
    
    mat = r1 - r0
      where
        r0 = fromBlocks [[(1><n)[0,0..]],[diagRect 0 diag0 n n]]
        r1 = fromBlocks [                [diagRect 0 diag1 n n], [(1><n)[0,0..]]]
        
        (diag0, diag1) = (subVector 1 n d0, subVector 0 n d1)
          where
            d0 = (dxs + dys*dy0dxs)/deltas
            d1 = (dxs + dys*dy1dxs)/deltas
            dy0dxs = join [fromList [0], dydxs]
            dy1dxs = join [dydxs, fromList [0]]
    
    zeros = (n><1)[0,0..]
    eye = ident n
    frontBunchingParam = 2.0
    diff = (fromBlocks [[zeros, eye]]) - (fromBlocks [[eye*(1+frontBunchingParam/(fromIntegral n)), zeros]])

    mat2 = diff <> mat
    rs = diff <> deltas
