{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Foil( Foil(..)
                 , toFoil
                 , panelizeNaca4
                 , loadFoil
                 , getUIUCFoil
                 ) where

import System.Directory(doesFileExist)
import Numeric.LinearAlgebra
import Data.Tuple.Utils(fst3)
import Network.HTTP(simpleHTTP, getRequest, getResponseBody)

import qualified HFoil.Naca4 as Naca4

type Panels a = [(a,a)]

data Elements a = SingleElement (Panels a)
                | MultiElement [Panels a]

instance Show (Elements a) where
  show (SingleElement x) = "{SingleElement: " ++ show (length x) ++ " nodes}"
  show (MultiElement x) = "{MultiElement: " ++ show (length x) ++ " elements, " ++
                          show (map length x) ++ " nodes == "++show (sum (map length x))++" total nodes}"

data Foil a = Foil { pNodes :: (Vector a, Vector a)
                   , pLengths :: Vector a
                   , pAngles :: Vector a
                   , pMidpoints :: (Vector a, Vector a)
                   , pTangents :: (Vector a, Vector a)
                   , pNormals :: (Vector a, Vector a)
                   , pUnitNormals :: (Vector a, Vector a)
                   , pName :: String
                   }

toFoil :: (Num (Vector a), RealFloat a, Container Vector a) =>
          String -> [(a, a)] -> Foil a
toFoil name xynodes = Foil { pNodes = (xNodes, yNodes)
                           , pLengths = lengths
                           , pAngles = zipVectorWith atan2 yTangents xTangents
                           , pMidpoints = (xMids, yMids)
                           , pTangents = (xTangents, yTangents)
                           , pNormals = (xNormals, yNormals)
                           , pUnitNormals = (xUnitNormals, yUnitNormals)
                           , pName = name
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


getUIUCFoil :: Read a => String -> IO (Elements a)
getUIUCFoil name = do
  let file = "http://www.ae.illinois.edu/m-selig/ads/coord/" ++ name ++ ".dat"
  dl <- simpleHTTP (getRequest file) >>= getResponseBody
  return (parseRawFoil dl)

parseRawFoil :: Read a => String -> Elements a
parseRawFoil raw
  -- bad data
  | any (\x -> 2 /= length x) (concat groupsOfLines) = error $ "parseRawFoil fail, bad foil data?" ++ show groupsOfLines
  -- single element
  | length elements == 1 = SingleElement (head elements)
  -- multi element
  | length elements > 1 = MultiElement elements
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
  
loadFoil :: FilePath -> IO (Maybe (Elements Double))
loadFoil filename = do
  exists <- doesFileExist filename
  if exists
    then do rawData <- readFile filename
            return (Just (parseRawFoil rawData))
    else do putStrLn $ "ERROR: file \"" ++ filename ++ "\" couldn't be found"
            return Nothing


panelizeNaca4 :: (Enum a, Floating (Vector a), RealFloat a, Field a) =>
                Naca4.Naca4 a -> Int -> Foil a
panelizeNaca4 foil nPanels = toFoil (Naca4.naca4_name foil) $ [(1,0)]++reverse lower++[(0,0)]++upper++[(1,0)]
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
