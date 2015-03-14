{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Foil
       ( Foil(..)
       , Element(..)
       , panelizeNaca4
       , loadFoil
       , getUIUCFoil
       ) where

import System.Directory ( doesFileExist )
import Numeric.LinearAlgebra hiding ( Element )
import Network.HTTP ( simpleHTTP, getRequest, getResponseBody )
import Foreign.Storable ( Storable )
import Text.Read ( readMaybe )
import Text.Parsec ( Parsec, ParsecT, Stream, (<?>), parse
                   , option, char, try, lookAhead, manyTill, digit, endOfLine
                   , many1, optional, anyChar, sepBy1, many, oneOf, skipMany )

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


-- make sure the nodes aren't reversed
toElement :: (Num (Vector a), RealFloat a, Container Vector a)
             => [(a, a)] -> Element a
toElement nodes
  | diff < 0 = el
  | otherwise = toElement' (reverse nodes)
  where
    el = toElement' nodes
    angles = toList (fAngles el)
    diff = sum $ zipWith (-) angles (drop 1 angles)

toElement' :: (Num (Vector a), RealFloat a, Container Vector a)
              => [(a, a)] -> Element a
toElement' xynodes =
  Element { fNodes = (xNodes, yNodes)
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


-- why isn't this standard???
poorMansStrip :: String -> String
poorMansStrip str = reverse $ dropWhile (== ' ') $ reverse $ dropWhile (== ' ') str

getUIUCFoil :: String -> IO (Either String (Foil Double))
getUIUCFoil name' = do
  let name = poorMansStrip name'
      file = "http://m-selig.ae.illinois.edu/ads/coord/" ++ name ++ ".dat"
  dl <- simpleHTTP (getRequest file) >>= getResponseBody
  return (parseRawFoil dl name)

loadFoil :: FilePath -> IO (Either String (Foil Double))
loadFoil filename' = do
  let filename = poorMansStrip filename'
  exists <- doesFileExist filename
  if exists
    then do rawData <- readFile filename
            return (parseRawFoil rawData filename) -- use filename as name
    else return (Left ("file \"" ++ filename ++ "\" couldn't be found"))

fst3 :: (a,b,c) -> a
fst3 (x,_,_) = x

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

    xs = vjoin [fromList [0],              xcs, fromList[1]]
    ys = vjoin [fromList [0], mapVector yt xcs, fromList[0]]
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
            dy0dxs = vjoin [fromList [0], dydxs]
            dy1dxs = vjoin [dydxs, fromList [0]]

    zeros = (n><1)[0,0..]
    eye = ident n
    frontBunchingParam = 2.0
    diff = (fromBlocks [[zeros, eye]]) - (fromBlocks [[eye*(1+frontBunchingParam/(fromIntegral n)), zeros]])

    mat2 = diff <> mat
    rs = diff <> deltas

parseRawFoil :: String -> String -> Either String (Foil Double)
parseRawFoil raw name = foil
  where
    foil = case parse foilP name raw of
      Left pe -> Left (raw ++ "\nError parsing the above data: " ++ (show pe))
      Right (_, els) -> Right (Foil (map toElement els) name)


headerP :: Parsec String () String
headerP = manyTill anyChar (try (lookAhead elementsP))

foilP :: Parsec String () (String, [[(Double, Double)]])
foilP = do
  s <- headerP <?> "header"
  els <- elementsP <?> "elements"
  return (s, els)

elementP :: Parsec String () [(Double, Double)]
elementP = many2 coordP
  where
    many2 :: (Stream s m t) => ParsecT s u m a -> ParsecT s u m [a]
    many2 p = do
      x0 <- p
      x1 <- p
      xs <- many p
      return (x0:x1:xs)


elementsP :: Parsec String () [[(Double, Double)]]
elementsP = sepBy1 elementP endOfLine

doubleP :: Parsec String () Double
doubleP = do
  neg <- option ' ' (char '-')
  leading <- option "0" (many1 digit)
  dot' <- char '.'
  trailing <- option "0" (many1 digit)
  let doublish = neg : leading ++ dot' : trailing
  case readMaybe doublish of
    Just x -> return x
    Nothing -> error $ "failed to read this supposed double: " ++ show doublish

coordP :: Parsec String () (Double, Double)
coordP = do
  let space' = oneOf [' ', '\t']
      spaces' = skipMany space' <?> "space-like"
  spaces'      <?> "leading whitespace"
  x <- doubleP <?> "first coordinate"
  spaces'      <?> "middle whitespace"
  y <- doubleP <?> "second coordinate"
  spaces'      <?> "trailing whitespace"
  _ <- optional endOfLine <?> "end of line"
  return (x,y)
