{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}

module HFoil.Flow( FlowSol(..)
                 , solveFlow
                 ) where

import Data.Packed.ST(newMatrix, writeMatrix, freezeMatrix)
import Control.Monad(forM_)
import Control.Monad.ST(runST)
import Numeric.LinearAlgebra hiding(i)
import Foreign.Storable

import HFoil.Foil

data FlowSol a = FlowSol { solFoil :: Foil a
                         , solVs :: Vector a
                         , solCps :: Vector a
                         , solVorticities :: [a]
                         , solAlpha :: a
                         , solForces :: (Vector a, Vector a)
                         , solCl :: a
                         , solCd :: a
                         , solCenterPressure :: (a,a)
                         }

solveFlow :: (Num (Vector a), RealFloat a, Field a) => Foil a -> a -> FlowSol a
solveFlow foil@(Foil elements _) alpha =
  FlowSol { solFoil = foil
          , solVs = vs
          , solCps = cps
          , solVorticities = vorticities
          , solAlpha = alpha
          , solForces = (xForces, yForces)
          , solCl = cl
          , solCd = cd
          , solCenterPressure = (xCp, yCp)
          }
  where
    -- surface velocities
    (vs,vorticities) = ( (mapVector (\q -> cos(q - alpha)) angles) + (mV <> qsGamma)
                       , map (qsGamma @>) [n-(length elements),n-1]
                       )
      where
        (mA, mV, b) = getAVb geometries angles kuttaIndices alpha
        qsGamma = flatten $ linearSolve mA b
        n = dim angles
    -- pressure coefficients
    cps = 1 - vs*vs
    
    -- forces and force coefficients
    xForces = -cps*(join $ map (fst . fNormals) elements)
    yForces = -cps*(join $ map (snd . fNormals) elements)
    xf = sumElements xForces
    yf = sumElements yForces
        
    cd =  xf*(cos alpha) + yf*(sin alpha)
    cl = -xf*(sin alpha) + yf*(cos alpha)

    -- centers of pressure
    (xCp, yCp) = ( (sumElements (xs*yForces)) / (sumElements yForces)
                 , (sumElements (ys*xForces)) / (sumElements xForces)
                 )
      where
        (xs,ys) = (\(x,y) -> (join x, join y)) $ unzip $ map fMidpoints elements

    kuttaIndices = ki 0 (map dim angles')
      where
        ki _ [] = []
        ki n0 (n:ns) = (n0,n0+n-1):(ki (n0+n) ns)

    angles = join angles'
    angles' = map fAngles elements

    -- (mSL,mCL,mSB,mCB)
    geometries = setupGeometries angles (xms, yms) (xns, yns) (xnps, ynps)
      where
        joinTuples (x,y) = (join x, join y)
        (xms,yms)   = joinTuples $ unzip $ map fMidpoints elements
        (xns,yns)   = joinTuples $ unzip $ map fInits     elements
        (xnps,ynps) = joinTuples $ unzip $ map fTails     elements


-- make influence matrix A and also the matrix V where V*[sources; vortex] == tangential speeds
--getAVb :: (Num (Vector a), RealFloat a, Container Vector a) =>
--          Foil a -> a -> (Matrix a, Matrix a)
getAVb :: (Floating a, Num (Vector a), Container Vector a, Storable a) =>
          (Matrix a, Matrix a, Matrix a, Matrix a)
          -> Vector a
          -> [(Int, Int)]
          -> a
          -> (Matrix a, Matrix a, Matrix a)
getAVb (mSL,mCL,mSB,mCB) angles kuttaIndices alpha =
  ( scale (1/(2*pi)) $ fromBlocks [[          mAij, fromColumns vsAin]
                                  ,[fromRows vsAnj, fromLists ann]]
  , scale (1/(2*pi)) $ fromBlocks [[ mSB - mCL
                                   , fromColumns vsTin]]
  , asColumn $ join $ (mapVector (\q -> sin $ q - alpha) angles):
                      (map (\(n0,nf) -> fromList [-cos((angles @> n0) - alpha) - cos((angles @> nf) - alpha)]) kuttaIndices)
  ) -- (for middle entry mV, mSL + mCB == mAij)
  where
    n = rows mSL
    
    -- sources influence matrix Aij
    mAij = mSL + mCB
    
    -- vortices influence vectors Ains
    vsAin = map (\(n0,nf) -> fromList $ map sumElements $ toRows $ subMatrix (0,n0) (n, nf-n0+1) m) kuttaIndices
      where
        m = (mCL - mSB)
    
    -- vortices tangential velocity outputs
    vsTin = map (\(n0,nf) -> fromList $ map sumElements $ toRows $ subMatrix (0,n0) (n, nf-n0+1) m) kuttaIndices
      where
        m = (mSL + mCB)

    -- sources kutta condition influence vector Anj
    vsAnj = map (\(n0,nf) -> (getRow n0 mSB) - (getRow n0 mCL) + (getRow nf mSB) - (getRow nf mCL)) kuttaIndices

    -- vortices kutta condition influence scalars ann
    ann = (flip map) kuttaIndices $ \(ni0,nif) -> (flip map) kuttaIndices $ \(nj0,njf) ->
      sumElements $ (getSubRow ni0 (nj0,njf) mSL) + (getSubRow ni0 (nj0,njf) mCB) + (getSubRow nif (nj0,njf) mSL) + (getSubRow nif (nj0,njf) mCB)

    getRow i mat = flatten $ subMatrix (i,0) (1,n) mat
    getSubRow i (j0,jf) mat = flatten $ subMatrix (i,j0) (1,jf-j0+1) mat


setupGeometries :: (RealFloat t, Storable t) =>
                   Vector t
                   -> (Vector t, Vector t)
                   -> (Vector t, Vector t)
                   -> (Vector t, Vector t)
                   -> (Matrix t, Matrix t, Matrix t, Matrix t)
setupGeometries angles (xms, yms) (xns, yns) (xnps, ynps) = runST $ do
  let n = dim angles
  mSL' <- newMatrix 0 n n
  mCL' <- newMatrix 0 n n
  mSB' <- newMatrix 0 n n
  mCB' <- newMatrix 0 n n
  
  forM_ [0..n-1] $ \i -> forM_ [0..n-1] $ \j -> do
    let qi = angles @> i
        qj = angles @> j
        
        xmi  = xms @> i
        ymi  = yms @> i
        xnj  = xns @> j
        ynj  = yns @> j
        xnjp = xnps @> j
        ynjp = ynps @> j
        
        s = sin (qi - qj)
        c = cos (qi - qj)
        
        l 
          | i == j    = 0
          | otherwise = log(r1/r0)
          where
            r1 = distance (xmi, ymi) (xnjp, ynjp)
            r0 = distance (xmi, ymi) (xnj,  ynj )
            distance (x1,y1) (x0,y0) = sqrt $ dx*dx + dy*dy
              where
                dx = x1 - x0
                dy = y1 - y0

        b 
          | i == j = pi
          | otherwise = atan2 (dyjp*dxj - dxjp*dyj)
                              (dxjp*dxj + dyjp*dyj)
          where
            dyj  = ymi - ynj
            dxj  = xmi - xnj
            dyjp = ymi - ynjp
            dxjp = xmi - xnjp
        
    _ <- writeMatrix mSL' i j (s*l)
    _ <- writeMatrix mCL' i j (c*l)
    _ <- writeMatrix mSB' i j (s*b)
    writeMatrix      mCB' i j (c*b)
  sl' <- freezeMatrix mSL'
  cl' <- freezeMatrix mCL'
  sb' <- freezeMatrix mSB'
  cb' <- freezeMatrix mCB'
  return (sl',cl',sb',cb')
