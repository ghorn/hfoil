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
                         , solStrengths :: Vector a
                         , solAlpha :: a
                         , solForces :: (Vector a, Vector a)
                         , solCl :: a
                         , solCd :: a
                         , solCenterPressure :: (a,a)
                         }

solveFlow :: (Num (Vector a), RealFloat a, Field a) => Foil a -> a -> FlowSol a
solveFlow foil alpha = FlowSol { solFoil = foil
                               , solVs = vs
                               , solCps = cps
                               , solStrengths = qsGamma
                               , solAlpha = alpha
                               , solForces = (xForces, yForces)
                               , solCl = cl
                               , solCd = cd
                               , solCenterPressure = (xCp, yCp)
                               }
  where
--    kuttaIndices = pKuttaIndices foil
    
    -- surface velocities
    (vs,qsGamma) = ( (mapVector (\q -> cos(q - alpha)) (fAngles foil)) + (mV <> qsGamma')
                   , qsGamma'
                   )
      where
        (mA, mV) = getAV foil
        b = getB foil alpha
        qsGamma' = flatten $ linearSolve mA b
    
    -- pressure coefficients
    cps = 1 - vs*vs
    
    -- forces and force coefficients
    xForces = -cps*(fst $ fNormals foil)
    yForces = -cps*(snd $ fNormals foil)
    xf = sumElements xForces
    yf = sumElements yForces
        
    cd =  xf*(cos alpha) + yf*(sin alpha)
    cl = -xf*(sin alpha) + yf*(cos alpha)

    -- centers of pressure
    (xCp, yCp) = ( (sumElements (xs*yForces)) / (sumElements yForces)
                 , (sumElements (ys*xForces)) / (sumElements xForces)
                 )
      where
        (xs,ys) = fMidpoints foil


getB :: (Floating a, Storable a) => Foil a -> a -> Matrix a
getB panels alpha = asColumn $ join [ mapVector (\q -> sin $ q - alpha) angles
                                    , fromList [-cos((angles @> 0) - alpha) - cos((angles @> (n-1)) - alpha)]
                                    ]
  where
    angles = fAngles panels
    n = (dim (fst (fNodes panels)))-1    

-- make influence matrix A and also the matrix V where V*[sources; vortex] == tangential speeds
getAV :: (Num (Vector a), RealFloat a, Container Vector a) =>
         Foil a -> (Matrix a, Matrix a)
getAV panels = ( scale (1/(2*pi)) $ fromBlocks [[       mAij,   asColumn vAin]
                                               ,[ asRow vAnj, fromLists [[ann]]]]
               , scale (1/(2*pi)) $ fromBlocks [ [mSB - mCL, asColumn $ fromList $ map sumElements (toRows $ mAij)]]
               ) -- (for last entry, mSL + mCB == mAij)
  where
    -- sources influence matrix Aij
    mAij = mSL + mCB
    
    -- vortices influence vector Ain
    vAin = fromList $ map sumElements $ toRows (mCL - mSB)
        
    -- sources kutta condition influence vector Anj
    vAnj = (getRow 0 mSB) - (getRow 0 mCL) + (getRow (n-1) mSB) - (getRow (n-1) mCL)

    -- vortices kutta condition influence scalar ann
    ann = sumElements $ (getRow 0 mSL) + (getRow 0 mCB) + (getRow (n-1) mSL) + (getRow (n-1) mCB)

    n = (dim (fst (fNodes panels)))-1
    
    getRow i mat = flatten $ subMatrix (i,0) (1,n) mat
--  getCol j mat = flatten $ subMatrix (0,j) (n,1) mat
      
    (mSL,mCL,mSB,mCB) = runST $ do
      mSL' <- newMatrix 0 n n
      mCL' <- newMatrix 0 n n
      mSB' <- newMatrix 0 n n
      mCB' <- newMatrix 0 n n
      
      forM_ [0..n-1] $ \i -> forM_ [0..n-1] $ \j -> do
        let qi = fAngles panels @> i
            qj = fAngles panels @> j
            s = sin (qi - qj)
            c = cos (qi - qj)
            l = lnrr panels (i,j)
            b = beta panels (i,j)
        _ <- writeMatrix mSL' i j (s*l)
        _ <- writeMatrix mCL' i j (c*l)
        _ <- writeMatrix mSB' i j (s*b)
        writeMatrix      mCB' i j (c*b)
      sl' <- freezeMatrix mSL'
      cl' <- freezeMatrix mCL'
      sb' <- freezeMatrix mSB'
      cb' <- freezeMatrix mCB'
      return (sl',cl',sb',cb')

lnrr :: (Floating a, Storable a) =>
        Foil a -> (Int, Int) -> a
lnrr panels (i,j)
  | i == j = 0
  | otherwise = log(r1/r0)
  where
    (xms, yms) = fMidpoints panels
    (xns, yns) = fNodes panels
    r1 = distance (xms @> i, yms @> i) ( xns @> (j+1), yns @> (j+1))
    r0 = distance (xms @> i, yms @> i) ( xns @>  j   , yns @>  j   )
    distance (x1,y1) (x0,y0) = sqrt $ dx*dx + dy*dy
      where
        dx = x1 - x0
        dy = y1 - y0

beta :: (RealFloat a, Storable a) =>
        Foil a -> (Int, Int) -> a
beta panels (i,j)
  | i == j    = pi
  | otherwise = atan2 (dyjp*dxj - dxjp*dyj)
                      (dxjp*dxj + dyjp*dyj)
  where
    dyj  = (yms @> i) - (yns @>  j   )
    dxj  = (xms @> i) - (xns @>  j   )
    dyjp = (yms @> i) - (yns @> (j+1))
    dxjp = (xms @> i) - (xns @> (j+1))
    (xms, yms) = fMidpoints panels
    (xns, yns) = fNodes panels
