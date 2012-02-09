{-# OPTIONS_GHC -Wall #-}

module Numeric.HFoil.Naca4( Naca4(..)
                          , coords
                          , yt
                          , dyt
                          , naca4
                          ) where

naca4 :: (Read a, Fractional a) => String -> Naca4 a
naca4 name@(m_:p_:t0:t1:[]) = Naca4 m p t ("NACA "++name)
  where
    m = 0.01 * read [m_]
    p = 0.1  * read [p_]
    t = 0.01 * read (t0:[t1])
naca4 _ = error "not a 4 digit airfoil"

-- m: max camber in hundredths of chord
-- p: position of max camber in tenths of chord
-- t: max thickness in hundredths of chord
data Naca4 a = Naca4 { naca4_m :: a
                     , naca4_p :: a
                     , naca4_t :: a
                     , naca4_name :: String
                     } deriving Show

--  xc: x/chord
yc :: (Ord a, Fractional a) => Naca4 a -> a -> a
yc (Naca4 {naca4_p = p, naca4_m = m}) xc
  | xc < 0    = error "xc < 0"
  | xc <= p   = m/(p*p)*(2*p*xc - xc*xc)
  | xc <= 1   = m/((1-p)*(1-p))*((1-2*p) + 2*p*xc - xc*xc)
  | otherwise = error "xc > 1"

dyc :: (Ord a, Fractional a) => Naca4 a -> a -> a
dyc (Naca4 {naca4_p = p, naca4_m = m}) xc
  | xc < 0    = error "xc < 0"
  | xc <= p   = m/(p*p)*(2*p - 2*xc)
  | xc <= 1   = m/((1-p)*(1-p))*(2*p - 2*xc)
  | otherwise = error "xc > 1"

yt :: (Ord a, Floating a) => Naca4 a -> a -> a
yt (Naca4 {naca4_t = t}) xc
  | xc < 0 = error "xc < 0"
  | xc > 1 = error $ "xc > 1"
  | otherwise = 5*t*(0.2969*sqrt(xc) - 0.1260*(xc) - 0.3537*(xc)**2 + 0.2843*(xc)**3 - 0.1015*(xc)**4)

dyt :: (Ord a, Floating a) => Naca4 a -> a -> a
dyt (Naca4 {naca4_t = t}) xc
  | xc < 0 = error "xc < 0"
  | xc > 1 = error "xc > 1"
  | otherwise = 5*t*(0.5*0.2969/sqrt(xc) - 0.1260 - 2*0.3537*(xc) + 3*0.2843*(xc)**2 - 4*0.1015*(xc)**3)

coords :: (Ord a, Floating a) => Naca4 a -> a -> ((a,a),(a,a))
coords foil xc 
  | naca4_m foil == 0 = ((xc,yt_), (xc,-yt_))
  | otherwise         = ((xu,yu ), (xl, yl ))
  where
    yt_ = yt foil xc
    yc_ = yc foil xc
    
    yu = yc_ + yt_ * (cos theta)
    yl = yc_ - yt_ * (cos theta)

    xu = xc  - yt_ * (sin theta)
    xl = xc  + yt_ * (sin theta)

    theta = atan $ dyc foil xc
