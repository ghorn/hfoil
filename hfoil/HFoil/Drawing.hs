{-# OPTIONS_GHC -Wall #-}

module HFoil.Drawing( toPic
                    , draw
                    ) where

import Graphics.Gloss

xSize, ySize :: Int
xSize = 800
ySize = 500

toPic :: Real a => Color -> [(a,a)] -> Picture
toPic col coords = scale (0.8*(fromIntegral xSize)) (0.8*(fromIntegral xSize))
                   $ translate (-0.5) 0
                   $ color col
                   $ line $ map (\(x,y) -> (realToFrac x, realToFrac y)) coords

draw :: [Picture] -> IO ()
draw pics = do
  display 
    (InWindow
     "hfoil"             -- window title
     (xSize, ySize)      -- window size
     (10, 650))          -- window position
    black                -- background color
    (pictures pics)      -- picture to display

--  let line = plot_lines_values ^= [[ (xc, yt (naca4 "0012") xc)
--                                   | xc <- [0,0.01..0.99::Double]]]
--             $ plot_lines_title ^= "naca 0012"
--             $ defaultPlotLines
--  
--      chart = layout1_title ^= "naca yo"
--              $ layout1_plots ^= [Left (toPlot line)]
--              $ defaultLayout1
--  
--  renderableToWindow (toRenderable chart) 640 480
--  _ <- renderableToPNGFile (toRenderable chart) 640 480 "mDiv_vs_tc.png"
--  return ()
