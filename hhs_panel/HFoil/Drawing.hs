{-# OPTIONS_GHC -Wall #-}

module HFoil.Drawing( drawNaca
                    ) where

import Graphics.Gloss

import HFoil.Naca4
import HFoil.Panels

drawNaca :: IO ()
drawNaca = do
  let picture = scale (0.8*(fromIntegral xSize)) (0.8*(fromIntegral xSize))
                $ translate (-0.5) 0
                $ color white
                $ line naca0012Coords
      
      naca0012Coords = map (\(x,y) -> (realToFrac x, realToFrac y)) $ toNodes (naca4 "0012" :: Naca4 Double) 200
      xSize = 400
      ySize = 150
  display 
    (InWindow
     "Hello World"       -- window title
     (xSize, ySize)          -- window size
     (10, 710))          -- window position
    black                        -- background color
    picture               -- picture to display
--    (pictures (map picture [map fst naca0012Coords, map snd naca0012Coords])) -- picture to display

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
