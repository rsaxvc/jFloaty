# jFloaty

JPEG -> RGBFFF decoder with reduced rounding

## Demo

This is a picture of my friend Russ on a JetSki, but horribly underexposed.

![](demo/input.jpg "Russ on a jet-ski in Table Rock Lake")

Lets fix it with 'auto input level' feature in GIMP.
But when we do so, the JPEG decoder has rounding error.
jFloaty is a patch to STB-Image that decodes JPEGs direct to float-32.

|          |      GIMP Alone | jFloaty+GIMP  | Notes |
|:--------:|:----------------------------------------:|:-------------------------------------------------------:|:-:|
|Image     | ![](demo/input.autolvl.png "RGB888-Img") | ![](demo/output.autolvl.png "RGBFFF-Img")               | There are more colors available. The banding is improved. |
|Histogram | ![](demo/input.hist.jpg "RGB888-Hist")   | ![](demo/output.hist.jpg "RGBFFF-Hist")                 | Evidence of more colors |
|Zoomed    | ![](demo/input.zoom.png   "RGB888-Zoom") | ![](demo/output.zoom.png "RGBFFF-Zoom")                 | Zoomed detail of face |
