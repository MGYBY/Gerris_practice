Some scripts and tricks useful for p.p.
1. Gmsh file for generating .stl files for solids.
2. Adding time labels in GfsView:
```
Label {
  r = 0 g = 0 b = 0
  shading = Constant
  maxlevel = 7
  font_size = 25
  raster_font = 0
  line_width = 1
} {
  x = 1.65 y = 0.825 z = 0
  label = "t = %2.2f s"
  symbol = 0
}
```
The maximum font size in GfsView is 10. To increase the max font size, modify the following code files:
**view/interface.c**, 
**view/mangled_interface.c**, 
**view/gfsview.glade**
