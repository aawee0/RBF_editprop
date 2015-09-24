# RBF_editprop

This piece of python code performs edit propagation on an image. Inputs are:
1) img.bmp - the original image
2) chg.bmp - same image with user-specified strokes (e.g. of different color) and black strokes as constraints. 
Pencil or any drawing tool with sharp border is preferential to make the strokes.

Both *.bmp files should be in the same directory as the main.py.
Names of input files are hardcoded and need no mention in the command line when launching.

Outputs are res1.png (plain edit propagation) and res2.png (texture-aware edit prop using Gabor filters).
