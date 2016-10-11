A short description of the range data
-----------------------------------------------------------------

*.mat file: is the image of Z-coordinates(floating-point image). 
*.ppm file: the object's picture

the Z value is the measured range of the object(the distance from a x-y 
rectangular table). Values are 0 and negative floating point values.
	* 0: means background
	* negative: the smaller, the further away from the viewer.

You need to choose your gaussian parameters based on z_depth characteristics.

