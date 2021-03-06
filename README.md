# Raytracing

[FINAL REPORT](./rapport.pdf)
___
Final animation (render time : 5h30 for 120 images, 512x512, 128 rays/pixel) :

![](anim.gif)

Meshes (need to be UV mapped for textures) :
- https://free3d.com/3d-model/small-satellite-308237.html
- https://free3d.com/3d-model/saturn-v1--741827.html
- https://free3d.com/3d-model/rocket-ship-v1--579030.html

## Getting started
Clone this repo and

```
cd raytracing
make
```

This will generate a `build` directory and the `main` executable.

Run it :
```
src/main
```

A file `image.png` is generated. 

Feel free to change parameters in the `src/main.cpp` file (don't forget to `make` again).

To clean everything : 
```
make clean
```

## Troubleshooting
- If the compiler `g++-11` isn't suppported on your machine, change it in the makefile. It is the only one allowing me to run natively openMP on a MacBook M1. 

Else, open an issue !