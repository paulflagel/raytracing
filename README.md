# Raytracing

Clone this repo and

```
cd raytracing
make
```

This will generate a `build` directory and the `raytracer` executable.

Run it :
```
./raytracer
```

A file `image.png` is generated. 

Feel free to change parameters in the `src/main.cpp` file (don't forget to `make` again).

To clean everything : 
```
make clean
```

### Troubleshooting
- If the compiler `g++-11` isn't suppported on your machine, change it in the makefile.
- Depending on your machine's specs, the number of threads used for parallelization may vary (8 here). I had to hardcode it in `src/RandomHelper/RandomHelper.c++` to create a thread-safe helper function. You can also change it.

Else, open an issue.