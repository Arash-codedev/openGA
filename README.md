# openGA
A free C++ Genetic Algorithm library

<br>

### Download

Download link:

https://github.com/Arash-codedev/openGA/archive/v1.0.5.zip

### Documentation

Documentation is available at:

https://github.com/Arash-codedev/openGA/blob/master/openGA.pdf

If you have found this library useful for your work, please cite the following paper:

[openGA, a C++ Genetic Algorithm library](https://www.researchgate.net/publication/320944800_openGA_a_C_Genetic_Algorithm_library)


```
@inproceedings{mohammadi2017openga,
  title={OpenGA, a C++ Genetic Algorithm Library},
  author={Mohammadi, Arash and Asadi, Houshyar and Mohamed, Shady and Nelson, Kyle and Nahavandi, Saeid},
  booktitle={Systems, Man, and Cybernetics (SMC), 2017 IEEE International Conference on},
  pages={2051--2056},
  year={2017},
  organization={IEEE}
}
```

### About 

This library is easy to use. Examples are provided in the document. The library is placed in a single file for comfort of the user. The standard C++ libraries are sufficent for this library and you do not need to install additional library for compiling.

### Demo

Youtube video:

[![Youtube video](https://img.youtube.com/vi/8T2Teo_Lwrc/0.jpg)](https://www.youtube.com/watch?v=8T2Teo_Lwrc)

### Supported platforms
- Linux
- Windows
- OSX

### Automatic code generator 

For the comfort of the client, openGA assist provides an interface to generate an automatic code. OpenGA assist is located at *assist/index.html*.

![openGA assist](https://user-images.githubusercontent.com/11730626/47605121-f276a980-da4d-11e8-9716-42ee1c27faee.png)

### FAQ 

**Why does the crossover function generate a single offspring instead of two?**
This library calls crossover function two times instead of calling it once to create two offspring. This is to reduce the burden from the client programmer.

**What is the shrink scale?**
The mutation jumps should reduce over generations to increase the accuracy of solution. Shrink scale is a generation dependent coefficient to scale the resolution.

**After optimization, how to access the population?**
To have an access to the last generation, you can use an object field called `last_generation`. A particular usage example can be
```
best_solution = ga_obj.last_generation.chromosomes[ga_obj.last_generation.best_chromosome_index];
```
In the current implementation of openGA, only the middle cost of the last generation is stored. The reason is that keeping the entire generations can lead to a huge memory usage and eventually crashing the application after an out of memory error for particular applications which keep huge amount of data.
If you are keen to keep the full information trace of all generations, you can perform it manually via the report function.

**Does multi-threading always improve the performance?**
Multi-threading itself has an overhead. Therefore, it is helpful when the evaluation process is heavy enough. If the evaluation process is so fast, parallel programming is not necessary.

**Instead of a function, can I set a class method as a GA operator or a cost function?**
Yes, C++ supports it. It can be achieved by using C++ standard library. For this purpose, you need to use `std::bind` and `std::placeholders`. For detailed implementation information search for C++ function template bind.
