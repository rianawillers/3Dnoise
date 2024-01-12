C code compiled in Visual Studio 2022

Use last flag: uint for unsigned long data / float for float data

Commandline execution from folder 3Dnoise:
x64\Debug\3Dnoise.exe ..\..\data\noise.ulong 21 120 7 150 uint > ..\out\noise_ulong.out
x64\Debug\3Dnoise.exe ..\..\data\noise.float 65 145 7 150 float > ..\out\noise_float.out

From within Visual Studio 2022:
 * Property Pages
   * Configuration Properties
     * Debugging
       * Command Arguments: ../../data/noise.float 65 145 7 150 float > ../out/noise_float.out
	   
NB: For screen output, do not pipe to out file

----------------------------------------------------------------------------------------------------

Python notebook requires the pyradi package