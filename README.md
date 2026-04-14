ROMS-KVIS Model
We modify the ROMS model to enable it to read the 3D visible solar radiation attenuation coefficient profile (KVIS) to replace the default Jerlov water type. The main modification is in lmd_swfrac.F, and the calculation of infrared radiation refers to [Lee et al., 2005].

If you want to use this ROMS-KVIS model, please follow these steps:
1. Add the varinfo.dat entries to your program.
2. Define #define DEKVIS in your head file.
3. Make sure lmd_swfrac is used.
4. Make sure there is a variable named KVIS in your 3D initial files.
5. Make sure you have swrad input (and ensure it is not zero throughout the simulation).

A diagnostic test is included in lmd_swfrac.F:
IF (i == 3 .AND. j == 3 .AND. k == 10) THEN
    PRINT *, 'Z(i,j): ', Z(i,j)
    PRINT *, 'KVIS_3d', KVIS_3d(i,j)
    PRINT *, 'lat', GRID(ng)%latr(i,j)
    PRINT *, 'swdk', swdk(i,j)
END IF


If this point is wet and the model has at least 10 layers, you can see the KVIS value at this location printed in the terminal (standard output) while running the model. Please note that ROMS indexing starts from 0 while MATLAB or EXCEL starts from 1, so the horizontal location indices may differ by one.
