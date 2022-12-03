# Lior and Nave's 3D Animations Assignment 1 (Decimation)

In this git repository there are two branches:
* Master branch = decimation with quadratic error calculation.
* Naive branch = decimation with tutorial 703 shortest edge and midpoint.

### MASTER BRANCH
***********************************
The assignment is coded in the folder [Assignment1](https://github.com/ThatGuyVanquish/3DAnimation_AS1/tree/master/tutorial/Assignment1)
#    

Our implementation added an object:
[MeshSimplification.cpp](https://github.com/ThatGuyVanquish/3DAnimation_AS1/blob/e17c7deb3582e7adf966487180c78e2ed671f91d/tutorial/Assignment1/MeshSimplification.cpp)

Mesh Simplification object has fields to store the mesh it is assigned to (so we could have multiple MeshSimplification objects, each responsible for a single Mesh object) and to calculate the amount of decimation levels.
### 
The object constructor calls a function [*createDecimatedMesh*](https://github.com/ThatGuyVanquish/3DAnimation_AS1/blob/e17c7deb3582e7adf966487180c78e2ed671f91d/tutorial/Assignment1/MeshSimplification.cpp#L152) which when initialized immediately calculates the mesh's necessary MeshData structs for each level of decimation - which is 10% of the previous level of the decimated mesh.

We calculate the cost and placement based on the [*calculateCostAndPos*](https://github.com/ThatGuyVanquish/3DAnimation_AS1/blob/e17c7deb3582e7adf966487180c78e2ed671f91d/tutorial/Assignment1/MeshSimplification.cpp#L170) lambda within the createDecimatedMesh method, using the equations given in the article and class PDF.
# 
We also changed [BasicScene.cpp](https://github.com/ThatGuyVanquish/3DAnimation_AS1/blob/master/tutorial/Assignment1/BasicScene.cpp) to encorporate our Mesh Simplification object.
###  

Object choice:
1) **objFiles** vector of string filenames to obj/off mesh files
2) **objIndex** int index within objFiles
3) **decimations** int amount of decimated mesh levels to compute
4) **recalcQsRate** int number of iterations before recalculating the Q matrices
