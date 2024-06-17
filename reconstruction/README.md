# Implicit Shape Reconstruction
Name: Svitlana Morkva

## Reports
### 1 - Visualization of the 'cat.off' point cloud with constrained points
| model name  | view 01                                                 | view 02                                                  |
| :---------: |---------------------------------------------------------|----------------------------------------------------------|
| cat         | <img align="center" src="./res/Cat1_1.png" width="300"> | <img align="center"  src="./res/Cat1_2.png" width="300"> |

### 2 - Grid with nodes colored w.r.t. the implicit function (using the non-axis-aligned grid as described in Section 2.3) 

Here I also added a parameter to have an ability to change from non- to axis-aligned grid and back.

| model name  | view 01                                                   | view 02                                                    |
| :---------: |-----------------------------------------------------------|------------------------------------------------------------|
| cat         | <img align="center" src="./res/Cat2_1.png" width="300">   | <img align="center"  src="./res/Cat2_2.png" width="300">   |
| luigi      | <img align="center" src="./res/Luigi2_1.png" width="300"> | <img align="center"  src="./res/Luigi2_2.png" width="300"> |

### 3 - Reconstructed surfaces
**Please also save the reconstructed shape (.off) in the results folder**

|                          sphere                          |                           cat                           |
|:--------------------------------------------------------:|:-------------------------------------------------------:| 
| <img align="center" src="./res/Sphere3.png" width="300"> | <img align="center"  src="./res/Cat3.png" width="300">  |
|                          luigi                           |                          horse                          |
| <img align="center" src="./res/Luigi3.png" width="300">  | <img align="center"  src="./res/Horse.png" width="300"> |
|                          hound                           |                                                         |
| <img align="center" src="./res/Hound3.png" width="300">  |                                                         |


**Please summarize your observations of what happens when you change the following parameters. Please feel free to add screenshots (by creating a separate table) to illustrate your explanation.**

| params                   | Your Observation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | 
| :---------------------:  |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| grid resolution          | When we increase the grid resolution - we increase the number of cells in each dimension, and so the number of points with the value of the implicit function (if the point is inside or outside of the figure). With increasing this value our reconstructed mesh will be smoother and less cubic, because we have more points on the surface. If the grid size is too small we also can lose some details, because the points of the mesh are too far away from the points on the grid. But the time consumption also will grow with the increasing of the grid size. |
| Wendland function radius | This parameter is responsible for determining how close should point of the mesh be to influence the value of the grid (and how much they will influence it). If this value is too big, points inside and outside of the surface will just cancel each other, and for many grid points we won't receive valid results, the reconstructed object will lose its details and will be much smoother. If the value is too small - the figure will have more edges and there can appear some holes in the places that don't have enough points.                               |
| polynomial degree        | The degree og the polynomial change the type of function that is used for the MLS. With increasing of the degree the figure by itself is smoother, but we get many flyaways.                                                                                                                                                                                                                                                                                                                                                                                            |

| params                   | Small                                                    | Big                                                      |
| :---------------------:  |----------------------------------------------------------|----------------------------------------------------------|
| grid resolution          | <img align="center" src="./res/GridS.png" width="300">   | <img align="center" src="./res/GridB.png" width="300">   |
| Wendland function radius | <img align="center" src="./res/RadS.png" width="300">    | <img align="center" src="./res/RadB.png" width="300">    |
| polynomial degree        | <img align="center" src="./res/DegreeS.png" width="300"> | <img align="center" src="./res/DegreeB.png" width="300"> |


**Please compare the computation time for constraints computation and MLS using brute force and a spatial index. Use hound.off for the comparison and use the same parameters for the timings with and without use of spatial index (please report the parameters that you used).**

Parameters: grid resolution = 50; wendland function radius = 0.05; polynomial degree = 0.

| step                    | brute force | spatial index |
| :---------------------: |:-----------:|:-------------:|
| constraints             |   3495 ms   |     32 ms     |
| MLS                     |  12955 ms   |    2303 ms    |



### 4 - Theory Question 1

**Prove that the normal of a surface defined by an implicit function is proportional to its gradient.**

[Link to the pdf with the proof](./results/Prove1.pdf)

### 5 - Theory Question 2

**Compute the closed-form gradient of the MLS approximation.**

I, unfortunately, did not finish this task. 

### 6 - Normal-based v.s. point-based reconstruction ("hound.off")

As you can see, the results received using the normal-based approach are much more realistic even when we use small grid resolution. Especially good this algorithm work with edges and small parts of the figure, this approach does not cut them off as the point-based one does.
Another advantage of this algorithm is that we are using much fewer points and there is no need to preprocess and add additional constraints.

| method       |                        view 01                         |                        view 02                         |
| :---------:  |:------------------------------------------------------:|:------------------------------------------------------:| 
| point-based  |  <img align="center" src="./res/PB1.png" width="300">  | <img align="center"  src="./res/PB2.png" width="300">  | 
| normal-based |  <img align="center" src="./res/NB1.png" width="300">  | <img align="center"  src="./res/NB2.png" width="300">  |


### 7 - MLS v.s. Screened Poisson Reconstruction v.s. RIMLS

**No implementation required, you can use [Meshlab](https://www.meshlab.net) to perform the comparisons.**

I tried to change the parameters of all methods, so their grid sizes will be as close as possible for more adequate analysis.

My realization of MLS has similar results to RIMLS, but their reconstructions are much smoother, without weird curves in some places. 
At the same moment, they are more precise with smaller detail, what you can see around the hat, cheeks, legs, and arms of Luigi. This algorithm has these advantages because it deals with outliers much better than a regular MLS.
On contrary, the Screened Poisson Reconstruction seems to have a completely different approach, the sizes of its grid differentiate, they are bigger on flat regions and smaller in curvier parts, which, I think, is a good advantage from the efficiency perspective as we don't need so many cells to represent shape.

| model names  |                            MLS                             |                             Possion                             |                             RIMLS                             | 
| :---------:  |:----------------------------------------------------------:|:---------------------------------------------------------------:|:-------------------------------------------------------------:|
| cat          |  <img align="center" src="./res/Cat_MLS.png" width="300">  |  <img align="center"  src="./res/Cat_Poisson.png" width="300">  |  <img align="center"  src="./res/Cat_RIMLS.png" width="300">  |
| luigi        | <img align="center" src="./res/Luigi_MLS.png" width="300"> | <img align="center"  src="./res/Luigi_Poisson.png" width="300"> | <img align="center"  src="./res/Luigi_RIMLS.png" width="300"> |

### 8 - PCA normals v.s. provided normals (luigi.off)

Parameters added: # of k-nearest neighbours, if automatic normals flipping is on, if first normal needed to be flipped (for automatic mode).

In this task I also implemented automatic normals flipping. I added the second table to show how it works in compare to the given normals. 

|     model names     |                                                                                  PCA normals                                                                                                                       |                      Provided normals                       | 
|:-------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------:|
| luigi (constraints) |                                                                                                                          <img align="center"  src="./res/Luigi8c_2.png" width="300">                                                                                                                           | <img align="center"  src="./res/Luigi8c_1.png" width="300"> |
| luigi (reconstruct) |                                                                                                                           <img align="center"  src="./res/Luigi8_2.png" width="300">                                                                                                                           | <img align="center"  src="./res/Luigi8_1.png" width="300">  |

In my opinion, I got pretty decent results with this method, but of course there are some inaccuracies, especially you can see it near the hat of the figure.<br/> 
These inaccuracies usually happen because of noise, sharp edges or if we don't have enough points in the region. 
It can also be hard to estimate a suitable # for k-nearest neighbours and even harder to estimate consistent normal directions without having the reference, as you can see in the next table.


| model names |                         PCA normals                         |                      Provided normals                       | 
|:-----------:|:-----------------------------------------------------------:|:-----------------------------------------------------------:|
|   sphere    | <img align="center"  src="./res/Sphere8_2.png" width="300"> | <img align="center"  src="./res/Sphere8_1.png" width="300"> |
|     cat     |  <img align="center"  src="./res/Cat8_2.png" width="300">   |  <img align="center"  src="./res/Cat8_1.png" width="300">   |
|    luigi    | <img align="center"  src="./res/Luigi8a_2.png" width="300"> | <img align="center"  src="./res/Luigi8a_1.png" width="300"> |
