# Shape Deformation

Name: Svitlana Morkva

## Report
### 1.1 - 1.4 Multiresolution Mesh Editing
|              model              | S                                                           | B                                                            | B'                                                           | S'                                                            |
|:-------------------------------:|-------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|---------------------------------------------------------------|
|            woody-lo             | <img align="center" src="./res/Woody_lo_S.png" width="300"> | <img align="center"  src="./res/Woody_lo_B.png" width="300"> | <img align="center" src="./res/Woody_lo_B2.png" width="300"> | <img align="center"  src="./res/Woody_lo_S2.png" width="300"> |
|            woody-hi             | <img align="center" src="./res/Woody_hi_S.png" width="300"> | <img align="center"  src="./res/Woody_hi_B.png" width="300"> | <img align="center" src="./res/Woody_hi_B2.png" width="300"> | <img align="center"  src="./res/Woody_hi_S2.png" width="300"> |
|   hand (translation+rotation)   | <img align="center" src="./res/Hand_S.png" width="300">     | <img align="center"  src="./res/Hand_B.png" width="300">     | <img align="center" src="./res/Hand_B2.png" width="300">     | <img align="center"  src="./res/Hand_S2.png" width="300">     |
| cylinder (translation+rotation) | <img align="center" src="./res/Cylinder_S.png" width="300"> | <img align="center"  src="./res/Cylinder_B.png" width="300"> | <img align="center" src="./res/Cylinder_B2.png" width="300"> | <img align="center"  src="./res/Cylinder_S2.png" width="300"> |

### 1.5 Real time mesh editing

| model | S' - real time                                          | model | S' - real time                                          |
| :-----------:  |---------------------------------------------------------|-----  |---------------------------------------------------------|
| bar            | <img align="center"  src="./res/Bar.gif" width="300">   | camel_head     | <img align="center"  src="./res/Camel.gif" width="300"> |
| bumpy_plane    | <img align="center"  src="./res/Plane.gif" width="400"> | cactus         | <img align="center"  src="./res/Cactus.gif" width="300"> |


### 2 Deformation transfer
| model | High-freq detail transfer                                  | Deformation transfer                                        |
| :-----------:  |------------------------------------------------------------|-------------------------------------------------------------|
| bar            | <img align="center" src="./res/HF_Bar.png" width="300">    | <img align="center"  src="./res/DT_Bar.png" width="300">    |
| bumpy_plane    | <img align="center" src="./res/HF_Plane.png" width="300">  | <img align="center"  src="./res/DT_Plane.png" width="300">  |
| camel_head     | <img align="center" src="./res/HF_Camel.png" width="300">  | <img align="center"  src="./res/DT_Camel.png" width="300">  |
| cactus         | <img align="center" src="./res/HF_Cactus.png" width="300"> | <img align="center"  src="./res/DT_Cactus.png" width="300"> |


#### Discussion of differences

As you can see from the results above, almost in every aspect Deformation transfer gives better results than High-freq detail transfer.
This algorithm produces much less distortion that can obviously be observed in bumps of "bumpy_plane" (especially on the edges), on the stem of "cactus" or on the "bar".
Deformation transfer doesn't lead to self-intersections as easily (tries it best to avoid them), but they can be easily observed in ears for "camel_head" or also in bumps of "bumpy_plane" 
for High-freq detail transfer. The only possible downside is that this updated algorithm doesn't preserve volumes of objects in the process of
minimization distortion, but I suppose that in most cases this is expected behaviour that produces the most natural result.
