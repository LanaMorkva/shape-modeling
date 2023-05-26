# Assignment 5

Name: Svitlana Morkva

Legi-Nr: 22-944-284

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Tasks

1) **Multiresolution Mesh Editing**: Provide screenshots for 4 different deformed meshes (woody-lo, woody-hi, hand, and cylinder). For each example, provide a rendering of S, B, B' and S'. Please make sure to include rotation deformation in addition to translation for the hand and cylinder meshes. (questions 1.1 - 1.4 of the assignment sheet)

2) **Real time mesh editing**: Provide animated gifs or short videos for 4 different deformed meshes (bar, bumpy_plane, camel_head, and cactus) showing that your algorithm can perform in real time. Please follow the more detailed instructions above the output table for this section. (question 1.5 of the assignment sheet)

3) **Deformation transfer**: Discuss and show the differences to the results obtained with the high-frequency detail transfer from part 1.4. on 4 different meshes (bar, bumpy_plane, camel_head, and cactus). (part 2 of the assignment sheet)



## Report
### 1.1 - 1.4 Multiresolution Mesh Editing
|              model              | S                                                           | B                                                            | B'                                                           | S'                                                            |
|:-------------------------------:|-------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|---------------------------------------------------------------|
|            woody-lo             | <img align="center" src="./res/Woody_lo_S.png" width="300"> | <img align="center"  src="./res/Woody_lo_B.png" width="300"> | <img align="center" src="./res/Woody_lo_B2.png" width="300"> | <img align="center"  src="./res/Woody_lo_S2.png" width="300"> |
|            woody-hi             | <img align="center" src="./res/Woody_hi_S.png" width="300"> | <img align="center"  src="./res/Woody_hi_B.png" width="300"> | <img align="center" src="./res/Woody_hi_B2.png" width="300"> | <img align="center"  src="./res/Woody_hi_S2.png" width="300"> |
|   hand (translation+rotation)   | <img align="center" src="./res/Hand_S.png" width="300">     | <img align="center"  src="./res/Hand_B.png" width="300">     | <img align="center" src="./res/Hand_B2.png" width="300">     | <img align="center"  src="./res/Hand_S2.png" width="300">     |
| cylinder (translation+rotation) | <img align="center" src="./res/Cylinder_S.png" width="300"> | <img align="center"  src="./res/Cylinder_B.png" width="300"> | <img align="center" src="./res/Cylinder_B2.png" width="300"> | <img align="center"  src="./res/Cylinder_S2.png" width="300"> |

### 1.5 Real time mesh editing

Show real time mesh editing using animated gifs or short videos. *Max 15 seconds per gif/video*. **You should at least showcase the following**: move multiple handles sequentially, show a rotation for at least 1 handle per mesh, add a new handle after you have already deformed using another handle, and include a large deformation where you move 1 handle far away and back again (without letting go in between). 

> You might find these links helpful for screen recording: [MacOS](https://support.apple.com/en-us/HT208721), [ubuntu](https://askubuntu.com/questions/4428/how-can-i-record-my-screen), [windows](https://support.microsoft.com/en-us/office/record-the-screen-d70508e8-25a3-4b97-b78a-a467b5372e21). Feel free to use any software for the recordings - please make sure to trim/crop the final video properly.

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
