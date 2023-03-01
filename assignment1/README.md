# Assignment 1
Name: Svitlana Morkva

Legi-Nr: 22-944-284


## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Tasks
- Add a text dump of the content of the two data structures for the provided mesh “plane.off”.

- Show three screenshots of the 'fandisk.off' model using 'per-face shading', 'per-vertex shading' and 'per-corner shading'. Describe the difference between them.

- Show screenshots of the provided meshes with each connected component colored differently. Show the number of connected components and the size of each component (measured in number of faces) for all the provided models.

- Show screenshots of the subdivided meshes.


## Reports
### text dump of the data structure for "plane.off"
| Vertices-to-Face        | Vertices-to-Vertex     |
|:------------------------|:-----------------------|
| V 0: 5 21               | V 0: 11 13 22          |
| V 1: 31                 | V 1: 9 12              |
| V 2: 10 26              | V 2: 10 15 23          |
| V 3: 0                  | V 3: 14 16             |
| V 4: 14 27 30           | V 4: 9 10 20 24        |
| V 5: 13 23 29           | V 5: 11 12 19 24       |
| V 6: 1 4 17             | V 6: 13 14 18 21       |
| V 7: 2 8 18             | V 7: 15 16 17 21       |
| V 8: 6 9 12 19 22 25    | V 8: 17 18 19 20 22 23 |
| V 9: 15 30 31           | V 9: 1 4 12 24         |
| V 10: 11 26 27          | V 10: 2 4 20 23        |
| V 11: 7 21 23           | V 11: 0 5 19 22        |
| V 12: 15 29 31          | V 12: 1 5 9 24         |
| V 13: 4 5 20            | V 13: 0 6 18 22        |
| V 14: 0 1 16            | V 14: 3 6 16 21        |
| V 15: 8 10 24           | V 15: 2 7 17 23        |
| V 16: 0 2 16            | V 16: 3 7 14 21        |
| V 17: 3 8 9 18 19 24    | V 17: 7 8 15 18 21 23  |
| V 18: 3 4 6 17 19 20    | V 18: 6 8 13 17 21 22  |
| V 19: 7 12 13 22 23 28  | V 19: 5 8 11 20 22 24  |
| V 20: 11 12 14 25 27 28 | V 20: 4 8 10 19 23 24  |
| V 21: 1 2 3 16 17 18    | V 21: 6 7 14 16 17 18  |
| V 22: 5 6 7 20 21 22    | V 22: 0 8 11 13 18 19  |
| V 23: 9 10 11 24 25 26  | V 23: 2 8 10 15 17 20  |
| V 24: 13 14 15 28 29 30 | V 24: 4 5 9 12 19 20   |



### Show three screenshots of the 'fandisk.off' model using different shading. Make sure you disable the wireframe and that the screenshots clearly show the differences between the different shading types.
| model name  | per-face shading                                      | per-vertex shading                                       | per-corner shading                                       |
| :---------: |-------------------------------------------------------|----------------------------------------------------------|----------------------------------------------------------|
| fandisk     | <img align="center" src="./res/Face.png" width="300"> | <img align="center"  src="./res/Vertex.png" width="300"> | <img align="center"  src="./res/Corner.png" width="300"> |

#### Briefly describe the difference between the different shading types.
The Per-Face Shading is the easiest one, where all vertices of the face are colored the same, so the color is constant in the face. In this method, we only calculate the normals of the face. <br>
In the Per-Vertex Shading, we calculate normals for each vertex as an average of normals of the surrounding faces, because of that there is going to be a change in color inside faces.<br>
The Per-Corner Shading is a smart combination of the previous 2 methods. At first, we calculate a normal for each face corner so that each vertex can have 2 possibly different normals.
Then we compare these normals and if the angle between them is bigger than some threshold - it's a sharp corner, we leave it untouched, otherwise, we average normals.
Like that we can preserve the sharp edges and smoothness of a mesh simultaneously. I received the best results with an angle equal to 80.

### Assign different colors to each connected component
| model name   | your results                                             | no. of components |            no. of faces per component            |
| :----------: |----------------------------------------------------------|:-----------------:|:------------------------------------------------:|
|bumpy_cube    | <img align="center"  src="./res/Bumpy.png" width="300">  |         1         |                       2496                       |
|bunny         | <img align="center"  src="./res/Bunny.png" width="300">  |         1         |                      27864                       |
|coffeecup     | <img align="center"  src="./res/Coffee.png" width="300"> |         2         |                    3360/2304                     |
|honda         | <img align="center"  src="./res/Car.png" width="300">    |        11         | 90/192/192/13216/704/1088/1088/1088/1088/736/736 |



### Screenshots of subdivided meshes. Make sure you enable the wireframe overlay.
| model name | original shape                                               | subdivided shape                                              |
| :--------: |--------------------------------------------------------------|---------------------------------------------------------------| 
| plane      | <img align="center"  src="./res/Plane_org.png" width="300">  | <img align="center"  src="./res/Plane_subd.png" width="300">  |
| sphere     | <img align="center"  src="./res/Sphere_org.png" width="300"> | <img align="center"  src="./res/Sphere_subd.png" width="300"> |
| bunny      | <img align="center"  src="./res/Bunny_org.png" width="300">  | <img align="center"  src="./res/Bunny_subd.png" width="300">  |
| bumpy_cube | <img align="center"  src="./res/Bumpy_org.png" width="300">  | <img align="center"  src="./res/Bumpy_subd.png" width="300">  |




