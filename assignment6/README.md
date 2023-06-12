# Assignment 6: Skinning & Skeletal Animation

Name: Svitlana Morkva

Legi-Nr: 22-944-284


## Required results
Edit this 'README.md' file to report all your results. Use the ```/res``` folder to store your results.

### Tasks

1. Read Sec. 1 carefully to get familar with the data format, problem formulation, and mathematical background.
2. (Task 2) two theoretical discussions 
3. (Task 3) visualize the animation of the input skeleton of the hand shape from two types of rotations (rotation matrices and quaternion_mode)
4. (Task 4) compute harmonic skinning weights on selected handles
5. (Task 5) per-vertex LBS + rotation/translation + Lerp
6. (Task 6) per-vertex LBS + dual quaternion + Nlerp
7. (Task 7) per-face LBS + quaternion + Slerp + Poisson Stitching
8. (Task 8.1) context-aware per-vertex LBS
9. (optional Task 8.2) context-aware per-face LBS
 
### Important Note
1. We do not provide template code for this assignment - feel free to use previous template code if you want
2. You can use any libigl functions (and other open-source code, but please add a reference in this case)
3. You are allowed to use your previous code (for example, you will find the Poisson Stitching technqiue quite similar to Deformation Transfer that you have implemented in Assignment 5; and also the provided handle selection in A5 might inspire you to design your handle selection tool in Task 4).
4. You are allowed to modify this report freely (but please try to stick to some table format of orangizing figures to avoid a 20-page long report)
5. Please feel free to try other skeletal animation/deformation data you can find online if you find the provided animation is not cool enough (for example [here](https://www.mixamo.com/#/), but note that they might be in different data format than the provide ones).
6. Please try to keep your code clean and self-explained (with necessary comments), since we will grade this assignment based on your code as well (~~please don't feel upset when you get bad animations, you can always blame our poor data preparation~~).

## Reports

### Task 2: Rotation Representation discussion
#### Task 2.1. compare different rotation representations

| Representions        |                                                                         Short Description                                                                         |                                                                                         pros                                                                                          |                                                                                                             cons                                                                                                             |
| :------------------: |:-----------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| rotation matrix      |                           Matrix 3x3 with columns of unit length, mutually orthogonal, with determinant=1. All elements are parameters.                           | - Simple representation <br/> - Rotation is a linear function of the parameters.<br/> - Trivial to compute the rotation and partial derivatives <br/> - Can use linear optimizations. |                                         - Need to add 6 non-linear constraints for the matrix to still be in SO(3) <br/> - Redundant information (only 3 DOF but uses 9 parameters)                                          |
| euler angles         |             3 parameters that represent the rotation around one of the coordinate axes that uses 3 nonlinear functions to compute rotation matrices.              |                        - Simple representation <br/> - Rotations and derivatives are easy to compute. <br/> - Provides user-friendly interface for animators.                         | - "Gimbal lock" - rotation loses one DOF in some cases <br/> - Poor interpolation of rotation <br/> - Need to choose the rotation sequence <br/> - Unsuitable for many applications as can't naturally represent most joints |
| axis angle           |                                           Rotation is represented by an axis vector (usually normalized) and an angle.                                            |                                                                - Compact representation <br/> - Intuitive to work with                                                                |                                                    - Can't interpolate linearly <br/> - Numerical errors <br/> - Infinite number of pairs that refer to the same rotation                                                    |
| quaternion_mode          | Represented by 4-dimensional vector, where first 3 components encode (with some construnction form, not straightforward) a unit axis and 4th - a rotation radians |                - Parial derivatives are linearly independent - no gimbal lock <br/> - Supports nice interpolation<br/> - Computationally efficient to apply rotations                 |                                             - Need to re-normalize after every integration step and impose additional constraints to maintain unit length<br/> - Less intuitive                                              |

#### Task 2.2. Theoretical question for dual quaternion_mode

|                                                                                                       Euler angles -> rotation  matrix                                                                                                        |  rotation matrix -> quaternion  |    quaternion + translation -> dual quaternion   |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------: |:-----------------------------------------------: |
| At first, we need to construct 3 matrices for the corresponding angle. They will encode the rotation around each of the axes. Then we can multiply these matrices in some specified order (12 cases) and it will produce the rotation matrix. | xxxxxx                          |xxxxxx                                            | 

### Task 3: animation of the skeleton
| from rotaton matrix  |  from quaternion_mode   |
| :------------------: |:------------------: |
| <img align="center"  src="./res/placeholder.gif" width="300">  |<img align="center"  src="./res/placeholder.gif" width="300">  |

### Task 4: computing harmonic skinning weights on selected handles
#### Task 4.1. handle selection
| shape name           |  joint 1            |  joint 2            |  joint 3            |
| :------------------: |:------------------: |:------------------: |:------------------: |
| hand | <img align="center"  src="./res/placeholder.png" width="300">  | <img align="center"  src="./res/placeholder.png" width="300"> | <img align="center"  src="./res/placeholder.png" width="300"> |


#### Task 4.2. skinning weights visualization
| shape name           |  joint 1            |  joint 2            |  joint 3            |
| :------------------: |:------------------: |:------------------: |:------------------: |
| hand | <img align="center"  src="./res/placeholder.png" width="300">  | <img align="center"  src="./res/placeholder.png" width="300"> | <img align="center"  src="./res/placeholder.png" width="300"> |

### Task 5/6/7: skeletal animation 
| Task 5: per-vertex + rotation + Lerp   | Task 6: per-vertex + quaternion + Nlerp      | Task 7: per-face + quaternion + Slerp  |
| :---------:                            |        :---------:                           |       :---------:                      |
|<img align="center"  src="./res/placeholder.gif" width="300">  |<img align="center"  src="./res/placeholder.gif" width="300">  |  <img align="center"  src="./res/placeholder.gif" width="300">  |


Your comments (how these setups different from each other, and please compare the visual effects)

| Task 5: per-vertex + rotation + Lerp   | Task 6: per-vertex + quaternion + Nlerp      | Task 7: per-face + quaternion + Slerp  |
| :---------:                            |        :---------:                           |       :---------:                      |
| xxxxxx            |xxxxxx            |xxxxxx            |


### Task 8.1: context-aware per-vertex LBS
#### Task 8.1.1 visualize the unposed example shapes
| shape name           |  pose 1             |   pose 2            |   pose 3            |
| :------------------: |:------------------: |:------------------: |:------------------: |
| human | <img align="center"  src="./res/placeholder.png" width="300">  | <img align="center"  src="./res/placeholder.png" width="300"> | <img align="center"  src="./res/placeholder.png" width="300"> |

#### Task 8.1.2 skeletal animition using context

| without context   | with context     | 
| :---------:                            |        :---------:                           |  
|<img align="center"  src="./res/placeholder.gif" width="300">  |<img align="center"  src="./res/placeholder.gif" width="300">  |  
