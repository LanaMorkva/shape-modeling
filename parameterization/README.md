# Mesh Parameterization

Name: Svitlana Morkva

## Reports
### (mandatory-1) parameterization and checkerboard texture models
#### cathead
| Method            | checkerboard textured models                              | Parameterization                                             |
| :--------------:  |-----------------------------------------------------------|--------------------------------------------------------------|
| Uniform (fixed)   | <img align="center" src="./res/Cat1.png" width="400">     | <img align="center"  src="./res/UVCat1.png" width="300">     |
| Cotangent (fixed) | <img align="center" src="./res/Cat2.png" width="400">     | <img align="center"  src="./res/UVCat2.png" width="300">     |
| LSCM (fixed)      | <img align="center" src="./res/Cat3.png" width="400">     | <img align="center"  src="./res/UVCat3.png" width="300">     |
| ARAP (fixed)      | <img align="center" src="./res/Cat4.png" width="400">     | <img align="center"  src="./res/UVCat4.png" width="300">     |
| LSCM (free)       | <img align="center" src="./res/Cat3Free.png" width="400"> | <img align="center"  src="./res/UVCat3Free.png" width="300"> |
| ARAP (free)       | <img align="center" src="./res/Cat4Free.png" width="400"> | <img align="center"  src="./res/UVCat4Free.png" width="300"> |

#### hemisphere
| Method            | checkerboard textured models                               | Parameterization                                              |
| :--------------:  |------------------------------------------------------------|---------------------------------------------------------------|
| Uniform (fixed)   | <img align="center" src="./res/Hemi1.png" width="400">     | <img align="center"  src="./res/UVHemi1.png" width="300">     |
| Cotangent (fixed) | <img align="center" src="./res/Hemi2.png" width="400">     | <img align="center"  src="./res/UVHemi2.png" width="300">     |
| LSCM (fixed)      | <img align="center" src="./res/Hemi3.png" width="400">     | <img align="center"  src="./res/UVHemi3.png" width="300">     |
| ARAP (fixed)      | <img align="center" src="./res/Hemi4.png" width="400">     | <img align="center"  src="./res/UVHemi4.png" width="300">     |
| LSCM (free)       | <img align="center" src="./res/Hemi3Free.png" width="400"> | <img align="center"  src="./res/UVHemi3Free.png" width="300"> |
| ARAP (free)       | <img align="center" src="./res/Hemi4Free.png" width="400"> | <img align="center"  src="./res/UVHemi4Free.png" width="300"> |

#### hemisphere_non_convex_boundary
| Method            | checkerboard textured models                                | Parameterization                                               |
| :--------------:  |-------------------------------------------------------------|----------------------------------------------------------------|
| Uniform (fixed)   | <img align="center" src="./res/2Hemi1.png" width="400">     | <img align="center"  src="./res/UV2Hemi1.png" width="300">     |
| Cotangent (fixed) | <img align="center" src="./res/2Hemi2.png" width="400">     | <img align="center"  src="./res/UV2Hemi2.png" width="300">     |
| LSCM (fixed)      | <img align="center" src="./res/2Hemi3.png" width="400">     | <img align="center"  src="./res/UV2Hemi3.png" width="300">     |
| ARAP (fixed)      | <img align="center" src="./res/2Hemi4.png" width="400">     | <img align="center"  src="./res/UV2Hemi4.png" width="300">     |
| LSCM (free)       | <img align="center" src="./res/2Hemi3Free.png" width="400"> | <img align="center"  src="./res/UV2Hemi3Free.png" width="300"> |
| ARAP (free)       | <img align="center" src="./res/2Hemi4Free.png" width="400"> | <img align="center"  src="./res/UV2Hemi4Free.png" width="300"> |

#### Octo_cut2
| Method            | checkerboard textured models                               | Parameterization                                              |
| :--------------:  |------------------------------------------------------------|---------------------------------------------------------------|
| Uniform (fixed)   | <img align="center" src="./res/Octo1.png" width="500">     | <img align="center"  src="./res/UVOcto1.png" width="400">     |
| Cotangent (fixed) | <img align="center" src="./res/Octo2.png" width="500">     | <img align="center"  src="./res/UVOcto2.png" width="400">     |
| LSCM (fixed)      | <img align="center" src="./res/Octo3.png" width="500">     | <img align="center"  src="./res/UVOcto3.png" width="400">     |
| ARAP (fixed)      | <img align="center" src="./res/Octo4.png" width="500">     | <img align="center"  src="./res/UVOcto4.png" width="400">     |
| LSCM (free)       | <img align="center" src="./res/Octo3Free.png" width="500"> | <img align="center"  src="./res/UVOcto3Free.png" width="600"> |
| ARAP (free)       | <img align="center" src="./res/Octo4Free.png" width="500"> | <img align="center"  src="./res/UVOcto4Free.png" width="600"> |


### (mandatory-2) distortion visualization
#### cathead
| mtd \ metric | Conformal (angle)                                             | Authalic (area)                                             | Isometric  (length)                                         |
|:------------:|---------------------------------------------------------------|-------------------------------------------------------------|-------------------------------------------------------------|
| LSCM (free)  | <img align="center" src="./res/ConformFree3.png" width="300"> | <img align="center"  src="./res/AuthFree3.png" width="300"> | <img align="center"  src="./res/IsomFree3.png" width="300"> |
| ARAP (free)  | <img align="center" src="./res/ConformFree4.png" width="300"> | <img align="center"  src="./res/AuthFree4.png" width="300"> | <img align="center"  src="./res/IsomFree4.png" width="300"> |


#### hemisphere
| mtd \ metric | Conformal (angle)                                          | Authalic (area)                                             | Isometric  (length)                                         |
|:------------:|------------------------------------------------------------|-------------------------------------------------------------|-------------------------------------------------------------|
| LSCM (free)  | <img align="center" src="./res/ConfHemi3.png" width="300"> | <img align="center"  src="./res/AuthHemi3.png" width="300"> | <img align="center"  src="./res/IsomHemi3.png" width="300"> |
| ARAP (free)  | <img align="center" src="./res/ConfHemi.png" width="300">  | <img align="center"  src="./res/AuthHemi.png" width="300">  | <img align="center"  src="./res/IsomHemi.png" width="300">  |



