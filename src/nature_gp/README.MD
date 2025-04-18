# Nature Geometry Processor Module

## Sub Modules

### XYZ
#### Contains code for handling 3D coordinate information

| Method                            | Implemented? | Notes                                      |
|----------------------------------|--------------|--------------------------------------------|
| cross_square_magnitude           | ❌           | Not implemented                            |
| set_linear_form02                | ❌           | Not implemented                            |
| set_linear_form12                | ✅           | Implemented                                |
| set_linear_form22                | ✅           | Implemented                                |
| set_linear_form23                | ✅           | Implemented                                |
| set_linear_form33                | ✅           | Implemented                                |
| set_linear_form34                | ✅           | Implemented                                |
| cross_crossed                    | ✅           | Implemented                                |
| cross_cross                      | ✅           | Implemented                                |
| dot_cross                        | ✅           | Implemented                                |
| crossed                          | ✅           | Implemented                                |
| subtracted                       | ✅           | Implemented                                |
| multiplied                       | ✅           | Implemented                                |
| divided                          | ✅           | Implemented                                |
| reversed                         | ✅           | Implemented                                |
| reverse                          | ✅           | Implemented                                |
| added                            | ✅           | Implemented                                |
| coords                           | ✅           | Implemented                                |
| change_coord                     | ✅           | Implements index validation                |
| coord_by_index                   | ✅           | Implements index validation                |
| set_coord_by_index               | ✅           | Implements index validation                |
| from                             | ✅           | Implemented                                |
| new                              | ✅           | Implemented                                |
| zero                             | ✅           | Implemented                                |
| multiply_xyz                     | ✅           | Implemented                                |
| x                                | ✅           | Implemented                                |
| y                                | ✅           | Implemented                                |
| z                                | ✅           | Implemented                                |
| set_x                            | ✅           | Implemented                                |
| set_y                            | ✅           | Implemented                                |
| set_z                            | ✅           | Implemented                                |
| set_coords                       | ✅           | Implemented                                |
| add                              | ✅           | Implemented                                |
| subtract                         | ✅           | Implemented                                |
| multiply                         | ✅           | Implemented                                |
| divide                           | ✅           | Implemented                                |
| dot                              | ✅           | Implemented                                |
| cross                            | ✅           | Implemented                                |
| cross_magnitude                  | ❌           | Not implemented                            |
| normalize                        | ✅           | Implements `NGP::resolution()` check       |
| normalized                       | ✅           | Implements `NGP::resolution()` check       |
| modulus                          | ✅           | Implemented                                |
| square_modulus                   | ❌           | Only `todo!()`                             |
| is_equal                         | ✅           | Uses `tolerance` and checks all components |
