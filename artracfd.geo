#------------------------------------------------------------------------------
#                                                                             -
#                   Internal Geometry Configuration                           -
#                                                                             -
# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen - 
#   plane; X is horizontal from west to east; Y is vertical from south to     -
#   north; Z axis is perpendicular to the screen and points from front to     -
#   back; The origin locates at the west-south-front corner of the            -
#   computational domain;                                                     - 
# - All coordinate values and velocity values here are relative values, they  -
#   should be normalized by the reference length and velocity respectively.   -                                       -
# - Please double check your input!                                           -
#                                                                             -
#------------------------------------------------------------------------------
geometry begin
#------------------------------------------------------------------------------
#
#                   >> Total Number of Objects <<
#
#------------------------------------------------------------------------------
count begin
1            # total number of objects (integer) 
count end
#------------------------------------------------------------------------------
#
#                    >> Geometry Information <<
#
#------------------------------------------------------------------------------
circle begin
3, 3, 0, 0.5, 0, 0, 0   # x, y, z, radius, u, v, w
circle end
#------------------------------------------------------------------------------
geometry end
#------------------------------------------------------------------------------
/* a good practice: end file with a newline */

