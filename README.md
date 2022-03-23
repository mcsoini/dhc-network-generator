# dhc-network-generator
Matlab scripts generating optimal network layouts for district heating and cooling based on OSM building and street data.

1. For a certain area (e.g. 1kmx1km pixel), load OSM data, i.e. building and street polygons
2. Make sure the street data forms a proper network, i.e. all intersections where we could lay a pipe are properly connected. 
3. Define street connection points by finding the closest point on a street for each building.
4. Simplify the street graph (e.g. remove all curves). The graph should ultimately only contain the street connection points as vertices and weighted edges corresponding to street distance. This is accomplished by iteratively simplifying the adjacency matrix.
5. Find the minimum spanning tree of the resulting graph. This corresponds to the minimum pipe-length network connecting all buildings.

Status: *Archived*
