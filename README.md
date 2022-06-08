The code computes the Delaunay Traingulation and Voronoi Diagram out of a set of input points given in convex position, in linear time.

Before run the code make sure a txt file (like "ExampleInput.txt") of x and y coordinates of points in convex position is included in the project folder.
To run the code use the command "python3 -i delaunay_voronoi.py < ExampleInput.txt" to plot both dealaunay and vornoi graph.
To run only one of those use either "python3 -i delaunay_voronoi.py < ExampleInput.txt delaunay" or "python3 -i delaunay_voronoi.py < ExampleInput.txt voronoi" to plot only one of the two graphs.
The code plots the desired graph, prints the DCEL infos and writes an output txt file with the half edges vertices of the computed delaunay triangulation.

These are the output plots:

![delauny](https://user-images.githubusercontent.com/77103965/172165893-95ea6378-d912-4bc0-8810-2784e76a1394.png)
![voronoi_result](https://user-images.githubusercontent.com/77103965/172165904-902d7a7a-7c34-4923-a9f1-c1cf45269e78.png)
