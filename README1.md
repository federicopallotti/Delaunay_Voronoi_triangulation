Before run the code make sure a txt file (like "ExampleInput.txt") of x and y coordinates of points in convex position is included in the project folder.
To run the code use the command "python3 -i delaunay_voronoi.py < ExampleInput.txt" to plot both dealaunay and vornoi graph.
To run only one of those use either "python3 -i delaunay_voronoi.py < ExampleInput.txt delaunay" or "python3 -i delaunay_voronoi.py < ExampleInput.txt voronoi" to plot only one of the two graphs.
The code plots the desired graph, prints the DCEL infos and writes an output txt file with the half edges vertices of the computed delaunay triangulation.
