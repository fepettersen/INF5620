'''
Consider a domain Omega = [0, 2] divided into the three elements [0, 1], [1, 1.2], and
[1.2, 2], with two nodes in each element (P1 elements). Suggest two different
element numberings and global node numberings for this mesh and set up the
corresponding nodes and elements lists in each case.
'''
#we can number the elements as follows

nodes = [0,1,1.2,2]
elements = [[0,1],[1,2],[2,3]]

#giving us
nodes[elements[1][0]] = nodes[1] = 1

'''Thereafter, subdivide the element [1.2, 2] into two new equal-sized elements.
Add the new node and the two new elements to each of the nodes and elements
lists.
'''
# Dividing element #2 in two new equally sized elements yields
nodes = [0,1,1.2,1.6,2]
elements = [[0,1],[1,2],[2,3],[3,4]]

# HPL: You can minimize the modifications by doing nodes.append(1.6) and
# elements[-1] = [2, 4]
# elements.append([4, 3])
# I have changed the text in the exercise to encourage this type of
# modifications.
