import numpy as np

class Point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
        pass
    def __repr__(self):
        string = ""
        string+= "{:21.14E}, {:21.14E}, {:21.14E}\n".format(self.x, self.y, self.z)
        return string

class Element:
    def __init__(self, points):
        self.points = points 
        pass

filename = "mesh/shape.vtk"

try:
    with open(filename) as vtk_file:
        text = vtk_file.read()
    pass
except:
    print("error: unable to open mesh file!")
    exit(1)
    
print(f"loading {filename}..., done!")
print("extracting mesh file elements...")

text = text.splitlines()

Points_list = []
Elements_list = []
Entities_list = []

flag_Points = False
flag_Cells  = False
flag_Entities  = False

for line in text:
    if line[0:6]=="POINTS":
        flag_Points = True
        continue
    
    if flag_Points:
        if line=="":
            flag_Points = False
            continue
        new_line = line.split(" ")
        Points_list.append(Point(float(new_line[0]), float(new_line[1]), float(new_line[2])))
        pass
    
    if line[0:5]=="CELLS":
        flag_Cells = True
        continue

    if flag_Cells:
        if line=="":
            flag_Cells = False
            continue
        new_line = line.split(" ")
        points = []
        for word in new_line:
            points.append(int(word))
            continue
        Elements_list.append(points)
        pass
    
    if line[0:6]=="LOOKUP":
        flag_Entities = True
        continue
    
    if flag_Entities:
        if line=="":
            flag_Entities = False
            continue
        new_line = line
        Entities_list.append(int(new_line))
        pass

    continue

file_0D = open("mesh/mesh/elements_0d.txt", "w")
file_1D = open("mesh/mesh/elements_1d.txt", "w")
file_2D = open("mesh/mesh/elements_2d.txt", "w")
file_3D = open("mesh/mesh/elements_3d.txt", "w")

counter_0D = 0
counter_1D = 0
counter_2D = 0
counter_3D = 0

counter_element = 0
for element in Elements_list:
    if len(element)==2:
        file_0D.write("{:21.14E} ".format(Points_list[element[1]].x))
        file_0D.write("{:21.14E} ".format(Points_list[element[1]].y))
        file_0D.write("{:21.14E} ".format(Points_list[element[1]].z))
        if len(Entities_list)>0:
            file_0D.write("{:d} ".format(Entities_list[counter_element]))
            pass
        file_0D.write("\n")
        counter_0D+=1
        counter_element+=1
        pass
    if len(element)==3:
        file_1D.write("{:21.14E} ".format(Points_list[element[1]].x))
        file_1D.write("{:21.14E} ".format(Points_list[element[1]].y))
        file_1D.write("{:21.14E} ".format(Points_list[element[1]].z))
        file_1D.write("{:21.14E} ".format(Points_list[element[2]].x))
        file_1D.write("{:21.14E} ".format(Points_list[element[2]].y))
        file_1D.write("{:21.14E} ".format(Points_list[element[2]].z))
        if len(Entities_list)>0:
            file_1D.write("{:d} ".format(Entities_list[counter_element]))
        else:
            file_1D.write("{:d} ".format(-1))
            pass
        file_1D.write("\n")
        counter_1D+=1
        counter_element+=1
        pass
    if len(element)==4:
        file_2D.write("{:21.14E} ".format(Points_list[element[1]].x))
        file_2D.write("{:21.14E} ".format(Points_list[element[1]].y))
        file_2D.write("{:21.14E} ".format(Points_list[element[1]].z))
        file_2D.write("{:21.14E} ".format(Points_list[element[2]].x))
        file_2D.write("{:21.14E} ".format(Points_list[element[2]].y))
        file_2D.write("{:21.14E} ".format(Points_list[element[2]].z))
        file_2D.write("{:21.14E} ".format(Points_list[element[3]].x))
        file_2D.write("{:21.14E} ".format(Points_list[element[3]].y))
        file_2D.write("{:21.14E} ".format(Points_list[element[3]].z))
        if len(Entities_list)>0:
            file_2D.write("{:d} ".format(Entities_list[counter_element]))
        else:
            file_2D.write("{:d} ".format(-1))
            pass
        file_2D.write("\n")
        counter_2D+=1
        counter_element+=1
        pass
    if len(element)==5:
        file_3D.write("{:21.14E} ".format(Points_list[element[1]].x))
        file_3D.write("{:21.14E} ".format(Points_list[element[1]].y))
        file_3D.write("{:21.14E} ".format(Points_list[element[1]].z))
        file_3D.write("{:21.14E} ".format(Points_list[element[2]].x))
        file_3D.write("{:21.14E} ".format(Points_list[element[2]].y))
        file_3D.write("{:21.14E} ".format(Points_list[element[2]].z))
        file_3D.write("{:21.14E} ".format(Points_list[element[3]].x))
        file_3D.write("{:21.14E} ".format(Points_list[element[3]].y))
        file_3D.write("{:21.14E} ".format(Points_list[element[3]].z))
        file_3D.write("{:21.14E} ".format(Points_list[element[4]].x))
        file_3D.write("{:21.14E} ".format(Points_list[element[4]].y))
        file_3D.write("{:21.14E} ".format(Points_list[element[4]].z))
        if len(Entities_list)>0:
            file_3D.write("{:d} ".format(Entities_list[counter_element]))
        else:
            file_3D.write("{:d} ".format(-1))
            pass
        file_3D.write("\n")
        counter_3D+=1
        counter_element+=1
        pass
    continue

file_0D.close()
file_1D.close()
file_2D.close()
file_3D.close()

print(f"total number of 0d elements: {counter_0D}")
print(f"total number of 1d elements: {counter_1D}")
print(f"total number of 2d elements: {counter_2D}")
print(f"total number of 3d elements: {counter_3D}")

if counter_0D+counter_1D+counter_2D+counter_3D==0:
    print("error: empty mesh file!")
    exit(1)

file = open("mesh/mesh/info.txt", "w")
file.write(f"{counter_0D} {counter_1D} {counter_2D} {counter_3D}\n")
if len(Entities_list)>0:
    file.write("1")
else:
    file.write("0")
    pass
file.close()

print("extracting mesh file elements was successful!")

exit(0)