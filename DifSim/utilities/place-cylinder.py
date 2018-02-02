#!/usr/bin/env python
import random
import math
import vtk

meanRadius = 28.640487/2.0

def placeCylinder(x, y, cylinderList):
    for cylinder in cylinderList:
        [cx, cy] = cylinder
        if math.sqrt((cx - x)*(cx - x) + (cy - y)*(cy - y)) < 2.0*radius:
            return False
    cylinderList.append([x, y])
    return True

cylinderList = []
cylinderList.append([0.0, 0.0, 0.0, meanRadius])
cylinderCollection = vtk.vtkAppendPolyData()

cylinder = vtk.vtkCylinderSource()
cylinder.SetRadius(meanRadius)
cylinder.SetHeight(2*meanRadius)
cylinder.SetResolution(1000)

rotate = vtk.vtkTransform()
rotate.RotateX(90.0)

polydataTransform = vtk.vtkTransformPolyDataFilter()
polydataTransform.SetTransform(rotate)
polydataTransform.SetInput(cylinder.GetOutput())
cylinderZ = polydataTransform.GetOutput()

cylinderWriter = vtk.vtkPolyDataWriter()
cylinderWriter.SetFileName("single-cylinder.vtk")
cylinderWriter.SetInput(cylinderZ)
cylinderWriter.Write()

f = open("single-cylinder.mdl", 'w')
f.write("structure POLYGON_LIST {\n  VERTEX_LIST {\n")
tFilter = vtk.vtkTriangleFilter()
tFilter.SetInput(cylinderZ)
cylinderTriangulation = tFilter.GetOutput()
cylinderTriangulation.Update()
cylinderPointData = cylinderTriangulation.GetPoints()

for i in range(cylinderPointData.GetNumberOfPoints()):
    d = cylinderPointData.GetPoint(i)
    x, y, z = d
    f.write("    [ %f, %f, %f ]\n" % (d[0], d[1], d[2]))

f.write("  }\n")
f.write("  ELEMENT_CONNECTIONS {\n")
for i in range(cylinderTriangulation.GetNumberOfPolys()):
    f.write("    [ %d" % (cylinderTriangulation.GetCell(i).GetPointId(0)))
    for j in range(cylinderTriangulation.GetCell(i).GetNumberOfPoints() - 1):
        f.write(", %d" % (cylinderTriangulation.GetCell(i).GetPointId(j+1)))
    f.write(" ]\n")
f.write("  }\n  FULLY_CLOSED = YES\n  ADD_EFFECTOR {\n")
f.write("    STATE = water_transport.T\n")
f.write("    DENSITY = 1\n")
f.write("    ELEMENT = ALL_ELEMENTS\n")
f.write("    POLE_ORIENTATION = POSITIVE_BACK\n  }\n}\n")
f.close()
cylinderCollection.AddInput(cylinderZ)
cylinderGlyphs = vtk.vtkGlyph3D()
cylinderGlyphs.SetSource(cylinder.GetOutput())
cylinderGlyphs.SetScaleModeToScaleByScalar()

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

data = vtk.vtkFloatArray()
data.SetName("data")
cylinderCenters = vtk.vtkPoints()
cylinderCenters.SetNumberOfPoints(len(cylinderList))
cylinderPoints = vtk.vtkPolyVertex()
cylinderPoints.GetPointIds().SetNumberOfIds(len(cylinderList))
i = 0
for cylinderObject in cylinderList:
    print "[%f, %f, %f, %f]" % (cylinderObject[0], cylinderObject[1], cylinderObject[2], cylinderObject[3])
    cylinderCenters.InsertPoint(i, cylinderObject[0], cylinderObject[1], cylinderObject[2])
    data.InsertNextValue(cylinderObject[3])
    cylinderPoints.GetPointIds().SetId(i, i)
    i += 1

cylinderDistribution = vtk.vtkUnstructuredGrid()
cylinderDistribution.Allocate(1,1)
cylinderDistribution.InsertNextCell(cylinderPoints.GetCellType(),
                                   cylinderPoints.GetPointIds())
cylinderDistribution.SetPoints(cylinderCenters)
cylinderDistribution.GetPointData().AddArray(data)
cylinderDistribution.GetPointData().SetActiveScalars("data")

cylinderGlyphs.SetInput(cylinderDistribution)


cylinderMapper = vtk.vtkPolyDataMapper()
cylinderActor = vtk.vtkActor()
cylinderMapper.SetInput(cylinderCollection.GetOutput())
cylinderActor.SetMapper(cylinderMapper)
ren.AddActor(cylinderActor)

outlineBox = vtk.vtkOutlineSource()
voxelHalfWidth = 15.0
outlineBox.SetBounds(-voxelHalfWidth, voxelHalfWidth, -voxelHalfWidth, voxelHalfWidth, -voxelHalfWidth, voxelHalfWidth)
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInput(outlineBox.GetOutput())
outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
ren.AddActor(outlineActor)

ren.SetBackground(0.1, 0.2, 0.4)

iren.Initialize()
renWin.Render()
iren.Start()
# equation for an ellipsoid

# (x/a)^2 + (y/b)^2 + (z/c)^2 = 1.0
