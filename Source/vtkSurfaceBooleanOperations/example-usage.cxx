
#include "vtkSurfaceBooleanOperations.h"

#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>

#include <vtkSmartPointer.h>
#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


int main(int argc, char ** argv) {

  vtkNew(vtkSphereSource,source1);
  source1->SetCenter(0.0, 0.0, 0.0);
  source1->SetRadius(1.0);

  vtkNew(vtkSphereSource,source2);
  source2->SetCenter(0.5, 0.0, 0.0);
  source2->SetRadius(1.0);


  vtkNew(vtkSurfaceBooleanOperations,booleanOperator);
  booleanOperator->AddInputConnection(source1->GetOutputPort());
  booleanOperator->AddInputConnection(source2->GetOutputPort());

  vtkNew(vtkPolyDataWriter,writer);
  writer->SetInputConnection(booleanOperator->GetOutputPort());

  booleanOperator->SetModeToUnion();
  writer->SetFileName("Union.vtk");
  writer->Update();

  booleanOperator->SetModeToIntersection();
  writer->SetFileName("Intersection.vtk");
  writer->Update();

  booleanOperator->SetModeToDifference();
  writer->SetFileName("Difference.vtk");
  writer->Update();

  return EXIT_SUCCESS;
}
