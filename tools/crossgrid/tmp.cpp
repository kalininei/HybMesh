  // Create a cell array to store the polygon in
  vtkSmartPointer<vtkCellArray> aCellArray =
    vtkSmartPointer<vtkCellArray>::New();
 
  // Define a polygonal hole with a clockwise polygon
  vtkSmartPointer<vtkPolygon> aPolygon =
   vtkSmartPointer<vtkPolygon>::New();
 
  aPolygon->GetPointIds()->InsertNextId(22);
  aPolygon->GetPointIds()->InsertNextId(23);
  aPolygon->GetPointIds()->InsertNextId(24);
  aPolygon->GetPointIds()->InsertNextId(25);
  aPolygon->GetPointIds()->InsertNextId(35);
  aPolygon->GetPointIds()->InsertNextId(45);
  aPolygon->GetPointIds()->InsertNextId(44);
  aPolygon->GetPointIds()->InsertNextId(43);
  aPolygon->GetPointIds()->InsertNextId(42);
  aPolygon->GetPointIds()->InsertNextId(32);
 
  aCellArray->InsertNextCell(aPolygon);
 
  // Create a polydata to store the boundary. The points must be the
  // same as the points we will triangulate.
  vtkSmartPointer<vtkPolyData> boundary =
   vtkSmartPointer<vtkPolyData>::New();
  boundary->SetPoints(aPolyData->GetPoints());
  boundary->SetPolys(aCellArray);
 
  // Triangulate the grid points
  vtkSmartPointer<vtkDelaunay2D> delaunay =
   vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
  delaunay->SetInput(aPolyData);
  delaunay->SetSource(boundary);
#else
  delaunay->SetInputData(aPolyData);
  delaunay->SetSourceData(boundary);
#endif
  delaunay->Update();
 
  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> meshMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  meshMapper->SetInputConnection(delaunay->GetOutputPort());
 
  vtkSmartPointer<vtkActor> meshActor =
    vtkSmartPointer<vtkActor>::New();
  meshActor->SetMapper(meshMapper);
  //meshActor->GetProperty()->SetEdgeColor(0,0,1); // Why aren't the edges aren't visible unless we set the representation to wireframe?
  //meshActor->GetProperty()->SetInterpolationToFlat();
  meshActor->GetProperty()->SetRepresentationToWireframe();
 
  vtkSmartPointer<vtkPolyDataMapper> boundaryMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  boundaryMapper->SetInputConnection(boundary->GetProducerPort());
#else
  boundaryMapper->SetInputData(boundary);
#endif
 
  vtkSmartPointer<vtkActor> boundaryActor =
    vtkSmartPointer<vtkActor>::New();
  boundaryActor->SetMapper(boundaryMapper);
  boundaryActor->GetProperty()->SetColor(1,0,0);
 
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // Add the actor to the scene
  renderer->AddActor(meshActor);
  renderer->AddActor(boundaryActor);
  renderer->SetBackground(.3, .6, .3); // Background color green
 
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
