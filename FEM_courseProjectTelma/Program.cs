using CourseProject;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "grid.txt";
string telmaGridTxt = "test.txt";
GridConverter gc = new GridConverter(telmaGridTxt, spaceGridPath);
gc.ConvertGrid();
//Grid grid = new(spaceGridPath);
//FEM fem = new(grid);
//fem.Compute();
////fem.PrintSolution();



////Console.WriteLine();
//Point2D points = new Point2D(2, 2.8);
//double res = fem.XValueAtPoint(points);
//Console.WriteLine($"Значение функции в точке ({points.X};{points.Y}) равно {res}.");

//res = fem.YValueAtPoint(points);
//Console.WriteLine($"Значение функции в точке ({points.X};{points.Y}) при t = {0.5} равно {res}.");