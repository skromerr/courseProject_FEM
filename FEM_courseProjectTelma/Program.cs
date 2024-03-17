using CourseProject;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "grid.txt";
string telmaGridTxt = "test.txt";
GridConverter gc = new GridConverter(telmaGridTxt, spaceGridPath);
double radius = 1.3;
double height = 0.1;
double width = 0.1;
double current = 1e10;
gc.MakeGrid(radius, height, width, current);
Grid grid = new(spaceGridPath);
FEM fem = new(grid);
fem.Compute();
////fem.PrintSolution();



////Console.WriteLine();
PointRZ points = new PointRZ(0.7, 0.23);
Console.WriteLine($"Вектор B в точке ({points.R};{points.Z}) равен ({fem.BrAtPoint(points)}, 0, {fem.BzAtPoint(points)}).");
double res = fem.AbsBAtPoint(points);
Console.WriteLine($"Значение B в точке ({points.R};{points.Z}) равно {res}.");

res = fem.AbsHAtPoint(points);
Console.WriteLine($"Значение H в точке ({points.R};{points.Z}) равно {res}.");

points = new PointRZ(1.17, -0.41);
Console.WriteLine($"Вектор B в точке ({points.R};{points.Z}) равен ({fem.BrAtPoint(points)}, 0, {fem.BzAtPoint(points)}).");
res = fem.AbsBAtPoint(points);
Console.WriteLine($"Значение B в точке ({points.R};{points.Z}) равно {res}.");

res = fem.AbsHAtPoint(points);
Console.WriteLine($"Значение H в точке ({points.R};{points.Z}) равно {res}.");