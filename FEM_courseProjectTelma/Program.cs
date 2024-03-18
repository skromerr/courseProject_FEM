using CourseProject;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "grid.txt";
string telmaGridTxt = "test.txt";
GridConverter gc = new GridConverter(telmaGridTxt, spaceGridPath);
double radius = 1;
double height = 0.1;
double width = 0.1;
double current = 1e7;
int xDroblenie = 8;
double bak = 10;
gc.MakeGrid(radius, height, width, current, xDroblenie, bak);
Grid grid = new(spaceGridPath);
FEM fem = new(grid);
fem.Compute();
////fem.PrintSolution();



////Console.WriteLine();
PointRZ point = new PointRZ(0.7, 0.23);
Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
double res = fem.AbsBAtPoint(point);
Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
res = fem.AbsHAtPoint(point);
Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

point = new PointRZ(0.7, -0.23);
Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
res = fem.AbsBAtPoint(point);
Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
res = fem.AbsHAtPoint(point);
Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

point = new PointRZ(1.17, -0.41);
Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
res = fem.AbsBAtPoint(point);
Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
res = fem.AbsHAtPoint(point);
Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

point = new PointRZ(1.46, 0.0);
Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
res = fem.AbsBAtPoint(point);
Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
res = fem.AbsHAtPoint(point);
Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");