using CourseProject;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "grid.txt";
string timeGridPath = "timeGrid.txt";
Grid grid = new(spaceGridPath, timeGridPath);
FEM fem = new(grid, new Test1(), true);
fem.Compute();
//fem.PrintSolution();

fem.PrintResearchSolution();


//Console.WriteLine();
PointRZ points = new PointRZ(2, 2.8);
double res = fem.ValueAtPoint(points);
Console.WriteLine($"Значение функции в точке ({points.R};{points.Z}) при t = {grid.Time[^1]} равно {res}.");
res = fem.test.U(points, grid.Time[^1]);
Console.WriteLine($"Точное значение функции в точке ({points.R};{points.Z}) при t = {grid.Time[^1]} равно {res}.");

res = fem.ValueAtPointAndHalfTime(points);
Console.WriteLine($"Значение функции в точке ({points.R};{points.Z}) при t = {0.5} равно {res}.");
res = fem.test.U(points, 0.5);
Console.WriteLine($"Точное значение функции в точке ({points.R};{points.Z}) при t = {0.5} равно {res}.");
Console.WriteLine();