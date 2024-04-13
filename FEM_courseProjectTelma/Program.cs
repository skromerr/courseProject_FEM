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
string resultPath = $"solution{xDroblenie}.txt";

gc.MakeGrid(radius, height, width, current, xDroblenie, bak);
Grid grid = new(spaceGridPath);
using (StreamWriter sw = new("log.txt", true))
{
    Console.WriteLine($"Сетка построена. Количество элементов - {grid.Elements.Length}. Количество узлов - {grid.Nodes.Count}.");
    sw.WriteLine($"Сетка построена. Количество элементов - {grid.Elements.Length}. Количество узлов - {grid.Nodes.Count}.");

    FEM fem = new(grid);
    fem.SetSlaeParametres(50000, 1e-15);

    //fem.Compute();
    //fem.PrintSolution(resultPath);

    fem.ReadSolutionFromFile(resultPath);

    Console.WriteLine("Впишите координаты точки, в которой хотите просмотреть значение. Формат: double double");

    string data;
    while ((data = Console.ReadLine()) != "0")
    {
        Console.WriteLine("-----------------------------------------------------------------------------------------------------");
        sw.WriteLine("-----------------------------------------------------------------------------------------------------");
        double[] values = data.Split(" ").Where(str => str != "").Select(Convert.ToDouble).ToArray();
        PointRZ point = new PointRZ(values[0], values[1]);
        double res = fem.AphiAtPoint(point);
        Console.WriteLine($"Значение компоненты Aphi в точке ({point.R};{point.Z}) равно {res:E4}.");
        sw.WriteLine($"Значение компоненты Aphi в точке ({point.R};{point.Z}) равно {res:E4}.");
        double res2 = fem.BzAtPoint(point);
        res = fem.BrAtPoint(point);
        Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({res:E4}, 0, {res2:E4}).");
        sw.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({res:E4}, 0, {res2:E4}).");
        res = fem.AbsBAtPoint(point);
        Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
        sw.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
        res = fem.AbsHAtPoint(point);
        Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");
        sw.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");
        Console.WriteLine("-----------------------------------------------------------------------------------------------------");
        sw.WriteLine("-----------------------------------------------------------------------------------------------------");
    }

}
//Console.WriteLine("-----------------------------------------------------------------------------------------------------");
//point = new PointRZ(0.7, -0.23);
//Console.WriteLine($"Значение компоненты Aphi в точке ({point.R};{point.Z}) равно {fem.AphiAtPoint(point)}.");
//Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
//res = fem.AbsBAtPoint(point);
//Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
//res = fem.AbsHAtPoint(point);
//Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

//Console.WriteLine("------------------------------------------------------------------------------------------------------");
//point = new PointRZ(1.17, -0.41);
//Console.WriteLine($"Значение компоненты Aphi в точке ({point.R};{point.Z}) равно {fem.AphiAtPoint(point)}.");
//Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
//res = fem.AbsBAtPoint(point);
//Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
//res = fem.AbsHAtPoint(point);
//Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

//Console.WriteLine("-------------------------------------------------------------------------------------------------------");
//point = new PointRZ(1.46, 0.0);
//Console.WriteLine($"Значение компоненты Aphi в точке ({point.R};{point.Z}) равно {fem.AphiAtPoint(point)}.");
//Console.WriteLine($"Вектор B в точке ({point.R};{point.Z}) равен ({fem.BrAtPoint(point)}, 0, {fem.BzAtPoint(point)}).");
//res = fem.AbsBAtPoint(point);
//Console.WriteLine($"Значение B в точке ({point.R};{point.Z}) равно {res}.");
//res = fem.AbsHAtPoint(point);
//Console.WriteLine($"Значение H в точке ({point.R};{point.Z}) равно {res}.");

//Console.WriteLine(GaussTriangle(ed, [new(0.0, 0.0), new(1.0, 2.0), new(1.0, 0.0)]));


//double ed(PointRZ point) => point.R * point.Z * point.R * point.Z;

//double GaussTriangle(Func<PointRZ, double> func, PointRZ[] vertices)
//{

//    const double x1a = 0.873821971016996;
//    const double x1b = 0.063089014491502;
//    const double x2a = 0.501426509658179;
//    const double x2b = 0.249286745170910;
//    const double x3a = 0.636502499121399;
//    const double x3b = 0.310352451033785;
//    const double x3c = 0.053145049844816;
//    const double w1 = 0.050844906370207;
//    const double w2 = 0.116786275726379;
//    const double w3 = 0.082851075618374;
//    double[] p1 = { x1a, x1b, x1b, x2a, x2b, x2b, x3a, x3b, x3a, x3c, x3b, x3c };
//    double[] p2 = { x1b, x1a, x1b, x2b, x2a, x2b, x3b, x3a, x3c, x3a, x3c, x3b };
//    double[] w = { w1, w1, w1, w2, w2, w2, w3, w3, w3, w3, w3, w3 };
//    double res = 0;

//    for (int i = 0; i < w.Length; i++)
//    {
//        PointRZ point = new();
//        point.R = (1 - p1[i] - p2[i]) * vertices[0].R + p1[i] * vertices[1].R + p2[i] * vertices[2].R;
//        point.Z = (1 - p1[i] - p2[i]) * vertices[0].Z + p1[i] * vertices[1].Z + p2[i] * vertices[2].Z;

//        res += func(point) * w[i];
//    }

//    return res /** 0.5*/;
//}