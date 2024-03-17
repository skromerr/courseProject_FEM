namespace CourseProject;

public class FEM
{
    private Grid grid;
    private SparseMatrix globalMatrix;
    private Vector globalVectorX;
    private Vector globalVectorY;
    private Vector solutionX;
    private Vector solutionY;
    private Solver slae;
    private Matrix alphas;
    private Vector localVector;
    private Matrix stiffnessMatrix;
    private Matrix massMatrix;
    private Point2D[] vertices;
    private FirstCondition[] firstConditions;
    private SecondCondition[] secondConditions;
    private string condPath => "conditions.txt";

    public FEM(Grid grid)
    {
        this.grid = grid;
        alphas = new(3);
        stiffnessMatrix = new(6);
        massMatrix = new(6);
        localVector = new(6);
        slae = new Solver(1000, 1e-16);

        vertices = new Point2D[3];
        solutionX = new Vector(grid.Nodes.Count);
        solutionY = new Vector(grid.Nodes.Count);

        SetBoundaryConditions();
    }

    private void SetBoundaryConditions()
    {
        int[] cond = File.ReadAllText(condPath).Split(' ').Select(value => Convert.ToInt32(value)).ToArray();
        int i;
        HashSet<int> fstCond = new HashSet<int>();
        HashSet<(int, int, int[])> sndCond = new HashSet<(int, int, int[])>();
        for (i = 0; i < 4; i++)
        {
            foreach (var edge in grid.Boundaries[i])
                if (cond[i] == 1)
                    foreach (var node in edge.Item2)
                        fstCond.Add(node);
                else
                    sndCond.Add((edge.Item1, i, edge.Item2));
        }

        firstConditions = new FirstCondition[fstCond.Count];
        i = 0;
        foreach (var node in fstCond)
            firstConditions[i++] = new(grid.Nodes[node], node);

        secondConditions = new SecondCondition[sndCond.Count];
        i = 0;
        foreach (var edge in sndCond)
            secondConditions[i++] = new(edge.Item1, edge.Item2, edge.Item3);
    }

    public void Compute()
    {

        BuildPortrait();


        AssemblyGlobalMatrix();
        //AccountSecondConditions();
        AccountFirstConditions();

        slae.SetSLAE(globalVectorX, globalMatrix);
        slae.CGM();
        Vector.Copy(slae.solution, solutionX);

        slae.SetSLAE(globalVectorY, globalMatrix);
        slae.CGM();
        Vector.Copy(slae.solution, solutionY);
    }

    public void AccountFirstConditions()
    {
        foreach (var fc in firstConditions)
        {
            globalMatrix.Di[fc.NodeNumber] = 1;
            globalVectorX[fc.NodeNumber] = 0;
            globalVectorY[fc.NodeNumber] = 0;
            for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
            {
                globalVectorX[globalMatrix.Jg[i]] -= 0;
                globalMatrix.Gg[i] = 0;
            }
            for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
            {
                for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
                {
                    if (globalMatrix.Jg[j] == fc.NodeNumber)
                    {
                        globalVectorX[i] -= 0;
                        globalMatrix.Gg[j] = 0;
                    }
                }
            }
        }
    }

    //public void AccountSecondConditions(int itime)
    //{
    //    foreach (var sc in secondConditions)
    //    {
    //        localVector.Fill(0);
    //        int ielem = sc.ElemNumber;

    //        vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
    //        vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
    //        vertices[2] = grid.Nodes[grid.Elements[ielem][2]];
    //        CalcuclateAlphas();

    //        PointRZ start = new(grid.Nodes[sc.Edge[0]].R, grid.Nodes[sc.Edge[0]].Z);
    //        PointRZ end = new(grid.Nodes[sc.Edge[2]].R, grid.Nodes[sc.Edge[2]].Z);
    //        double lambda = grid.Materials[grid.Elements[ielem][^1]].Lambda;
    //        int localNum;

    //        for (int i = 0; i < 3; i++)
    //        {
    //            localNum = GetLocalNumber(ielem, sc.Edge[i]);
    //            localVector[localNum] = lambda * GaussEdge(localNum, sc.EdgeType, start, end);
    //        }

    //        AddElementToVector(ielem);
    //    }
    //}

    //private int GetLocalNumber(int ielem, int globalNumber)
    //{
    //    for (int i = 0; i < 6; i++)
    //    {
    //        if (grid.Elements[ielem][i] == globalNumber) return i;
    //    }

    //    throw new Exception($"У элемента {ielem} нет узла с номером {globalNumber}.");
    //}

    public void BuildPortrait()
    {
        var list = new HashSet<int>[grid.Nodes.Count].Select(_ => new HashSet<int>()).ToList();
        foreach (var element in grid.Elements)
        {
            foreach (var pos in element.Nodes)
            {
                foreach (var node in element.Nodes)
                {
                    if (pos > node)
                    {
                        list[pos].Add(node);
                    }
                }
            }
        }

        int count = list.Sum(childList => childList.Count);

        globalMatrix = new(grid.Nodes.Count, count);
        globalVectorX = new(grid.Nodes.Count);
        globalVectorY = new(grid.Nodes.Count);

        globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Count; i++)
            globalMatrix.Ig[i + 1] = globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childList in list)
        {
            foreach (var value in childList.Order())
            {
                globalMatrix.Jg[k++] = value;
            }
        }
    }

    private void AddElement(int i, int j, double value)
    {
        if (i == j)
        {
            globalMatrix.Di[i] += value;
            return;
        }

        for (int icol = globalMatrix.Ig[i]; icol < globalMatrix.Ig[i + 1]; icol++)
        {
            if (globalMatrix.Jg[icol] == j)
            {
                globalMatrix.Gg[icol] += value;
                return;
            }
        }
    }

    private void AssemblyGlobalMatrix()
    {

        globalMatrix.Clear();
        globalVectorX.Fill(0);
        globalVectorY.Fill(0);

        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            AssemblyLocalMatrixes();

            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    AddElement(grid.Elements[ielem][i], grid.Elements[ielem][j], stiffnessMatrix[i, j]);

            AssemblyLocallVector(ielem);

            AddElementToVector(ielem);

            stiffnessMatrix.Clear();
            massMatrix.Clear();
            localVector.Fill(0);
        }
    }

    void AssemblyLocallVector(int ielem)
    {
        for (int i = 0; i < 6; i++)
        {
            double psi(Point2D point)
                => GetPsi(point, i);
            localVector[i] = GaussTriangle(psi);
        }

        localVector = Math.Abs(DeterminantD()) * localVector;
    }

    private void AddElementToVector(int ielem)
    {
        double J = grid.Materials[grid.Elements[ielem].MaterialNumber].J;
        double Jx = J;
        double Jy = J;
        for (int i = 0; i < 6; i++)
        {
            globalVectorX[grid.Elements[ielem][i]] += Jx * localVector[i];
            globalVectorY[grid.Elements[ielem][i]] += Jy * localVector[i];
        }
    }

    private double GaussTriangle(Func<Point2D, double> func)
    {

        const double x1a = 0.873821971016996;
        const double x1b = 0.063089014491502;
        const double x2a = 0.501426509658179;
        const double x2b = 0.249286745170910;
        const double x3a = 0.636502499121399;
        const double x3b = 0.310352451033785;
        const double x3c = 0.053145049844816;
        const double w1 = 0.050844906370207;
        const double w2 = 0.116786275726379;
        const double w3 = 0.082851075618374;
        double[] p1 = { x1a, x1b, x1b, x2a, x2b, x2b, x3a, x3b, x3a, x3c, x3b, x3c };
        double[] p2 = { x1b, x1a, x1b, x2b, x2a, x2b, x3b, x3a, x3c, x3a, x3c, x3b };
        double[] w = { w1, w1, w1, w2, w2, w2, w3, w3, w3, w3, w3, w3 };
        double res = 0;

        for (int i = 0; i < w.Length; i++)
        {
            Point2D point = new();
            point.X = (1 - p1[i] - p2[i]) * vertices[0].X + p1[i] * vertices[1].X + p2[i] * vertices[2].X;
            point.Y = (1 - p1[i] - p2[i]) * vertices[0].Y + p1[i] * vertices[1].Y + p2[i] * vertices[2].Y;

            res += func(point) * w[i];
        }

        return res * 0.5;
    }

    //private double GaussEdge(int numberPsi, Func<PointRZ, double, int, double> theta, int funcNum, PointRZ fstPoint, PointRZ sndPoint)
    //{
    //    double[] p = { 0.0,
    //                 1.0 / 3.0 * Math.Sqrt(5 - 2 * Math.Sqrt(10.0 / 7.0)),
    //                 -1.0 / 3.0 * Math.Sqrt(5 - 2 * Math.Sqrt(10.0 / 7.0)),
    //                 1.0 / 3.0 * Math.Sqrt(5 + 2 * Math.Sqrt(10.0 / 7.0)),
    //                 -1.0 / 3.0 * Math.Sqrt(5 + 2 * Math.Sqrt(10.0 / 7.0))};

    //    double[] w = { 128.0 / 225.0,
    //                        (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
    //                        (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
    //                        (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0,
    //                        (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0 };
    //    double res = 0;

    //    double lengthEdge = Math.Sqrt((fstPoint.R - sndPoint.R) * (fstPoint.R - sndPoint.R) +
    //                          (fstPoint.Z - sndPoint.Z) * (fstPoint.Z - sndPoint.Z));

    //    for (int i = 0; i < w.Length; i++)
    //    {
    //        PointRZ point = new();
    //        point.R = (sndPoint.R - fstPoint.R) * (1 + p[i]) / 2.0 + fstPoint.R;
    //        point.Z = (sndPoint.Z - fstPoint.Z) * (1 + p[i]) / 2.0 + fstPoint.Z;

    //        res += theta(point, t, funcNum) * basis(point, numberPsi) * w[i] * point.R;
    //    }

    //    return lengthEdge * res / 2.0;
    //}

    private double DeterminantD()
         => (vertices[1].X - vertices[0].X) * (vertices[2].Y - vertices[0].Y) -
            (vertices[2].X - vertices[0].X) * (vertices[1].Y - vertices[0].Y);

    private void CalcuclateAlphas()
    {
        double dD = DeterminantD();

        alphas[0, 0] = (vertices[1].X * vertices[2].Y - vertices[2].X * vertices[1].Y) / dD;
        alphas[0, 1] = (vertices[1].Y - vertices[2].Y) / dD;
        alphas[0, 2] = (vertices[2].X - vertices[1].X) / dD;

        alphas[1, 0] = (vertices[2].X * vertices[0].Y - vertices[0].X * vertices[2].Y) / dD;
        alphas[1, 1] = (vertices[2].Y - vertices[0].Y) / dD;
        alphas[1, 2] = (vertices[0].X - vertices[2].X) / dD;

        alphas[2, 0] = (vertices[0].X * vertices[1].Y - vertices[1].X * vertices[0].Y) / dD;
        alphas[2, 1] = (vertices[0].Y - vertices[1].Y) / dD;
        alphas[2, 2] = (vertices[1].X - vertices[0].X) / dD;
    }

    private void AssemblyLocalMatrixes()
    {
        double dD = Math.Abs(DeterminantD());
        CalcuclateAlphas();
        for (int i = 0; i < stiffnessMatrix.Size; i++)
            for (int j = 0; j <= i; j++)
            {
                double gradGradFunc(Point2D point)
                    => GetDpsiDx(point, i) * GetDpsiDx(point, j) + GetDpsiDy(point, i) * GetDpsiDy(point, j);
                stiffnessMatrix[i, j] = stiffnessMatrix[j, i] = GaussTriangle(gradGradFunc);
            }
        stiffnessMatrix = dD * stiffnessMatrix;
    }

    private (double, double, double) getL(Point2D point)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.X + alphas[0, 2] * point.Y;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.X + alphas[1, 2] * point.Y;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.X + alphas[2, 2] * point.Y;

        return (l1, l2, l3);
    }

    private double GetPsi(Point2D point, int numPsi)
    {
        (var l1, var l2, var l3) = getL(point);

        switch (numPsi)
        {
            case 0:
                return l1 * (2 * l1 - 1);
            case 1:
                return l2 * (2 * l2 - 1);
            case 2:
                return l3 * (2 * l3 - 1);
            case 3:
                return 4 * l1 * l2;
            case 4:
                return 4 * l2 * l3;
            case 5:
                return 4 * l3 * l1;
            default:
                return 0;
        }
    }

    private double GetDpsiDx(Point2D point, int numPsi)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.X + alphas[0, 2] * point.Y;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.X + alphas[1, 2] * point.Y;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.X + alphas[2, 2] * point.Y;

        switch (numPsi)
        {
            case 0:
                return 4 * alphas[0, 1] * l1 - alphas[0, 1];
            case 1:
                return 4 * alphas[1, 1] * l2 - alphas[1, 1];
            case 2:
                return 4 * alphas[2, 1] * l3 - alphas[2, 1];
            case 3:
                return 4 * (alphas[0, 1] * l2 + alphas[1, 1] * l1);
            case 4:
                return 4 * (alphas[1, 1] * l3 + alphas[2, 1] * l2);
            case 5:
                return 4 * (alphas[2, 1] * l1 + alphas[0, 1] * l3);
            default:
                return 0;
        }
    }

    private double GetDpsiDy(Point2D point, int numPsi)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.X + alphas[0, 2] * point.Y;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.X + alphas[1, 2] * point.Y;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.X + alphas[2, 2] * point.Y;

        switch (numPsi)
        {
            case 0:
                return 4 * alphas[0, 2] * l1 - alphas[0, 2];
            case 1:
                return 4 * alphas[1, 2] * l2 - alphas[1, 2];
            case 2:
                return 4 * alphas[2, 2] * l3 - alphas[2, 2];
            case 3:
                return 4 * (alphas[0, 2] * l2 + alphas[1, 2] * l1);
            case 4:
                return 4 * (alphas[1, 2] * l3 + alphas[2, 2] * l2);
            case 5:
                return 4 * (alphas[2, 2] * l1 + alphas[0, 2] * l3);
            default:
                return 0;
        }
    }

    private int FindElement(Point2D point)
    {
        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            double s01 = Math.Abs((vertices[1].X - vertices[0].X) * (point.Y - vertices[0].Y) -
                     (point.X - vertices[0].X) * (vertices[1].Y - vertices[0].Y));

            double s12 = Math.Abs((vertices[2].X - vertices[1].X) * (point.Y - vertices[1].Y) -
                     (point.X - vertices[1].X) * (vertices[2].Y - vertices[1].Y));

            double s20 = Math.Abs((vertices[0].X - vertices[2].X) * (point.Y - vertices[2].Y) -
                     (point.X - vertices[2].X) * (vertices[0].Y - vertices[2].Y));

            double dD = Math.Abs(DeterminantD());

            if (Math.Abs(dD - (s01 + s12 + s20)) < 1e-10)
                return ielem;
        }
        return -1;
    }

    public double XValueAtPoint(Point2D point)
    {
        double res = 0;

        int ielem = FindElement(point);

        if (ielem != -1)
        {
            CalcuclateAlphas();
            for (int i = 0; i < 6; i++)
            {
                res += solutionX[grid.Elements[ielem][i]] * GetPsi(point, i);
            }
        }
        return res;
    }

    public double YValueAtPoint(Point2D point)
    {
        double res = 0;

        int ielem = FindElement(point);

        if (ielem != -1)
        {
            CalcuclateAlphas();
            for (int i = 0; i < 6; i++)
            {
                res += solutionY[grid.Elements[ielem][i]] * GetPsi(point, i);
            }
        }
        return res;
    }
}
