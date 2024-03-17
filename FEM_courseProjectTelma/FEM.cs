using System.Numerics;

namespace CourseProject;

public class FEM
{
    private Grid grid;
    private SparseMatrix globalMatrix;
    private Vector globalVector;
    private Vector solution;
    private Solver slae;
    private Matrix alphas;
    private Vector localVector;
    private Matrix stiffnessMatrix;
    private Matrix massMatrix;
    private PointRZ[] vertices;
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

        vertices = new PointRZ[3];
        solution = new Vector(grid.Nodes.Count);
    }

    public void Compute()
    {

        BuildPortrait();


        AssemblyGlobalMatrix();
        //AccountSecondConditions();
        AccountFirstConditions();

        slae.SetSLAE(globalVector, globalMatrix);
        slae.CGM();
        Vector.Copy(slae.solution, solution);

    }

    public void AccountFirstConditions()
    {
        foreach (var fc in grid.Boundary)
        {
            globalMatrix.Di[fc.NodeNumber] = 1;
            globalVector[fc.NodeNumber] = 0;
            for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
            {
                globalVector[globalMatrix.Jg[i]] -= 0;
                globalMatrix.Gg[i] = 0;
            }
            for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
            {
                for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
                {
                    if (globalMatrix.Jg[j] == fc.NodeNumber)
                    {
                        globalVector[i] -= 0;
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
        globalVector = new(grid.Nodes.Count);

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
        globalVector.Fill(0);

        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            AssemblyLocalMatrixes();

            stiffnessMatrix = 1 / PhysicsConstants.VacuumPermeability * (stiffnessMatrix + massMatrix);

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
            double psi(PointRZ point)
                => GetPsi(point, i);
            localVector[i] = GaussTriangle(psi);
        }

        localVector = Math.Abs(DeterminantD()) * localVector;
    }

    private void AddElementToVector(int ielem)
    {
        double J = grid.Elements[ielem].Current;
        for (int i = 0; i < 6; i++)
        {
            globalVector[grid.Elements[ielem][i]] += J * localVector[i];
        }
    }

    private double GaussTriangle(Func<PointRZ, double> func)
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
            PointRZ point = new();
            point.R = (1 - p1[i] - p2[i]) * vertices[0].R + p1[i] * vertices[1].R + p2[i] * vertices[2].R;
            point.Z = (1 - p1[i] - p2[i]) * vertices[0].Z + p1[i] * vertices[1].Z + p2[i] * vertices[2].Z;

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
         => (vertices[1].R - vertices[0].R) * (vertices[2].Z - vertices[0].Z) -
            (vertices[2].R - vertices[0].R) * (vertices[1].Z - vertices[0].Z);

    private void CalcuclateAlphas()
    {
        double dD = DeterminantD();

        alphas[0, 0] = (vertices[1].R * vertices[2].Z - vertices[2].R * vertices[1].Z) / dD;
        alphas[0, 1] = (vertices[1].Z - vertices[2].Z) / dD;
        alphas[0, 2] = (vertices[2].R - vertices[1].R) / dD;

        alphas[1, 0] = (vertices[2].R * vertices[0].Z - vertices[0].R * vertices[2].Z) / dD;
        alphas[1, 1] = (vertices[2].Z - vertices[0].Z) / dD;
        alphas[1, 2] = (vertices[0].R - vertices[2].R) / dD;

        alphas[2, 0] = (vertices[0].R * vertices[1].Z - vertices[1].R * vertices[0].Z) / dD;
        alphas[2, 1] = (vertices[0].Z - vertices[1].Z) / dD;
        alphas[2, 2] = (vertices[1].R - vertices[0].R) / dD;
    }

    private void AssemblyLocalMatrixes()
    {
        double dD = Math.Abs(DeterminantD());
        CalcuclateAlphas();
        for (int i = 0; i < stiffnessMatrix.Size; i++)
            for (int j = 0; j <= i; j++)
            {
                double gradGradFunc(PointRZ point)
                    => (GetDpsiDx(point, i) * GetDpsiDx(point, j) + GetDpsiDy(point, i) * GetDpsiDy(point, j)) * point.R;
                stiffnessMatrix[i, j] = stiffnessMatrix[j, i] = GaussTriangle(gradGradFunc);

                double massFunc(PointRZ point)
                    => GetPsi(point, i) * GetPsi(point, j) / point.R;
            }
        stiffnessMatrix = dD * stiffnessMatrix;
    }

    private (double, double, double) getL(PointRZ point)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.R + alphas[0, 2] * point.Z;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.R + alphas[1, 2] * point.Z;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.R + alphas[2, 2] * point.Z;

        return (l1, l2, l3);
    }

    private double GetPsi(PointRZ point, int numPsi)
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

    private double GetDpsiDx(PointRZ point, int numPsi)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.R + alphas[0, 2] * point.Z;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.R + alphas[1, 2] * point.Z;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.R + alphas[2, 2] * point.Z;

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

    private double GetDpsiDy(PointRZ point, int numPsi)
    {
        double l1 = alphas[0, 0] + alphas[0, 1] * point.R + alphas[0, 2] * point.Z;
        double l2 = alphas[1, 0] + alphas[1, 1] * point.R + alphas[1, 2] * point.Z;
        double l3 = alphas[2, 0] + alphas[2, 1] * point.R + alphas[2, 2] * point.Z;

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

    private int FindElement(PointRZ point)
    {
        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
            vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
            vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

            double s01 = Math.Abs((vertices[1].R - vertices[0].R) * (point.Z - vertices[0].Z) -
                     (point.R - vertices[0].R) * (vertices[1].Z - vertices[0].Z));

            double s12 = Math.Abs((vertices[2].R - vertices[1].R) * (point.Z - vertices[1].Z) -
                     (point.R - vertices[1].R) * (vertices[2].Z - vertices[1].Z));

            double s20 = Math.Abs((vertices[0].R - vertices[2].R) * (point.Z - vertices[2].Z) -
                     (point.R - vertices[2].R) * (vertices[0].Z - vertices[2].Z));

            double dD = Math.Abs(DeterminantD());

            if (Math.Abs(dD - (s01 + s12 + s20)) < 1e-10)
                return ielem;
        }
        return -1;
    }

    public double AphiAtPoint(PointRZ point)
    {
        double res = 0;

        int ielem = FindElement(point);

        if (ielem != -1)
        {
            CalcuclateAlphas();
            for (int i = 0; i < 6; i++)
            {
                res += solution[grid.Elements[ielem][i]] * GetPsi(point, i);
            }
        }
        return res;
    }
    public double BzAtPoint(PointRZ point)
    {
        double hx = 1e-7;

        return ((point + new PointRZ(hx, 0)).R * AphiAtPoint(point + new PointRZ(hx, 0)) - (point - new PointRZ(hx, 0)).R * AphiAtPoint(point - new PointRZ(hx, 0))) / (2.0 * hx * point.R);
    }

    public double BrAtPoint(PointRZ point)
    {
        double hy = 1e-7;

        return -(AphiAtPoint(point + new PointRZ(0, hy)) - AphiAtPoint(point - new PointRZ(0, hy))) / (2.0 * hy);
    }

    public double AbsBAtPoint(PointRZ point)
    => Math.Sqrt(Math.Pow(BrAtPoint(point), 2) + Math.Pow(BzAtPoint(point), 2));

    public double AbsHAtPoint(PointRZ point)
        => AbsBAtPoint(point) * PhysicsConstants.VacuumPermeability;
}
