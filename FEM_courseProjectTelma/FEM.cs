namespace CourseProject;

public class FEM
{
   private Grid grid;
   private SparseMatrix globalMatrix;
   private Vector globalVector;
   private Vector[] layers;
   private Solver slae;
   private Matrix alphas;
   private Vector localVector;
   private Matrix stiffnessMatrix;
   private Matrix massMatrix;
   private PointRZ[] vertices;
   private FirstCondition[] firstConditions;
   private SecondCondition[] secondConditions;
   private Vector halfTimeSolution;
   public ITest test;
   private bool isPhysical;
   private string condPath => "conditions.txt";

   public FEM(Grid grid, ITest test, bool isPhysical)
   {
      this.grid = grid;
      this.isPhysical = isPhysical;
      alphas = new(3);
      stiffnessMatrix = new(6);
      massMatrix = new(6);
      localVector = new(6);
      slae = new Solver(1000, 1e-16);
      this.test = test;

      vertices = new PointRZ[3];
      layers = new Vector[2].Select(_ => new Vector(grid.Nodes.Count)).ToArray();
      halfTimeSolution = new Vector(grid.Nodes.Count);

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
      try
      {
         ArgumentNullException.ThrowIfNull(test, "Set the test.");
      }
      catch(Exception ex)
      {
         Console.WriteLine(ex.Message);
         throw;
      }

      BuildPortrait();
      PrepareStartConditions();

      for (int itime = 2; itime < grid.Time.Length; itime++)
      {
         AssemblyGlobalMatrix(itime);
         AccountSecondConditions(itime);
         AccountFirstConditions(itime);
         
         slae.SetSLAE(globalVector, globalMatrix);
         slae.CGM();

         Vector.Copy(layers[1], layers[0]);
         Vector.Copy(slae.solution, layers[1]);

         PrintLayerResult(itime);
      }
   }

   public void PrepareStartConditions()
   {
      // первый слой
      for (int i = 0; i < layers[0].Length; i++)
      {
         layers[0][i] = test.U(grid.Nodes[i], grid.Time[0]);
            
      }

      // второй слой
      if (isPhysical)
      {
         AssemblyGlobalMatrix(1);
         AccountSecondConditions(1);
         AccountFirstConditions(1);

         slae.SetSLAE(globalVector, globalMatrix);
         slae.CGM();
         Vector.Copy(slae.solution, layers[1]);

      }

      else
      {
         for (int i = 0; i < layers[1].Length; i++)
         {
            layers[1][i] = test.U(grid.Nodes[i], grid.Time[1]);
         }

         slae.solution = layers[1];
      }

      PrintLayerResult(1);
   }

   public void AccountFirstConditions(int itime)
   {
      foreach (var fc in firstConditions)
      {
         globalMatrix.Di[fc.NodeNumber] = 1;
         globalVector[fc.NodeNumber] = test.U(fc.point, grid.Time[itime]);
         for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
         {
            globalVector[globalMatrix.Jg[i]] -= test.U(fc.point, grid.Time[itime]) * globalMatrix.Gg[i];
            globalMatrix.Gg[i] = 0;
         }
         for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
         {
            for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
            {
               if (globalMatrix.Jg[j] == fc.NodeNumber)
               {
                  globalVector[i] -= test.U(fc.point, grid.Time[itime]) * globalMatrix.Gg[j];
                  globalMatrix.Gg[j] = 0;
               }
            }
         }
      }
   }

   public void AccountSecondConditions(int itime)
   {
      foreach (var sc in secondConditions)
      {
         localVector.Fill(0);
         int ielem = sc.ElemNumber;

         vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
         vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
         vertices[2] = grid.Nodes[grid.Elements[ielem][2]];
         CalcuclateAlphas();

         PointRZ start = new(grid.Nodes[sc.Edge[0]].R, grid.Nodes[sc.Edge[0]].Z);
         PointRZ end = new(grid.Nodes[sc.Edge[2]].R, grid.Nodes[sc.Edge[2]].Z);
         double lambda = grid.Materials[grid.Elements[ielem][^1]].Lambda;
         int localNum;

         for (int i = 0; i < 3; i++)
         {
            localNum = GetLocalNumber(ielem, sc.Edge[i]);
            localVector[localNum] = lambda * GaussEdge(localNum, test.Theta, sc.EdgeType, grid.Time[itime], start, end);
         }

         AddElementToVector(ielem);
      }
   }

   private int GetLocalNumber(int ielem, int globalNumber)
   {
      for (int i = 0; i < 6; i++)
      {
         if (grid.Elements[ielem][i] == globalNumber) return i;
      }

      throw new Exception($"У элемента {ielem} нет узла с номером {globalNumber}.");
   }

   public void BuildPortrait()
   {
      var list = new HashSet<int>[grid.Nodes.Count].Select(_ => new HashSet<int>()).ToList();
      foreach (var element in grid.Elements)
      {
         foreach (var pos in element)
         {
            foreach (var node in element)
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

   private void AssemblyGlobalMatrix(int itime)
   {
      double t01 = grid.Time[itime] - grid.Time[itime - 1];
      double t02 = itime > 1 ? grid.Time[itime] - grid.Time[itime - 2] : 1;
      double t12 = itime > 1 ? grid.Time[itime - 1] - grid.Time[itime - 2] : 1;

      globalMatrix.Clear();
      globalVector.Fill(0);

      for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
      {
         vertices[0] = grid.Nodes[grid.Elements[ielem][0]];
         vertices[1] = grid.Nodes[grid.Elements[ielem][1]];
         vertices[2] = grid.Nodes[grid.Elements[ielem][2]];

         AssemblyLocalMatrixes();

         massMatrix = grid.Materials[grid.Elements[ielem][^1]].Sigma * massMatrix;
         stiffnessMatrix = grid.Materials[grid.Elements[ielem][^1]].Lambda * stiffnessMatrix + (1.0 / t01 + (itime > 1 ? 1.0 / t02 : 0)) * massMatrix;

         for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
               AddElement(grid.Elements[ielem][i], grid.Elements[ielem][j], stiffnessMatrix[i, j]);

         AssemblyLocallVector(ielem, itime, t01, t02, t12);
         
         AddElementToVector(ielem);

         stiffnessMatrix.Clear();
         massMatrix.Clear();
         localVector.Fill(0);
      }
   }

   void AssemblyLocallVector(int ielem, int itime, double t01, double t02, double t12)
   {
      double lambda = grid.Materials[grid.Elements[ielem][^1]].Lambda;
      double sigma = grid.Materials[grid.Elements[ielem][^1]].Sigma;

      for (int i = 0; i < 6; i++)
      {
         localVector[i] = GaussTriangle(i, grid.Time[itime], lambda, sigma);
      }

      localVector = Math.Abs(DeterminantD()) * localVector;

      Vector qj1 = new(6);
      Vector qj2 = new(6);
      if (itime == 1)
      {
         for (int i = 0; i < qj1.Length; i++)
         {
            for (int j = 0; j < qj1.Length; j++)
               qj1[i] += massMatrix[i, j] * layers[0][grid.Elements[ielem][j]];
         }

         localVector += 1.0 / t01 * qj1;
      }
      else
      {
         for (int i = 0; i < qj1.Length; i++)
         {
            for (int j = 0; j < qj1.Length; j++)
            {
               qj2[i] += massMatrix[i, j] * layers[0][grid.Elements[ielem][j]];

               qj1[i] += massMatrix[i, j] * layers[1][grid.Elements[ielem][j]];
            }
         }

         localVector += t02 / (t12 * t01) * qj1 - t01 / (t02 * t12) * qj2;
      }
   }

   private void AddElementToVector(int ielem)
   {
      for (int i = 0; i < 6; i++)
      {
         globalVector[grid.Elements[ielem][i]] += localVector[i];
      }
   }

   private double GaussTriangle(int numberPsi, double t, double lambda, double sigma)
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

         res += test.F(point, t, sigma, lambda) * basis(point, numberPsi) * w[i] * 0.5 * point.R;
      }

      return res;
   }

   private double GaussEdge(int numberPsi, Func<PointRZ, double,int, double> theta, int funcNum, double t, PointRZ fstPoint, PointRZ sndPoint)
   {
      double[] p = { 0.0,
                     1.0 / 3.0 * Math.Sqrt(5 - 2 * Math.Sqrt(10.0 / 7.0)),
                     -1.0 / 3.0 * Math.Sqrt(5 - 2 * Math.Sqrt(10.0 / 7.0)),
                     1.0 / 3.0 * Math.Sqrt(5 + 2 * Math.Sqrt(10.0 / 7.0)),
                     -1.0 / 3.0 * Math.Sqrt(5 + 2 * Math.Sqrt(10.0 / 7.0))};

      double[] w = { 128.0 / 225.0,
                            (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
                            (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
                            (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0,
                            (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0 };
      double res = 0;

      double lengthEdge = Math.Sqrt((fstPoint.R - sndPoint.R) * (fstPoint.R - sndPoint.R) +
                            (fstPoint.Z - sndPoint.Z) * (fstPoint.Z - sndPoint.Z));

      for (int i = 0; i < w.Length; i++)
      {
         PointRZ point = new();
         point.R = (sndPoint.R - fstPoint.R) * (1 + p[i]) / 2.0 + fstPoint.R;
         point.Z = (sndPoint.Z - fstPoint.Z) * (1 + p[i]) / 2.0 + fstPoint.Z;                                           

         res += theta(point, t, funcNum) * basis(point, numberPsi) * w[i] * point.R;
      }

      return lengthEdge * res / 2.0;
   }

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

      stiffnessMatrix[0, 0] = dD *  (alphas[0, 1] * alphas[0, 1] + alphas[0, 2] * alphas[0, 2])
         * (3 * vertices[0].R + vertices[1].R + vertices[2].R) / 10;


      stiffnessMatrix[1, 0] = -dD *  (alphas[0, 1] * alphas[1, 1] + alphas[0, 2] * alphas[1, 2])
         * (2 * vertices[0].R + 2 * vertices[1].R + vertices[2].R) / 30;
      stiffnessMatrix[0, 1] = stiffnessMatrix[1, 0];


      stiffnessMatrix[1, 1] = dD *  (alphas[1, 1] * alphas[1, 1] + alphas[1, 2] * alphas[1, 2])
               * (vertices[0].R + 3 * vertices[1].R + vertices[2].R) / 10;


      stiffnessMatrix[2, 0] = -dD *  (alphas[0, 1] * alphas[2, 1] + alphas[0, 2] * alphas[2, 2])
         * (2 * vertices[0].R + vertices[1].R + 2 * vertices[2].R) / 30;
      stiffnessMatrix[0, 2] = stiffnessMatrix[2, 0];


      stiffnessMatrix[2, 1] = -dD *  (alphas[1, 1] * alphas[2, 1] + alphas[1, 2] * alphas[2, 2])
               * (vertices[0].R + 2 * vertices[1].R + 2 * vertices[2].R) / 30;
      stiffnessMatrix[1, 2] = stiffnessMatrix[2, 1];


      stiffnessMatrix[2, 2] = dD *  (alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2])
         * (vertices[0].R + vertices[1].R + 3 * vertices[2].R) / 10;


      stiffnessMatrix[3, 0] = dD *
         ((3 * alphas[0, 1] * alphas[0, 1] + 14 * alphas[0, 1] * alphas[1, 1] +
         3 * alphas[0, 2] * alphas[0, 2] + 14 * alphas[0, 2] * alphas[1, 2]) * vertices[0].R +
         (-2 * alphas[0, 1] * alphas[0, 1] + 3 * alphas[0, 1] * alphas[1, 1] -
         2 * alphas[0, 2] * alphas[0, 2] + 3 * alphas[0, 2] * alphas[1, 2]) * vertices[1].R +
         (-alphas[0, 1] * alphas[0, 1] + 3 * alphas[0, 1] * alphas[1, 1] -
         alphas[0, 2] * alphas[0, 2] + 3 * alphas[0, 2] * alphas[1, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[0, 3] = stiffnessMatrix[3, 0];


      stiffnessMatrix[3, 1] = dD *
         ((3 * alphas[0, 1] * alphas[1, 1] + 3 * alphas[0, 2] * alphas[1, 2] -
          2 * alphas[1, 1] * alphas[1, 1] - 2 * alphas[1, 2] * alphas[1, 2]) * vertices[0].R +
         (14 * alphas[0, 1] * alphas[1, 1] + 14 * alphas[0, 2] * alphas[1, 2] +
         3 * alphas[1, 1] * alphas[1, 1] + 3 * alphas[1, 2] * alphas[1, 2]) * vertices[1].R +
         (3 * alphas[0, 1] * alphas[1, 1] + 3 * alphas[0, 2] * alphas[1, 2] -
         alphas[1, 1] * alphas[1, 1] - alphas[1, 2] * alphas[1, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[1, 3] = stiffnessMatrix[3, 1];


      stiffnessMatrix[3, 2] = -dD *
         ((alphas[0, 1] * alphas[2, 1] + alphas[0, 2] * alphas[2, 2] +
         2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[0].R +
         (2 * alphas[0, 1] * alphas[2, 1] + 2 * alphas[0, 2] * alphas[2, 2] +
         alphas[1, 1] * alphas[2, 1] + alphas[1, 2] * alphas[2, 2]) * vertices[1].R +
         (-3 * alphas[0, 1] * alphas[2, 1] - 3 * alphas[0, 2] * alphas[2, 2] -
         3 * alphas[1, 1] * alphas[2, 1] - 3 * alphas[1, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[2, 3] = stiffnessMatrix[3, 2];


      stiffnessMatrix[3, 3] = dD *  4 *
               ((alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[1, 1] +
               alphas[0, 2] * alphas[0, 2] + 2 * alphas[0, 2] * alphas[1, 2] + 3 *
               alphas[1, 1] * alphas[1, 1] + 3 * alphas[1, 2] * alphas[1, 2]) * vertices[0].R +
               (3 * alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[1, 1] +
               3 * alphas[0, 2] * alphas[0, 2] + 2 * alphas[0, 2] * alphas[1, 2] +
               alphas[1, 1] * alphas[1, 1] + alphas[1, 2] * alphas[1, 2]) * vertices[1].R +
               (alphas[0, 1] * alphas[0, 1] + alphas[0, 1] * alphas[1, 1] +
               alphas[0, 2] * alphas[0, 2] + alphas[0, 2] * alphas[1, 2] +
               alphas[1, 1] * alphas[1, 1] + alphas[1, 2] * alphas[1, 2]) * vertices[2].R) / 15;


      stiffnessMatrix[4, 0] = dD *
               ((3 * alphas[0, 1] * alphas[1, 1] + 3 * alphas[0, 1] * alphas[2, 1] +
               3 * alphas[0, 2] * alphas[1, 2] + 3 * alphas[0, 2] * alphas[2, 2]) * vertices[0].R +
               (-alphas[0, 1] * alphas[1, 1] - 2 * alphas[0, 1] * alphas[2, 1] -
               alphas[0, 2] * alphas[1, 2] - 2 * alphas[0, 2] * alphas[2, 2]) * vertices[1].R +
               (-2 * alphas[0, 1] * alphas[1, 1] - alphas[0, 1] * alphas[2, 1] -
               2 * alphas[0, 2] * alphas[1, 2] - alphas[0, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[0, 4] = stiffnessMatrix[4, 0];


      stiffnessMatrix[4, 1] = -dD *
               ((alphas[1, 1] * alphas[1, 1] - 3 * alphas[1, 1] * alphas[2, 1] +
               alphas[1, 2] * alphas[1, 2] - 3 * alphas[1, 2] * alphas[2, 2]) * vertices[0].R +
               (-3 * alphas[1, 1] * alphas[1, 1] - 14 * alphas[1, 1] * alphas[2, 1] -
               3 * alphas[1, 2] * alphas[1, 2] - 14 * alphas[1, 2] * alphas[2, 2]) * vertices[1].R +
               (2 * alphas[1, 1] * alphas[1, 1] - 3 * alphas[1, 1] * alphas[2, 1] +
               2 * alphas[1, 2] * alphas[1, 2] - 3 * alphas[1, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[1, 4] = stiffnessMatrix[4, 1];


      stiffnessMatrix[4, 2] = dD *
               ((3 * alphas[1, 1] * alphas[2, 1] + 3 * alphas[1, 2] * alphas[2, 2] -
               alphas[2, 1] * alphas[2, 1] - alphas[2, 2] * alphas[2, 2]) * vertices[0].R +
               (3 * alphas[1, 1] * alphas[2, 1] + 3 * alphas[1, 2] * alphas[2, 2] -
               2 * alphas[2, 1] * alphas[2, 1] - 2 * alphas[2, 2] * alphas[2, 2]) * vertices[1].R +
               (14 * alphas[1, 1] * alphas[2, 1] + 14 * alphas[1, 2] * alphas[2, 2] +
               3 * alphas[2, 1] * alphas[2, 1] + 3 * alphas[2, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[2, 4] = stiffnessMatrix[4, 2];


      stiffnessMatrix[4, 3] = dD *  2 *
               ((alphas[0, 1] * alphas[1, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[1, 1] + 2 * alphas[1, 1] * alphas[2, 1]
               + 2 * alphas[1, 2] * alphas[1, 2] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[0].R +
               (2 * alphas[0, 1] * alphas[1, 1] + 6 * alphas[0, 1] * alphas[2, 1] +
               2 * alphas[0, 2] * alphas[1, 2] + 6 * alphas[0, 2] * alphas[2, 2] +
               alphas[1, 1] * alphas[1, 1] + 2 * alphas[1, 1] * alphas[2, 1]
               + alphas[1, 2] * alphas[1, 2] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[1].R +
               (2 * alphas[0, 1] * alphas[1, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               2 * alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[1, 1] + alphas[1, 1] * alphas[2, 1] +
               2 * alphas[1, 2] * alphas[1, 2] + alphas[1, 2] * alphas[2, 2]) * vertices[2].R) / 15;
      stiffnessMatrix[3, 4] = stiffnessMatrix[4, 3];


      stiffnessMatrix[4, 4] = dD *  4 *
               ((alphas[1, 1] * alphas[1, 1] + alphas[1, 1] * alphas[2, 1] +
               alphas[1, 2] * alphas[1, 2] + alphas[1, 2] * alphas[2, 2] +
               alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2]) * vertices[0].R +
               (alphas[1, 1] * alphas[1, 1] + 2 * alphas[1, 1] * alphas[2, 1] +
               alphas[1, 2] * alphas[1, 2] + 2 * alphas[1, 2] * alphas[2, 2] +
               3 * alphas[2, 1] * alphas[2, 1] + 3 * alphas[2, 2] * alphas[2, 2]) * vertices[1].R +
               (3 * alphas[1, 1] * alphas[1, 1] + 2 * alphas[1, 1] * alphas[2, 1] +
               3 * alphas[1, 2] * alphas[1, 2] + 2 * alphas[1, 2] * alphas[2, 2] +
               alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2]) * vertices[2].R) / 15;


      stiffnessMatrix[5, 0] = dD *
               ((3 * alphas[0, 1] * alphas[0, 1] + 14 * alphas[0, 1] * alphas[2, 1] +
               3 * alphas[0, 2] * alphas[0, 2] + 14 * alphas[0, 2] * alphas[2, 2]) * vertices[0].R +
               (-alphas[0, 1] * alphas[0, 1] + 3 * alphas[0, 1] * alphas[2, 1] -
               alphas[0, 2] * alphas[0, 2] + 3 * alphas[0, 2] * alphas[2, 2]) * vertices[1].R +
               (-2 * alphas[0, 1] * alphas[0, 1] + 3 * alphas[0, 1] * alphas[2, 1] -
               2 * alphas[0, 2] * alphas[0, 2] + 3 * alphas[0, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[0, 5] = stiffnessMatrix[5, 0];


      stiffnessMatrix[5, 1] = -dD *
               ((alphas[0, 1] * alphas[1, 1] + alphas[0, 2] * alphas[1, 2] +
               2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[0].R +
               (-3 * alphas[0, 1] * alphas[1, 1] - 3 * alphas[0, 2] * alphas[1, 2] -
               3 * alphas[1, 1] * alphas[2, 1] - 3 * alphas[1, 2] * alphas[2, 2]) * vertices[1].R +
               (2 * alphas[0, 1] * alphas[1, 1] + 2 * alphas[0, 2] * alphas[1, 2] +
               alphas[1, 1] * alphas[2, 1] + alphas[1, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[1, 5] = stiffnessMatrix[5, 1];


      stiffnessMatrix[5, 2] = dD *
         ((3 * alphas[0, 1] * alphas[2, 1] + 3 * alphas[0, 2] * alphas[2, 2] -
         2 * alphas[2, 1] * alphas[2, 1] - 2 * alphas[2, 2] * alphas[2, 2]) * vertices[0].R +
         (3 * alphas[0, 1] * alphas[2, 1] + 3 * alphas[0, 2] * alphas[2, 2] -
         alphas[2, 1] * alphas[2, 1] - alphas[2, 2] * alphas[2, 2]) * vertices[1].R +
         (14 * alphas[0, 1] * alphas[2, 1] + 14 * alphas[0, 2] * alphas[2, 2] +
         3 * alphas[2, 1] * alphas[2, 1] + 3 * alphas[2, 2] * alphas[2, 2]) * vertices[2].R) / 30;
      stiffnessMatrix[2, 5] = stiffnessMatrix[5, 2];


      stiffnessMatrix[5, 3] = dD *  2 *
               ((alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[1, 1] +
               2 * alphas[0, 1] * alphas[2, 1] + alphas[0, 2] * alphas[0, 2] +
               2 * alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               6 * alphas[1, 1] * alphas[2, 1] + 6 * alphas[1, 2] * alphas[2, 2]) * vertices[0].R +
               (2 * alphas[0, 1] * alphas[0, 1] + alphas[0, 1] * alphas[1, 1] +
               2 * alphas[0, 1] * alphas[2, 1] + 2 * alphas[0, 2] * alphas[0, 2] +
               alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[1].R +
               (2 * alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[1, 1] +
               alphas[0, 1] * alphas[2, 1] + 2 * alphas[0, 2] * alphas[0, 2] +
               2 * alphas[0, 2] * alphas[1, 2] + alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2]) * vertices[2].R) / 15;
      stiffnessMatrix[3, 5] = stiffnessMatrix[5, 3];


      stiffnessMatrix[5, 4] = dD *  2 *
               ((2 * alphas[0, 1] * alphas[1, 1] + alphas[0, 1] * alphas[2, 1] +
               2 * alphas[0, 2] * alphas[1, 2] + alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2] +
               2 * alphas[2, 1] * alphas[2, 1] + 2 * alphas[2, 2] * alphas[2, 2]) * vertices[0].R +
               (2 * alphas[0, 1] * alphas[1, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               2 * alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               alphas[1, 1] * alphas[2, 1] + alphas[1, 2] * alphas[2, 2] +
               2 * alphas[2, 1] * alphas[2, 1] + 2 * alphas[2, 2] * alphas[2, 2]) * vertices[1].R +
               (6 * alphas[0, 1] * alphas[1, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               6 * alphas[0, 2] * alphas[1, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               2 * alphas[1, 1] * alphas[2, 1] + 2 * alphas[1, 2] * alphas[2, 2] +
               alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2]) * vertices[2].R) / 15;
      stiffnessMatrix[4, 5] = stiffnessMatrix[5, 4];

      stiffnessMatrix[5, 5] = dD *  4 *
               ((alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               alphas[0, 2] * alphas[0, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               3 * alphas[2, 1] * alphas[2, 1] + 3 * alphas[2, 2] * alphas[2, 2]) * vertices[0].R +
               (alphas[0, 1] * alphas[0, 1] + alphas[0, 1] * alphas[2, 1] +
               alphas[0, 2] * alphas[0, 2] + alphas[0, 2] * alphas[2, 2] +
               alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2]) * vertices[1].R +
               (3 * alphas[0, 1] * alphas[0, 1] + 2 * alphas[0, 1] * alphas[2, 1] +
               3 * alphas[0, 2] * alphas[0, 2] + 2 * alphas[0, 2] * alphas[2, 2] +
               alphas[2, 1] * alphas[2, 1] + alphas[2, 2] * alphas[2, 2]) * vertices[2].R) / 15;


      massMatrix[0, 0] = dD * (5 * vertices[0].R + vertices[1].R + vertices[2].R) / 420;


      massMatrix[1, 0] = -dD * (4 * vertices[0].R + 4 * vertices[1].R - vertices[2].R) / 2520;
      massMatrix[0, 1] = massMatrix[1, 0];


      massMatrix[1, 1] = dD * (vertices[0].R + 5 * vertices[1].R + vertices[2].R) / 420;


      massMatrix[2, 0] = -dD * (4 * vertices[0].R - vertices[1].R + 4 * vertices[2].R) / 2520;
      massMatrix[0, 2] = massMatrix[2, 0];


      massMatrix[2, 1] = dD * (vertices[0].R - 4 * vertices[1].R - 4 * vertices[2].R) / 2520;
      massMatrix[1, 2] = massMatrix[2, 1];

      massMatrix[2, 2] = dD * (vertices[0].R + vertices[1].R + 5 * vertices[2].R) / 420;


      massMatrix[3, 0] = dD * (3 * vertices[0].R - 2 * vertices[1].R - vertices[2].R) / 630;
      massMatrix[0, 3] = massMatrix[3, 0];


      massMatrix[3, 1] = -dD * (2 * vertices[0].R - 3 * vertices[1].R + vertices[2].R) / 630;
      massMatrix[1, 3] = massMatrix[3, 1];


      massMatrix[3, 2] = -dD * (3 * vertices[0].R + 3 * vertices[1].R + vertices[2].R) / 630;
      massMatrix[2, 3] = massMatrix[3, 2];


      massMatrix[3, 3] = dD * 4 * (3 * vertices[0].R + 3 * vertices[1].R + vertices[2].R) / 315;


      massMatrix[4, 0] = -dD * (vertices[0].R + 3 * vertices[1].R + 3 * vertices[2].R) / 630;
      massMatrix[0, 4] = massMatrix[4, 0];


      massMatrix[4, 1] = -dD * (vertices[0].R - 3 * vertices[1].R + 2 * vertices[2].R) / 630;
      massMatrix[1, 4] = massMatrix[4, 1];


      massMatrix[4, 2] = -dD * (vertices[0].R + 2 * vertices[1].R - 3 * vertices[2].R) / 630;
      massMatrix[2, 4] = massMatrix[4, 2];


      massMatrix[4, 3] = dD * 2 * (2 * vertices[0].R + 3 * vertices[1].R + 2 * vertices[2].R) / 315;
      massMatrix[3, 4] = massMatrix[4, 3];


      massMatrix[4, 4] = dD * 4 * (vertices[0].R + 3 * vertices[1].R + 3 * vertices[2].R) / 315;


      massMatrix[5, 0] = dD * (3 * vertices[0].R - vertices[1].R - 2 * vertices[2].R) / 630;
      massMatrix[0, 5] = massMatrix[5, 0];


      massMatrix[5, 1] = -dD * (3 * vertices[0].R + vertices[1].R + 3 * vertices[2].R) / 630;
      massMatrix[1, 5] = massMatrix[5, 1];


      massMatrix[5, 2] = -dD * (2 * vertices[0].R + vertices[1].R - 3 * vertices[2].R) / 630;
      massMatrix[2, 5] = massMatrix[5, 2];


      massMatrix[5, 3] = dD * 2 * (3 * vertices[0].R + 2 * vertices[1].R + 2 * vertices[2].R) / 315;
      massMatrix[3, 5] = massMatrix[5, 3];


      massMatrix[5, 4] = dD * 2 * (2 * vertices[0].R + 2 * vertices[1].R + 3 * vertices[2].R) / 315;
      massMatrix[4, 5] = massMatrix[5, 4];


      massMatrix[5, 5] = dD * 4 * (3 * vertices[0].R + vertices[1].R + 3 * vertices[2].R) / 315;
   }

   private double basis(PointRZ point, int numPsi)
   {
      double l1 = alphas[0, 0] + alphas[0, 1] * point.R + alphas[0, 2] * point.Z;
      double l2 = alphas[1, 0] + alphas[1, 1] * point.R + alphas[1, 2] * point.Z;
      double l3 = alphas[2, 0] + alphas[2, 1] * point.R + alphas[2, 2] * point.Z;
      
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
   public double ValueAtPoint(PointRZ point)
   {
      double res = 0;

      int ielem = FindElement(point);

      if (ielem != -1)
      {
         CalcuclateAlphas();
         for (int i = 0; i < 6; i++)
         {
            res += layers[1][grid.Elements[ielem][i]] * basis(point, i);
         }
      }
      return res;
   }

   public double ValueAtPointAndHalfTime(PointRZ point)
   {
      double res = 0;


      int ielem = FindElement(point);

      if (ielem != -1)
      {
         CalcuclateAlphas();
         for (int i = 0; i < 6; i++)
         {
            res += halfTimeSolution[grid.Elements[ielem][i]] * basis(point, i);
         }
      }
      return res;
   }

   private void PrintLayerResult(int itime)
   {
      Vector error;
      Vector exact = new(grid.Nodes.Count);
      double residual;
      bool appEnd = true;
      if (itime == 1) appEnd = false;

      using StreamWriter sw1 = new("results/solution" + test.GetType() + ".txt", appEnd),
         sw2 = new("results/layerRes" + test.GetType() + ".csv", appEnd),
         sw3 = new("results/currentRes.txt", appEnd);

      for (int i = 0; i < exact.Length; i++)
      {
         exact[i] = test.U(grid.Nodes[i], grid.Time[itime]);
      }

      error = exact - slae.solution;

      residual = Math.Sqrt(error * error / error.Length);

      sw3.Write($"{grid.Time[itime]} ");
      for (int i = 0; i < slae.solution.Length; i++)
      {
         sw1.Write($"{slae.solution[i]} ");
         sw3.Write($"{slae.solution[i]} ");
      }
      sw1.WriteLine();
      sw3.WriteLine();

      if (Math.Abs(grid.Time[itime] - 0.5) < 1e-12) Vector.Copy(slae.solution, halfTimeSolution);

      if (itime == 1)
         sw2.WriteLine("t; Погрешность");
      sw2.WriteLine($"{grid.Time[itime]:0.0000}; {residual:0.00E+0}");
      Console.WriteLine($"{grid.Time[itime]:0.0000}; {residual:0.00E+0}");
   }

   public void PrintSolution()
   {
      
      Console.WriteLine("Численное решение");
      for (int i = 0; i < layers[1].Length; i++)
      {
         Console.WriteLine($"{layers[1][i]}");
      }
      Vector exactSolution = new(grid.Nodes.Count);
      Console.WriteLine("Точное решение");
      for (int i = 0; i < exactSolution.Length; i++)
      {
         exactSolution[i] = test.U(grid.Nodes[i], grid.Time[^1]);
         Console.WriteLine($"{exactSolution[i]}");
      }
      Console.WriteLine("Погрешность");
      Vector inaccuracy = layers[1] - exactSolution;
      //int internalCount = 0;
      for (int i = 0; i < inaccuracy.Length; i++)
      {
         //if (Math.Abs(inaccuracy[i]) > 1e-10) internalCount++;
         Console.WriteLine($"{inaccuracy[i]}");
      }
      Console.WriteLine("Относительная погрешность");
      Console.WriteLine($"{Math.Sqrt(inaccuracy * inaccuracy / inaccuracy.Length)}");
      //Console.WriteLine($"{Math.Sqrt(inaccuracy * inaccuracy / internalCount)}");
   }

   public void PrintResearchSolution()
   {
      Console.WriteLine("t = 0.5");
      Console.WriteLine($"{halfTimeSolution[3]}\n{halfTimeSolution[16]}\n" +
         $"{halfTimeSolution[5]}\n{halfTimeSolution[8]}");
      Console.WriteLine($"---------\n{test.U(grid.Nodes[3], 0.5)}\n{test.U(grid.Nodes[16], 0.5)}\n" +
         $"{test.U(grid.Nodes[5], 0.5)}\n{test.U(grid.Nodes[8], 0.5)}");
      Console.WriteLine("t = 1");
      Console.WriteLine($"{layers[1][3]}\n{layers[1][16]}\n" +
         $"{layers[1][5]}\n{layers[1][8]}");
      Console.WriteLine($"---------\n{test.U(grid.Nodes[3], 1)}\n{test.U(grid.Nodes[16], 1)}\n" +
   $"{test.U(grid.Nodes[5], 1)}\n{test.U(grid.Nodes[8], 1)}");
   }
}
