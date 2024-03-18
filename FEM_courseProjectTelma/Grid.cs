using System.Collections;
using System.Drawing;
using System.Reflection;
using System.Runtime.CompilerServices;

namespace CourseProject;

public class FiniteElement
{
    public int[] Nodes { get; set; } = new int[6];

    public double Current { get; set; }

    public int this[int i]
    {
        get => Nodes[i]; set => Nodes[i] = value;
    }
}

public class Grid
{
   public List<PointRZ> Nodes { get; set; }
   public FiniteElement[] Elements { get; }
    public HashSet<FirstCondition> Boundary { get; init; }
    public double R { get => _rPoints is null ? 0.001 : _rPoints[1]; }
    public double r0 { get => _rPoints is null ? 0.001 : _rPoints[0]; }

    private double[] _rPoints;
    private int[] _rSteps;
    private double[] _rCoefs;
    private double[] _rCur;

    private double[] _zPoints;
    private int[] _zSteps;
    private double[] _zCoefs;
    private double[] _zCur;

    public Grid(string spaceGridPath)
    {
        List<int[]> boundaries = new List<int[]>();
        string[] data;
        int rSplits, zSplits;

        // чтение данных из файла сетки
        using (StreamReader sr = new(spaceGridPath))
        {
            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            rSplits = Convert.ToInt32(data[1]);

            _rPoints = new double[rSplits + 1];
            _rPoints[0] = Convert.ToDouble(data[0]);
            _rSteps = new int[rSplits];
            _rCoefs = new double[rSplits];
            _rCur = new double[rSplits];

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rPoints[i + 1] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rSteps[i] = Convert.ToInt32(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rCoefs[i] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rCur[i] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            zSplits = Convert.ToInt32(data[1]);

            _zPoints = new double[zSplits + 1];
            _zPoints[0] = Convert.ToDouble(data[0]);
            _zSteps = new int[zSplits];
            _zCoefs = new double[zSplits];
            _zCur = new double[zSplits];

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zPoints[i + 1] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zSteps[i] = Convert.ToInt32(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zCoefs[i] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _zCur[i] = Convert.ToDouble(data[i]);
            }
        }

        // Генерация сетки и глобальная нумерация узлов
        var rUniq = new double[_rSteps.Sum() + 1];
        var zUniq = new double[_zSteps.Sum() + 1];

        rUniq[0] = _rPoints[0];
        int ridx = 1;
        for (int i = 0; i < rSplits; i++)
        {
            double sumCoef = 0;
            for (int j = 0; j < _rSteps[i]; j++)
            {
                sumCoef += Math.Pow(_rCoefs[i], j);
            }

            double r_step = (_rPoints[i + 1] - _rPoints[i]) / sumCoef;

            for (int j = 0; j < _rSteps[i] - 1; j++)
            {
                rUniq[ridx] = rUniq[ridx - 1] + r_step;
                r_step *= _rCoefs[i];
                ridx++;
            }

            rUniq[ridx++] = _rPoints[i + 1];
        }

        zUniq[0] = _zPoints[0];
        int zidx = 1;
        for (int i = 0; i < zSplits; i++)
        {
            double sumCoef = 0;
            for (int j = 0; j < _zSteps[i]; j++)
            {
                sumCoef += Math.Pow(_zCoefs[i], j);
            }

            double z_step = (_zPoints[i + 1] - _zPoints[i]) / sumCoef;

            for (int j = 0; j < _zSteps[i] - 1; j++)
            {
                zUniq[zidx] = zUniq[zidx - 1] + z_step;
                z_step *= _zCoefs[i];
                zidx++;
            }

            zUniq[zidx++] = _zPoints[i + 1];
        }

        Nodes = new List<PointRZ>();

        for (int i = 0; i < zUniq.Length; i++)
            for (int j = 0; j < rUniq.Length; j++)
                Nodes.Add(new(rUniq[j], zUniq[i]));

        Console.WriteLine(Nodes.Count);
        // генерация элементов
        Elements = new FiniteElement[2 * (rUniq.Length - 1) * (zUniq.Length - 1)];

        int elidx = 0;
        zidx = 0;
        int kidx = 0;
        for (int i = 0; i < _zSteps.Length; i++)
        {
            for (int j = 0; j < _zSteps[i]; j++)
            {
                for (int k = 0; k < _rSteps.Length; k++)
                {
                    for (int m = 0; m < _rSteps[k]; m++)
                    {
                        Elements[elidx] = new();
                        Elements[elidx].Nodes[0] = (zidx + 1) * rUniq.Length + kidx;
                        Elements[elidx].Nodes[1] = (zidx + 1) * rUniq.Length + kidx + 1;
                        Elements[elidx].Nodes[2] = zidx * rUniq.Length + kidx;
                        Elements[elidx].Current = _rCur[k] * _zCur[i];
                        elidx++;
                        Elements[elidx] = new();
                        Elements[elidx].Nodes[0] = (zidx + 1) * rUniq.Length + kidx + 1;
                        Elements[elidx].Nodes[1] = zidx * rUniq.Length + kidx;
                        Elements[elidx].Nodes[2] = zidx * rUniq.Length + kidx + 1;
                        Elements[elidx].Current = _rCur[k] * _zCur[i];
                        elidx++;
                        kidx++;
                    }
                }
                kidx = 0;
                zidx++;
            }
        }

        // указываем границы с краевыми условиями
        Boundary = new HashSet<FirstCondition>();

        // нижняя граница
        for (int i = 1; i < rUniq.Length; i++)
            boundaries.Add([Nodes.Count - i, 0, Nodes.Count - (i + 1)]);

        // боковые границы
        for (int i = 0; i < zUniq.Length - 1; i++)
        {
            //boundaries.Add([rUniq.Length * i, 0, rUniq.Length * (i + 1)]);  // левая
            boundaries.Add([rUniq.Length * (i + 1) - 1, 0, rUniq.Length * (i + 2) - 1]);  // правая
        }

        // верхняя граница
        for (int i = 0; i < rUniq.Length - 1; i++)
            boundaries.Add([i, 0, i + 1]);


        NumberNodes(boundaries);
    }

   private void NumberNodes(List<int[]> boudnaries)
   {
      Dictionary<(int,int), (int,int)> hashtable = new Dictionary<(int, int), (int, int)>();
      (int, int) pointInfo;
      int num = Nodes.Count;
      (int, int)[] key = new (int,int)[3];

      for (int ielem = 0;  ielem < Elements.Length; ielem++)
      {
         key[0] = (Math.Min(Elements[ielem][0], Elements[ielem][1]), Math.Max(Elements[ielem][0], Elements[ielem][1]));
         key[1] = (Math.Min(Elements[ielem][1], Elements[ielem][2]), Math.Max(Elements[ielem][1], Elements[ielem][2]));
         key[2] = (Math.Min(Elements[ielem][0], Elements[ielem][2]), Math.Max(Elements[ielem][0], Elements[ielem][2]));

         for (int i = 0; i < 3; i++)
         {
            if (hashtable.ContainsKey(key[i]))
            {
               pointInfo = hashtable[key[i]];
               Elements[ielem][i + 3] = pointInfo.Item1;
            }
            else
            {
               Elements[ielem][i + 3] = num;
               pointInfo = (num, ielem);
               hashtable.Add(key[i], pointInfo);
               Nodes.Add((Nodes[key[i].Item1] + Nodes[key[i].Item2]) / 2);
               num++;
            }
         }
      }

      foreach (var edge in boudnaries)
      {
         key[0] = (Math.Min(edge[0], edge[2]), Math.Max(edge[0], edge[2]));

         if (hashtable.ContainsKey(key[0]))
         {
            pointInfo = hashtable[key[0]];
            edge[1] = pointInfo.Item1;
                foreach (var node in edge)
                    Boundary.Add(new(Nodes[node], node));
         }
         else
         {
            throw new Exception("Такая грань не была описана в элементах");
         }
      }
   }
}

public static class PhysicsConstants
{
    public const double VacuumPermeability = 4.0 * Math.PI * 1E-07;
}
