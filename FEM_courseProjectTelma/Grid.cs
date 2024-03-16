﻿using System.Collections;
using System.Drawing;

namespace CourseProject;

public struct Material
{
   public double Lambda { get; init; }
   public double Sigma { get; init; }
   public string Name { get; init; }

   public Material(double lambda, double sigma, string name)
   {
      Lambda = lambda;
      Sigma = sigma;
      Name = name;
   }
}

public class Grid
{
   private double tStart;
   private double tEnd;
   private int tSteps;
   private double tStep;

   public List<PointRZ> Nodes { get; set; }
   public double[] Time { get; }
   public int[][] Elements { get; }
   public Material[] Materials { get; }
   public HashSet<(int, int[])>[] Boundaries { get; set; }

   public Grid(string spaceGridPath, string timeGridPath)
   {
      (int[], string)[] boundaries;
      using (StreamReader sr = new (spaceGridPath))
      {
         string[] data;
         data = sr.ReadLine()!.Split().ToArray();
         Materials = new Material[Convert.ToInt32(data[0])];

         for (int i = 0; i < Materials.Length; i++)
         {
            data = sr.ReadLine()!.Split(" ").ToArray();
            Materials[i] = new Material(Convert.ToDouble(data[0]), Convert.ToDouble(data[1]), data[2]);
         }

         data = sr.ReadLine()!.Split().ToArray();
         var nodesCount = Convert.ToInt32(data[0]);
         Nodes = new List<PointRZ>();
         for (int i = 0; i < nodesCount; i++)
         {
            data = sr.ReadLine()!.Split(" ").ToArray();
            Nodes.Add(new PointRZ(Convert.ToDouble(data[0]), Convert.ToDouble(data[1])));
         }

         data = sr.ReadLine()!.Split().ToArray();
         Elements = new int[Convert.ToInt32(data[0])].Select(_ => new int[7]).ToArray();
         for (int i = 0; i < Elements.Length; i++)
         {
            data = sr.ReadLine()!.Split(" ").ToArray();
            Elements[i][0] = Convert.ToInt32(data[0]);
            Elements[i][1] = Convert.ToInt32(data[1]);
            Elements[i][2] = Convert.ToInt32(data[2]);
            Elements[i][6] = Convert.ToInt32(data[3]);
         }

         data = sr.ReadLine()!.Split().ToArray();
         boundaries = new (int[], string)[int.Parse(data[0])];
         for (int i = 0; i < boundaries.Length; i++)
         {
            boundaries[i].Item1 = new int[3];
            data = sr.ReadLine()!.Split(" ").ToArray();
            boundaries[i].Item1[0] = Convert.ToInt32(data[0]);
            boundaries[i].Item1[2] = Convert.ToInt32(data[1]);
            boundaries[i].Item2 = data[2];
         }
      }

      using (StreamReader sr = new(timeGridPath))
      {
         string[] data;

         data = sr.ReadLine()!.Split(" ").ToArray();
         tStart = Convert.ToDouble(data[0]);
         tEnd = Convert.ToDouble(data[1]);
         tSteps = Convert.ToInt32(data[2]);
      }

      tStep = (tEnd - tStart) / tSteps;
      Time = new double[tSteps + 1];

      Time[0] = tStart;
      Time[^1] = tEnd;

      for (int i = 1; i < tSteps; i++) 
      {
         Time[i] = tStart + i * tStep;
      }

      Boundaries = new HashSet<(int, int[])>[4].Select(_ => new HashSet<(int, int[])>()).ToArray();

      NumberNodes(boundaries);
   }

   private void NumberNodes((int[], string)[] boudnaries)
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
         key[0] = (Math.Min(edge.Item1[0], edge.Item1[2]), Math.Max(edge.Item1[0], edge.Item1[2]));

         if (hashtable.ContainsKey(key[0]))
         {
            pointInfo = hashtable[key[0]];
            edge.Item1[1] = pointInfo.Item1;
            Boundaries[getType(edge.Item2)].Add((pointInfo.Item2, edge.Item1));
         }
         else
         {
            throw new Exception("Такая грань не была описана в элементах");
         }
      }

      int getType(string type)
      {
         return type switch
         {
            "Bottom" => 0,
            "Right" => 1,
            "Top" => 2,
            "Left" => 3,
            _ => 0,
         };
      }
   }
}