using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProject;

public class GridConverter
{
    public string PathIn { get; init; }
    public string PathOut { get; init; }

    public GridConverter(string pathIn, string pathOut)
    {
        PathIn = pathIn;
        PathOut = pathOut;
    }

    public void ConvertGrid()
    {
        string[] data;
        Point2D[] coords;
        using StreamReader sr = new StreamReader(PathIn);

        sr.ReadLine();
        data = sr.ReadLine().Split(" ");
        var numCoords = Convert.ToInt32(data[0]);
        coords = new Point2D[numCoords];

        for (int i = 0; i < numCoords; i++)
        {
            data = sr.ReadLine().Split("\t").Where(str => str != "").ToArray();
            coords[i] = new Point2D(Convert.ToDouble(data[0]), Convert.ToDouble(data[1]));
        }

        List<int[]> Elements = new List<int[]>();
        List<int[]> Boundaries = new List<int[]>();

        data = sr.ReadLine().Split(" ");
        var numObj = Convert.ToInt32(data[0]);
        for (int i = 0; i < numObj; i++)
        {
            data = sr.ReadLine().Split(" ").Where(str => str != "").ToArray();
            if (data[0] == "Triangle")
            {
                Elements.Add([Convert.ToInt32(data[5]), Convert.ToInt32(data[6]), Convert.ToInt32(data[7]), Convert.ToInt32(data[3])]);
            }
            else if (data[0] == "Segment")
            {
                Boundaries.Add([Convert.ToInt32(data[5]), Convert.ToInt32(data[6]), Convert.ToInt32(data[3])]);
            }
        }

        data = sr.ReadLine().Split(" ");
        var numMaterials = Convert.ToInt32(data[0]);
        List<string> MaterialsName = new List<string>();
        Dictionary<int, int> MaterialNumber = new Dictionary<int, int>();

        for (int i = 0; i < numMaterials;i++)
        {
            data = sr.ReadLine().Split(" ").Where(str => str != "").ToArray();
            MaterialsName.Add(data[1]);
            MaterialNumber.Add(Convert.ToInt32(data[0]), i);
        }

        data = sr.ReadLine().Split(" ");
        var numBoundaryTypes = Convert.ToInt32(data[0]);
        Dictionary<int, string> BoundTypes = new Dictionary<int, string>();
        for (int i = 0; i < numBoundaryTypes;i++)
        {
            data = sr.ReadLine().Split(" ").Where(str => str != "").ToArray();

            BoundTypes.Add(i, data[1]);
        }
        sr.Close();
        using StreamWriter sw = new StreamWriter(PathOut);

        sw.WriteLine(numMaterials);
        for (int i = 0; i < numMaterials; i++)
        {
            sw.WriteLine("ВпишиJ " + MaterialsName[i]);
        }

        sw.WriteLine(numCoords);
        for (int i = 0; i < numCoords; i++)
            sw.WriteLine($"{coords[i].X} {coords[i].Y}");

        sw.WriteLine(Elements.Count);
        for (int i = 0; i < Elements.Count; i++)
            sw.WriteLine($"{Elements[i][0]} {Elements[i][1]} {Elements[i][2]} {MaterialNumber[Elements[i][3]]}");

        sw.WriteLine(Boundaries.Count);
        for (int i = 0; i < Boundaries.Count; i++)
            sw.WriteLine($"{Boundaries[i][0]} {Boundaries[i][1]} {BoundTypes[Boundaries[i][2]]}");
        sw.Close();
    }
}
