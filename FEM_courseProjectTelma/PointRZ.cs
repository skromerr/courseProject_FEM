namespace CourseProject;

public class PointRZ
{
   public double R { get; set; }
   public double Z { get; set; }

   public PointRZ()
   {
      R = 0;
      Z = 0;
   }
   public PointRZ(double r, double z)
   {
      R = r;
      Z = z;
   }

   public static PointRZ operator +(PointRZ a, PointRZ b) => new(a.R + b.R, a.Z + b.Z);

   public static PointRZ operator -(PointRZ a, PointRZ b) => new(a.R - b.R, a.Z - b.Z);

   public static PointRZ operator *(double coef, PointRZ a) => new(coef * a.R, coef * a.Z);

   public static PointRZ operator /(PointRZ a, double coef) => new(a.R / coef, a.Z / coef);

   public override string ToString() => R.ToString() + ' ' + Z.ToString();
}
