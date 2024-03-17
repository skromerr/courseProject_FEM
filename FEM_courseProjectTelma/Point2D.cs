namespace CourseProject;

public class Point2D
{
   public double X { get; set; }
   public double Y { get; set; }

   public Point2D()
   {
      X = 0;
      Y = 0;
   }
   public Point2D(double x, double y)
   {
      X = x;
      Y = y;
   }

   public static Point2D operator +(Point2D a, Point2D b) => new(a.X + b.X, a.Y + b.Y);

   public static Point2D operator -(Point2D a, Point2D b) => new(a.X - b.X, a.Y - b.Y);

   public static Point2D operator *(double coef, Point2D a) => new(coef * a.X, coef * a.Y);

   public static Point2D operator /(Point2D a, double coef) => new(a.X / coef, a.Y / coef);

   public override string ToString() => X.ToString() + ' ' + Y.ToString();
}
