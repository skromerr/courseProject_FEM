namespace CourseProject;

public interface ITest
{
   public double U(PointRZ point, double t);

   public double F(PointRZ point, double t, double sigma, double lambda);

   public double Theta(PointRZ point, double t, int funcNum);
}

public class Test1 : ITest
{
   public double U(PointRZ point, double t)
      => point.R * point.Z * t;

   public double F(PointRZ point, double t, double sigma, double lambda)
      => sigma * point.R * point.Z - lambda * t * point.Z / point.R;

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return -point.R * t;
         case 1:
            return point.Z * t;
         case 2:
            return point.R * t;
         case 3:
            return -point.Z * t;
         default:
            return 0;
      }
   }
}

public class Test2 : ITest
{
   public double U(PointRZ point, double t)
      => point.R * point.Z * t * t;

   public double F(PointRZ point, double t, double sigma, double lambda)
      => 2 * sigma * t * point.R * point.Z - lambda * t * t * point.Z / point.R;

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return -point.R * t * t;
         case 1:
            return point.Z * t * t;
         case 2:
            return point.R * t * t;
         case 3:
            return -point.Z * t * t;
         default:
            return 0;
      }
   }
}

public class Test3 : ITest
{
   public double U(PointRZ point, double t)
      => point.R * point.Z * t * t * t ;

   public double F(PointRZ point, double t, double sigma, double lambda)
      => 3 * sigma * t * t * point.R * point.Z - lambda * t * t * t * point.Z / point.R;

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return -point.R * t * t * t;
         case 1:
            return point.Z * t * t * t;
         case 2:
            return point.R * t * t * t;
         case 3:
            return -point.Z * t * t * t;
         default:
            return 0;
      }
   }
}

public class Test4 : ITest
{
   public double U(PointRZ point, double t)
      => Math.Pow(t, 3);

   public double F(PointRZ point, double t, double sigma, double lambda)
      => 3 * sigma * t * t;

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return 0;
         case 1:
            return 0;
         case 2:
            return 0;
         case 3:
            return 0;
         default:
            return 0;
      }
   }
}

public class Test5 : ITest
{
   public double U(PointRZ point, double t)
      => Math.Sin(t);

   public double F(PointRZ point, double t, double sigma, double lambda)
      => sigma * Math.Cos(t);

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return 0;
         case 1:
            return 0;
         case 2:
            return 0;
         case 3:
            return 0;
         default:
            return 0;
      }
   }
}

public class Test6 : ITest
{
   public double U(PointRZ point, double t)
      => Math.Pow(point.R, 3);

   public double F(PointRZ point, double t, double sigma, double lambda)
      => -lambda * 9 * point.R;

   public double Theta(PointRZ point, double t, int funcNum)
   {
      switch (funcNum)  // 0 - bottom, 1 - right, 2 - top, 3 - left
      {
         case 0:
            return -3 * point.R * point.R;
         case 1:
            return 0;
         case 2:
            return 3 * point.R * point.R;
         case 3:
            return 0;
         default:
            return 0;
      }
   }
}