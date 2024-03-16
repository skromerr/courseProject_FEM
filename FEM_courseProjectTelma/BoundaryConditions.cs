namespace CourseProject;

public class FirstCondition
{
   public PointRZ point { get; }
   public int NodeNumber { get; }

   public FirstCondition(PointRZ node, int nodeNumber)
   {
      point = node;
      NodeNumber = nodeNumber;
   }
}

public class SecondCondition
{
   public int[] Edge { get; }
   public int ElemNumber { get; }
   public int EdgeType { get; }   // 0 - bottom, 1 - right
                                       // 2 - top, 3 - left

   public SecondCondition(int elemNumber, int edgeType, int[] edge)
   {
      ElemNumber = elemNumber;
      EdgeType = edgeType;
      Edge = edge;
   }
}
