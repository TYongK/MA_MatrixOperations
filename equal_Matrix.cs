using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(DataTree<double> Tree_Matrix_A, DataTree<double> Tree_Matrix_B, double epsilon, ref object A)
  {


    A = MatrixAreEqual_Grasshopper(Tree_Matrix_A, Tree_Matrix_B, epsilon);


  }

  // <Custom additional code> 
  //returns "string" or "double[][]"
  static dynamic MatrixAreEqual_Original(double[][] matrixA, double[][] matrixB, double epsilon)
  {
    // true if all values in matrixA == values in matrixB
    int aRows = matrixA.Length;
    int aCols = matrixA[0].Length;

    int bRows = matrixB.Length;
    int bCols = matrixB[0].Length;

    if (aRows != bRows || aCols != bCols)
    {
      return "Non-conformable matrices";
    }

    for (int i = 0; i < aRows; ++i) // each row of A and B
    {
      for (int j = 0; j < aCols; ++j) // each col of A and B
      {
        //if (matrixA[i][j] != matrixB[i][j])

        if (Math.Abs(matrixA[i][j] - matrixB[i][j]) > epsilon)
        {
          return false;
        }
      }
    }
    return true;
  }

  //returns "string" or "DataTree<double>"
  static dynamic MatrixAreEqual_Grasshopper(DataTree<double> matrixA, DataTree<double> matrixB, double epsilon)
  {
    // true if all values in matrixA == values in matrixB
    IList<List<double>> matrixA_Branches = matrixA.Branches;
    int aRows = matrixA_Branches.Count; int aCols = matrixA_Branches[0].Count;

    IList<List<double>> matrixB_Branches = matrixB.Branches;
    int bRows = matrixB_Branches.Count; int bCols = matrixB_Branches[0].Count;

    if (aRows != bRows || aCols != bCols)
    {
      return "Non-conformable matrices";
    }

    for (int i = 0; i < aRows; ++i) // each row of A and B
    {
      for (int j = 0; j < aCols; ++j) // each col of A and B
      {
        //if (matrixA[i][j] != matrixB[i][j])
        if (Math.Abs(matrixA_Branches[i][j] - matrixB_Branches[i][j]) > epsilon)
        {
          return false;
        }
      }
    }

    return true;
  }
  // Resource: https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
  // Resource: https://www.storyofmathematics.com/equivalent-matrices
  // Editor: Taeyong Kim
  // </Custom additional code> 
}
