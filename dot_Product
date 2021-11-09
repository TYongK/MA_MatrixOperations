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
  private void RunScript(DataTree<double> Tree_axb_Matrix, DataTree<double> Tree_bxa_Matrix, ref object A)
  {

    if(MatrixProduct_Grasshopper(Tree_axb_Matrix, Tree_bxa_Matrix).GetType() == typeof(string))
    {
      A = "Non-conformable matrices";
    }
    else
    {
      double[][] m = MatrixProduct_Grasshopper(Tree_axb_Matrix, Tree_bxa_Matrix);

      //////Define Tree//////
      DataTree < double > doubleTree = new DataTree<double>();

      // Raws = Branches
      for(int i = 0; i < m.Length;i++)
      {
        GH_Path pth = new GH_Path(i);

        // Cols = Index
        for (int j = 0; j < m[i].Length;j++)
        {
          doubleTree.Add(m[i][j], pth);
        }
      }
      //////Define Tree End//////

      A = doubleTree;
    }


  }

  // <Custom additional code> 
  static double[][] MatrixCreate(int rows, int cols)
  {
    double[][] result = new double[rows][];

    for (int i = 0; i < rows; ++i)
    {
      result[i] = new double[cols];
    }
    return result;
  }

  //returns "string" or "double[][]"
  static dynamic MatrixProduct_Original(double[][] matrixA, double[][] matrixB)
  {
    int aRows = matrixA.Length; int aCols = matrixA[0].Length;
    int bRows = matrixB.Length; int bCols = matrixB[0].Length;
    if (aCols != bRows)
    {
      return "Non-conformable matrices";
    }

    double[][] result = MatrixCreate(aRows, bCols);

    for (int i = 0; i < aRows; ++i) // each row of A
      for (int j = 0; j < bCols; ++j) // each col of B
        for (int k = 0; k < aCols; ++k) // could use k less-than bRows
          result[i][j] += matrixA[i][k] * matrixB[k][j];

    return result;
  }

  //returns "string" or "DataTree<double>"
  static dynamic MatrixProduct_Grasshopper(DataTree<double> matrixA, DataTree<double> matrixB)
  {
    IList<List<double>> matrixA_Branches = matrixA.Branches;
    int aRows = matrixA_Branches.Count; int aCols = matrixA_Branches[0].Count;

    IList<List<double>> matrixB_Branches = matrixB.Branches;
    int bRows = matrixB_Branches.Count; int bCols = matrixB_Branches[0].Count;

    if (aCols != bRows)
    {
      return "Non-conformable matrices";
    }

    double[][] result = MatrixCreate(aRows, bCols);

    for (int i = 0; i < aRows; ++i) // each row of A
    {
      for (int j = 0; j < bCols; ++j) // each col of B
      {
        for (int k = 0; k < aCols; ++k) // could use k less-than bRows
        {
          result[i][j] += matrixA_Branches[i][k] * matrixB_Branches[k][j];
        }
      }
    }

    return result;
  }
  // Resource#1: https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
  // Resource#2: https://en.wikipedia.org/wiki/Matrix_multiplication
  // Editor: Taeyong Kim.
  // </Custom additional code> 
}
