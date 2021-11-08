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
  private void RunScript(DataTree<double> Tree_axb_Matrix, List<double> Array_1xb_Vector, ref object A)
  {

    if(MatrixVectorProduct_Grasshopper(Tree_axb_Matrix, Array_1xb_Vector.ToArray()).GetType() == typeof(string))
    {
      A = "Non-conformable matrices";
    }
    else
    {
      A = MatrixVectorProduct_Grasshopper(Tree_axb_Matrix, Array_1xb_Vector.ToArray());
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
  static dynamic MatrixVectorProduct_Original(double[][] matrix, double[] vector)
  {
    // result of multiplying an n x m matrix by a m x 1
    // column vector (yielding an n x 1 column vector)
    int mRows = matrix.Length; int mCols = matrix[0].Length;
    int vRows = vector.Length;
    if (mCols != vRows)
    {
      return "Non-conformable matrices";
    }
    double[] result = new double[mRows];
    {
      for (int i = 0; i < mRows; ++i)
      {
        for (int j = 0; j < mCols; ++j)
        {
          result[i] += matrix[i][j] * vector[j];
        }
      }
    }
    return result;
  }

  //returns "string" or "DataTree<double>"
  static dynamic MatrixVectorProduct_Grasshopper(DataTree<double> matrix, double[] vector)
  {
    // result of multiplying an n x m matrix by a m x 1
    // column vector (yielding an n x 1 column vector)
    IList<List<double>> matrix_Branches = matrix.Branches;
    int mRows = matrix_Branches.Count; int mCols = matrix_Branches[0].Count;

    int vRows = vector.Length;


    if (mCols != vRows)
    {
      return "Non-conformable matrices";
    }
    double[] result = new double[mRows];
    {
      for (int i = 0; i < mRows; ++i)
      {
        for (int j = 0; j < mCols; ++j)
        {
          result[i] += matrix_Branches[i][j] * vector[j];
        }
      }
    }
    return result;
  }

  // Resource#1: https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
  // Resource#2: https://mathinsight.org/matrix_vector_multiplication
  // Editor: Taeyong Kim.
  // </Custom additional code> 
}
